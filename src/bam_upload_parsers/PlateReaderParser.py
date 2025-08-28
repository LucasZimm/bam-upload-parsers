# -----------------------------------------------------------------------------
# File:        parsers.py
# Original Author:  Bastian Ruehle
# Maintainer:       Lucas Zimmermann
# Email:            lucas.zimmermann@bam.de
# Version:          2.0.0 not sure
# Created:          2024 (by Bastian Ruehle and Ingo Breßler)
# Modified:         2025 (by Lucas Zimmermann)
#
# Copyright (c) 2024, Bastian Ruehle,
# Federal Institute for Materials Research and Testing (BAM)
# Copyright (c) 2025, BAM
#


import codecs
import json
import math
import os.path
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union

import h5py
import nmrglue as ng
import numpy as np

# from bam_masterdata.datamodel.object_types import
from bam_masterdata.logger import logger
from bam_masterdata.metadata.entities import CollectionType
from bam_masterdata.parsing import AbstractParser
from structlog._config import BoundLoggerLazyProxy

from .utils import metadata_to_instance


class PlateReaderParser(AbstractParser):
    def parse(
        self, files: list[str], collection: CollectionType, logger: BoundLoggerLazyProxy
    ) -> None:
        """
        Parse ASCII files exported from Spectramax measurements.

        Parameters
        ----------
        files: List
            The list of files that are being parsed.
        collection: CollectionType
            Collection for adding the Objects.
        logger:
            Logger for logging errors, warnings and notes.

        Returns
        -------
        None
        """
        for file in files:
            if file.endswith(".txt"):
                break
        else:
            return {}

        metadata: dict[str, str | float | int] = {}

        fc = [
            i.replace("\r", "").replace("\n", "").split("\t")
            for i in codecs.open(file, "r", "utf-16-le")
        ]
        tmp = 0
        x_data: (
            list[float]
            | list[str]
            | list[list[float]]
            | list[list[str]]
            | list[np.ndarray]
        ) = []
        y_data: list[float] | list[list[float]] | list[np.ndarray] = []
        legends = []
        titles = []

        for i, l in enumerate(fc):
            if isinstance(l, list) and len(l) > 0 and l[0] == "Plate:":
                tmp += 1

                if l[5] == "Absorbance":
                    offset = 0
                elif l[5] == "Fluorescence":
                    offset = 1

                metadata["PLATEREADER.READ_TYPE"] = l[4]  # Endpoint, Spectrum
                metadata["PLATEREADER.READ_MODE"] = l[
                    5
                ]  # Absorbance, Fluorescence, Luminsecence
                metadata["PLATEREADER.NUMBER_OF_WELLS"] = l[18 + offset]

                if l[4] == "Spectrum":
                    if l[5] == "Absorbance":
                        metadata["PLATEREADER.ABSORBANCE_START_WAVELENGTH"] = float(
                            l[11]
                        )
                        metadata["PLATEREADER.ABSORBANCE_END_WAVELENGTH"] = float(l[12])
                        metadata["PLATEREADER.ABSORBANCE_RESOLUTION"] = float(l[13])
                        titles.append("Absorbance Spectrum\n")
                    elif l[5] == "Fluorescence":
                        metadata["PLATEREADER.SPECTRUM_TYPE"] = l[23]
                        if l[23] == "Excitation Sweep":
                            metadata["PLATEREADER.EMISSION_WAVELENGTH"] = float(l[24])
                            metadata["PLATEREADER.EXCITATION_START_WAVELENGTH"] = float(
                                l[12]
                            )
                            metadata["PLATEREADER.EXCITATION_END_WAVELENGTH"] = float(
                                l[13]
                            )
                            metadata["PLATEREADER.EXCITATION_RESOLUTION"] = float(l[14])
                            titles.append(
                                f"Excitation Spectrum\nEm.: {float(l[24])} nm"
                            )
                        elif l[23] == "Emission Sweep":
                            metadata["PLATEREADER.EXCITATION_WAVELENGTH"] = float(l[24])
                            metadata["PLATEREADER.EMISSION_START_WAVELENGTH"] = float(
                                l[12]
                            )
                            metadata["PLATEREADER.EMISSION_END_WAVELENGTH"] = float(
                                l[13]
                            )
                            metadata["PLATEREADER.EMISSION_RESOLUTION"] = float(l[14])
                            titles.append(f"Emission Spectrum\nEx.: {float(l[24])} nm")
                    j = i + 1
                    tmp_x = []
                    tmp_y = []
                    legends.append(
                        [fc[j][k] for k in range(2, len(fc[j])) if fc[j + 1][k] != ""]
                    )
                    j += 1
                    while fc[j][0] != "":
                        tmp_x.append(float(fc[j][0]))
                        tmp_y.append(
                            np.array(
                                [float(k) for k in fc[j][2:] if k != ""],
                                dtype="float32",
                            )
                        )
                        j += 1
                    x_data.append(np.array(tmp_x))
                    y_data.append(np.array(tmp_y))
                elif l[4] == "Endpoint":
                    if l[5] == "Absorbance":
                        metadata["PLATEREADER.ABSORBANCE_WAVELENGTH"] = float(l[15])
                        titles.append(
                            f"Absorbance Measurement\nAbs.: {float(l[15])} nm"
                        )
                    elif l[5] == "Fluorescence":
                        metadata["PLATEREADER.EMISSION_WAVELENGTH"] = float(l[16])
                        metadata["PLATEREADER.EXCITATION_WAVELENGTH"] = float(l[20])
                        titles.append(
                            f"Fluorescence Measurement\nEx.: {float(l[20])} nm, Em.: {float(l[16])} nm"
                        )
                    j = i + 1
                    x_data.append(
                        [
                            fc[j][k]
                            for k in range(2, len(fc[j]))
                            if fc[j + 1][k] != "" and fc[j + 1][k] != "#SAT"
                        ]
                    )
                    j += 1
                    y_data.append(
                        np.array(
                            [float(k) for k in fc[j][2:] if k != "" and k != "#SAT"],
                            dtype="float32",
                        )
                    )
                    legends.append([])
        # TODO find platereader class or add in bam masterdata
        # platereader = metadata_to_instance(metadata, PlateReader())
        # collection.add(platereader)
