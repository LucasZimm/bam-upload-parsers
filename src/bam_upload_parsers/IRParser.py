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
from bam_masterdata.datamodel.object_types import Dls, Ftir, Nmr, Sem, Tem
from bam_masterdata.logger import logger
from bam_masterdata.metadata.entities import CollectionType
from bam_masterdata.parsing import AbstractParser
from structlog._config import BoundLoggerLazyProxy

from .utils import metadata_to_instance


class IRParser(AbstractParser):
    def parse(
        self,
        files: list[str],
        collection: CollectionType,
        logger,
    ) -> None:
        """
        Parse ASCII files exported from IR measurements.

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
        filename_spa_list = []
        filename_csv_list = []
        for i in files:
            if i.lower().endswith(".spa"):
                filename_spa_list.append(i)
            elif i.lower().endswith(".csv"):
                filename_csv_list.append(i)

        if len(filename_spa_list) == 0:
            return
        for filename_spa in filename_spa_list:
            metadata: dict[str, str | float | int] = {}
            fc = [i for i in open(filename_spa, "rb")]
            for i, l in enumerate(fc):
                if b"Background gemessen am " in l:
                    metadata["START_DATE"] = datetime.strftime(
                        datetime.strptime(
                            l[l.find(b" am  ") + 8 : l.find(b" (GMT")].decode(),
                            "%b %d %H:%M:%S %Y",
                        ),
                        "%Y-%m-%d %H:%M:%S",
                    )
                    tmp = (
                        fc[i + 2][fc[i + 2].find(b":\t ") + 3 : fc[i + 2].find(b"\r\n")]
                        .decode()
                        .replace(",", ".")
                        .split(" ")
                    )
                    metadata["FTIR.RESOLUTION"] = int(float(tmp[0]))
                    metadata["FTIR.START_WAVENUMBER"] = float(tmp[2])
                    metadata["FTIR.END_WAVENUMBER"] = float(tmp[4])
                    metadata["FTIR.INSTRUMENT"] = fc[i + 3][
                        fc[i + 3].find(b": ") + 2 : fc[i + 3].find(b"\r\n")
                    ].decode()
                    metadata["FTIR.ACCESSORY"] = "FTIR_ACCESSORY_GOLDEN_GATE"
            if isinstance(metadata["START_DATE"], str):
                metadata["START_DATE"] = datetime.strptime(
                    metadata["START_DATE"], "%Y-%m-%d %H:%M:%S"
                )
            ftir = metadata_to_instance(metadata, Ftir())
            collection.add(ftir)
