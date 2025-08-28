# -----------------------------------------------------------------------------
# File:        NMRParser.py
# Original Author:  Bastian Ruehle
# Maintainer:       Bastian Ruehle, Jose M. Pizarro, Lucas Zimmermann
# Email:            bastian.ruehle@bam.de, jose.pizarro-blanco@bam.de, lucas.zimmermann@bam.de
# Created:          2024 (by Bastian Ruehle and Ingo Breßler)
# Modified:         2025 (by Lucas Zimmermann)
#
# Copyright (c) 2024, Bastian Ruehle,
# Federal Institute for Materials Research and Testing (BAM)
# Copyright (c) 2025, BAMResearch


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


class NMRParser(AbstractParser):
    def parse(
        self, files: list[str], collection: CollectionType, logger: BoundLoggerLazyProxy
    ) -> None:
        """
        Parse jcamp dx NMR files.

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
            if not (
                file.endswith(".jdx")
                or file.endswith(".jx")
                or file.endswith(".jcamp")
                or file.endswith(".fid")
            ):
                continue

            PROPERTY_TYPE_CODE_DICT_OI_XPULSE = {
                # 'TITLE': '$NAME',
                "LONGDATE": "START_DATE",
                "SPECTROMETERDATASYSTEM": "NMR.INSTRUMENT",
                ".OBSERVENUCLEUS": "NMR.NUCLEUS_DIRECT",
                "$INDIRECTNUCLEUS": "NMR.NUCLEUS_INDIRECT",
                ".SOLVENTNAME": "NMR.SOLVENT",
                ".OBSERVEFREQUENCY": "NMR.FREQUENCY",
                ".EXPERIMENT": "NMR.EXPERIMENT",
                "$NS": "NMR.SCANS",
                # '': 'NMR.START_CHEMICAL_SHIFT',
                # '': 'NMR.END_CHEMICAL_SHIFT',
                # '': 'NMR.IS_QNMR',
                # '': 'NMR.PULSE_ANGLE',
                "$RELAXATIONDELAY": "NMR.INTERPULSE_DELAY",
                ".ACQUISITIONTIME": "NMR.ACQUISITION_TIME",
            }

            PROPERTY_TYPE_CODE_DICT_VARIAN = {
                # 'samplename': '$NAME',
                "time_run": "START_DATE",
                "systemname_": "NMR.INSTRUMENT",
                "tn": "NMR.NUCLEUS_DIRECT",
                # '$INDIRECTNUCLEUS': 'NMR.NUCLEUS_INDIRECT',
                "solvent": "NMR.SOLVENT",
                "sfrq": "NMR.FREQUENCY",
                "apptype": "NMR.EXPERIMENT",
                "nt": "NMR.SCANS",
                # '': 'NMR.START_CHEMICAL_SHIFT',  # obsSW
                # '': 'NMR.END_CHEMICAL_SHIFT',  # obsSW
                # '': 'NMR.IS_QNMR',
                # '': 'NMR.PULSE_ANGLE', # pw/pw90 * 90
                "d1": "NMR.INTERPULSE_DELAY",
                "ACQtime": "NMR.ACQUISITION_TIME",
            }
            metadata: dict[str, str | float | int] = {}

            if (
                file.lower().endswith(".jdx")
                or file.lower().endswith(".dx")
                or file.lower().endswith(".jcamp")
            ):
                dic, raw_data = ng.jcampdx.read(file)
                # create complex array
                data = np.empty((raw_data[0].shape[0],), dtype="complex128")
                data.real = raw_data[0][:]
                data.imag = raw_data[1][:]
                # process
                data = data[:8192]  # reduce size to make lp algorithm feasible
                data = ng.proc_lp.lp(
                    data=data, pred=8, mode="b", order=15, append="before"
                )
                data = ng.proc_base.zf_auto(data)
                data = ng.proc_base.fft(data)
                data = ng.process.proc_autophase.autops(
                    data=data, fn="acme", disp=False
                )  # disp=False turns off convergence messages
                data = ng.proc_base.di(data)
                # data = ng.process.proc_bl.baseline_corrector(data)
                data -= data.min()
                data /= data.max() / 100
                # determine the ppm scale
                udic = ng.jcampdx.guess_udic(dic, data)
                udic[0]["car"] = 0
                udic[0]["sw"] = float(dic["$SWEEPWIDTH"][0])
                uc = ng.fileiobase.uc_from_udic(udic)
                ppm_scale = uc.ppm_scale()
                thresh = 5 * (np.max(data[-100:]) - np.mean(data[-100:])) + np.mean(
                    data[-100:]
                )  # Guess as 5 times the noise in the last 100 points
                peaks = ng.peakpick.pick(data, thresh)
                peak_heights = sorted(data[[int(i[0]) for i in peaks]])
                offset = 7.29 - ppm_scale[int(peaks[0][0])]
                ppm_scale += offset
            else:
                dic, data = ng.varian.read(file)
                # process
                if len(data.shape) > 1:
                    data = np.mean(data, axis=0)
                data = ng.proc_base.zf_auto(data)
                data = ng.proc_base.fft(data)
                data = ng.process.proc_autophase.autops(
                    data=data, fn="acme", disp=False
                )  # disp=False turns off convergence messages
                data = ng.proc_base.di(data)
                # data = ng.process.proc_bl.baseline_corrector(data)
                data -= data.min()
                data /= data.max() / 100
                # determine the ppm scale
                sw = dic["procpar"]["obsSW"]["values"][0]
                sw = sw[sw.find("(") + 1 : sw.find(")")].split(",")
                ppm_scale = np.linspace(float(sw[0]), float(sw[1]), len(data))
                thresh = 5 * (np.max(data[-100:]) - np.mean(data[-100:])) + np.mean(
                    data[-100:]
                )  # Guess as 5 times the noise in the last 100 points
                peaks = ng.peakpick.pick(data, thresh)
                peak_heights = sorted(data[[int(i[0]) for i in peaks]])

            # Oxford Instruments X-Pulse
            if "SPECTROMETERDATASYSTEM" in dic.keys() and dic["SPECTROMETERDATASYSTEM"][
                0
            ].startswith("Oxford Instruments X-Pulse"):
                for k in PROPERTY_TYPE_CODE_DICT_OI_XPULSE.keys():
                    if k not in dic.keys():
                        continue
                    if k in (".OBSERVENUCLEUS", "$INDIRECTNUCLEUS"):
                        metadata[PROPERTY_TYPE_CODE_DICT_OI_XPULSE[k]] = dic[k][
                            0
                        ].replace("^", "")
                    elif k in (
                        ".OBSERVEFREQUENCY",
                        "$RELAXATIONDELAY",
                        ".ACQUISITIONTIME",
                    ):
                        metadata[PROPERTY_TYPE_CODE_DICT_OI_XPULSE[k]] = float(
                            dic[k][0]
                        )
                    elif k in ("$NS",):
                        metadata[PROPERTY_TYPE_CODE_DICT_OI_XPULSE[k]] = int(dic[k][0])
                    elif k in ("TITLE",):
                        metadata[PROPERTY_TYPE_CODE_DICT_OI_XPULSE[k]] = dic[k][0][
                            : dic[k][0].find("\n")
                        ]
                    elif k == "LONGDATE":
                        metadata[PROPERTY_TYPE_CODE_DICT_OI_XPULSE[k]] = (
                            datetime.strftime(
                                datetime.strptime(
                                    dic[k][0].replace("T", " "), "%d/%m/%Y %H:%M:%S"
                                ),
                                "%Y-%m-%d %H:%M:%S",
                            )
                        )
                    else:
                        metadata[PROPERTY_TYPE_CODE_DICT_OI_XPULSE[k]] = dic[k][0]

                if "$SWEEPWIDTH" in dic.keys() and ".OBSERVEFREQUENCY" in dic.keys():
                    sweep = round(
                        float(dic["$SWEEPWIDTH"][0])
                        / float(dic[".OBSERVEFREQUENCY"][0])
                        / 2,
                        2,
                    )
                    metadata["NMR.START_CHEMICAL_SHIFT"] = -sweep
                    metadata["NMR.END_CHEMICAL_SHIFT"] = sweep
            # Varian
            elif "procpar" in dic.keys() and dic["procpar"]["systemname_"]["values"][
                0
            ].startswith("NMR500-vnmrs500"):
                dic = dic["procpar"]
                for k in PROPERTY_TYPE_CODE_DICT_VARIAN.keys():
                    if k not in dic.keys():
                        continue
                    if k in ("sfrq", "d1", "ACQtime"):
                        metadata[PROPERTY_TYPE_CODE_DICT_VARIAN[k]] = float(
                            dic[k]["values"][0]
                        )
                    elif k in ("nt",):
                        metadata[PROPERTY_TYPE_CODE_DICT_VARIAN[k]] = int(
                            dic[k]["values"][0]
                        )
                    elif k == "time_run":
                        metadata[PROPERTY_TYPE_CODE_DICT_VARIAN[k]] = datetime.strftime(
                            datetime.strptime(
                                dic[k]["values"][0].replace("T", " "), "%Y%m%d %H%M%S"
                            ),
                            "%Y-%m-%d %H:%M:%S",
                        )
                    else:
                        metadata[PROPERTY_TYPE_CODE_DICT_VARIAN[k]] = dic[k]["values"][
                            0
                        ]

                if "obsSW" in dic.keys():
                    metadata["NMR.START_CHEMICAL_SHIFT"] = sw[1]
                    metadata["NMR.END_CHEMICAL_SHIFT"] = sw[0]
                if "pw" in dic.keys() and "pw90" in dic.keys():
                    metadata["NMR.PULSE_ANGLE"] = round(
                        float(dic["pw"]["values"][0])
                        / float(dic["pw90"]["values"][0])
                        * 90,
                        2,
                    )

            # Check Controlled Vocabularies
            if "NMR.EXPERIMENT" in metadata.keys():
                if (
                    "1D EXPERIMENT" in metadata["NMR.EXPERIMENT"].upper()
                    or "STD1D" in metadata["NMR.EXPERIMENT"].upper()
                ):
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_STANDARD"
                elif "HSQC" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_HSQC"
                elif "HMQC" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_HMQC"
                elif "HMBC" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_HMBC"
                elif "COSY" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_COSY"
                elif "TOCSY" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_TOCSY"
                elif "NOESY" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_NOESY"
                elif "ROESY" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_ROESY"
                elif "DEPT" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_DEPT"
                elif "INVERSION" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_DEPT"
                elif "INVERSION" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_IR"
                elif "ECHO" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_SPIN_ECHO"
                elif "MAS" in metadata["NMR.EXPERIMENT"].upper():
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_CP_MAS"
                else:
                    metadata["NOTES"] = f"Experiment: {metadata['NMR.EXPERIMENT']}\n"
                    metadata["NMR.EXPERIMENT"] = "NMR_EXP_OTHER"

            if "NMR.SOLVENT" in metadata.keys():
                if (
                    "CDCl3" in metadata["NMR.SOLVENT"].upper()
                    or "chloroform" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_CDCL3"
                elif (
                    "CD2Cl2" in metadata["NMR.SOLVENT"].upper()
                    or "dcm" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_CD2CL2"
                elif (
                    "D2O" in metadata["NMR.SOLVENT"].upper()
                    or "water" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_D2O"
                elif (
                    "CD3OD" in metadata["NMR.SOLVENT"].upper()
                    or "methanol" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_CD3OD"
                elif "dmso" in metadata["NMR.SOLVENT"].lower():
                    metadata["NMR.SOLVENT"] = "NMR_SOL_DMSO_D6"
                elif "aceton" in metadata["NMR.SOLVENT"].lower():
                    metadata["NMR.SOLVENT"] = "NMR_SOL_ACETONE_D6"
                elif (
                    "CD3CN" in metadata["NMR.SOLVENT"].upper()
                    or "acetonitril" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_CD3CN"
                elif (
                    "THF" in metadata["NMR.SOLVENT"].upper()
                    or "tetrahydrofuran" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_THF-D8"
                elif "tolu" in metadata["NMR.SOLVENT"].lower():
                    metadata["NMR.SOLVENT"] = "NMR_SOL_TOLUENE_D8"
                elif (
                    "C6D5Cl" in metadata["NMR.SOLVENT"].upper()
                    or "chlorobenz" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_C6D5Cl"
                elif (
                    "C6D6" in metadata["NMR.SOLVENT"].upper()
                    or "benz" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_C6D6"
                elif (
                    "TFE" in metadata["NMR.SOLVENT"].upper()
                    or "trifluorethanol" in metadata["NMR.SOLVENT"].lower()
                ):
                    metadata["NMR.SOLVENT"] = "NMR_SOL_TFE_D3"
                else:
                    metadata["NOTES"] = f"Solvent: {metadata['NMR.SOLVENT']}\n"
                    metadata["NMR.SOLVENT"] = "NMR_SOL_OTHER"

            if "NMR.NUCLEUS_DIRECT" in metadata.keys():
                if metadata["NMR.NUCLEUS_DIRECT"][0] not in (
                    str(i) for i in range(0, 10)
                ):
                    pos = min(
                        j
                        for j in (
                            metadata["NMR.NUCLEUS_DIRECT"].find(str(i))
                            for i in range(0, 10)
                        )
                        if j != -1
                    )
                    metadata["NMR.NUCLEUS_DIRECT"] = (
                        metadata["NMR.NUCLEUS_DIRECT"][pos:]
                        + metadata["NMR.NUCLEUS_DIRECT"][:pos]
                    )
                metadata["NMR.NUCLEUS_DIRECT"] = (
                    f"NMR_NUC_{metadata['NMR.NUCLEUS_DIRECT']}"
                )

            if "NMR.NUCLEUS_INDIRECT" in metadata.keys():
                if metadata["NMR.NUCLEUS_INDIRECT"][0] not in (
                    str(i) for i in range(0, 10)
                ):
                    pos = min(
                        j
                        for j in (
                            metadata["NMR.NUCLEUS_INDIRECT"].find(str(i))
                            for i in range(0, 10)
                        )
                        if j != -1
                    )
                    metadata["NMR.NUCLEUS_INDIRECT"] = (
                        metadata["NMR.NUCLEUS_INDIRECT"][pos:]
                        + metadata["NMR.NUCLEUS_INDIRECT"][:pos]
                    )
                metadata["NMR.NUCLEUS_INDIRECT"] = (
                    f"NMR_NUC_{metadata['NMR.NUCLEUS_INDIRECT']}"
                )
            if isinstance(metadata["START_DATE"], str):
                metadata["START_DATE"] = datetime.strptime(
                    metadata["START_DATE"], "%Y-%m-%d %H:%M:%S"
                )
            nmr = metadata_to_instance(metadata, Nmr())
            collection.add(nmr)
