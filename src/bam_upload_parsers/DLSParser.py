# -----------------------------------------------------------------------------
# File:        DLSParser.py
# Original Author:  Bastian Ruehle
# Maintainer:       Bastian Ruehle, Jose M. Pizarro, Lucas Zimmermann
# Email:            bastian.ruehle@bam.de, jose.pizarro-blanco@bam.de, lucas.zimmermann@bam.de
# Created:          2024 (by Bastian Ruehle and Ingo Breßler)
# Modified:         2025 (by Lucas Zimmermann)
#
# Copyright (c) 2024, Bastian Ruehle,
# Federal Institute for Materials Research and Testing (BAM)
# Copyright (c) 2025, BAMResearch


from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from bam_masterdata.datamodel.object_types import Dls
from bam_masterdata.logger import logger
from bam_masterdata.metadata.entities import CollectionType
from bam_masterdata.parsing import AbstractParser

from .utils import metadata_to_instance


class DLSParser(AbstractParser):
    def parse(self, files, collection, logger):
        """
        Parse ASCII files exported from DLS measurements.

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
        csv_files = [file for file in files if file.endswith(".csv")]
        if not csv_files:
            return
        # Mapping of extracted values to their field names in OpenBIS
        for file in csv_files:
            PROPERTY_TYPE_CODE_DICT = {
                # 'Sample Name': '$NAME',
                "Measurement Start Date And Time/Size Measurement Result": "START_DATE",
                "Name/Material Settings/Material/Sample Settings/Sample Settings/Size Measurement Result": "DLS.MATERIAL",
                "Name/Dispersant Settings/Dispersant/Sample Settings/Sample Settings/Size Measurement Result": "DLS.DISPERSANT",
                "Analysis Model/Size Analysis Settings/Size Analysis/Size Measurement Settings/Measurement Settings/Size Measurement Result": "DLS.ANALYSISMODEL",
                "Cell Name/Cell Description Settings/View Settings/Cell Settings/Cell/Sample Settings/Sample Settings/Size Measurement Result": "DLS.CELLDESCRIPTION",
                "Description/Cell Description Settings/View Settings/Cell Settings/Cell/Sample Settings/Sample Settings/Size Measurement Result": "DLS.CELLDESCRIPTION",
                "Fka Model/Zeta Fka Parameter Settings/Zeta F Ka Parameter Settings/Zeta Analysis Settings/Zeta Analysis/Zeta Measurement Settings/Measurement Settings/Zeta Measurement Result": "DLS.FKAMODEL",
            }
            metadata: dict[str, str | float | int] = {}
            attenuators: list[float] = []
            temperatures: list[float] = []
            zaverages: list[float] = []
            polydispersities: list[float] = []
            zeta_averages: list[float] = []
            numbers: list[list[float]] | np.ndarray = []
            volumes: list[list[float]] | np.ndarray = []
            intensities: list[list[float]] | np.ndarray = []
            sizes: list[list[float]] | np.ndarray = []
            zetas: list[list[float]] | np.ndarray | list[np.ndarray] = []
            pots: list[list[float]] | np.ndarray | list[np.ndarray] = []
            int_sizes: list[list[float]] | list[float] | tuple[float, ...] = []
            int_widths: list[list[float]] | list[float] | tuple[float, ...] = []
            int_areas: list[list[float]] | list[float] | tuple[float, ...] = []
            vol_sizes: list[list[float]] | list[float] | tuple[float, ...] = []
            vol_widths: list[list[float]] | list[float] | tuple[float, ...] = []
            vol_areas: list[list[float]] | list[float] | tuple[float, ...] = []
            num_sizes: list[list[float]] | list[float] | tuple[float, ...] = []
            num_widths: list[list[float]] | list[float] | tuple[float, ...] = []
            num_areas: list[list[float]] | list[float] | tuple[float, ...] = []
            zeta_potentials: list[list[float]] | list[float] | tuple[float, ...] = []
            zeta_widths: list[list[float]] | list[float] | tuple[float, ...] = []
            zeta_areas: list[list[float]] | list[float] | tuple[float, ...] = []
            intercepts: list[float] = []
            cumulants_errors: list[float] = []
            multimodal_errors: list[float] = []
            conductivities: list[float] = []
            voltages: list[float] = []

            runs = -1
            li = [
                int_areas,
                int_sizes,
                int_widths,
                vol_areas,
                vol_sizes,
                vol_widths,
                num_areas,
                num_sizes,
                num_widths,
                zeta_areas,
                zeta_potentials,
                zeta_widths,
            ]

            # Read file content
            fc = [
                i.replace("\r", "").replace("\n", "").split(";")
                for i in open(file, encoding="latin1")
            ]

            # Extract data (averaging in case of multiple entries)
            for i, c in enumerate(fc):
                if c[0] == "":
                    continue
                if c[0] in PROPERTY_TYPE_CODE_DICT.keys():
                    # Check if it is a String field (those have the full name in the PROPERTY_TYPE_CODE_DICT), and if so, just use the value directly
                    if (
                        c[0]
                        == "Measurement Start Date And Time/Size Measurement Result"
                    ):
                        metadata[PROPERTY_TYPE_CODE_DICT[c[0]]] = (
                            c[1].replace("T", " ").replace("Z", "").split(".")[0]
                        )
                    else:
                        metadata[PROPERTY_TYPE_CODE_DICT[c[0]]] = c[1]
                else:
                    # If it is a float field, append the results for averaging
                    if "Unclassified " in c[0]:
                        continue

                    if c[0].startswith("Sizes [nm] (Run "):
                        assert isinstance(sizes, list)
                        sizes.append([float(j) for j in c[1:]])
                    elif c[0].startswith("Intensity [%] (Run "):
                        assert isinstance(intensities, list)
                        intensities.append([float(j) for j in c[1:]])
                    elif c[0].startswith("Volume [%] (Run "):
                        assert isinstance(volumes, list)
                        volumes.append([float(j) for j in c[1:]])
                    elif c[0].startswith("Number [%] (Run "):
                        assert isinstance(numbers, list)
                        numbers.append([float(j) for j in c[1:]])
                    elif c[0].startswith("Zeta Potential [mV] (Run "):
                        assert isinstance(pots, list)
                        pots.append([float(j) for j in c[1:]])
                    elif c[0].startswith("Counts [kcps] (Run "):
                        assert isinstance(zetas, list)
                        zetas.append([float(j) for j in c[1:]])
                    elif c[0].startswith("Data for "):
                        runs += 1
                        for j in li:
                            j.append([])
                    elif (
                        c[0]
                        == "Attenuator/Actual Instrument Settings/Actual Instrument Settings/Size Measurement Result"
                    ):
                        attenuators.append(float(c[1]))
                    elif (
                        c[0]
                        == "Temperature (°C)/Actual Instrument Settings/Actual Instrument Settings/Size Measurement Result"
                    ):
                        temperatures.append(float(c[1]))
                    elif (
                        c[0]
                        == "Z-Average (nm)/Cumulants Result/Cumulants Result/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        zaverages.append(float(c[1]))
                    elif (
                        c[0]
                        == "Polydispersity Index (PI)/Cumulants Result/Cumulants Result/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        polydispersities.append(float(c[1]))
                    elif (
                        c[0]
                        == "Intercept/Cumulants Result/Cumulants Result/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        intercepts.append(float(c[1]))
                    elif (
                        c[0]
                        == "Fit Error/Cumulants Result/Cumulants Result/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        cumulants_errors.append(float(c[1]))
                    elif (
                        c[0]
                        == "Fit Error/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        multimodal_errors.append(float(c[1]))
                    elif (
                        c[0]
                        == "Zeta Potential (mV)/Zeta Analysis Result/Zeta Analysis Result/Zeta Measurement Result"
                    ):
                        zeta_averages.append(float(c[1]))
                    elif (
                        c[0]
                        == "Measured Voltage (V)/Zeta Analysis Result/Zeta Analysis Result/Zeta Measurement Result"
                    ):
                        voltages.append(float(c[1]))
                    elif (
                        c[0]
                        == "Conductivity (mS/cm)/Zeta Analysis Result/Zeta Analysis Result/Zeta Measurement Result"
                    ):
                        conductivities.append(float(c[1]))
                    elif c[0].startswith("Data for "):
                        runs += 1
                        for j in li:
                            j.append([])
                    elif (
                        c[0]
                        == "Mean/Size Peak/Particle Size Intensity Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        int_sizes[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Standard Deviation/Size Peak/Particle Size Intensity Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        int_widths[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Area/Size Peak/Particle Size Intensity Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        int_areas[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Mean/Size Peak/Particle Size Volume Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        vol_sizes[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Standard Deviation/Size Peak/Particle Size Volume Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        vol_widths[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Area/Size Peak/Particle Size Volume Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        vol_areas[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Mean/Size Peak/Particle Size Number Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        num_sizes[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Standard Deviation/Size Peak/Particle Size Number Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        num_widths[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Area/Size Peak/Particle Size Number Distribution Peaks ordered by area/Size Analysis Result/Size Analysis Result/Size Measurement Result"
                    ):
                        num_areas[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Mean/Size Peak/Zeta Peaks/Zeta Analysis Result/Zeta Analysis Result/Zeta Measurement Result"
                    ):
                        zeta_potentials[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Standard Deviation/Size Peak/Zeta Peaks/Zeta Analysis Result/Zeta Analysis Result/Zeta Measurement Result"
                    ):
                        zeta_widths[runs].append(float(c[1]))
                    elif (
                        c[0]
                        == "Area/Size Peak/Zeta Peaks/Zeta Analysis Result/Zeta Analysis Result/Zeta Measurement Result"
                    ):
                        zeta_areas[runs].append(float(c[1]))

            # Average peaks across runs
            all_items = [
                [int_areas, int_sizes, int_widths],
                [vol_areas, vol_sizes, vol_widths],
                [num_areas, num_sizes, num_widths],
                [zeta_areas, zeta_potentials, zeta_widths],
            ]
            # Fill rugged arrays with nans for converting to numpy, then average across runs, ignoring nans and using slice notation [:] to modify in-place
            for li in all_items:
                max_len = max([len(jt) for it in li for jt in it])
                k = 0
                for it in li:
                    for jt in it:
                        while len(jt) < max_len:
                            jt.append(np.nan)
                    it[:] = np.nanmean(it[:], axis=0)

            # Sort peaks by Area
            if len(int_areas) > 0:
                int_areas, int_sizes, int_widths = zip(
                    *sorted(list(zip(int_areas, int_sizes, int_widths)), reverse=True)
                )
            if len(vol_areas) > 0:
                vol_areas, vol_sizes, vol_widths = zip(
                    *sorted(list(zip(vol_areas, vol_sizes, vol_widths)), reverse=True)
                )
            if len(num_areas) > 0:
                num_areas, num_sizes, num_widths = zip(
                    *sorted(list(zip(num_areas, num_sizes, num_widths)), reverse=True)
                )
            if len(zeta_areas) > 0:
                zeta_areas, zeta_potentials, zeta_widths = zip(
                    *sorted(
                        list(zip(zeta_areas, zeta_potentials, zeta_widths)),
                        reverse=True,
                    )
                )

            # Calculate PD for individual Peaks
            int_pis = np.array(int_widths) / np.array(int_sizes)
            vol_pis = np.array(vol_widths) / np.array(vol_sizes)
            num_pis = np.array(num_widths) / np.array(num_sizes)

            # Average the numerical values for certain fields and add them to the metadata dict:
            metadata["DLS.ATTENUATOR"] = int(np.array(attenuators).mean())
            metadata["DLS.TEMPERATURE"] = round(np.array(temperatures).mean(), 2)
            metadata["DLS.ZAVG"] = round(np.array(zaverages).mean(), 2)
            metadata["DLS.PDI"] = round(np.array(polydispersities).mean(), 3)
            metadata["DLS.INTERCEPT"] = round(np.array(intercepts).mean(), 3)
            metadata["DLS.CUMULANTSFITERROR"] = round(
                np.array(cumulants_errors).mean(), 6
            )
            metadata["DLS.MULTIMODALFITERROR"] = round(
                np.array(multimodal_errors).mean(), 6
            )
            if len(zeta_areas) > 0:
                metadata["DLS.ZETA"] = round(np.array(zeta_averages).mean(), 2)
                metadata["DLS.VOLT"] = round(np.array(voltages).mean(), 2)
                metadata["DLS.COND"] = round(np.array(conductivities).mean(), 4)

            # For the different distribution types, report data for the first three peaks
            for i in range(0, min(3, len(int_sizes))):
                metadata[f"DLS.PK{i + 1}INT"] = round(int_sizes[i], 2)
                metadata[f"DLS.PK{i + 1}INTWIDTH"] = round(int_widths[i], 2)
                metadata[f"DLS.PK{i + 1}INTPD"] = round(int_pis[i], 2)
            for i in range(0, min(3, len(vol_sizes))):
                metadata[f"DLS.PK{i + 1}VOL"] = round(vol_sizes[i], 2)
                metadata[f"DLS.PK{i + 1}VOLWIDTH"] = round(vol_widths[i], 2)
                metadata[f"DLS.PK{i + 1}VOLPD"] = round(vol_pis[i], 2)
            for i in range(0, min(3, len(num_sizes))):
                metadata[f"DLS.PK{i + 1}NUM"] = round(num_sizes[i], 2)
                metadata[f"DLS.PK{i + 1}NUMWIDTH"] = round(num_widths[i], 2)
                metadata[f"DLS.PK{i + 1}NUMPD"] = round(num_pis[i], 2)
            for i in range(0, min(3, len(zeta_potentials))):
                metadata[f"DLS.PK{i + 1}ZETA"] = round(zeta_potentials[i], 2)
                metadata[f"DLS.PK{i + 1}ZETAWIDTH"] = round(zeta_widths[i], 2)

            dls = metadata_to_instance(metadata, Dls())
            collection.add(dls)
