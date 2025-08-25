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

from . import _dm3_lib_modified as dm3


# help functions for mapping metadata to masterdata object types
def metadata_to_instance(metadata: dict, instance):
    props = metadata_to_masterdata(metadata, instance)
    for k, v in props.items():
        setattr(instance, k, v)
    return instance


def metadata_to_masterdata(metadata: dict, object_instance):
    data = metadata
    object_prop_list = {
        prop_name: data_value
        for prop_name, prop_obj in object_instance._property_metadata.items()
        if (data_value := data.get(prop_obj.code)) is not None
    }
    print(object_prop_list)
    return object_prop_list


class DLSParser(AbstractParser):
    def parse(self, files, collection, logger):
        """
        Parse ASCII files exported from DLS measurements.

        Parameters
        ----------
        files: Union[list, tuple]
            The list of files that are being parsed.
        create_preview: bool = True
            If set to True, a preview image is created in the same folder as the parsed files. Default is True.

        Returns
        -------
        dict[str, Union[str, float, int]]
            The dictionary with the parsed metadata.
        """
        for file in files:
            if file.endswith(".csv"):
                break
        else:
            return
            # Mapping of extracted values to their field names in OpenBIS
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
                if c[0] == "Measurement Start Date And Time/Size Measurement Result":
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
                    list(zip(zeta_areas, zeta_potentials, zeta_widths)), reverse=True
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
        metadata["DLS.CUMULANTSFITERROR"] = round(np.array(cumulants_errors).mean(), 6)
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


class SEMParser(AbstractParser):
    def parse(self, files, collection, logger):
        """
        Parse SEM images saved in tif file format.

        Parameters
        ----------
        files: Union[list, tuple]
            The list of files that are being parsed.
        create_preview: bool = True
            If set to True, a preview image called "SEM_preview.png" is created in the same folder as the parsed files. Default is True.

        Returns
        -------
        dict[str, Union[str, float, int]]
            The dictionary with the parsed metadata.
        """

        metadata: dict[str, str | float | int | set] = {}
        metadata["AcquisitionDateTime"] = ""
        metadata["WorkingDistance"] = set()
        metadata["AccelerationVoltage"] = set()
        metadata["Detector"] = set()
        metadata["PixelSizeX"] = set()
        metadata["PixelSizeY"] = set()
        metadata["ImageSizeX"] = set()
        metadata["ImageSizeY"] = set()

        preview_image_list = {}
        preview_image_list["SEM"] = {}
        preview_image_list["TSEM"] = {}

        def extract_metadata_tif(
            file_path: str,
            metadata: dict[str, str | float | int | set],
            preview_image_list: dict[str, np.ndarray],
        ) -> tuple[dict[str, str | float | int | set], dict[str, np.ndarray]]:
            """
            Extract metadata from tif files.

            Parameters
            ----------
            file_path: str
                The path to the parsed file.
            metadata: dict[str, Union[str, float, int, set]]
                The metadata dict that contains the metadata already parsed from other images.
            preview_image_list: dict[str, np.ndarray]
                The preview image dict that contains the preview images already created from other images.

            Returns
            -------
            tuple[dict[str, Union[str, float, int, set]], dict[str, np.ndarray]]
                The updated metadata and preview image dicts.
            """

            # Read the file
            with open(file_path, "rb") as f:
                fc = f.read()

            # Parse tif header
            if fc[0] == 0x49 and fc[1] == 0x49:
                bo = "little"
            elif fc[0] == 0x4D and fc[1] == 0x4D:
                bo = "big"
            else:
                exit(-1)

            # Read IFD (only works for baseline tiff with one IFD)
            offset = int.from_bytes(fc[0x04:0x08], byteorder=bo)
            ifd_length = int.from_bytes(fc[offset : offset + 2], byteorder=bo)
            ifd = [fc[i : i + 12] for i in range(offset + 2, ifd_length * 12, 12)]
            start_metadata1 = 0
            start_metadata2 = 0
            end_metadata = 0
            for t in ifd:
                tag_id = int.from_bytes(t[0:2], byteorder=bo)
                value = int.from_bytes(t[8:], byteorder=bo)
                print(f"Tag: {hex(tag_id)}, Value: {value}")
                if (
                    int.from_bytes(t[0:2], byteorder=bo) == 0x0111
                ):  # Tag ID for TIFFTAG_STRIPOFFSETS, i.e., start of image data
                    end_metadata = int.from_bytes(t[8:], byteorder=bo)
                elif (
                    int.from_bytes(t[0:2], byteorder=bo) == 0x8546
                ):  # Tag ID for first metadata block
                    start_metadata1 = int.from_bytes(t[8:], byteorder=bo)
                elif (
                    int.from_bytes(t[0:2], byteorder=bo) == 0x8547
                ):  # Tag ID for second metadata block
                    start_metadata2 = int.from_bytes(t[8:], byteorder=bo)
                elif int.from_bytes(t[0:2], byteorder=bo) == 0xFFF1:
                    start_metadata1 = int.from_bytes(t[8:], byteorder=bo)

                elif int.from_bytes(t[0:2], byteorder=bo) == 0x100:
                    metadata["ImageSizeX"].add(str(value))
                elif int.from_bytes(t[0:2], byteorder=bo) == 0x101:
                    metadata["ImageSizeY"].add(str(value))

            if start_metadata1 == 0 or end_metadata == 0:
                exit(-1)

            # Parse metadata from header and store them in a dictionary
            header = fc[start_metadata1:end_metadata].decode("latin-1")
            header = [i.replace("\r", "") for i in header.split("\n")]
            mode = ""

            for i, t in enumerate(header):
                if t == "AP_WD":
                    metadata["WorkingDistance"].add(
                        header[i + 1].replace("WD = ", "").strip()
                    )
                elif t == "AP_PIXEL_SIZE":
                    tmp_pixel_size = header[i + 1].replace("Pixel Size = ", "").strip()
                    metadata["PixelSizeX"].add(tmp_pixel_size)
                    metadata["PixelSizeY"].add(tmp_pixel_size)
                elif t == "DP_DETECTOR_CHANNEL":
                    mode = header[i + 1].replace("Signal A = ", "").strip()
                    metadata["Detector"].add(mode)
                elif t == "AP_MANUALKV":
                    metadata["AccelerationVoltage"].add(
                        header[i + 1].replace("EHT Target = ", "").strip()
                    )
                elif t == "AP_DATE":
                    date = header[i + 1].replace("Date :", "")
                    time = "00:00:00"
                    for j, s in enumerate(header):
                        if s == "AP_TIME":
                            time = header[j + 1].replace("Time :", "")
                            break
                    metadata["AcquisitionDateTime"] = datetime.strptime(
                        date + " " + time, "%d %b %Y %H:%M:%S"
                    ).strftime("%Y-%m-%d %H:%M:%S")

            return metadata

        def get_metadata_for_openbis(
            metadata: dict[str, str | float | int | set],
        ) -> dict[str, str | float | int]:
            """
            Write the extracted metadata to a dictionary that can be used with OpenBIS.

            Parameters
            ----------
            metadata: dict[str, Union[str, float, int, set]]
                The metadata dict that contains the metadata already parsed from other images.

            Returns
            -------
            dict[str, Union[str, float, int]]
                The updated metadata dict.
            """

            PROPERTY_TYPE_CODE_DICT = {
                "WorkingDistance": "SEM.WORKINGDISTANCE",
                "Detector": "SEM.DETECTOR",
                "PixelSizeX": "SEM.PIXELSIZEX",
                "PixelSizeY": "SEM.PIXELSIZEY",
                "ImageSizeX": "SEM.IMAGESIZEX",
                "ImageSizeY": "SEM.IMAGESIZEY",
                "AccelerationVoltage": "SEM.ACCELERATIONVOLTAGE",
                "AcquisitionDateTime": "START_DATE",
            }

            md: dict[str, str | float | int] = {}
            for k in PROPERTY_TYPE_CODE_DICT.keys():
                if k in metadata.keys():
                    if isinstance(metadata[k], set):
                        tmp = list(metadata[k])
                        tmp = ";".join(tmp)
                    else:
                        tmp = metadata[k]
                    if len(tmp) > 0:
                        md[PROPERTY_TYPE_CODE_DICT[k]] = tmp

            return md

        object_ids = []
        metadatas = []
        for file in files:
            if file.endswith(".tif"):
                metadatas.append(
                    extract_metadata_tif(file, metadata, preview_image_list)
                )
        for metadata in metadatas:
            metadata = get_metadata_for_openbis(metadata)
            if metadata != {}:
                sem = metadata_to_instance(metadata, Sem())
                collection.add(sem)

        return object_ids


class TEMParser(AbstractParser):
    def parse(
        self, files: list[str], collection: CollectionType, logger: BoundLoggerLazyProxy
    ) -> None:
        """
        Parse TEM images saved in dm3, emd, or tif file format

        Parameters
        ----------
        files: Union[list, tuple]
            The list of files that are being parsed.
        create_preview: bool = True
            If set to True, a preview image called "TEM_preview.png" is created in the same folder as the parsed files. Default is True.

        Returns
        -------
        dict[str, Union[str, float, int]]
            The dictionary with the parsed metadata.
        """

        metadata: dict[str, str | float | int | set] = {}
        metadata["Instrument"] = ""
        metadata["AcquisitionDateTime"] = ""
        metadata["OperatingMode"] = set()
        metadata["ProjectorMode"] = set()
        metadata["AccelerationVoltage"] = set()
        metadata["Detector"] = set()
        metadata["PixelSizeX"] = set()
        metadata["PixelSizeY"] = set()
        metadata["ImageSizeX"] = set()
        metadata["ImageSizeY"] = set()
        metadata["Magnification"] = set()
        metadata["CameraLength"] = set()
        metadata["SA_ApertureDiameter"] = set()
        metadata["SA_AperturePosX"] = set()
        metadata["SA_AperturePosY"] = set()
        metadata["UserID"] = set()
        metadata["SampleID"] = set()
        metadata["SampleDescription"] = set()
        metadata["SpotIndex"] = set()
        metadata["GunLensSetting"] = set()
        metadata["Aperture[C2].Name"] = set()
        metadata["Aperture[OBJ].Name"] = set()

        preview_image_list = {}
        preview_image_list["TEM;Imaging"] = {}
        preview_image_list["TEM;Diffraction"] = {}
        preview_image_list["STEM;Diffraction"] = {}

        def extract_metadata_tif(
            file_path: str,
            metadata: dict[str, str | float | int | set],
            preview_image_list: dict[str, np.ndarray],
        ) -> tuple[dict[str, str | float | int | set], dict[str, np.ndarray]]:
            """
            Extract metadata from tif files.

            Parameters
            ----------
            file_path: str
                The path to the parsed file.
            metadata: dict[str, Union[str, float, int, set]]
                The metadata dict that contains the metadata already parsed from other images.
            preview_image_list: dict[str, np.ndarray]
                The preview image dict that contains the preview images already created from other images.

            Returns
            -------
            tuple[dict[str, Union[str, float, int, set]], dict[str, np.ndarray]]
                The updated metadata and preview image dicts.
            """

            # xml namespaces
            ns = {
                "xmlns:nil": "http://schemas.fei.com/Metadata/v1/2013/07",
                "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            }

            # Read the file and search for xml metadata
            with open(file_path, "rb") as f:
                fc = f.read()

            fc = fc[
                fc.rfind(b"<?xml version=") : fc.rfind(b"</Metadata>")
                + len(b"</Metadata>")
            ].decode("utf-8")

            if fc == "":  # no xml data found
                return metadata

            # Parse the metadata and store them in a dictionary
            tree = ET.ElementTree(ET.fromstring(fc))
            root = tree.getroot()

            # metadata['Instrument'] = root.find('Instrument', ns).find('Manufacturer', ns).text
            # metadata['Instrument'] += ' ' + root.find('Instrument', ns).find('InstrumentModel', ns).text
            metadata["Instrument"] = (
                root.find("Instrument", ns).find("InstrumentModel", ns).text
            )
            metadata["AcquisitionDateTime"] = (
                root.find("Acquisition", ns)
                .find("AcquisitionDatetime", ns)
                .text.replace("T", " ")
            )
            metadata["AcquisitionDateTime"] = metadata["AcquisitionDateTime"][
                : metadata["AcquisitionDateTime"].rfind(".")
            ]
            tmp_operating_mode = root.find("Optics", ns).find("OperatingMode", ns).text
            metadata["OperatingMode"].add(tmp_operating_mode)
            tmp_projector_mode = root.find("Optics", ns).find("ProjectorMode", ns).text
            metadata["ProjectorMode"].add(tmp_projector_mode)
            metadata["AccelerationVoltage"].add(
                "{} keV".format(
                    int(root.find("Optics", ns).find("AccelerationVoltage", ns).text)
                    // 1000
                )
            )
            metadata["Detector"].add(
                root.find("BinaryResult", ns).find("Detector", ns).text
            )
            metadata["SpotIndex"].add(
                root.find("Optics", ns).find("SpotIndex", ns).text
            )
            metadata["GunLensSetting"].add(
                root.find("Optics", ns).find("GunLensSetting", ns).text
            )
            if (
                root.find("Core", ns) is not None
                and root.find("Core", ns).find("UserID", ns) is not None
            ):  # Exported from new software version (contains more fields)
                metadata["UserID"].add(root.find("Core", ns).find("UserID", ns).text)
            if (
                root.find("Sample", ns) is not None
                and root.find("Sample", ns).find("SampleID", ns) is not None
            ):  # Exported from new software version (contains more fields)
                metadata["SampleID"].add(
                    root.find("Sample", ns).find("SampleID", ns).text
                )
            if (
                root.find("Sample", ns) is not None
                and root.find("Sample", ns).find("SampleDescription", ns) is not None
            ):  # Exported from new software version (contains more fields)
                metadata["SampleDescription"].add(
                    root.find("Sample", ns).find("SampleDescription", ns).text
                )

            tmp_pixel_size_x = (
                root.find("BinaryResult", ns).find("PixelSize", ns).find("X", ns).text
            )
            tmp_pixel_size_x += " " + root.find("BinaryResult", ns).find(
                "PixelSize", ns
            ).find("X", ns).get("unit")
            metadata["PixelSizeX"].add(tmp_pixel_size_x)
            tmp_pixel_size_y = (
                root.find("BinaryResult", ns).find("PixelSize", ns).find("Y", ns).text
            )
            tmp_pixel_size_y += " " + root.find("BinaryResult", ns).find(
                "PixelSize", ns
            ).find("Y", ns).get("unit")
            metadata["PixelSizeY"].add(tmp_pixel_size_y)

            for p in root.find("CustomPropertyGroup", ns).find("CustomProperties", ns):
                if p.get("name") == "Aperture[C2].Name":
                    metadata["Aperture[C2].Name"].add(p.get("value"))
                if p.get("name") == "Aperture[OBJ].Name":
                    metadata["Aperture[OBJ].Name"].add(p.get("value"))

            if tmp_operating_mode == "TEM":
                metadata["ImageSizeX"].add(
                    root.find("BinaryResult", ns)
                    .find("ImageSize", ns)
                    .find("X", ns)
                    .text
                )
                metadata["ImageSizeY"].add(
                    root.find("BinaryResult", ns)
                    .find("ImageSize", ns)
                    .find("Y", ns)
                    .text
                )
                if tmp_projector_mode == "Imaging":  # Brightfield
                    metadata["Magnification"].add(
                        root.find("Optics", ns).find("NominalMagnification", ns).text
                    )
                else:  # ED
                    pass  # Magnification metadata field not present in ED
            else:  # STEM
                metadata["CameraLength"].add(
                    root.find("Optics", ns).find("CameraLength", ns).text
                )
                metadata["ImageSizeX"].add(
                    root.find("ScanSettings", ns)
                    .find("ScanSize", ns)
                    .find("Width", ns)
                    .text
                )
                metadata["ImageSizeY"].add(
                    root.find("ScanSettings", ns)
                    .find("ScanSize", ns)
                    .find("Height", ns)
                    .text
                )
                for p in root.find("CustomPropertyGroup", ns).find(
                    "CustomProperties", ns
                ):
                    if p.get("name") == "StemMagnification":
                        metadata["Magnification"].add(p.get("value"))
                        break

            # These xml fields are always present, but they are not really importrant for the other modes
            if tmp_projector_mode == "Diffraction":
                for i in root.iter("Aperture"):
                    if i.find("Name", ns).text == "SA":  # Selected Area Aperture
                        if (
                            i.find("Type", ns).text == "None"
                        ):  # Check if Aperture is actually used
                            break
                        metadata["CameraLength"].add(
                            root.find("Optics", ns).find("CameraLength", ns).text
                        )
                        metadata["SA_ApertureDiameter"].add(i.find("Diameter", ns).text)
                        metadata["SA_AperturePosX"].add(
                            i.find("PositionOffset", ns).find("X", ns).text
                        )
                        metadata["SA_AperturePosY"].add(
                            i.find("PositionOffset", ns).find("Y", ns).text
                        )
                    else:
                        continue

            return metadata, preview_image_list

        def extract_metadata_emd(
            file_path: str,
            metadata: dict[str, str | float | int | set],
            preview_image_list: dict[str, np.ndarray],
        ):
            """
            Extract metadata from emd files.

            Parameters
            ----------
            file_path: str
                The path to the parsed file.
            metadata: dict[str, Union[str, float, int, set]]
                The metadata dict that contains the metadata already parsed from other images.
            preview_image_list: dict[str, np.ndarray]
                The preview image dict that contains the preview images already created from other images.

            Returns
            -------
            dict[str, Union[str, float, int, set]], dict[str, np.ndarray]
                The updated metadata and preview image dicts.
            """

            with h5py.File(file_path, "r") as f:
                # Check if there is an image in the emd file (not the case for edx spectra)
                if f.get("Data/Image") is None:
                    return metadata, preview_image_list
                # Take the first image (currently we don't have multiple images per file)
                image_id = [im for im in f.get("Data/Image")][0]

                # Parse imagedata into numpy array
                img = np.squeeze(f.get(f"Data/Image/{image_id}/Data"))

                # Parse the metadata from a bytearray into a string
                all_metadata = (
                    bytearray(np.squeeze(f.get(f"Data/Image/{image_id}/Metadata")))
                    .decode("utf-8")
                    .replace("\x00", "")
                )

            if all_metadata == "":  # no json data found
                return metadata

            # Convert to json
            all_metadata = json.loads(all_metadata)

            # metadata['Instrument'] = all_metadata['Instrument']['Manufacturer']
            # metadata['Instrument'] += ' ' + all_metadata['Instrument']['InstrumentModel']
            metadata["Instrument"] = all_metadata["Instrument"]["InstrumentModel"]
            metadata["AcquisitionDateTime"] = all_metadata["Acquisition"][
                "AcquisitionStartDatetime"
            ]["DateTime"]
            metadata["AcquisitionDateTime"] = datetime.utcfromtimestamp(
                int(metadata["AcquisitionDateTime"])
            ).strftime("%Y-%m-%d %H:%M:%S")
            tmp_operating_mode = all_metadata["Optics"]["OperatingMode"]
            tmp_projector_mode = all_metadata["Optics"]["ProjectorMode"]
            metadata["AccelerationVoltage"].add(
                "{} keV".format(
                    int(all_metadata["Optics"]["AccelerationVoltage"]) // 1000
                )
            )
            metadata["Detector"].add(all_metadata["BinaryResult"]["Detector"])
            tmp_pixel_size_x = all_metadata["BinaryResult"]["PixelSize"]["width"]
            tmp_pixel_size_x += " " + all_metadata["BinaryResult"]["PixelUnitX"]
            metadata["PixelSizeX"].add(tmp_pixel_size_x)
            tmp_pixel_size_y = all_metadata["BinaryResult"]["PixelSize"]["height"]
            tmp_pixel_size_y += " " + all_metadata["BinaryResult"]["PixelUnitY"]
            metadata["PixelSizeY"].add(tmp_pixel_size_y)
            metadata["SpotIndex"].add(all_metadata["Optics"]["SpotIndex"])
            metadata["GunLensSetting"].add(all_metadata["Optics"]["GunLensSetting"])
            metadata["Aperture[C2].Name"].add(
                all_metadata["CustomProperties"]["Aperture[C2].Name"]["value"]
            )
            metadata["Aperture[OBJ].Name"].add(
                all_metadata["CustomProperties"]["Aperture[OBJ].Name"]["value"]
            )

            if (
                "Core" in all_metadata and "UserId" in all_metadata["Core"]
            ):  # Exported from new software version (contains more fields)
                metadata["UserID"].add(all_metadata["Core"]["UserId"])
            if (
                "Sample" in all_metadata and "SampleId" in all_metadata["Sample"]
            ):  # Exported from new software version (contains more fields)
                metadata["SampleID"].add(all_metadata["Sample"]["SampleId"])
            if (
                "Sample" in all_metadata
                and "SampleDescription" in all_metadata["Sample"]
            ):  # Exported from new software version (contains more fields)
                metadata["SampleDescription"].add(
                    all_metadata["Sample"]["SampleDescription"]
                )

            if tmp_operating_mode == "1":
                tmp_operating_mode = "TEM"
                metadata["OperatingMode"].add("TEM")
                metadata["ImageSizeX"].add(
                    all_metadata["BinaryResult"]["ImageSize"]["width"]
                )
                metadata["ImageSizeY"].add(
                    all_metadata["BinaryResult"]["ImageSize"]["height"]
                )
            elif tmp_operating_mode == "2":
                tmp_operating_mode = "STEM"
                metadata["OperatingMode"].add("STEM")
                metadata["ImageSizeX"].add(all_metadata["Scan"]["ScanSize"]["width"])
                metadata["ImageSizeY"].add(all_metadata["Scan"]["ScanSize"]["height"])

            if tmp_projector_mode == "1":
                tmp_projector_mode = "Diffraction"
                metadata["ProjectorMode"].add("Diffraction")
                metadata["CameraLength"].add(all_metadata["Optics"]["CameraLength"])
                if tmp_operating_mode == "STEM":
                    metadata["Magnification"].add(
                        all_metadata["CustomProperties"]["StemMagnification"]["value"]
                    )
            elif tmp_projector_mode == "2":
                tmp_projector_mode = "Imaging"
                metadata["ProjectorMode"].add("Imaging")
                metadata["Magnification"].add(
                    all_metadata["Optics"]["NominalMagnification"]
                )

            # These xml fields are always present, but they are not really important for the other modes
            if tmp_projector_mode == "Diffraction":
                for i in all_metadata["Optics"]["Apertures"].keys():
                    j = all_metadata["Optics"]["Apertures"][i]
                    if j["Name"] == "SA" and j["Type"] != "None":
                        metadata["SA_ApertureDiameter"].add(j["Diameter"])
                        metadata["SA_AperturePosX"].add(j["PositionOffset"]["x"])
                        metadata["SA_AperturePosY"].add(j["PositionOffset"]["y"])

            # Extract preview image
            preview_img = np.asarray(img, dtype="float32")
            preview_image_list[f"{tmp_operating_mode};{tmp_projector_mode}"][
                tmp_pixel_size_x
            ] = preview_img  # Overwrite previous preview images

            return metadata, preview_image_list

        def extract_metadata_dm3(
            file_path: str,
            metadata: dict[str, str | float | int | set],
            preview_image_list: dict[str, np.ndarray],
        ):
            """
            Extract metadata from dm3 files.

            Parameters
            ----------
            file_path: str
                The path to the parsed file.
            metadata: dict[str, Union[str, float, int, set]]
                The metadata dict that contains the metadata already parsed from other images.
            preview_image_list: dict[str, np.ndarray]
                The preview image dict that contains the preview images already created from other images.

            Returns
            -------
            dict[str, Union[str, float, int, set]], dict[str, np.ndarray]
                The updated metadata and preview image dicts.
            """

            # TODO: Test with STEM and SAED Images
            f = dm3.DM3(file_path)
            allMetadata = f.info

            if allMetadata == "":  # no metadata found
                return metadata

            acq_date = allMetadata["acq_date"].decode("utf-8").split(".")

            metadata["Instrument"] = allMetadata["micro"].decode("utf-8")
            metadata["AcquisitionDateTime"] = "{}-{}-{} {}".format(
                acq_date[2],
                acq_date[1],
                acq_date[0],
                allMetadata["acq_time"].decode("utf-8"),
            )
            tmp_operating_mode = allMetadata["OperatingMode"].decode("utf-8")
            metadata["OperatingMode"].add(tmp_operating_mode)
            tmp_projector_mode = allMetadata["ProjectorMode"].decode("utf-8")
            tmp_projector_mode = tmp_projector_mode[0] + tmp_projector_mode[1:].lower()
            metadata["ProjectorMode"].add(tmp_projector_mode)
            metadata["AccelerationVoltage"].add(
                "{} keV".format(int(float(allMetadata["hv"].decode("utf-8")) // 1000))
            )
            metadata["Detector"].add(allMetadata["Detector"].decode("utf-8"))
            tmp_pixel_size_x = allMetadata["ScaleX"].decode("utf-8")
            tmp_pixel_size_x += " " + allMetadata["UnitX"].decode("utf8")
            metadata["PixelSizeX"].add(tmp_pixel_size_x)
            tmp_pixel_size_y = allMetadata["ScaleY"].decode("utf-8")
            tmp_pixel_size_y += " " + allMetadata["UnitY"].decode("utf8")
            metadata["PixelSizeY"].add(tmp_pixel_size_y)
            metadata["Magnification"].add(allMetadata["mag"].decode("utf8"))

            if tmp_operating_mode == "Diffraction" and tmp_projector_mode != "STEM":
                metadata["CameraLength"].add(
                    allMetadata["CameraLength"].decode("utf-8")
                )
                # SA-Aperture setting does not seem to be logged on this instrument or in this file format

            # Extract preview image
            preview_img = np.asarray(f.Image, dtype="float32").copy()
            preview_image_list[f"{tmp_operating_mode};{tmp_projector_mode}"][
                tmp_pixel_size_x
            ] = preview_img  # Overwrite previous preview images

            return metadata, preview_image_list

        def get_metadata_for_openbis(
            metadata: dict[str, str | float | int | set],
        ) -> dict[str, str | float | int]:
            """
            Write the extracted metadata to a dictionary that can be used with OpenBIS.

            Parameters
            ----------
            metadata: dict[str, Union[str, float, int, set]]
                The metadata dict that contains the metadata already parsed from other images.

            Returns
            -------
            dict[str, Union[str, float, int]]
                The updated metadata dict.
            """

            PROPERTY_TYPE_CODE_DICT = {
                "AcquisitionDateTime": "START_DATE",
                "Instrument": "TEM.INSTRUMENT",
                "OperatingMode": "TEM.OPERATINGMODE",
                "ProjectorMode": "TEM.PROJECTORMODE",
                "AccelerationVoltage": "TEM.ACCELERATIONVOLTAGE",
                "Detector": "TEM.DETECTOR",
                "PixelSizeX": "TEM.PIXELSIZEX",
                "PixelSizeY": "TEM.PIXELSIZEY",
                "ImageSizeX": "TEM.IMAGESIZEX",
                "ImageSizeY": "TEM.IMAGESIZEY",
                "Magnification": "TEM.MAGNIFICATION",
                "CameraLength": "TEM.CAMERALENGTH",
                "SA_ApertureDiameter": "TEM.SAED_APERTUREDIAMETER",
                "SA_AperturePosX": "TEM.SAED_APERTUREPOSX",
                "SA_AperturePosY": "TEM.SAED_APERTUREPOSY",
                "UserID": "TEM.OPERATOR",
                "SampleID": "TEM.SHORT_NAME",
                "SampleDescription": "TEM.SAMPLE_DESCRIPTION",
                "SpotIndex": "TEM.SPOT_INDEX",
                "GunLensSetting": "TEM.GUN_LENS_SETTING",
                "Aperture[C2].Name": "TEM.C2_APERTURE_NAME",
                "Aperture[OBJ].Name": "TEM.OBJ_APERTURE_NAME",
            }

            md = {}
            for k in PROPERTY_TYPE_CODE_DICT.keys():
                metadata
                if k in metadata.keys():
                    if isinstance(metadata[k], set):
                        tmp = list(metadata[k])
                        tmp.sort()
                        tmp = ";".join(tmp)
                    else:
                        tmp = metadata[k]
                    if len(tmp) > 0:
                        md[PROPERTY_TYPE_CODE_DICT[k]] = tmp

            return md

        for file in files:
            if file.lower().endswith(".tif") or file.lower().endswith(".tiff"):
                metadata, preview_image_list = extract_metadata_tif(
                    file, metadata, preview_image_list
                )
            elif file.lower().endswith(".emd"):
                metadata, preview_image_list = extract_metadata_emd(
                    file, metadata, preview_image_list
                )
            elif file.lower().endswith(".dm3"):
                metadata, preview_image_list = extract_metadata_dm3(
                    file, metadata, preview_image_list
                )

        if len(metadata) == 0:
            return {}

        object_ids = []
        metadatas = []
        for file in files:
            if file.endswith(".tif"):
                metadatas.append(
                    extract_metadata_tif(file, metadata, preview_image_list)
                )
        for metadata[0] in metadatas:
            metadata = get_metadata_for_openbis(metadata)
            if metadata != {}:
                sem = metadata_to_instance(metadata, Tem())
                collection.add(sem)

        return object_ids


class NMRParser(AbstractParser):
    def parse(
        self, files: list[str], collection: CollectionType, logger: BoundLoggerLazyProxy
    ) -> None:
        """
        Parse jcamp dx NMR files.

        Parameters
        ----------
        files: Union[list, tuple]
            The list of files that are being parsed.
        create_preview: bool = True
            If set to True, a preview image is created in the same folder as the parsed files. Default is True.

        Returns
        -------
        dict[str, Union[str, float, int]]
            The dictionary with the parsed metadata.
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

            nmr = metadata_to_instance(metadata, Nmr())
            collection.add(nmr)


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
        files: Union[list, tuple]
            The list of files that are being parsed.
        create_preview: bool = True
            If set to True, a preview image is created in the same folder as the parsed files. Default is True.

        Returns
        -------
        dict[str, Union[str, float, int]]
            The dictionary with the parsed metadata.
        """

        for i in files:
            filename_spa_list = []
            filename_csv_list = []
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

            ftir = metadata_to_instance(metadata, Ftir())
            collection.add(ftir)


class PlateReaderParser(AbstractParser):
    def parse(
        self, files: list[str], collection: CollectionType, logger: BoundLoggerLazyProxy
    ) -> None:
        """
        Parse ASCII files exported from Spectramax measurements.

        Parameters
        ----------
        files: Union[list, tuple]
            The list of files that are being parsed.
        create_preview: bool = True
            If set to True, a preview image is created in the same folder as the parsed files. Default is True.

        Returns
        -------
        dict[str, Union[str, float, int]]
            The dictionary with the parsed metadata.
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
