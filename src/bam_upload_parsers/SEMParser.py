# -----------------------------------------------------------------------------
# File:        SEMParser.py
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
from bam_masterdata.datamodel.object_types import Sem
from bam_masterdata.logger import logger
from bam_masterdata.parsing import AbstractParser
from structlog._config import BoundLoggerLazyProxy

from .utils import metadata_to_instance


class SEMParser(AbstractParser):
    def parse(self, files, collection, logger):
        """
        Parse SEM images saved in tif file format.

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

        metadata: dict[str, str | float | int | set] = {}
        metadata["Name"] = "Name"
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
                # print(f"Tag: {hex(tag_id)}, Value: {value}")
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
            # print("md:", metadata)
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
                "Name": "$NAME",
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
                        if PROPERTY_TYPE_CODE_DICT[k] == "START_DATE":
                            md[PROPERTY_TYPE_CODE_DICT[k]] = datetime.strptime(
                                tmp, "%Y-%m-%d %H:%M:%S"
                            )
                        else:
                            md[PROPERTY_TYPE_CODE_DICT[k]] = tmp
            # print("md:", md)
            return md

        metadatas = []
        for file in files:
            if file.endswith(".tif"):
                metadatas.append(
                    extract_metadata_tif(file, metadata, preview_image_list)
                )
        for metadata in metadatas:
            metadata = get_metadata_for_openbis(metadata)
            if metadata != {}:
                sem = metadata_to_instance(
                    metadata,
                    Sem(),
                )
                collection.add(sem)
