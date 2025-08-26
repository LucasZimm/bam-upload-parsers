from .DLSParser import DLSParser
from .IRParser import IRParser
from .NMRParser import NMRParser
from .PlateReaderParser import PlateReaderParser
from .SEMParser import SEMParser
from .TEMParser import TEMParser

# Add more metadata if needed
sem_parser_entry_point = {
    "name": "SEMParser",
    "description": "Parse SEM images saved in tif file format.",
    "parser_class": SEMParser,
}
plate_reader_parser_entry_point = {
    "name": "PlateReaderParser",
    "description": "Parse ASCII files exported from Spectramax measurements.",
    "parser_class": PlateReaderParser,
}
ir_parser_entry_point = {
    "name": "IRParser",
    "description": "Parse ASCII files exported from IR measurements.",
    "parser_class": IRParser,
}
nmr_parser_entry_point = {
    "name": "NMRParser",
    "description": "Parse jcamp dx NMR files.",
    "parser_class": NMRParser,
}
tem_parser_entry_point = {
    "name": "TEMParser",
    "description": "Parse TEM images saved in dm3, emd, or tif file Format",
    "parser_class": TEMParser,
}

dls_parser_entry_point = {
    "name": "DLSParser",
    "description": "Parse ASCII files exported from DLS measurements.",
    "parser_class": DLSParser,
}
