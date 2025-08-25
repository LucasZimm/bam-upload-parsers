from bam_masterdata.datamodel.object_types import ExperimentalStep
from bam_masterdata.logger import logger
from bam_masterdata.metadata.entities import CollectionType


class TestMasterdataParserExample:
    def test_parse(self, parser):
        collection = CollectionType()
        parser.parse([], collection, logger)

        assert len(collection.attached_objects) == 0
        id1 = collection.add(ExperimentalStep(name="test"))
        parser.parse([], collection, logger)
        assert len(collection.attached_objects) == 1
        id2 = collection.add(ExperimentalStep(name="test2"))
        collection.add_relationship(id1, id2)
        assert len(collection.relationships) == 1
