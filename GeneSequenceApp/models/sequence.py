import logging
import unittest
import tempfile
import os
import datetime

from sqlalchemy import Column, Integer, String, ForeignKey

import GeneSequenceApp
from GeneSequenceApp.db import session, Base
from GeneSequenceApp.models.genbank_entry import GenBankEntry
from GeneSequenceApp.datetime_helper import datetime_within_range

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

class Sequence(Base):
    __tablename__ = "sequences"
    GenBank_ID = Column(String, ForeignKey("genbank_entries.GenBank_ID"), primary_key = True)
    Sequence_ID = Column(Integer, primary_key = True) 
    Sequence = Column(String)
    Date_Added = Column(Integer)

    def date_added(self):
        return datetime.datetime.fromtimestamp(self.Date_Added)

    @classmethod
    def add(cls, genbank_id, sequence):
        previous_entries = session().query(Sequence).filter(Sequence.GenBank_ID == genbank_id).all()
        
        if len(previous_entries):
            next_sequence_id = max([entry.Sequence_ID for entry in previous_entries]) + 1
        else:
            next_sequence_id = 0

        new_seq = cls(
            GenBank_ID=genbank_id,
            Sequence_ID=next_sequence_id,
            Sequence=sequence,
            Date_Added=int(datetime.datetime.now().timestamp()))

        new_seq.insert()
        return new_seq 

    def insert(self):
        session().add(self)
        session().commit()

class SequenceTests(unittest.TestCase):

    @staticmethod
    def fill_db_with_test_data():
        pass

    def setUp(self):
        self.db_fd, GeneSequenceApp.app.config['DATABASE'] = tempfile.mkstemp()
        self.app = GeneSequenceApp.app.test_client()
        GeneSequenceApp.db.create_if_not_exists(True)
        GenBankEntry.fill_db_with_test_data()
        self.fill_db_with_test_data()

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(GeneSequenceApp.app.config['DATABASE'])

    def test_AddingSequenceForExistingEntryProducesNoError(self):

        with GeneSequenceApp.app.test_request_context("/"):
            GeneSequenceApp.app.preprocess_request()

            seq = Sequence.add("422900783", "ACGTACGTACGTACGT")
            self.assertEqual("422900783", seq.GenBank_ID)
            self.assertEqual(0, seq.Sequence_ID)
            self.assertEqual("ACGTACGTACGTACGT", seq.Sequence)
            self.assertTrue(
                datetime_within_range(seq.date_added(), datetime.datetime.now(), datetime.timedelta(seconds=1))
            )

if __name__ == "__main__":

    unittest.main() 