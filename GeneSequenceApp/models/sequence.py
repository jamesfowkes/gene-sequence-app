import logging
import unittest
import tempfile
import os
import datetime
import pickle

from sqlalchemy import Column, Integer, String, ForeignKey

import GeneSequenceApp
from GeneSequenceApp.db import session, Base
from GeneSequenceApp.models.genbank_entry import GenBankEntry
from GeneSequenceApp.datetime_helper import datetime_within_range

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

class Sequence(Base):
    __tablename__ = "sequences"
    GenBank_ID = Column(String, ForeignKey("genbank_entries.GenBank_ID"), primary_key = True)
    Sequence_ID = Column(Integer, primary_key = True) 
    Sequence = Column(String)
    Date_Added = Column(Integer)
    Desc = Column(String)

    def date_added(self):
        return datetime.datetime.fromtimestamp(self.Date_Added)

    @staticmethod
    def get_next_sequence_id(genbank_id):
        previous_entries = session().query(Sequence).filter(Sequence.GenBank_ID == genbank_id).all()
        
        if len(previous_entries):
            next_sequence_id = max([entry.Sequence_ID for entry in previous_entries]) + 1
        else:
            next_sequence_id = 0

        return next_sequence_id

    @classmethod
    def create_from_seqrecord(cls, seq_record):
        genbank_id = seq_record.annotations['gi']
        seq = str(seq_record.seq)

        get_module_logger().info(
            "Adding Sequence with genbank ID {} and sequence {}".format(genbank_id, cls.abbreviate(seq)))
        
        new_seq = cls(
            GenBank_ID=genbank_id,
            Sequence_ID=cls.get_next_sequence_id(genbank_id),
            Sequence=seq,
            Date_Added=int(datetime.datetime.now().timestamp()))

        new_seq.insert()
        return new_seq 

    def insert(self):
        session().add(self)
        session().commit()

    def length(self):
        return len(self.Sequence)

    def abbreviated(self):
        return self.abbreviate(self.Sequence)

    @classmethod
    def get_for_genbank_id(cls, genbank_id):
        return session().query(Sequence).filter(Sequence.GenBank_ID == genbank_id).all()

    @staticmethod
    def abbreviate(seq):
        if len(seq) <= 30:
            return seq
        else:
            return seq[0:10] + "..." + seq[-10:]

class SequenceTests(unittest.TestCase):

    @staticmethod
    def fill_db_with_test_data():
        Sequence.add("62638183", "AAAAAAAAAAAAAAAA")
        Sequence.add("62638183", "CCCCCCCCCCCCCCCC")

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

    def test_AddingSequenceForNonExistentEntryProducesError(self):
        with GeneSequenceApp.app.test_request_context("/"):
            GeneSequenceApp.app.preprocess_request()

            with self.assertRaises(Exception):
                seq = Sequence.add("NOTAGENBANKID", "ACGTACGTACGTACGT")
    
    def test_GetForGenBankIDReturnsCorrectSequences(self):
        with GeneSequenceApp.app.test_request_context("/"):
                GeneSequenceApp.app.preprocess_request()

                sequences = Sequence.get_for_genbank_id("62638183")
                self.assertEqual(2, len(sequences))
                self.assertTrue("AAAAAAAAAAAAAAAA" in [seq.Sequence for seq in sequences])
                self.assertTrue("CCCCCCCCCCCCCCCC" in [seq.Sequence for seq in sequences])

    def test_AbbreviatedReturnsFullSequenceForLessThanEqualTo30Chars(self):
        self.assertEqual("", Sequence.abbreviate(""))
        self.assertEqual("A", Sequence.abbreviate("A"))
        self.assertEqual("AB", Sequence.abbreviate("AB"))
        self.assertEqual("ABC", Sequence.abbreviate("ABC"))
        self.assertEqual("ABCDEFGHIJKLMNOPQRSTUVWXYZABC", Sequence.abbreviate("ABCDEFGHIJKLMNOPQRSTUVWXYZABC"))
        self.assertEqual("ABCDEFGHIJKLMNOPQRSTUVWXYZABCD", Sequence.abbreviate("ABCDEFGHIJKLMNOPQRSTUVWXYZABCD"))

    def test_AbbreviateReturnsPartialFullSequenceForGreaterThan30Chars(self):
        self.assertEqual("ABCDEFGHIJ...VWXYZABCDE", Sequence.abbreviate("ABCDEFGHIJKLMNOPQRSTUVWXYZABCDE"))

    def test_CreateFromSeqRecord(self):
    
        with GeneSequenceApp.app.test_request_context("/"):
            GeneSequenceApp.app.preprocess_request()
            
            with open(__location__ + "/testfiles/422900759.pickle", 'rb') as f:
                record = pickle.load(f)

            GenBankEntry.create_from_seqrecord(record, "")
            entry = Sequence.create_from_seqrecord(record)
            self.assertEqual(entry.Sequence, str(record.seq))

if __name__ == "__main__":

    unittest.main() 