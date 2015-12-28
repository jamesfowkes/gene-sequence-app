import unittest
import tempfile
import logging
import os
import pickle
import datetime

from sqlalchemy import Column, Integer, String

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import GeneSequenceApp
from GeneSequenceApp.db import session, Base
from GeneSequenceApp.datetime_helper import datetime_within_range

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

def split_into_genus_and_species(organism):
    return organism.split(" ")

class GenBankEntry(Base):
    __tablename__ = "genbank_entries"
    GenBank_ID = Column(String, primary_key = True)
    Accession_Number = Column(String)
    Desc = Column(String)
    Genus = Column(String)
    Species = Column(String)
    Silk_Type = Column(String)
    Date_Added = Column(Integer)
    
    @classmethod
    def all(self):
        return session().query(GenBankEntry).all()

    @classmethod
    def get_by_genus_species(cls, genus = None, species = None):
        query = session().query(GenBankEntry)
        if genus is not None:
            query = query.filter(GenBankEntry.Genus == genus)
        if species is not None:
            query = query.filter(GenBankEntry.Species == species)
        return query.all()    

    @classmethod
    def get_by_silk_type(cls, silk_type = None):
        query = session().query(GenBankEntry)
        if silk_type is not None:
            query = query.filter(GenBankEntry.Silk_Type == silk_type)
        return query.all()

    def date_added(self):
        return datetime.datetime.fromtimestamp(self.Date_Added)

    @classmethod
    def create_from_seqrecord(cls, seq_record, silk_type):
        (genus, species) = split_into_genus_and_species(seq_record.annotations['organism'])
        
        new_entry = cls(
            GenBank_ID=seq_record.annotations['gi'],
            Accession_Number=seq_record.id,
            Desc=seq_record.description,
            Genus=genus,
            Species=species,
            Silk_Type=silk_type,
            Date_Added=int(datetime.datetime.now().timestamp()))
        
        new_entry.insert()
        return new_entry

    def insert(self):
        session().add(self)
        session().commit()

    @staticmethod    
    def fill_db_with_test_data():
        
        GenBankEntry(GenBank_ID="62638183", Accession_Number="AY994149.1", Desc="Latrodectus hesperus egg case silk protein-1 (ECP-1) mRNA, complete cds", Genus="Latrodectus", Species="hesperus", Silk_Type="ECP", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="399932052", Accession_Number="JX262192", Desc="Latrodectus hesperus clone 2525 aggregate gland silk factor 2 mRNA, complete cds", Genus="Latrodectus", Species="hesperus", Silk_Type="AcSp", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="422900783", Accession_Number="JX978182", Desc="Latrodectus geometricus clone LgSD7 aciniform spidroin 1 (AcSp1) gene, partial cds", Genus="Latrodectus", Species="geometricus", Silk_Type="AcSp", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="1263288", Accession_Number="U47856.1", Desc="Araneus diadematus fibroin-4 mRNA, partial cds", Genus="Araneus", Species="diadematus", Silk_Type="Fib-4", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="392997885", Accession_Number="JX102566.1", Desc="Megahexura fulva fibroin 1 (fib1) mRNA, partial cds", Genus="Megahexura", Species="fulva", Silk_Type="Fib-1", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="170672094", Accession_Number="EU394445.1", Desc="Latrodectus hesperus minor ampullate spidroin 1-like protein mRNA, partial cds", Genus="Latrodectus", Species="hesperus", Silk_Type="MaSp1", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="89113991", Accession_Number="DQ399324", Desc="Deinopis spinosa clone DS28 MiSp mRNA, partial cds", Genus="Deinopis", Species="spinosa", Silk_Type="MiSp", Date_Added=0).insert()
        GenBankEntry(GenBank_ID="89114011", Accession_Number="DQ399334.1", Desc="Uloborus diversus clone US101 MaSp2 mRNA, partial cds", Genus="Uloborus", Species="diversus", Silk_Type="MaSp2", Date_Added=0).insert()


class GenBankEntryTests(unittest.TestCase):

    def setUp(self):
        self.db_fd, GeneSequenceApp.app.config['DATABASE'] = tempfile.mkstemp()
        self.app = GeneSequenceApp.app.test_client()
        GeneSequenceApp.db.create_if_not_exists(True)
        GenBankEntry.fill_db_with_test_data()
        
    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(GeneSequenceApp.app.config['DATABASE'])

    def test_FilterByGenusSpecies(self):
        with GeneSequenceApp.app.test_request_context("/"):
            GeneSequenceApp.app.preprocess_request()
            results = GenBankEntry.get_by_genus_species(genus = "Latrodectus")
            self.assertEqual(4, len(results))
            results = GenBankEntry.get_by_genus_species(species = "hesperus")
            self.assertEqual(3, len(results))
            results = GenBankEntry.get_by_genus_species(genus = "Latrodectus", species = "hesperus")
            self.assertEqual(3, len(results))

    def test_FilterBySilkType(self):
        with GeneSequenceApp.app.test_request_context("/"):
            GeneSequenceApp.app.preprocess_request()
            results = GenBankEntry.get_by_silk_type(silk_type = "AcSp")
            self.assertEqual(2, len(results))

    def test_CreateFromSeqRecord(self):
    
        
        with GeneSequenceApp.app.test_request_context("/"):
            GeneSequenceApp.app.preprocess_request()
            
            with open(__location__ + "/testfiles/422900759.pickle", 'rb') as f:
                record = pickle.load(f)

            entry = GenBankEntry.create_from_seqrecord(record, "AcSp1")
            self.assertEqual("422900759", entry.GenBank_ID)
            self.assertEqual("JX978171.1", entry.Accession_Number)
            self.assertEqual("Latrodectus hesperus clone 28K13 aciniform spidroin 1 (AcSp1) gene, complete cds.", entry.Desc)
            self.assertEqual("Latrodectus", entry.Genus)
            self.assertEqual("hesperus", entry.Species)
            self.assertEqual("AcSp1", entry.Silk_Type)
            self.assertTrue(
                datetime_within_range(entry.date_added(), datetime.datetime.now(), datetime.timedelta(seconds=1))
            )


if __name__ == "__main__":

    unittest.main() 