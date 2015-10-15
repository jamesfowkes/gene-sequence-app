from GeneSequenceApp.db import session, Base
import logging
from sqlalchemy import Column, Integer, String, ForeignKey

def get_module_logger():
	return logging.getLogger(__name__)
	
get_module_logger().setLevel(logging.INFO)

class Sequence(Base):
	__tablename__ = "sequences"
	GenBank_ID = Column(String, ForeignKey("genbank_entries.GenBank_ID"), primary_key = True)
	Sequence_ID = Column(Integer, primary_key = True) 
	Sequence = Column(String)
	Date_Added = Column(Integer)
