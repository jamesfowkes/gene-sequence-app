from GeneSequenceApp.db import session, Base
import logging
from sqlalchemy import Column, Integer, String, ForeignKey

def get_module_logger():
	return logging.getLogger(__name__)
	
get_module_logger().setLevel(logging.INFO)

class Note(Base):
	__tablename__ = "notes"
	GenBank_ID = Column(String, ForeignKey("sequences.GenBank_ID"), primary_key = True)
	Sequence_ID = Column(Integer, ForeignKey("sequences.Sequence_ID"), primary_key = True) 
	Note_ID = Column(Integer, primary_key = True)
	Note = Column(String)
	Date_Added = Column(Integer)
	