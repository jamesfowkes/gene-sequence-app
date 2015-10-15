from GeneSequenceApp.db import session, Base
import logging
from sqlalchemy import Column, Integer, String, ForeignKey

def get_module_logger():
	return logging.getLogger(__name__)
	
get_module_logger().setLevel(logging.INFO)

class Type(Base):
	__tablename__ = "types"
	GenBank_ID = Column(String, ForeignKey("sequences.GenBank_ID"), primary_key = True)
	Sequence_ID = Column(Integer, ForeignKey("sequences.Sequence_ID"), primary_key = True) 
	Type_ID = Column(Integer, primary_key = True)
	Type = Column(String)
	Date_Added = Column(Integer)
	