import logging

from GeneSequenceApp.biopython.entrez_fasta import entrez_get_fasta_seqrecord
from GeneSequenceApp.models.genbank_entry import GenBankEntry
from GeneSequenceApp.models.sequence import Sequence

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

def try_add_from_genbank(id, silk_type):

	get_module_logger().info("Calling entrez for ID {}".format(id))
	fasta = entrez_get_fasta_seqrecord(id)
	success = False

	if fasta is not None:
		get_module_logger().info("Data for ID {} found. Trying to add to database.".format(id))
		try:
			get_module_logger().info("Creating GenBankEntry.")
			GenBankEntry.create_from_seqrecord(fasta, silk_type)
			get_module_logger().info("Creating Sequence.")
			Sequence.create_from_seqrecord(fasta)
			success = True
		except:
			raise
	else:
		get_module_logger().info("Data for ID {} not found.")

	return success
