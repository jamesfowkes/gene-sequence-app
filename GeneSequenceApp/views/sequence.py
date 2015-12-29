import logging

from flask import render_template
from flask import request, redirect, url_for

from GeneSequenceApp import app
from GeneSequenceApp.views.sequence_entry_form import SequenceEntryForm
from GeneSequenceApp.models.sequence import Sequence
from GeneSequenceApp.models.genbank_entry import GenBankEntry

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

@app.route('/sequences/<genbank_id>')
def render_sequences_view(genbank_id):
	if GenBankEntry.exists(genbank_id):
		sequences = Sequence.get_for_genbank_id(genbank_id)
		return render_template(
			"sequence.view.template.html", title = "Sequences", genbank_id=genbank_id, sequences=sequences)
	else:
		return render_template(
			"sequence.view.template.html", title = "Sequences", genbank_id=genbank_id, sequences=sequences)
	
@app.route('/sequences/add', methods = ("GET" ,"POST"))
def render_sequence_add():
	form_to_validate = SequenceEntryForm(GenBankEntry.all())

	if form_to_validate.validate_on_submit():
		try:
			Sequence.add(form.genbank_id.data, form.sequence.data)
		except:
			return render_template("sequence.view.template.html", title = "Sequences")
	else:
		return render_template("sequence.form.template.html", title = "Add new sequence", form=form_to_validate)
	