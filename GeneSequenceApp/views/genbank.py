import logging

from flask import render_template
from flask import request, redirect, url_for

from GeneSequenceApp import app

from GeneSequenceApp.views.genbank_entry_form import GenBankEntryForm
from GeneSequenceApp.controllers import genbank
from GeneSequenceApp.models.genbank_entry import GenBankEntry

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

@app.route('/genbank')
def render_genbank_view():
	entries = GenBankEntry.all()

	get_module_logger().info("Rendering table for {} GenBank entries.".format(len(entries)))

	return render_template("genbank.view.template.html", title = "GenBank Interface", genbank_entries = entries)
	
@app.route('/genbank/add', methods = ("GET" ,"POST"))
def render_genbank_add():
	form_to_validate = GenBankEntryForm()
	genbank_data = {}
	if form_to_validate.validate_on_submit():
		if genbank.try_add_from_genbank(form_to_validate.genbank_id.data, form_to_validate.silk_type.data):
			return redirect(url_for('render_genbank_view'))
		else:
			return render_template("genbank.view.template.html", title = "GenBank Interface", genbank_data = "Failed!")
	else:
		return render_template("genbank.form.template.html", title = "GenBank Interface", form = form_to_validate)
