import logging

from flask import render_template, flash, request, redirect, url_for

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
	
	if form_to_validate.validate_on_submit():
		get_module_logger().info("GenBank add form validated.")
		if genbank.try_add_from_genbank(form_to_validate.genbank_id.data, form_to_validate.silk_type.data):
			flash("GenBank ID {} successfully added!".format(form_to_validate.genbank_id.data), 'alert-success')
			return redirect(url_for('render_genbank_view'))
		else:
			flash("Failed to add GenBank ID {}!".format(form_to_validate.genbank_id.data), 'bad')
			return redirect(url_for('render_genbank_view'))
	else:
		return render_template("genbank.form.template.html", title = "GenBank Interface", form = form_to_validate)
