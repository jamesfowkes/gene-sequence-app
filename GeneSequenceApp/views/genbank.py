from flask import render_template
from flask import request
from GeneSequenceApp import app

from genbank_entry_form import GenbankEntryForm
@app.route('/genbank')
def render_genbank_view():
	form = GenbankEntryForm()
	return render_template("genbank.template.html", title = "Genbank Interface", form = form)

@app.route('/genbank/add', methods = ("GET" ,"POST"))
def render_genbank_add():
	form_to_validate = GenbankEntryForm()
	genbank_data = {}
	if form_to_validate.validate_on_submit():
		genbank_data["ID"] = form_to_validate.genbank_ID.data
		return render_template("genbank.template.html", title = "Genbank Interface", genbank_data = genbank_data)