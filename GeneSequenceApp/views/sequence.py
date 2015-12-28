import logging

from flask import render_template
from flask import request, redirect, url_for

from GeneSequenceApp import app

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

@app.route('/sequences/<genbank_id>')
def render_sequences_view(genbank_id):
	return render_template("sequence.view.template.html", title = "Sequences")
	
@app.route('/sequences/add', methods = ("GET" ,"POST"))
def render_genbank_add():
	form_to_validate = SequenceEntryForm()
	