from flask import render_template
from flask import request
from GeneSequenceApp import app

@app.route('/')
def hello_world():
	return render_template("index.template.html", title = "a nice title for your page", heading1="some other text")
	
@app.route('/complement/')
def complement_seq():
	seq = request.args.get("single_input")
	return render_template("index.template.html", title = "Complement sequence", heading1= "Complement sequence", sequence=seq)
	
@app.route('/complement/input')
def request_seq():
	return render_template("single_input.template.html", title = "Input sequence", prompt = "Enter sequence into available field", submit_url = "complement_seq" )
	
