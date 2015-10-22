from flask_wtf import Form
from wtforms import StringField, validators

class GenbankEntryForm(Form):
	genbank_ID = StringField("Genbank ID", validators = [validators.Required()])
	