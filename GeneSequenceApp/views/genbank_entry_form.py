from flask_wtf import Form
import wtforms

from GeneSequenceApp.models.silk_types import get_silk_types

class GenBankEntryForm(Form):
	genbank_id = wtforms.StringField("GenBank ID", validators = [wtforms.validators.Required()])

	silk_type_choices = list(zip(get_silk_types(), get_silk_types()))

	silk_type = wtforms.SelectField("Silk Type", choices = silk_type_choices)

	submit = wtforms.SubmitField("Add")	