from flask_wtf import Form
import wtforms

from GeneSequenceApp.models.genbank_entry import GenBankEntry

class __SequenceEntryForm(Form):

	genbank_id = wtforms.SelectField("Genbank Entry")
	sequence = wtforms.TextField("Sequence")
	submit = wtforms.SubmitField("Add")

def SequenceEntryForm(entries):
	genbank_id_choice_texts = ["{} ({})".format(entry.Desc, entry.GenBank_ID) for entry in entries]
	genbank_id_choice_values = [entry.GenBank_ID for entry in entries]

	genbank_id_choices = list(zip(genbank_id_choice_texts, genbank_id_choice_values))

	form = __SequenceEntryForm()
	form.genbank_id.choices = genbank_id_choices

	return form