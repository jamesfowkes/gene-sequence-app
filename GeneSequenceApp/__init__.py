import os

from flask import Flask
from flask_bootstrap import Bootstrap

from GeneSequenceApp.secret import set_secret_key

app = Flask(__name__)

set_secret_key(app)

Bootstrap(app)

def try_load_config_from_envvar(var):
	if var not in os.environ:
		return

	try:
		app.config.from_envvar(var)
	except:
		pass

	try:
		path = "GeneSequenceApp.config.{}".format(os.environ[var])
		app.config.from_object(path)
	except:
		raise

app.config.from_object("GeneSequenceApp.config.default.Config")
try_load_config_from_envvar("GENESEQUENCEAPP_CONFIG")

import GeneSequenceApp.views.content
import GeneSequenceApp.views.genbank
import GeneSequenceApp.views.genbank_entry_form

import GeneSequenceApp.models.sequence
import GeneSequenceApp.models.genbank_entry
import GeneSequenceApp.models.note
import GeneSequenceApp.models.type
