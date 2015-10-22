from flask import Flask
app=Flask(__name__)

app.config.from_object("GeneSequenceApp.default")
app.config.from_envvar("GENESEQUENCEAPP_CONFIG_PATH")

import GeneSequenceApp.views.content
import GeneSequenceApp.views.genbank

import GeneSequenceApp.models.sequence
import GeneSequenceApp.models.genbank_entry
import GeneSequenceApp.models.note
import GeneSequenceApp.models.type
