from flask import Flask
app=Flask(__name__)

app.config.from_object("GeneSequenceApp.default")
#app.config.from_envvar("GENESEQUENCEAPP_CFG")

import GeneSequenceApp.views.content

import GeneSequenceApp.controllers.sequence
import GeneSequenceApp.controllers.genbank_entry
import GeneSequenceApp.controllers.note
import GeneSequenceApp.controllers.type
