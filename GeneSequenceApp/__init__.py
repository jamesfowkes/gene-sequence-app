from flask import Flask
app=Flask(__name__)

import GeneSequenceApp.views.content

import GeneSequenceApp.controllers.sequence
import GeneSequenceApp.controllers.genbank_entry
import GeneSequenceApp.controllers.note
import GeneSequenceApp.controllers.type
