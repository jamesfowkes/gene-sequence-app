import logging

from GeneSequenceApp import app, db

if __name__ == "__main__":
	logging.basicConfig(level=logging.INFO)

	db.create_if_not_exists()

	app.run()