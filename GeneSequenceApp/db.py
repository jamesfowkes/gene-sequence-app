from GeneSequenceApp import app

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy_utils import database_exists, create_database

import logging

Base = declarative_base()

def get_module_logger():
    return logging.getLogger(__name__)
    
get_module_logger().setLevel(logging.INFO)

_session = None
_engine = None

def create_if_not_exists(force_creation=False):
    connect_db()
    
    if not database_exists(_engine.url) or force_creation:
        create_database(_engine.url)
        Base.metadata.create_all(bind=_engine)

    _session.close()

def session():
    return _session

@app.before_request
def before_request():
    connect_db()

@app.teardown_request
def teardown_request(exception):
    if _session is not None:
        _session.close()

def connect_db():
    global _session
    global _engine
    _engine = create_engine("sqlite:///" + app.config['DATABASE'])
    Session = sessionmaker(bind=_engine)
    _session = Session()