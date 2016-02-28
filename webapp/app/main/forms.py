from flask.ext.wtf import Form
from wtforms import (SubmitField, BooleanField,
                        StringField, HiddenField, TextAreaField)
from wtforms.validators import Optional
# from wtforms import ValidationError, widgets

class valChainForm(Form):
    randd = BooleanField('Research and Development',
            validators=[Optional()])
    preclin = BooleanField('Pre-clinical research',
            validators=[Optional()])
    clintrials = BooleanField('Clinical trials',
            validators=[Optional()])
    manu = BooleanField('Manufacturing',
            validators=[Optional()])
    sale = BooleanField('Marketing and sales',
            validators=[Optional()])
    firstname = StringField()
    lastname = StringField()
    dname = HiddenField()
    submit = SubmitField('Submit')

class AuthorProfile(Form):
    lastname = StringField('Lastname')
    forename = StringField('Firstname')
    initials = StringField('Initials')
    identifiers = StringField('Indentifiers')
    affil = TextAreaField('Affiliation')
    dname = HiddenField()
    submit = SubmitField('Submit')
