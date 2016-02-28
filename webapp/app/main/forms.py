from flask.ext.wtf import Form
from wtforms import (SubmitField, BooleanField,
                        StringField, HiddenField, TextAreaField,
                        SelectField)
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

class Pub(Form):
    title = TextAreaField('Article title')
    abstract = TextAreaField('Abstract')
    category = SelectField('Category', choices=[
            ("stem cell", "Stem Cell"),
            ("vaccine", "Vaccine"),
            ("therapeutic peptide", "Therapeutic peptide"),
            ("therapeutic protein", "Therapeutic protein"),
        ])
    pubdate = StringField('Publication date')
    coauthor = TextAreaField('Coauthors')
    author_firstname = StringField('Author firstname')
    author_lastname = StringField('Author lastname')
    pii = StringField('pii')
    doi = StringField('doi')
    keywords = StringField('Keywords')
    dname = HiddenField()
    submit = SubmitField('Submit')
