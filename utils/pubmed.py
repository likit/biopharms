from Bio import Entrez
from datetime import datetime
import dataentry

EMAIL = 'preeyano@msu.edu'
Entrez.email = EMAIL

category = 'vaccine'
query = "vaccine"
year = 2015
search_term = "%s AND Thailand[AFFL] %d[PDAT]" % (query, year)

handle = Entrez.egquery(term=search_term)
record = Entrez.read(handle)
for row in record['eGQueryResult']:
    if row['DbName'] == 'pubmed':
        print(row['DbName'], row['Count'])

handle = Entrez.esearch(db="pubmed", term=search_term)
record = Entrez.read(handle)

n=0
for article_id in record['IdList']:
    pub_data = {}
    authors = []
    handle = Entrez.efetch(db='pubmed', id=article_id,
                rettype="medline", retmode='xml')
    records = Entrez.read(handle)

    print('Getting article Ids...')
    # Get article IDs
    for id in records[0]['PubmedData']['ArticleIdList']:
        if 'IdType' in id.attributes:
            id_type = id.attributes['IdType']
            pub_data[id_type] = str(id)

    print('Getting article Info...')
    # Get some article info
    article_date = records[0]['MedlineCitation']\
                                ['Article'].get('ArticleDate', '')
    if article_date:
        pub_data['ArticleDate'] = datetime(int(article_date[0]['Year']),
                                            int(article_date[0]['Month']),
                                            int(article_date[0]['Day']))
    else:
        pub_data['ArticleDate'] = ''
    pub_data['ArticleTitle'] = unicode(records[0]['MedlineCitation']\
                                ['Article'].get('ArticleTitle', ''))
    abstracts = records[0]['MedlineCitation']\
                                ['Article'].get('Abstract', '')
    # Concat abstract text
    if abstracts:
        abstract_text = []
        for ab in abstracts['AbstractText']:
            abstract_text.append(ab)
        pub_data['Abstract'] = ''.join(abstract_text)
    else:
        pub_data['Abstract'] = ''

    pub_data['JournalMedlineTA'] = unicode(records[0]['MedlineCitation']\
                                ['MedlineJournalInfo'].get('MedlineTA', ''))
    pub_data['JournalCountry'] = unicode(records[0]['MedlineCitation']\
                                ['MedlineJournalInfo'].get('Country', ''))

    # Get authors
    print('Getting authors Info...')
    for author in records[0]['MedlineCitation']['Article']['AuthorList']:
        author_data = {}
        author_data['LastName'] = author.get('LastName', '')
        author_data['ForeName'] = author.get('ForeName', '')
        author_data['Initials'] = author.get('Initials', '')
        author_data['Identifier'] = author.get('Identifier', '')

        author_data['Affiliation'] = []
        for affl in author['AffiliationInfo']:
            author_data['Affiliation'].append(affl.get('Affiliation', ''))

        # print(author['LastName'], author['Initials'], author['ForeName'])
        authors.append(author_data)

    pub_data['BiopharmCategory'] = category

    keywords = []
    if records[0]['MedlineCitation']['KeywordList']:
        for kwlist in records[0]['MedlineCitation']['KeywordList']:
            for kw in kwlist:
                keywords.append(kw)
        pub_data['Keywords'] = ','.join(keywords)
        break

    if n ==5 : break
    n += 1

print(pub_data['Keywords'])
print(pub_data['BiopharmCategory'])
# Add data to the database
# dataentry.add_pubmed(pub_data, authors)
