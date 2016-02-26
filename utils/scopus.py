# coding:utf-8
import os
import requests
import xml.etree.ElementTree as et
import json
import sys
import time
import pymongo
import re
from datetime import datetime
from optparse import OptionParser
from dataentry import add_scopus

parser = OptionParser()
parser.add_option('-y', '--year', dest='pub_year', help='Year', type='string')
parser.add_option('-q', '--query', dest='query', help='Query', type='string')
parser.add_option('-g', '--group', dest='group', help='Group', type='string')
parser.add_option('-n', '--number', dest='number', help='Number', type='int')
parser.add_option('-p', '--item-per-page', dest='perpage',
                    help='Items per page', type='int')
parser.add_option('--affil', dest='affil', help='Affiliation', type='string')

(options, args) = parser.parse_args()

# db = pymongo.Connection()['data']

api_key = os.environ['SCOPUS_APIKEY']

'''User defined variables.'''
item_per_page = options.perpage or 200

#TODO: if no pub_year specified, fallback to the current year
pub_year = options.pub_year

outfile_prefix = options.query
if options.affil:
    query = 'tak(%s) AND pub-date IS %s AND AFFIL(%s)' \
                            % (options.query, pub_year, options.affil)
else:
    query = 'tak(%s) AND pub-date IS %s' \
                            % (options.query, pub_year)
sleeptime = 5

url = 'http://api.elsevier.com/content/search/index:SCIDIR'
params = {'apiKey': api_key, 'query': query}
apikey = {'apiKey': api_key}

print >> sys.stderr, 'Downloading papers published in %s' % pub_year
print >> sys.stderr, 'Query is "%s"' % (query)

r = requests.get(url, params=params)
total_results = int(r.json()['search-results']['opensearch:totalResults'])
print >> sys.stderr, 'Total articles found = %d' % total_results
max_articles = options.number or total_results

page = 0
done = False
for start in range(0,total_results+1, item_per_page):
    print >> sys.stderr, \
            'Waiting %d sec to download from page %d... (%d articles/page)' \
                                        % (sleeptime, page+1, item_per_page)
    time.sleep(sleeptime)
    # op = open('%s_page_%d_%d' % (outfile_prefix, page+1, pub_year), 'w')
    url = 'http://api.elsevier.com/content/search/index:SCIDIR'
    params = {'apiKey': api_key,
                'query': query,
                'start': start,
                'count': item_per_page}

    articles = requests.get(url, params=params).json()['search-results']
    for n, entry in enumerate(articles['entry'], start=1):

        print >> sys.stderr, '%d) %s..%s' \
                % (n, entry['dc:title'][:80], entry['dc:creator'][:30])

        article_url = 'http://api.elsevier.com/content/article/%s' % (entry['dc:identifier'])
        article = requests.get(article_url, params=apikey)
        content = et.fromstring(article.text.encode('utf-8'))
        article_dict= {}
        for child in content.getchildren():
            for i in child:
                tag = re.match('(\{([^}]+)\})(.*)',
                        i.tag.encode('utf-8')).groups()[-1]
                if tag not in article_dict:
                    if i.text:
                        article_dict[tag] = [i.text.encode('utf-8')]
                    else:
                        article_dict[tag] = [None]
                else:
                    if i.text:
                        article_dict[tag].append(i.text.encode('utf-8'))
                    else:
                        article_dict[tag].append(None)

        article_dict['group'] = options.group
        pub_data = {}
        pub_data['pii'] = article_dict['pii'][0] or ''
        pub_data['doi'] = article_dict['doi'][0] or ''
        pub_data['ArticleDate'] = article_dict['coverDate'][0] or ''
        pub_data['Keywords'] = ','.join(article_dict['subject'])
        pub_data['ArticleTitle'] = article_dict['title']
        pub_data['Abstract'] = article_dict['description']

        authors = []
        for au in article_dict['creator']:
            lastname, forename = au.split(',')
            forename = forename.lstrip()
            initials = forename[0].upper() + lastname[0].upper()
            authors.append(
                    {
                        'ForeName': forename,
                        'LastName': lastname,
                        'Initials': initials,
                        'Affiliation': ''
                        }
                    )

        # if pub_data:
        #     for key, value in pub_data.iteritems():
        #         print('{} \n\t {}'.format(key, value))
        for au in authors:
            print('{} {} ({})'.format(au['ForeName'],
                                    au['LastName'],
                                    au['Initials']))

        add_scopus(pub_data, authors)
        break
        if n >= max_articles:
            done = True
            break

        # existing_article = db.scopus.find_one({'identifier': article_dict['identifier']})
        # if not existing_article:
        #     db.scopus.insert(article_dict, safe=True)
        # else:
        #     print >> sys.stderr, 'The article already exists.'
    if done:
        print >> sys.stderr, 'Done.'
        break
    print
    page += 1
