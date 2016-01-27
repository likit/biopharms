import pymongo
import copy
import requests

db = pymongo.Connection()['data']

def get_author_data(article):
    for auth in article['creator']:
        url = 'http://api.elsevier.com/content/search/author'

        lastname, firstname = auth.split(',')
        lastname = lastname.lstrip().rstrip()
        firstname = firstname.lstrip().rstrip().replace('.', '')
        initials = '.'.join(firstname) + '.'
        print auth, ':', firstname, lastname, initials
        params = {'apikey': '871232b0f825c9b5f38f8833dc0d8691',
                     'query': 'AUTHLASTNAME(%s)' % lastname}
        resp = requests.get(url, params=params).json()
        if 'search-results' in resp:
            for au in resp['search-results']['entry']:
                #if au['affiliation-current']['affiliation-country'] == 'Thailand':
                if 'given-name' not in au['preferred-name']:
                    continue

                if 'subject-area' in au:
                    for subjarea in au['subject-area']:
                        subjarea['subject'] = subjarea.pop('$')

                if ((lastname == au['preferred-name']['surname'] and
                            initials == au['preferred-name']['initials']) or
                            (firstname == au['preferred-name']['given-name'] and
                                lastname == au['preferred-name']['surname'])):
                    print au['preferred-name']['given-name'], au['preferred-name']['surname']
                    if 'affiliation-current' in au:
                        print au['affiliation-current']['affiliation-id'], au['affiliation-current']['affiliation-country']
                    if 'initials' in au:
                        print au['preferred-name']['initials']
                    author = db.authors.find_one({'dc:identifier': au['dc:identifier']})
                    if author == None:
                        au['pubs'] = [article['doi'][0]]
                        if 'affiliation-current' in au:
                            au['affiliation'] = [au['affiliation-current']]
                        db.authors.insert(au, safe=True)
                        print '\tdata added'
                    else:
                        print '\tauthor exists'
                        if article['doi'][0] not in author['pubs']:
                            print '\tadding new article'
                            newauthor = copy.deepcopy(author)
                            newauthor['pubs'].append(article['doi'][0])
                            if ('affiliation-current' in au) and ('affiliation' in newauthor):
                                newauthor['affiliation'].append(au['affiliation-current'])
                            db.authors.update({'dc:identifier': au['dc:identifier']},
                                    newauthor, safe=True)
                            print '\tdata updated'


for n, article in enumerate(db.scopus.find({'group': 'vaccine'})):
    print n, article['title'][0], article['doi'][0]
    get_author_data(article)
