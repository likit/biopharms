from datetime import datetime
from . import main
from .. import graph
from collections import defaultdict

from flask import (Flask, render_template, session,
                    request, redirect, url_for, flash, jsonify)


ALLLABELS = ['PUBMED', 'SCOPUS']

@main.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')


@main.route('/main_page', methods=['GET', 'POST'])
def main_page():
    labels = graph.cypher.execute("MATCH (n) return labels(n);")
    select_labels = set()
    for lab in labels:
        for l in lab[0]: select_labels.add(l)
    return render_template('main.html',
            alllabels=ALLLABELS, select_labels=list(select_labels))


@main.route('/_get_people_list')
def get_people_list():
    cypher_command = \
            "MATCH (n:AUTHOR) return n;"
    people = graph.cypher.execute(cypher_command)
    alist = []
    for p in people:
        firstname = p.n.properties.get('ForeName', '')
        lastname = p.n.properties.get('LastName', '')
        initials = p.n.properties.get('Initials', '')
        try:
            affl = p.n.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        alist.append([firstname, lastname, initials, affl, 0])

    return jsonify(data=alist)


@main.route('/_get_coauthor_list', methods=['GET', 'POST'])
def get_coauthor_list():
    firstname, lastname = request.args.get('fullname').split('|')
    cypher_command = \
            "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(c)<-[g:COAUTHOR]-(f:AUTHOR) return f;" % (firstname, lastname)
    people = graph.cypher.execute(cypher_command)
    alist = []
    for p in people:
        firstname = p.f.properties.get('ForeName', '')
        lastname = p.f.properties.get('LastName', '')
        initials = p.f.properties.get('Initials', '')
        try:
            affl = p.f.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        alist.append([firstname, lastname, initials, affl, 0])

    return jsonify(data=alist)


@main.route('/view_person')
def view_person():
    fullname = request.args.get('fullname')
    firstname, lastname = fullname.split('|')
    cypher_command = \
            "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'}) return n;" \
            % (firstname, lastname)
    author = graph.cypher.execute(cypher_command)
    cypher_command = \
            "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(f:ARTICLE) return f;" \
            % (firstname, lastname)
    pubs = graph.cypher.execute(cypher_command)
    kw=set()
    biopharmcat = set()
    pubyears = {}
    for y in range(2005,2016):
        pubyears[y] = 0

    for p in pubs:
        if p.f.properties['Keywords']:
            for word in (p.f.properties['Keywords'].split(',')):
                kw.add(word.lower())
        if p.f.properties['BiopharmCategory']:
            biopharmcat.add(p.f.properties['BiopharmCategory'])
        if p.f.properties['ArticleDate']:
            pubdate = datetime.strptime(p.f.properties['ArticleDate'],
                    '%Y-%m-%dT%H:%M:%S')
            pubyears[pubdate.year] += 1

    pubyear_data = [{
            'values': [{'key': '2005', 'y': pubyears[2005]},
                        {'key': '2006', 'y': pubyears[2006]},
                        {'key': '2007', 'y': pubyears[2007]},
                        {'key': '2008', 'y': pubyears[2008]},
                        {'key': '2009', 'y': pubyears[2009]},
                        {'key': '2010', 'y': pubyears[2010]},
                        {'key': '2011', 'y': pubyears[2011]},
                        {'key': '2012', 'y': pubyears[2012]},
                        {'key': '2013', 'y': pubyears[2013]},
                        {'key': '2014', 'y': pubyears[2014]},
                        {'key': '2015', 'y': pubyears[2015]},
                        ],
            'key': 'Year',
            }];

    print(pubyear_data)
    try:
        affl=author.one.properties['Affiliation'][0]
    except:
        affl = 'Affiliation Not Available'
    initials = author.one.properties.get('Initials', '')

    return render_template('person.html',
            fullname=fullname,
            firstname=firstname,
            lastname=lastname,
            initials=initials,
            affl=affl,
            keywords=kw,
            labels=author.one.labels,
            categories=biopharmcat,
            pubyear_data=pubyear_data,
            )


@main.route('/_get_pub_list', methods=['GET', 'POST'])
def get_pub_list():
    firstname, lastname = request.args.get('fullname').split('|')
    cypher_command = \
            "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(c) return c;" \
            % (firstname, lastname)
    people = graph.cypher.execute(cypher_command)
    alist = []
    for p in people:
        title = p.c.properties['ArticleTitle']
        abstract = p.c.properties['Abstract']
        pubdate = p.c.properties['ArticleDate']
        category = p.c.properties['BiopharmCategory']
        keywords = p.c.properties['Keywords']
        alist.append([title, abstract, pubdate, keywords, category])

    return jsonify(data=alist)


@main.route('/pub_summary')
def pub_summary():
    ctg = request.args.get('category', 'all')  # category

    def get_data(p, category_data, keyword_data, kw):
        if p.n.properties['ArticleDate']:
            pubdate = datetime.strptime(p.n.properties['ArticleDate'],
                    '%Y-%m-%dT%H:%M:%S')

            category_data[(p.n.properties['BiopharmCategory'])]\
                    [pubdate.year] += 1

            if p.n.properties['Keywords']:
                for word in (p.n.properties['Keywords'].split(',')):
                    kw[word.lower()] += 1
                    keyword_data[word][pubdate.year] += 1

    kw = defaultdict(int)
    category_data = defaultdict(lambda: defaultdict(int))
    keyword_data = defaultdict(lambda: defaultdict(int))

    cypher_command = \
            "MATCH (n:ARTICLE) return n;"
    pubs = graph.cypher.execute(cypher_command)

    for p in pubs:
        if ctg == 'all':
            get_data(p, category_data, keyword_data, kw)
        else:
            if p.n.properties['BiopharmCategory'] == ctg:
                get_data(p, category_data, keyword_data, kw)

    keyword_year_sum = []
    for k,v in sorted(kw.iteritems(), key=lambda (k,v): v, reverse=True)[:10]:
        v = keyword_data[k]
        keyword_year_data = []
        for kk,vv in v.iteritems():
            keyword_year_data.append({'key': kk, 'y': vv})
        if keyword_year_data:
            keyword_year_sum.append({'values': [x for x in keyword_year_data], 'key': k})

    category_year_sum = []
    for k,v in category_data.iteritems():
        category_year_data = []
        for kk,vv in v.iteritems():
            category_year_data.append({'key': kk, 'y': vv})
        if category_year_data:
            category_year_sum.append({'values':
                [x for x in category_year_data], 'key': k})

    return render_template('pubsummary.html',
            category_year_sum=category_year_sum,
            keyword_year_sum=keyword_year_sum,
            category=ctg,
            category_data=category_data,
            )
