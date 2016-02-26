import urllib
from datetime import datetime
from . import main
from .. import graph
from collections import defaultdict, namedtuple
from forms import valChainForm

from flask import (Flask, render_template, session,
                    request, redirect, url_for, flash, jsonify)


ALLLABELS = ['PUBMED', 'SCOPUS']
ALL_CATEGORIES = ['vaccine', 'stem cell',
                    'therapeutic antibody',
                    'therapeutic peptide']

@main.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')


@main.route('/top_rank', methods=['GET', 'POST'])
def top_rank():
    labels = graph.cypher.execute("MATCH (n) return labels(n);")
    select_labels = set()
    for lab in labels:
        for l in lab[0]: select_labels.add(l)
    category = request.args.get('category', 'all')  # category
    return render_template('toprank.html',
            alllabels=ALLLABELS, select_labels=list(select_labels),
            all_categories=ALL_CATEGORIES,
            category=category)

@main.route('/_get_toprank_list')
def get_toprank_list():
    # category = 'all'
    category = request.args.get('category', 'all')  # category
    alist = []
    if category == "all":
        command = "start n=node(*) match (n:AUTHOR)--(c:ARTICLE) return n, count(c) as connections order by connections desc limit 50;"
    else:
        command = 'start n=node(*) match (n:AUTHOR)--(c:ARTICLE {BiopharmCategory:"%s"}) return n, count(c) as connections order by connections desc limit 50;' % (category)
    results = graph.cypher.execute(command)

    for r in results:
        firstname = r.n.properties.get('ForeName', '')
        lastname = r.n.properties.get('LastName', '')
        initials = r.n.properties.get('Initials', '')
        pub_count = r.connections
        try:
            affl = r.n.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        if firstname.strip() and lastname.strip():
            # alist.append([firstname, lastname, initials, affl, 0])
            alist.append([firstname, lastname, initials, affl, 0, pub_count])

    return jsonify(data=alist)


@main.route('/top_connectors', methods=['GET', 'POST'])
def top_connectors():
    labels = graph.cypher.execute("MATCH (n) return labels(n);")
    select_labels = set()
    for lab in labels:
        for l in lab[0]: select_labels.add(l)
    category = request.args.get('category', 'all')  # category
    return render_template('top_connectors.html',
            alllabels=ALLLABELS, select_labels=list(select_labels),
            all_categories=ALL_CATEGORIES,
            category=category)

@main.route('/_get_top_connectors_list')
def get_top_connectors_list():
    # category = 'all'
    category = request.args.get('category', 'all')  # category
    alist = []
    if category == "all":
        command = \
                "start n=node(*) MATCH (n:AUTHOR)-[r:COAUTHOR]->(c:ARTICLE)<-[g:COAUTHOR]-(f:AUTHOR) return n, count(distinct(f)) as connections order by connections desc limit 50;"
    else:
        command = \
        "start n=node(*) MATCH (n:AUTHOR)-[r:COAUTHOR]->(c:ARTICLE {BiopharmCategory:'%s'})<-[g:COAUTHOR]-(f:AUTHOR) return n, count(distinct(f)) as connections order by connections desc limit 50;" % (category)
    results = graph.cypher.execute(command)

    if category == "all":
        command = "match (n:AUTHOR {ForeName:\"%s\", LastName:\"%s\"})--(c:ARTICLE) return n, count(c) as numpub;"
    else:
        command = "match (n:AUTHOR {ForeName:\"%s\", LastName:\"%s\"})--(c:ARTICLE {BiopharmCategory:\"%s\"}) return count(c) as numpub;"
    for r in results:
        firstname = r.n.properties.get('ForeName', '')
        lastname = r.n.properties.get('LastName', '')
        initials = r.n.properties.get('Initials', '')
        coauthor_count = r.connections

        if category == 'all':
            results = graph.cypher.execute(command % (firstname, lastname))
        else:
            results = graph.cypher.execute(command % (firstname, lastname, category))
        for c in results:
            numpub = c.numpub
        try:
            affl = r.n.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        if firstname.strip() and lastname.strip():
            # alist.append([firstname, lastname, initials, affl, 0])
            ratio = round(float(coauthor_count)/numpub, 2)
            alist.append([firstname, lastname, initials, affl, 0,
                            coauthor_count, numpub, ratio])

    return jsonify(data=alist)

@main.route('/main_page', methods=['GET', 'POST'])
def main_page():
    labels = graph.cypher.execute("MATCH (n) return labels(n);")
    select_labels = set()
    for lab in labels:
        for l in lab[0]: select_labels.add(l)
    category = request.args.get('category', 'all')  # category
    return render_template('main.html',
            alllabels=ALLLABELS, select_labels=list(select_labels),
            all_categories=ALL_CATEGORIES,
            category=category)

@main.route('/main_page_scopus', methods=['GET', 'POST'])
def main_page_scopus():
    labels = graph.cypher.execute("MATCH (n:SCOPUS) return labels(n);")
    select_labels = set()
    for lab in labels:
        for l in lab[0]: select_labels.add(l)
    category = request.args.get('category', 'all')  # category
    return render_template('main_scopus.html',
            alllabels=ALLLABELS, select_labels=list(select_labels),
            all_categories=ALL_CATEGORIES,
            category=category)


@main.route('/affil', methods=['GET', 'POST'])
def show_affil():
    affil_labels = request.args.get('labels').replace('::', ':')
    affil_name = request.args.get('name')
    affiliation = affil_labels+'|'+affil_name
    labels = graph.cypher.execute("MATCH (n) return labels(n);")
    select_labels = set()
    for lab in labels:
        for l in lab[0]: select_labels.add(l)
    category = request.args.get('category', 'all')  # category
    return render_template('affiliation.html',
            alllabels=ALLLABELS, select_labels=list(select_labels),
            all_categories=ALL_CATEGORIES,
            affiliation=affiliation,
            affil_labels=affil_labels,
            affil_name=affil_name,
            category=category)


@main.route('/_get_affiliation_list')
def get_affil_list():
    # category = 'all'
    print('Affiliation', request.args.get('affiliation'))
    affil_labels, affil_name = request.args.get('affiliation').split('|')
    category = request.args.get('category', 'all')  # category
    alist = []
    if category == "all":
        command = 'START n=node(*) MATCH (n:AUTHOR)-[:IN*..]->(m:%s {name:"%s"}) return distinct(n),m' % (affil_labels, affil_name)
    else:
        command = 'START n=node(*) MATCH (p:ARTICLE {BiopharmCategory:"%s"})<-[COAUTHOR]-(n:AUTHOR)-[:IN*..]->(m:%s {name:"%s"}) return distinct(n),m' % (category, affil_labels, affil_name)

    results = graph.cypher.execute(command)
    for r in results:
        firstname = r.n.properties.get('ForeName', '')
        lastname = r.n.properties.get('LastName', '')
        initials = r.n.properties.get('Initials', '')
        try:
            affl = r.n.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        if firstname.strip() and lastname.strip():
            alist.append([firstname, lastname, initials, affl, 0])

    return jsonify(data=alist)


@main.route('/_get_people_list')
def get_people_list():
    # category = 'all'
    print('Category', request.args.get('category'))
    category = request.args.get('category', 'all')
    if category == "all":
        cypher_command = "MATCH (n:AUTHOR) return distinct(n);"
    else:
        cypher_command = \
            'MATCH (n:AUTHOR)-[COAUTHOR]->(g:ARTICLE {BiopharmCategory: "%s"}) return distinct(n);' % category
    people = graph.cypher.execute(cypher_command)
    alist = []
    for p in people:
        firstname = p.n.properties.get('ForeName', '')
        lastname = p.n.properties.get('LastName', '')
        initials = p.n.properties.get('Initials', '')
        valchain = p.n.properties.get('manualValueChain', [])
        valchain = ','.join(valchain)
        try:
            affl = p.n.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        if firstname.strip() and lastname.strip():
            alist.append([firstname, lastname, initials, affl, valchain, 0])

    return jsonify(data=alist)

@main.route('/_get_people_list_scopus')
def get_people_list_scopus():
    # category = 'all'
    print('Category', request.args.get('category'))
    category = request.args.get('category', 'all')
    if category == "all":
        cypher_command = "MATCH (n:AUTHOR:SCOPUS) return distinct(n);"
    else:
        cypher_command = \
                'MATCH (n:AUTHOR:SCOPUS)-[COAUTHOR]->(g:ARTICLE {BiopharmCategory: "%s"}) return distinct(n);' % category
    people = graph.cypher.execute(cypher_command)
    alist = []
    for p in people:
        firstname = p.n.properties.get('ForeName', '')
        lastname = p.n.properties.get('LastName', '')
        initials = p.n.properties.get('Initials', '')
        valchain = p.n.properties.get('manualValueChain', [])
        valchain = ','.join(valchain)
        try:
            affl = p.n.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        if firstname.strip() and lastname.strip():
            alist.append([firstname, lastname, initials, affl, valchain, 0])

    return jsonify(data=alist)


@main.route('/_get_coauthor_list', methods=['GET', 'POST'])
def get_coauthor_list():
    firstname, lastname = request.args.get('fullname').split('|')
    cypher_command = \
            "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(c)<-[g:COAUTHOR]-(f:AUTHOR) return distinct(f);" % (firstname, lastname)
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

@main.route('/view_keyword')
def view_keyword():
    #TODO: add dname
    keyword = request.args.get('keyword')
    cypher_command = \
        "MATCH (n:KEYWORD {word:\"%s\"})-[]-(ARTICLE)-[]-(a:AUTHOR) return n, a;" \
        % ('HIV')
    authors = graph.cypher.execute(cypher_command)
    coauthor_nodes = []
    coauthor_edges = []
    coauthor_nodes.append(
            {'data': {
                'id': keyword,
                'name': keyword,
                'favShape': 'ellipse',
                'favColor': '#993399'
                }}
            )
    coauthor_dict = {}

    Coauthor = namedtuple('Coauthor', ['firstname', 'lastname'])

    alist = []

    for p in authors:
        firstname = p.a.properties.get('ForeName', '')
        lastname = p.a.properties.get('LastName', '')
        coauthor_key = '%s-%s' % (firstname, lastname)
        if coauthor_key not in coauthor_dict:
            coauthor_dict[coauthor_key] = Coauthor(firstname, lastname)
        initials = p.a.properties.get('Initials', '')
        try:
            affl = p.a.properties['Affiliation'][0]
        except:
            affl = 'Not Available'
        if firstname.strip() and lastname.strip():
            alist.append([firstname, lastname, initials, affl, 0])


    for k, p in coauthor_dict.iteritems():
        print('coauthor', p)
        node = {'data': {
            'id': p.lastname,
            'name': p.firstname + '\n' + p.lastname,
            'favShape': 'ellipse',
            'favColor': '#0066ff',
            'href': urllib.unquote("%s" %
                        url_for('main.view_person', fullname='%s|%s' %
                        (p.firstname, p.lastname)))
                }
            }
        edge = {'data': {
            'id': p.lastname + lastname,
            'source': p.lastname,
            'target': keyword,
            }}
        coauthor_nodes.append(node)
        coauthor_edges.append(edge)

    coauthor_graph = {'nodes': coauthor_nodes, 'edges': coauthor_edges}

    return render_template('keywords.html',
            keyword=keyword,
            coauthor_graph=coauthor_graph,
            alist=alist,
           )

@main.route('/valuechain', methods=['POST', 'GET'])
def value_chain():
    form = valChainForm(request.form)
    if form.validate_on_submit():
        if form.dname.data == 'scopus':
            cypher_command = \
                    "MATCH (n:AUTHOR:SCOPUS {ForeName:'%s', LastName:'%s'}) return n;" \
                    % (form.firstname.data, form.lastname.data)
        else:
            cypher_command = \
                    "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'}) return n;" \
                    % (form.firstname.data, form.lastname.data)
        authors = graph.cypher.execute(cypher_command)
        flash('Data added successfully.')
        for author in authors:
            print(author.n.properties['ForeName'],
                    author.n.properties['LastName'])

            print(author.n.properties.get('manualValueChain', None))
            manualValueChain = set()
            if form.randd.data:
                manualValueChain.add('randd')
            if form.preclin.data:
                manualValueChain.add('preclin')
            if form.clintrials.data:
                manualValueChain.add('clintrials')
            if form.manu.data:
                manualValueChain.add('manu')
            if form.sale.data:
                manualValueChain.add('sale')
            author.n.properties['manualValueChain'] = manualValueChain
            print(author.n.properties.get('manualValueChain', None))
            author.n.push()
        return redirect(url_for('main.view_person', dname=form.dname.data,
            fullname=('%s|%s' % (form.firstname.data, form.lastname.data))))


@main.route('/view_person', methods=['GET'])
def view_person():
    fullname = request.args.get('fullname')
    dname = request.args.get('dname')
    print(fullname, dname)
    firstname, lastname = fullname.split('|')
    if dname == 'scopus':
        cypher_command = \
                "MATCH (n:AUTHOR:SCOPUS {ForeName:'%s', LastName:'%s'}) return n;" \
                % (firstname, lastname)
    else:
        cypher_command = \
                "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'}) return n;" \
                % (firstname, lastname)
    author = graph.cypher.execute(cypher_command)

    # cypher_command = \
    #         "start n = node(*) match p=(n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:IN]->(m) return n,m ORDER BY r.year DESC LIMIT 5;" % (firstname, lastname)
    # curr_affil = graph.cypher.execute(cypher_command)
    # for c in curr_affil:
    #     print(c)
    # all_affils = []
    # for c in curr_affil:
    #     command = 'start n = node(*) match p=(n:%s {name:"%s"})-[r:IN*..3]->(m) return n,m' \
    #                     % (':'.join(c.m.labels), c.m.properties['name'])
    #     print(command)
    #     results = graph.cypher.execute(command)
    #     print(results)
    #     affil = set()
    #     for a in results:
    #         affil.add(('::'.join(a.m.labels), a.m.properties['name']))
    #     affil.add(('::'.join(a.n.labels), a.n.properties['name']))
    #     all_affils.append(affil)

    # for a in all_affils:
    #     print(a)

    if dname == 'scopus':
        cypher_command = \
            "MATCH (n:AUTHOR:SCOPUS {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(f:ARTICLE) return f;" \
            % (firstname, lastname)
    else:
        cypher_command = \
            "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(f:ARTICLE) return f;" \
            % (firstname, lastname)
    pubs = graph.cypher.execute(cypher_command)
    kw=set()
    biopharmcat = set()
    pubyears = {}
    #TODO: replace a hardcoded range
    for y in range(2005,2016):
        pubyears[y] = 0

    for p in pubs:
        if p.f.properties['Keywords']:
            for word in (p.f.properties['Keywords'].split(',')):
                kw.add(word)
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

    # print(pubyear_data)
    try:
        affl=author.one.properties['Affiliation'][0]
    except:
        affl = 'Affiliation Not Available'
    initials = author.one.properties.get('Initials', '')
    valuechain_keys = {'clintrials': 'Clinical trials',
            'preclin': 'Preclinical research',
            'manu': 'Manufacturing',
            'sale': 'Marketing and Sale',
            'randd': 'Research and Development'}
    valuechain = set()
    for v in author.one.properties.get('manualValueChain', set()):
        valuechain.add(valuechain_keys[v])

    if dname == 'scopus':
        cypher_command = \
                "MATCH (n:AUTHOR:SCOPUS {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(c)<-[g:COAUTHOR]-(f:AUTHOR) return f;" % (firstname, lastname)
    else:
        cypher_command = \
                "MATCH (n:AUTHOR {ForeName:'%s', LastName:'%s'})-[r:COAUTHOR]->(c)<-[g:COAUTHOR]-(f:AUTHOR) return f;" % (firstname, lastname)

    coauthors = graph.cypher.execute(cypher_command)
    coauthor_nodes = []
    coauthor_edges = []
    coauthor_nodes.append(
            {'data': {
                'id': lastname,
                'name': firstname + '\n' + lastname,
                'favShape': 'ellipse',
                'favColor': '#993399'
                }}
            )
    coauthor_dict = {}
    coauthor_count = defaultdict(int)
    coauth_id = 0
    Coauthor = namedtuple('Coauthor', ['firstname', 'lastname'])
    for p in coauthors:
        cofirstname = p.f.properties.get('ForeName', '')
        colastname = p.f.properties.get('LastName', '')
        coauthor_key = '%s-%s' % (cofirstname, colastname)
        coauthor_count[coauthor_key] += 1
        if coauthor_key not in coauthor_dict:
            coauthor_dict[coauthor_key] = Coauthor(cofirstname, colastname)

    top_coauthors = [(coauthor_dict[k], coauthor_count[k]) for k in
                        coauthor_dict]
    top_coauthors = sorted(top_coauthors,
                            key=lambda x: x[1], reverse=True)
    if len(top_coauthors) > 30:
        top_coauthors = top_coauthors[:30]

    for item in top_coauthors:
        p = item[0]
        node = {'data': {
            'id': p.lastname,
            'name': p.firstname + '\n' + p.lastname,
            'favShape': 'ellipse',
            'favColor': '#0066ff',
            'href': urllib.unquote("%s" %
                        url_for('main.view_person', dname=dname,
                            fullname='%s|%s' % (p.firstname, p.lastname)))
                }
            }
        edge = {'data': {
            'id': p.lastname + lastname,
            'source': lastname,
            'target': p.lastname,
            'weight': item[1],
            'label': str(item[1])
            }}
        coauthor_nodes.append(node)
        coauthor_edges.append(edge)

    coauthor_graph = {'nodes': coauthor_nodes, 'edges': coauthor_edges}
    valchain_form = valChainForm()
    # if valchain_form.validate_on_submit():
    #     print(valchain_form.data.randd)

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
            # affiliations=all_affils,
            coauthor_graph=coauthor_graph,
            form=valchain_form,
            valuechain=valuechain,
            dname=dname
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

    # TODO: A list of categories should be stored somewhere
    all_categories = set()  # store all categories from the database

    for p in pubs:
        if ctg == 'all':
            get_data(p, category_data, keyword_data, kw)
        else:
            if p.n.properties['BiopharmCategory'] == ctg:
                get_data(p, category_data, keyword_data, kw)

        all_categories.add(p.n.properties['BiopharmCategory'])

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

    for c in category_year_sum:
        print c
    return render_template('pubsummary.html',
            category_year_sum=category_year_sum,
            keyword_year_sum=keyword_year_sum,
            category=ctg,
            category_data=category_data,
            all_categories=sorted(all_categories),
            )

@main.route('/show_value_chain/', methods=['GET'])
def show_value_chain():
    ctg = request.args.get('category', 'all')  # category
    if ctg == 'all':
        cypher_command = \
            'MATCH (n:ARTICLE)-[:COAUTHOR]-(a:AUTHOR) where has(n.BiopharmCategory) return distinct(a);'
    else:
        cypher_command = \
            'MATCH (n:ARTICLE)-[:COAUTHOR]-(a:AUTHOR) where has(n.BiopharmCategory) and n.BiopharmCategory="%s" return distinct(a);' % ctg

    authors = graph.cypher.execute(cypher_command)
    value_chain = defaultdict(int)
    alist = defaultdict(list)
    for author in authors:
        if author.a.properties['manualValueChain']:
            for v in author.a.properties['manualValueChain']:
                value_chain[v] += 1
                firstname = author.a.properties.get('ForeName', '')
                lastname = author.a.properties.get('LastName', '')
                initials = author.a.properties.get('Initials', '')
                try:
                    affl = author.a.properties['Affiliation'][0]
                except:
                    affl = 'Not Available'
                alist[v].append([firstname, lastname, initials, affl, 0])

    value_chain_keys = {
            'preclin': 'Preclinical research',
            'clintrials': 'Clinical trials',
            'randd': 'Research and development',
            'manu': 'Manufacturing',
            'sale': 'Marketing and Sale'
            }
    value_chain_data = [{
            'values': [{'key': value_chain_keys[k],
                'value': value_chain[k]} for k in
                    ['randd', 'preclin', 'clintrials', 'manu', 'sale']],
            'key': 'Phase',
            }];
    print(value_chain_data)

    return render_template('value_chain.html',
            all_categories=sorted(ALL_CATEGORIES),
            value_chain_data=value_chain_data,
            category=ctg,
            alist=alist,
            )

@main.route('/about/')
def about():
    return render_template('about.html')
