'''Creates adn adds appropriate nodes and relationships for

a given data in Neo4j DB.

'''
import py2neo
from py2neo import Graph, Node, Relationship, Rev, Path
from affiliation import add_affil
from datetime import datetime

#TODO: replace hardcoded values with variables
NEO4J_PASSWORD = 'neo4j'
# py2neo.authenticate("188.166.235.1:7474", "neo4j", NEO4J_PASSWORD) 
# graph = Graph("http://188.166.235.1:7474/db/data")
py2neo.authenticate("127.0.0.1:7474", "neo4j", NEO4J_PASSWORD) 
graph = Graph("http://localhost:7474/db/data")


def create_keyword_nodes(pub):
    labels = ['KEYWORD']
    kw_nodes = []
    if pub.properties['Keywords']:
        for kw in pub.properties['Keywords'].split(','):
            node = Node.cast(labels, {'word': kw})
            kw_nodes.append(node)

        for i in range(len(kw_nodes)-1):
            path = Path(kw_nodes[i], 'RELATE', kw_nodes[i+1])
            graph.create(path)
        r = Relationship(pub, 'HAS', kw_nodes[0])
        graph.create(r)
    else:
        return

def add_pubmed(pub_data, authors):
    # Create publication
    print("\tCreating publication node")
    labels = ['PUBMED', 'ARTICLE']
    node = Node.cast(labels, pub_data)
    pub_node, = graph.create(node)
    # print("Node - ", pub_node)

    # add keyword graph
    create_keyword_nodes(pub_node)
    article_date = pub_data['ArticleDate']

    # Create authors
    print('\tSearching/creating author nodes..')
    labels = ['AUTHOR']
    for author in authors:
        print('\t\tSearching for %s, %s' % \
                (author['LastName'], author['ForeName']))
        cypher_command = \
            'MATCH (n:AUTHOR {LastName: "%s", ForeName: "%s"}) return n;' \
            % (author['LastName'], author['ForeName'])
        found_authors = graph.cypher.execute(cypher_command)
        firstnames = set()
        for auth in found_authors:
            firstnames.add(auth.n.properties['ForeName'].lower())

        if author['ForeName'].lower() not in firstnames:
            print('\t\tNot found.')
            node = Node.cast(labels, author)
            auth, = graph.create(node)
            r = Relationship(auth, 'COAUTHOR', pub_node, ambiguous=False)
            coauthor_path = graph.create(r)
            # add_affil(auth, graph, article_date)
        else:
            print('\t\tFound %d persons.' % len(found_authors))
            for auth in found_authors:
                r = Relationship(auth.n, 'COAUTHOR', pub_node, ambiguous=True)
                auth.n.properties['Affiliation'] += author['Affiliation']
                coauthor_path = graph.create(r)
                # assume that the latest affil is current
                # add_affil(auth.n, graph, article_date)


def add_scopus(pub_data, authors):
    # Create publication
    print("\tCreating publication node")
    labels = ['SCOPUS', 'ARTICLE']
    node = Node.cast(labels, pub_data)
    pub_node, = graph.create(node)
    # print("Node - ", pub_node)

    # add keyword graph
    create_keyword_nodes(pub_node)
    article_date = pub_data['ArticleDate']

    # Create authors
    print('\tSearching/creating author nodes..')
    labels = ['AUTHOR']
    for author in authors:
        print('\t\tSearching for %s, %s' % \
                (author['ForeName'], author['LastName']))
        cypher_command = \
            'MATCH (n:AUTHOR {LastName: "%s", ForeName: "%s"}) return n;' \
            % (author['LastName'], author['ForeName'])
        found_authors = graph.cypher.execute(cypher_command)
        firstnames = set()
        for auth in found_authors:
            firstnames.add(auth.n.properties['ForeName'].lower())

        if author['ForeName'].lower() not in firstnames:
            print('\t\tNot found.')
            node = Node.cast(labels, author)
            auth, = graph.create(node)
            r = Relationship(auth, 'COAUTHOR', pub_node, ambiguous=False)
            coauthor_path = graph.create(r)
            # add_affil(auth, graph, article_date)
        else:
            print('\t\tFound %d persons.' % len(found_authors))
            for auth in found_authors:
                r = Relationship(auth.n, 'COAUTHOR', pub_node, ambiguous=True)
                # auth.n.properties['Affiliation'] += author['Affiliation']
                auth.n.properties['Affiliation'] = []
                coauthor_path = graph.create(r)
                # assume that the latest affil is current
                # add_affil(auth.n, graph, article_date)


if __name__=='__main__':
    pass
    # socialnet = SocialNetwork()
    # create_node_with_label_properties_cast()
    # create_relationship_with_properties()
    # executeSimpleCypherQuery()
    # executeCypherQueryInTransaction()
    # create_paths()
