'''Creates adn adds appropriate nodes and relationships for

a given data in Neo4j DB.

'''
import py2neo
from py2neo import Graph, Node, Relationship, Rev, Path
from affiliation import add_affil
from datetime import datetime

NEO4J_PASSWORD = 'Intrinity0'
py2neo.authenticate("localhost:7474", "neo4j", NEO4J_PASSWORD) 
graph = Graph("http://localhost:7474/db/data")


def add_pubmed(pub_data, authors):
    # Create publication
    print("\tCreating publication node")
    labels = ['PUBMED', 'ARTICLE']
    node = Node.cast(labels, pub_data)
    pub_node, = graph.create(node)
    # print("Node - ", pub_node)

    article_date = datetime.strptime(pub_data['ArticleDate'],
                                                '%Y-%m-%dT%H:%M:%S')
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
        if len(found_authors) == 0:
            print('\t\tNot found.')
            node = Node.cast(labels, author)
            auth, = graph.create(node)
            r = Relationship(auth, 'COAUTHOR', pub_node, ambiguous=False)
            coauthor_path = graph.create(r)
            add_affil(auth, graph, article_date)
        elif len(found_authors) >= 1:
            print('\t\tFound %d persons.' % len(found_authors))
            for auth in found_authors:
                r = Relationship(auth.n, 'COAUTHOR', pub_node, ambiguous=True)
                auth.n.properties['Affiliation'] += author['Affiliation']
                coauthor_path = graph.create(r)
                add_affil(auth, graph, article_date)


if __name__=='__main__':
    pass
    # socialnet = SocialNetwork()
    # create_node_with_label_properties_cast()
    # create_relationship_with_properties()
    # executeSimpleCypherQuery()
    # executeCypherQueryInTransaction()
    # create_paths()
