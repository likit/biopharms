'''Creates adn adds appropriate nodes and relationships for

a given data in Neo4j DB.

'''
import py2neo
from py2neo import Graph, Node, Relationship, Rev, Path

NEO4J_PASSWORD = 'Intrinity0'
py2neo.authenticate("localhost:7474", "neo4j", NEO4J_PASSWORD) 
graph = Graph("http://localhost:7474/db/data")


def create_node_with_label_properties_cast():
    labels = ['FirstLabel', 'SecondLabel']

    properties = {'name': 'MyPython3'}
    node = Node.cast(labels, properties)
    result_node, = graph.create(node)
    print("Node - ", result_node)


def create_relationship_with_properties():
    amy = Node('FEMALE', name='Amy')
    kristine = Node('FEMALE', name='Kristine')
    sheryl = Node('FEMALE', name="Sheryl")
    kristine_amy = Relationship(kristine, 'FRIEND', amy, size=2005)
    amy_sheryl = Relationship(amy, Rev('FRIEND'), sheryl, since=2001)
    result_nodes = graph.create(kristine_amy, amy_sheryl)


def executeSimpleCypherQuery():
    results = graph.cypher.execute("MATCH (n) return n.name as name, labels(n) as labels")
    for index in range(len(results)):
        record = results[index]
        print("Printing Record = ", index, " - Name =",
                record.name, ", Labels = ", record.labels)
        print("End - execution of Simple Cypher Query")


def executeCypherQueryInTransaction():
    tx = graph.cypher.begin()
    tx.append("CREATE (n:Node1{name:'John'}) RETURN n")
    tx.append("CREATE (n:Node1{name:'Russel1'}) RETURN n")
    tx.append("CREATE (n:Node1{name:'Smith'}) RETURN n")
    results = tx.commit()
    for result in results:
        for record in result:
            print(record.n)
    print("End - execution of Cypher Query in Transaction")


def create_paths():
    bradley, matthew, lisa = Node(name="Bradley"), Node(name="Matthew"), \
                                Node(name="Lisa")
    path_1 = Path(bradley, "Knows", matthew, Rev("Knows"), lisa)
    graph.create(path_1)
    john, annie, ripley = Node(name="John"), Node(name="Annie"), Node(name="Ripley")
    path_2 = Path(john, "Knows", annie, "Knows", ripley)
    path_3 = path_1.append("Knows", path_2)
    result_path = graph.create(path_3)

    print('Print raw data')
    print('Nodes in the Path-1 = ', result_path[0].nodes)
    print('Relationships in the Path-1 = ', result_path[0].relationships)
    print('Print - all relationships')
    for rels in result_path[0].relationships:
        print(rels)


class SocialNetwork(object):
    def __init__(self):
        py2neo.authenticate("localhost:7474", "neo4j", NEO4J_PASSWORD) 
        self.graph = Graph("http://localhost:7474/db/data")
        self.people = self.create_people(self.graph)
        self.friend_path = self.create_friends(self.graph, self.people)
        self.movie = self.create_movies(self.graph)
        self.rate_movies(self.graph, self.movie, self.people)

    def create_people(self, graph):
        print('Creating people')
        bradley = Node('MALE', 'TEACHER', name='Bradley',
                surname='Green', age=24, country='US')
        matthew = Node('MALE', 'STUDENT', name='Matthew',
                surname='Cooper', age=36, country='US')
        lisa = Node('FEMALE', name='Lisa', surname='Adams',
                age=15, country='Canada')
        john = Node('MALE', name='John', surname='Godman',
                age=24, country='Mexico')
        annie = Node('FEMALE', name='Annie', surname='Behr',
                age=25, country='Canada')
        ripley = Node('MALE', name='Ripley',
                surname='Aniston', country='US')
        graph.create(bradley, matthew, lisa, john, annie, ripley)
        print('People Created')
        people = {'bradley':bradley, 'matthew':matthew, 'lisa':lisa, 'john':john,
                'annie':annie, 'ripley':ripley}
        return people


    def create_friends(self, graph, people):
        print('Creating relationships between people')
        path1 = \
                Path(people['bradley'], 'FRIEND', people['matthew'], 'FRIEND',
                        people['lisa'], 'FRIEND', people['john'])
        path2 = \
                path1.prepend(people['lisa'], Rev('FRIEND'))
        path3 = \
                Path(people['annie'], 'FRIEND', people['ripley'], 'FRIEND',
                        people['lisa'])
        path4 = \
                Path(people['bradley'], 'TEACHES', people['matthew'])
        friends_path = graph.create(path2, path3, path4)
        return friends_path


    def create_movies(self, graph):
        print("Creating Movies")
        firstBlood = Node('MOVIE', name='First Blood')
        avengers = Node('MOVIE', name="Avengers")
        matrix = Node('MOVIE', name='matrix')
        graph.create(firstBlood, avengers, matrix)
        print("Movies created")
        return {'firstBlood': firstBlood, 'avengers':avengers, 'matrix':matrix}


    def rate_movies(self, graph, movies, people):
        print("Start rating movies")
        matthew_firstBlood = \
                Relationship(people['matthew'], 'HAS_RATED',
                        movies['firstBlood'], ratings=4)
        john_firstBlood = \
                Relationship(people['john'], 'HAS_RATED',
                        movies['firstBlood'], ratings=4)
        annie_firstBlood = \
                Relationship(people['annie'], 'HAS_RATED',
                        movies['firstBlood'], ratings=4)
        ripley_firstBlood = \
                Relationship(people['ripley'], 'HAS_RATED',
                        movies['firstBlood'], ratings=4)
        lisa_avengers = \
                Relationship(people['lisa'], 'HAS_RATED',
                        movies['avengers'], ratings=5)
        matthew_avengers = \
                Relationship(people['matthew'], 'HAS_RATED',
                        movies['avengers'], ratings=4)
        annie_avengers = \
                Relationship(people['annie'], 'HAS_RATED',
                        movies['avengers'], ratings=4)

        movie_path = graph.create(matthew_firstBlood, john_firstBlood,
                annie_firstBlood, ripley_firstBlood, lisa_avengers,
                matthew_avengers, annie_avengers)
        print("Finished rating the movies")
        return movie_path


def add_pubmed(pub_data, authors):
    # Create publication
    print("\tCreating publication node")
    labels = ['PUBMED', 'ARTICLE']
    node = Node.cast(labels, pub_data)
    pub_node, = graph.create(node)
    # print("Node - ", pub_node)

    # Create authors
    print('\tSearching/creating author nodes..')
    labels = ['AUTHOR']
    for author in authors:
        print('\t\tSearching for %s, %s' % \
                (author['LastName'], author['ForeName']))
        cypher_command = \
            "MATCH (n:AUTHOR {LastName: '%s', ForeName: '%s'}) return n;" \
            % (author['LastName'], author['ForeName'])
        found_authors = graph.cypher.execute(cypher_command)
        if len(found_authors) == 0:
            print('\t\tNot found.')
            node = Node.cast(labels, author)
            auth_node, = graph.create(node)
            r = Relationship(auth_node, 'COAUTHOR', pub_node, ambiguous=False)
            coauthor_path = graph.create(r)
        elif len(found_authors) >= 1:
            print('\t\tFound %d persons.' % len(found_authors))
            for auth in found_authors:
                r = Relationship(auth.n, 'COAUTHOR', pub_node, ambiguous=True)
                auth.n.properties['Affiliation'] += author['Affiliation']
                coauthor_path = graph.create(r)


if __name__=='__main__':
    pass
    # socialnet = SocialNetwork()
    # create_node_with_label_properties_cast()
    # create_relationship_with_properties()
    # executeSimpleCypherQuery()
    # executeCypherQueryInTransaction()
    # create_paths()
