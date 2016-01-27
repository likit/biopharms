import py2neo
import datetime
from py2neo import Graph, Node, Relationship, Rev, Path
from nltk.chunk import ChunkParserI
from nltk.chunk.util import conlltags2tree
from nltk.corpus import gazetteers, stopwords
from nltk import word_tokenize, ne_chunk, pos_tag, RegexpParser

# CALAIS_API_KEY = "XsR52pdASdSaxj0k2sbtfngjaI9D8BVh"

class LocationChunker(ChunkParserI):
    def __init__(self):
        self.locations = set(gazetteers.words())
        self.lookahead = 0

        for loc in self.locations:
            nwords = loc.count(' ')

            if nwords > self.lookahead:
                self.lookahead = nwords

    def iob_locations(self, tagged_sent):
        i = 0
        l = len(tagged_sent)
        inside = False
        # print('iob_locations was called.')

        while i < l:
            word, tag = tagged_sent[i]
            j = i + 1
            k = j + self.lookahead
            nextwords, nexttags = [], []
            loc = False
            while j < k:
                if ' '.join([word] + nextwords) in self.locations:
                    if inside:
                        yield word, tag, 'I-LOCATION'
                    else:
                        yield word, tag, 'B-LOCATION'

                    for nword, ntag in zip(nextwords, nexttags):
                        yield nword, ntag, 'I-LOCATION'

                    loc, inside = True, True
                    i = j
                    break

                if j < l:
                    nextword, nexttag = tagged_sent[j]
                    nextwords.append(nextword)
                    nexttags.append(nexttag)
                    j += 1
                else:
                    break
            if not loc:
                inside = False
                i += 1
                yield word, tag, 'O'

    def parse(self, tagged_sent):
        iobs = self.iob_locations(tagged_sent)
        return conlltags2tree(iobs)

NEO4J_PASSWORD = 'Intrinity0'
py2neo.authenticate("localhost:7474", "neo4j", NEO4J_PASSWORD) 
graph = Graph("http://localhost:7474/db/data")

# eng_stopwords = stopwords.words('english')

institutions = {
        'university',
        'department',
        'ministry',
        'institute',
        'research',
        'faculty',
        'center',
        'school',
        'hospital',
        'division',
        }

grammar = """
PLACE: {<DT>?<NN.*>+<IN>?<NN.*>+}
"""

cp = RegexpParser(grammar)

def add_affil(auth, graph, article_date=datetime.datetime.utcnow()):
    '''Adds affiliations to a given author or the db.'''
    # cypher_command = 'MATCH (n:AUTHOR) return n;'
    # found_authors = graph.cypher.execute(cypher_command)
    # print(len(found_authors))
    # for auth in found_authors: break

    # auth = auth.n

    location = LocationChunker()
    for af in auth.properties['Affiliation']:
        for loc in af.split(';'):
            affset = []
            for loc_chunk in loc.split(','):
                loc_chunk_token = word_tokenize(loc_chunk)
                find_loc = True
                for n in loc_chunk_token:
                    if n.lower() in institutions:
                        aff_name = ' '.join(loc_chunk_token)
                        node = graph.merge_one(n.upper(),
                                'name', aff_name)
                        node.labels.add('AFFILIATION')
                        node.push()
                        affset.append(node)
                        find_loc = False
                if find_loc:
                    p = location.parse([(c, 'NN') for c in loc_chunk_token])
                    for subtree in p.subtrees():
                        if subtree.label() == 'LOCATION':
                            loc_name = \
                            ' '.join([x[0] for x in subtree.leaves()])
                            node = graph.merge_one('LOCATION',
                                    'name', loc_name)
                            node.labels.add('AFFILIATION')
                            node.push()
                            affset.append(node)
                            break  # may need to be removed
            if len(affset) > 1:
                for i in range(len(affset)-1):
                    path = Path(affset[i], 'IN', affset[i+1])
                    graph.create(path)

            if len(affset) > 0:
                rel = Relationship(auth, 'IN', affset[0], year=article_date,
                        ambiguous=False)
                graph.create(rel)

if __name__=='__main__':
    add_affil(graph,)
