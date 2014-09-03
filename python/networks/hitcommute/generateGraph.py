"""
    Program that generates a graph like representation for use 
    with the graphviz package. Format is .dot
"""

import os, string, sys
from networkx import *
import pylab

class GraphGenerator:
    def __init__(self, prefix, snap, ext):
        self.prefix = prefix
        self.snap = snap
        self.ext = ext
        self.sizelist = []
        self.edgesDutyCycle = {}
        
    def readAffMatrix(self, snap):
        fname = self.prefix + str(snap) + self.ext
        try:
            f = open(fname, 'r')
            for line in f.readlines():
                contents = line.strip().split()
                v1 = int(contents[0])
                v2 = int(contents[1])
                if self.edgesDutyCycle.has_key((v1, v2)):
                    self.edgesDutyCycle[(v1, v2)] = self.edgesDutyCycle[(v1, v2)] + 1
                else:
                    self.edgesDutyCycle[(v1, v2)] = 0.0
            self.sizelist.append(len(self.edgesDutyCycle))
        except IOError, (errno, stderr):
            print 'Could not open file %s %s' %(errno, stderr)
        else:
            f.close()
            
    def makeGraphStructure(self):
        self.G = XGraph()
        keys = self.edgesDutyCycle.keys()
        for key in keys:
            self.G.add_edge(key[0], key[1], self.edgesDutyCycle[key])
        (dist, path) = networkx.path.floyd_warshall_array(self.G)
        #print dist
        pylab.figure(figsize=(8, 11))
        pos = spring_layout(self.G)
        networkx.draw_networkx_nodes(self.G, pos, node_color='grey', alpha=0.5, node_size=300)
        networkx.draw_networkx_edges(self.G, pos, width=1)
        networkx.draw_networkx_labels(self.G, pos, font_size=9)
        # networkx.write_dot(self.G, 'dhfr.dot')
        pylab.xticks([])
        pylab.yticks([])
        pylab.savefig("dhfr.png")
        #xticks([])
        #ticks([])
        #pylab.show()
        
            
    def writeGraph(self, outFilename, format='dot'):
        print 'Generating graphviz compatible format'
        keys = self.edgesDutyCycle.keys()
        keystrings = []
        
        # REMOVE DUPLICATE KEYS
        for key_pair1 in keys:
            for key_pair2 in keys:
                if key_pair1[0] == key_pair2[1] and key_pair1[1] == key_pair2[0]:
                    keys.remove(key_pair2)
                    
        colorDict = {}
        for key in keys:
            if self.edgesDutyCycle[key] >= 0.0 and self.edgesDutyCycle[key] <= 0.1:
                colorDict[key] = 'violet'
            elif self.edgesDutyCycle[key] > 0.1 and self.edgesDutyCycle[key] <= 0.2:
                colorDict[key] = 'indigo'
            elif self.edgesDutyCycle[key] > 0.2 and self.edgesDutyCycle[key] <= 0.3:
                colorDict[key] = 'blue1'
            elif self.edgesDutyCycle[key] > 0.3 and self.edgesDutyCycle[key] <= 0.4:
                colorDict[key] = 'blue4'
            elif self.edgesDutyCycle[key] > 0.4 and self.edgesDutyCycle[key] <= 0.5:
                colorDict[key] = 'green'
            elif self.edgesDutyCycle[key] > 0.5 and self.edgesDutyCycle[key] <= 0.6:
                colorDict[key] = 'yellowgreen'
            elif self.edgesDutyCycle[key] > 0.6 and self.edgesDutyCycle[key] <= 0.7:
                colorDict[key] = 'yellow'
            elif self.edgesDutyCycle[key] > 0.7 and self.edgesDutyCycle[key] <= 0.8:
                colorDict[key] = 'orange'
            elif self.edgesDutyCycle[key] > 0.8 and self.edgesDutyCycle[key] <= 0.9:
                colorDict[key] = 'orangered'
            elif self.edgesDutyCycle[key] > 0.9 and self.edgesDutyCycle[key] <= 1.0:
                colorDict[key] = 'red'
    
        # PRINT GRAPH
        if format == 'dot':
            self.writeDotFile(outFilename, keys, colorDict)
        elif format == 'gml':
            self.writeGMLFile(outFilename, keys, colorDict)
        else:
            print 'Unknown format asked for'
            
    def writeGMLFile(self, outFilename, keys, colorDict):
        f = open(outFilename, 'w')
        print >> f, 'graph ['
        print >> f, '\tdirected 0'
        for i in range(1, 162):
            print >> f, '\tnode ['
            print >> f, '\t\t id ' + str(i)
            print >> f, '\t\t label \"' + str(i) + '\"'
            print >> f, '\t]'
            print >> f
        for key in keys:
            print >> f, '\tedge ['
            print >> f, '\t\t source ' + str(key[0])
            print >> f, '\t\t target ' + str(key[1])
            print >> f, '\t\t label ' + '\"' + str(self.edgesDutyCycle[key]) + '\"'
            print >> f, '\t]'
        print >> f, ']'
        f.close()
    
    def writeDotFile(self, outFilename, keys, colorDict):
        f = open(outFilename, 'w')
        self.printHeader(f)
        print >> f, 'graph ' + outFilename[0:-4] + '{'
        print >> f, '\tnode ['
        print >> f, '\t\tfontname = \"Arial\"'
        print >> f, '\t\tshape = \"circle\"'
        print >> f, '\t\twidth = \"0.50000\"'
        print >> f, '\t\theight = \"0.500000\"'
        print >> f, '\t\tcolor = \"black\"'
        print >> f, '\t]'
        
        for key in keys:
##            if (key[0] != key[1]+1) and (key[0] != key[1]-1) and (key[1] != key[0]+1) and (key[1] != key[0]-1) and (key[0] != key[1]+2)and (key[0] != key[1]-2) and (key[1] != key[0]+2) and (key[1] != key[0]-2):        
            print >> f, '\t \"' + str(key[0]) + '\" -- \"' + str(key[1]) #+ '\" [color =\"' + colorDict[key] + '\"];' 
        print >> f, '}'
        f.close()
#        g1 = pydot.graph_from_edges(keystrings)
#        g1.write_png('dhfr.png', prog='dot')
#        g1.write('dhfr.dot', prog='dot', format='raw')     
    
    def printHeader(self, f):
        print >> f, '/* Output generated by generateGraph.py. Program that '
        print >> f, ' * outputs a graph in the dot format for further processing.'
        print >> f, ' * Author : Arvind Ramanathan'
        print >> f, ' * Date: '
        print >> f, ' */'
        
    
    def driver(self, skip):
        for i in range(1, self.snap+1, skip):
            self.readAffMatrix(i)
        normFactor = self.snap / skip
        print normFactor
        for key in self.edgesDutyCycle.keys():
            self.edgesDutyCycle[key] = self.edgesDutyCycle[key]/normFactor
            
    def buildGraph(self):
        keys = self.edgesDutyCycle.keys()
        keystrings = []
        
        # REMOVE DUPLICATE KEYS
        for key_pair1 in keys:
            for key_pair2 in keys:
                if key_pair1[0] == key_pair2[1] and key_pair1[1] == key_pair2[0]:
                    keys.remove(key_pair2)
        
        # FOR EVERY KEY BUILD the adjacency representation
        graphList = []
        for i in range(0, 161):
            graphList.append([])
        for i in range(1, 162):
                for key in keys:
                    if key[0] == i:
                        graphList[i-1].append({key[1]:self.edgesDutyCycle[key]})
        graph = {}
        for i in range(0, 161):
            graph[i+1] = graphList[i]
        print graph
        return graphList


    
if __name__ == '__main__':
    path = 'contacts/'
    g = GraphGenerator(path + 'affinity_', 7000, '.txt')
    g.driver(10)
    # g.makeGraphStructure()
    g.writeGraph('dhfr.gml', 'gml')
    # print g.buildGraph()
    # d, p = g.bellman_ford(g.edgesDutyCycle, 161)
    #print d
    #print p

"""    
    from numpy import *
    import pylab
    from pylab import *
    pylab.plot(g.sizelist, color='blue', lw=2)
    pylab.show()
"""