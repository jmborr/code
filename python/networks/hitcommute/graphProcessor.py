"""
    Program to read and write an igraph graph for drawing images, etc
"""

import os, string, sys

import igraph
from networkx import *
import pylab

class graphProcessor:
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
                if self.edgesDutyCycle.has_key((v1-1, v2-1)):
                    self.edgesDutyCycle[(v1-1, v2-1)] = self.edgesDutyCycle[(v1-1, v2-1)] + 1
                else:
                    self.edgesDutyCycle[(v1-1, v2-1)] = 0.0
            self.sizelist.append(len(self.edgesDutyCycle))
            # print self.edgesDutyCycle.keys()
        except IOError, (errno, stderr):
            print 'Could not open file %s %s' %(errno, stderr)
        else:
            f.close()
            
    def constructGraph(self, fname):
        self.G = igraph.Graph.Read_GML(fname)
        layout = self.G.layout('kk')
        self.weightList = []
        for edge in self.G.get_edgelist():
            self.weightList.append(self.edgesDutyCycle[edge])
        self.MST = self.G.spanning_tree(self.weightList)
        for i in range(0, 161):
            if len(self.MST.neighbors(i)) == 1:
                print i,
        print
        for i in range(0, 161):
            if len(self.MST.neighbors(i)) == 2:
                print i,
        print
        for i in range(0, 161):
            if len(self.MST.neighbors(i)) == 3:
                print i,
        print
        for i in range(0, 161):
            if len(self.MST.neighbors(i)) == 4:
                print i,
        print
        for i in range(0, 161):
            if len(self.MST.neighbors(i)) == 5:
                print i,
        print
        
    
    def makeGraphStructure(self):
        self.G = XGraph()
        keys = self.edgesDutyCycle.keys()
        for key in keys:
            self.G.add_edge(key[0], key[1], self.edgesDutyCycle[key])
        #print dist
        pylab.figure(figsize=(12, 20))
        pos = spring_layout(self.G)
        networkx.draw_networkx_nodes(self.G, pos, node_color='grey', alpha=0.5, node_size=300)
        networkx.draw_networkx_edges(self.G, pos, width=1)
        networkx.draw_networkx_labels(self.G, pos, font_size=9)
        
        MST = XGraph()
        for edge in self.MST.get_edgelist():
            MST.add_edge(edge[0], edge[1])
        networkx.draw_networkx_nodes(MST, pos, node_color='r', alpha=1.0, node_size=300)
        networkx.draw_networkx_edges(MST, pos, width=8, edge_color='g', alpha=0.5)
        
        # networkx.write_dot(self.G, 'dhfr.dot')
        pylab.xticks([])
        pylab.yticks([])
        pylab.savefig("dhfr_all.png")
        #xticks([])
        #ticks([])
        #pylab.show()
           
    def driver(self, skip):
        for i in range(1, self.snap+1, skip):
            self.readAffMatrix(i)
        normFactor = self.snap / skip
        print normFactor
        for key in self.edgesDutyCycle.keys():
            self.edgesDutyCycle[key] = self.edgesDutyCycle[key]/normFactor
        
if __name__=='__main__':
    path = 'contacts/'
    gProcess = graphProcessor(path + '/affinity_' , 7000, '.txt')
    gProcess.driver(10)
    gProcess.constructGraph('dhfr.gml')
    #gProcess.makeGraphStructure()
        