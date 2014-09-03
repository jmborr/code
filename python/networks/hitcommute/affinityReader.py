"""
    Program to read the affinity matrices from a trajectory and 
    process it to obtain the Hitting/ Commute Time values
"""

import os, string, sys
from numpy import *
from scipy import sparse
from networkx import *
import pylab
from pylab import *


class affinityReader:
    def __init__(self, fprefix, nsnap, ext):
        self.fprefix = fprefix
        self.ext = ext
        self.nsnap = nsnap
        self.affinityMatrix = []
        
    def readMatrix(self, snap):
        fname = self.fprefix + str(snap) + self.ext
        print fname
        try:
            afm = []
            f = open(fname, 'r')
            for line in f.readlines():
                contents = line.strip().split()
                tmp = []
                for content in contents:
                    tmp.append(float(content))
                afm.append(tmp)
            self.affinityMatrix.append(afm)
        except IOError, (errno, stderr):
            print 'Could not open file: %s, %s' %(errno, stderr)
        else:
             f.close()
    
    def readSparseMatrix(self, snap):
        fname = self.fprefix + str(snap) + self.ext
        print fname
        try:
            f = open(fname, 'r')
            I = []
            J = []
            V = []
            for line in f.readlines():
                contents = line.strip().split()
                I.append(int(contents[0])-1)
                J.append(int(contents[1])-1) 
                V.append(float(contents[2]))
            Aij = zeros((max(I) + 1, max(I) + 1), float)
            for i in range(0, len(I)):
                Aij[I[i]][J[i]] = V[i]
            self.affinityMatrix.append(Aij)
        except IOError, (errno, stderr):
            print 'Could not open file: %s, %s' %(errno, stderr)
        else:
             f.close()   
            
             
    def driver(self, skip):
        for i in range(1, self.nsnap+1, skip):
            #self.readMatrix(i)
            self.readSparseMatrix(i)
            
    def getAffinityMatrix(self, i):
        Aij = self.affinityMatrix[i]
        return Aij
    
    def getAverageAffinityMatrix(self):
        sum = zeros((len(self.affinityMatrix[0]), len(self.affinityMatrix[0])))
        for afm in self.affinityMatrix:
            sum = sum + afm
        sum = sum/ len(self.affinityMatrix)
        f = open('avgAff_matrix.txt', 'w')
        for i in range(0, len(sum)):
            for j in range(0, len(sum[i])):
                print >> f, sum[i][j],
            print >> f
        f.close()
        return sum
    
    def makeGraph(self, avgAFM, target):
        G = XGraph()
        for i in range(0, len(avgAFM)):
            for j in range(i, len(avgAFM)):
                if(avgAFM[i][j] > 0.0):
                    G.add_edge(i, j, 1.0/avgAFM[i][j])
        #pylab.figure(figsize=(12, 20))
        #pos = spring_layout(G)
        #networkx.draw_networkx_nodes(G, pos, node_color='grey', alpha=0.5, node_size=300)
        #networkx.draw_networkx_edges(G, pos, width=1.0)
        #networkx.draw_networkx_labels(G, pos, font_size=9)
        
        for i in range(0, len(avgAFM)):
            #path = []
            #distance = 0.0
            (distance, path) = networkx.path.bidirectional_dijkstra(G, i, target)
            print path, len(path)
            #PGraph = Graph()
            #for i in range(0, len(path)-1):
             #   PGraph.add_edge(path[i], path[i+1])
            #networkx.draw_networkx_nodes(PGraph, pos, node_color='b', alpha=1.0, node_size=300)
            #networkx.draw_networkx_edges(PGraph, pos, width=8, edge_color='g', alpha=0.5)
            #pylab.xticks([])
            #pylab.yticks([])
            #pylab.savefig('dhfr_' + str(i) + '_' + str(target)+ '.png')
            
if __name__ == '__main__':
    from hitCommuteComputer import hitCommuteComputer
    path = 'contacts/'
    affReader = affinityReader( path + ('affinity_'), 7000, '.txt')
    affReader.driver(10)
    hitTimeArray = []
    comTimeArray = []
    avgAFM = affReader.getAverageAffinityMatrix()
    affReader.makeGraph(avgAFM, 160)
    
#    for i in range(0, len(affReader.affinityMatrix)):
#        Aij = affReader.getAffinityMatrix(i)
#        # print Aij.shape
#        hcCalculator = hitCommuteComputer(Aij)
#        [Hz, Cz] = hcCalculator.driver(Aij)
#        hitTimeArray.append(Hz)
#        comTimeArray.append(Cz)
#    avgHitTimes = average(array(hitTimeArray), 0)
#    print avgHitTimes.mean(1).tolist()
#    print avgHitTimes.std(1).tolist()
#    avgComTimes = average(array(comTimeArray), 0)
#    
#    from statComputer import StatComputer
#    sc = StatComputer(avgHitTimes)
#    sc.computeCorrelations(array(hitTimeArray), 2)
##    sc.computePDF(hitTimeArray[0], 41)
##    sc.writeToFile('1RX2-correl.txt')
#    sc.plotCorrelMatrix('1RX2-HitCommuteCorrel.png')
##    sc.plotHitTimes()