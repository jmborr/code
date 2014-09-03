"""
    Program to read and see how distance variations affect shortest paths in protein structures
"""

import os, string, sys
from numpy import *
from scipy import sparse
from networkx import *
import pylab
from pylab import *

class distanceReader:
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
                if V[i] > 0.0:
                    Aij[I[i]][J[i]] = 1
            self.affinityMatrix.append(Aij)
        except IOError, (errno, stderr):
            print 'Could not open file: %s, %s' %(errno, stderr)
        else:
             f.close()              
             
    def driver(self, skip):
        for i in range(1, self.nsnap+1, skip):
            self.readSparseMatrix(i)
            
    def getAverageAffinityMatrix(self):
        sum = zeros((len(self.affinityMatrix[0]), len(self.affinityMatrix[0])))
        for afm in self.affinityMatrix:
            sum = sum + afm
        sum = sum/ len(self.affinityMatrix)
        f = open('dist_0_1_matrix.txt', 'w')
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
                if(avgAFM[i][j] >= 1.0):
                    G.add_edge(i, j, 1.0)
 
        for i in range(0, len(avgAFM)):
            (distance, path) = networkx.path.bidirectional_dijkstra(G, i, target)
            print path, len(path)
            
if __name__ == '__main__':
    from hitCommuteComputer import hitCommuteComputer
    path = 'contacts/'
    affReader = distanceReader( path + ('affinity_'), 7000, '.txt')
    affReader.driver(10)
    hitTimeArray = []
    comTimeArray = []
    avgAFM = affReader.getAverageAffinityMatrix()
    affReader.makeGraph(avgAFM, 160)