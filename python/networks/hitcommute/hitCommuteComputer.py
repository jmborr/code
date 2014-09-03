"""
    Program that takes as input a set of input affinity matrices
    and outputs the Hitting/Commute times of these affinity 
    matrices generated over a course of time series data
    such as MD/ MC simulations.
    
    Affinity matrices are defined as:
    $A_{ij} = \frac{N_{ij}}{\sqrt{N_{i}.N_{j}}$
    
    where: 
    $N_{ij}$ refers to the number of atoms from residue i 
            in contact with residue j
    $N_{i}$ refers to the number of heavy atoms in residue i
    $N_{j}$ refers to the number of heavy atoms in residue j
"""
import os, string, sys
from numpy import *
from scipy import sparse, linalg

class hitCommuteComputer:
    def __init__(self, Aij):
        """
        Constructor
        
        @attention: works only on Affinity matrices defined above
        
        @type: An n x n matrix with Affinity entries computed from either a single structure or an ensemble of structures
        @param: $A_{ij}$ is a matrix
        
        @rtype: None
        @return: A new object of type hitCommuteComputer
        """
        self.size = 0
        self.size = len(Aij)
        assert self.size > 0, 'Bad argument: size: %d' %(self.size)
        self.hitTime = []
        self.comTime = []
        self.laplacian = []
        for i in range(0, self.size):
            self.hitTime.append([])
            self.comTime.append([])
            
    def computeLaplacian(self, Aij):
        """
        Optional method to access the Laplacian of a graph
        @rtype: A matrix of n x n, where n is the number of residues in the matrix
        """
        D = sum(Aij, 1)
        self.laplacian = spdiags(D, 0, len(D), len(D)) - A
        
        return self.laplacian
            
    def computeHitTimes(self, Aij):
        """
        Hit Time computer - main function to calculate the hit times of a graph
        @param Aij: Affinity matrix
        
        @rtype: A matrix of n x n where n is the number of residues in the matrix
        @return: Hit times computed from the graph
        """
        
        D = sum(Aij, 1)    #D is the degree vector
                
        # Compute stationary distribution
        st = D /sum(D)
        
        # Markov transition matrix, row normalized
        P = dot(diag(D**(-1)), Aij)
        W  = ones((len(st),1)) * st.T
        Z = linalg.pinv(eye(P.shape[0], P.shape[1]) - P + W)
        
        Hz = ones((len(st),1)) * diag(Z).T - Z 
        Hz = Hz / (ones((len(st),1)) * st.T)
        Hz = Hz.T
        
        self.hitTime = Hz
        return self.hitTime
    
    def computeCommuteTimes(self):
        """
        Commute Time computer - main function to calculate the commute times of a graph. 
        
        @attention: Commute times are just defined to be the sum of the hitTimes(i, j) + hitTimes(j, i)
        
        @param None: No parameters are expected, since the hitTimes have to be computed
        
        @rtype: A matrix of n x n where n is the number of residues in the matrix
        @return: Hit times computed from the graph
        """
        self.comTime = self.hitTime + self.hitTime.T
        return self.comTime

    """    
    def computeMeanOccupationTimes(self):
        N = []
        for i in range(0, len(self.hitTime)):
            N.append([])
            for j in range(0, len(self.hitTime)):
                N[i].append([])
                for k in range(0, len(self.hitTime)):
                    N[i][j].append([])

        for k in range(0, len(self.hitTime)):
            for i in range(0, len(self.hitTime)):
                for j in range(0, len(self.hitTime)):
                    N[k][i][j] = self.hitTime[k][i] - self.hitTime[k][j] + self.hitTime[j][j] - self.hitTime[i][j]
                    
        return array(N)
    """                
    
    def driver(self, Aij):
        """
        driver routine to compute both hit times and commute times
        @param Aij:  Affinity matrix
        
        @attention: meaningful if and only if Aij is >0 and <= n, where n is the number of residues
        
        @rtype: A list of two arrays Hz, Cz
        @return: hit Times (Hz) and commute times (Cz) of the graph computed
        """
        Hz = self.computeHitTimes(Aij)
        Cz = self.computeCommuteTimes()
        return [Hz, Cz]

      
#if __name__=='__main__':    
#    a = [[5., 3., 4., 2.], [3., 4., 2., 2.], [4., 1., 4., 2.], [2., 3., 2., 5.]]
#    Aij = array(a)
#    print Aij
#    # print Aij.shape
#    hccomputer = hitCommuteComputer(Aij)
#    Hz = hccomputer.computeHitTimes(Aij)
#    print Hz
#   
        