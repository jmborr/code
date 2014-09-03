"""
    Program that does the front end processing of data from 
    computing affinity matrices, statistics about hitting and
    commute times, etc.
"""

from numpy import *
import pylab
from pylab import *

class StatComputer:
    def __init__(self, hitTimeArray):
        try:
            # Process for recievers
            self.avgHitTimes_R = average(hitTimeArray, 1)
            self.stdHitTimes_R = std(hitTimeArray, 1)
            # Process for Broadcasters
            self.avgHitTimes_B = average(hitTimeArray, 0)
            self.stdHitTimes_B = std(hitTimeArray, 0)
        except IOError, message:
            print 'Encountered error in plotter: %s' %(message)
            sys.exit()
            
    def computePDF(self, hitTimeArray, resi):
        self.broadcastHitTimes = hitTimeArray[:][resi]
        print self.broadcastHitTimes.shape
        self.receivingHitTimes = hitTimeArray[resi][:]
        print self.receivingHitTimes.shape
            
    def computeCorrelations(self, hitTimeList, sndRcv = 1):
        # Compute average hitting times
        sz = len(hitTimeList) 
        self.hiList = hitTimeList[:].mean(sndRcv)
        
        #Compute the mean of this hitting times
        self.avgHitTimes = self.hiList.mean(0)
        self.hiList = self.hiList - self.avgHitTimes
        self.correlMat = dot(self.hiList.T, self.hiList)
        # Here we are done with unnormalized correlation matrix
        
        # for variances        
        self.indHitTimes_R = hitTimeList[:].mean(sndRcv)
        tmparray = average(self.indHitTimes_R, 0) ** 2
        self.indHitTimes_R = self.indHitTimes_R**2
        tmp2 = average(self.indHitTimes_R - tmparray, 0)
        
        # normalize the correlation matrix
        for i in range(0, len(tmp2)):
            for j in range(0, len(tmp2)):
                self.correlMat[i][j] = self.correlMat[i][j]/(sz * math.sqrt(tmp2[i] * tmp2[j]))
            
    def writeToFile(self, filename='correl.txt'):
        f = open(filename, 'w')
        for i in range(0, len(self.correlMat)):
            for j in range(0, len(self.correlMat)):
                print >> f, self.correlMat[i][j],
            print >> f
        
    def plotCorrelMatrix(self, filename='fig.png'):
        pylab.contourf(self.correlMat)
        pylab.colorbar()
        pylab.savefig(filename)
        pylab.show()
        
    def plotHitTimes(self):
        pylab.plot(self.avgHitTimes_R, color='blue', lw=2)
        pylab.show()
