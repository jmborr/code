"""
    Program to read a PDB file, parse it to obtain 
    an affinity matrix for the protein of interest.
"""
from numpy import *

class PDBReader:
    def __init__(self, filename):
        atmMap = {'C' : 12.00, 'N' : 14.00, 'O' : 16.00, 'H' : 1.01, 'P' : 30.0, 'S' : 32.00}
        self.x_Coords = []
        self.y_Coords = []
        self.z_Coords = []
        self.masses = []
        self.reslist = []
        try:
            f = open(filename, 'r')
            i = 0
            for line in f.readlines():
                if line[0:4] == 'ATOM' or line[0:6] =='HETATM':
                    
                    #Process residues
                    resnum = int(line[23:30])
                    self.reslist.append(resnum)
                    
                    # Process coordinates
                    x = float(line[30:37])
                    y = float(line[38:45])
                    z = float(line[46:53])
                    self.x_Coords.append(x)
                    self.y_Coords.append(y)
                    self.z_Coords.append(z)
                    
                    # Process residue numbers
                    atmName = line[77:78]
                    self.masses.append(atmMap[atmName])
                else:
                    i = i + 1
        except IOError, (errno, stderr):
            print 'Could not open file: %s %s' %(errno, stderr)
        finally:
            f.close()
            print 'No. of atoms read: %d' %(len(self.reslist))

                                                  
    def getCoords(self):
        return [self.x_Coords, self.y_Coords, self.z_Coords]
    
    def computeAffinityMatrix(self):
        natm = len(self.x_Coords)
        nres = 315#max(self.reslist) - min(self.reslist) + 1
        print 'No. of residues: %d' %(nres)
        self.affinityMatrix = zeros((nres, nres), float)
        
        # Compute where one residue begins and another ends
        pnt = [0] * nres
        tmpres = 0
        pnt[0] = 0
        for i in range(0, natm):
            if(self.reslist[i] > self.reslist[i-1]):
                tmpres = tmpres + 1
                pnt[tmpres] = i
        print pnt

        #Now that we have markers, compute affinity
        ni = 0
        nj = 0
        for i in range(0, nres):
            alist_i = []
            stp1 = 0
            if i == nres-1:
                stp1 = natm
            else:
                stp1 = pnt[i+1]
            for k in range(pnt[i], stp1):
                #if self.masses[k] > 2.0:
                alist_i.append(k)
            ni = len(alist_i)
            if ni == 0:
                ni == 1
            for j in range(0, nres):
                alist_j = []
                stp2 = 0
                if j == nres-1:
                    stp2 = natm
                else:
                    stp2 = pnt[j+1]
                for l in range(pnt[j], stp2):
                    #if self.masses[l] > 2.0:
                    alist_j.append(l)
                nj = len(alist_j)
                if ni == 0 or nj == 0:
                    ni = 1
                    nj = 1
                    
                if i == j:
                    self.affinityMatrix[i][j] = 0.0
                else:
                    self.affinityMatrix[i][j] = self.computeAffinity(alist_i, alist_j) / math.sqrt(ni * nj)
        return self.affinityMatrix

                
    def computeAffinity(self, alist_i, alist_j):
        nij = 0;
        for atom_i in alist_i:
            for atom_j in alist_j:
                tmp_x = self.x_Coords[atom_i] - self.x_Coords[atom_j]
                tmp_y = self.y_Coords[atom_i] - self.y_Coords[atom_j]
                tmp_z = self.z_Coords[atom_i] - self.z_Coords[atom_j]
                d = math.sqrt(tmp_x**2 + tmp_y**2 + tmp_z**2)
                if d <= 4.0:
                    nij = nij + 1;
        return float(nij)
                
if __name__== '__main__':
    filename = 'PDBs/1RX2.pdb'
    pdbread = PDBReader(filename)
    Aij = pdbread.computeAffinityMatrix()
    Aij = Aij[0:315, 0:315]
    from hitCommuteComputer import hitCommuteComputer
    hcc = hitCommuteComputer(Aij)
    [Hz, Cz] = hcc.driver(Aij)
    import pylab
    from pylab import *
    pylab.contourf(Cz[0:159, 0:159])
    pylab.colorbar()
    pylab.savefig(filename[0:-4] + '_comTimes.png')
    pylab.clf()
#    N = hcc.computeMeanOccupationTimes()
#    import pylab
#    from pylab import *
#    for i in range(0, len(N)):
#        pylab.contourf(N[i])
#        pylab.colorbar()
#        pylab.savefig('meanOcc_' + str(i+1) + '.png')
#        pylab.clf()
    print Hz.max()
    f = open(filename[0:-4] + '_hittingTimeStats.txt', 'w')
    HrMean = Hz.mean(1).tolist()
    HrStdv = Hz.std(1).tolist()
    HbMean = Hz.mean(0).tolist()
    for hitTime in HrMean:
        print >> f, hitTime
#    print >> f, HrStdv
    print Cz.max()
#    print Cz
        