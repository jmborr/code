import numpy



class atom:
    """store info for each atom"""

    def __init__(self,d={},l=''):
        if d:
            self.__dict__=d
        elif l:
            (i,t,j,m,x,y,z,s)=l.split()
            self.i=int(i); self.type=t; self.Number=int(j);
            self.mass=float(m)
            self.r=numpy.array([float(x),float(y),float(z)])
            self.Sigma=float(s)

    def print2seq(self):
        """print info in one line"""
        l=(self.i,self.type,self.Number,self.mass,
           self.r[0],self.r[0],self.r[2],
           self.Sigma)
        return '%3d %4s %3d %9.4f %9.4f %9.4f %4.1f'%l

class freq:
    """store info for each frequency"""

    def __init__(self,d={}):
        self.__dict__={'iMode':None, 'f':None, 'displ':numpy.array(None)}
        if d:
            self.__dict__=d    

    def print2seq(self,atoms):
        """print into into a sequence object"""
        buf='! Mode# Frequency\n'
        buf+='%5d %9.4f\n\n'%(self.iMode,self.f)
        buf+='! Displacements\n!                  X         Y         Z\n'
        for i in range(len(atoms)):
            at=atoms[i]; d=self.displ[i]
            l=(at.i,at.type,at.Number,d[0],d[1],d[2])
            buf+='%3d %4s %3d %9.4f %9.4f %9.4f\n'%l
        return buf
        
class aCLIMAX:
    """python object for aCLIMAX input file"""

    def __init__(self):
        self.atoms=[]
        self.freqs=[]

    def outputFreqFile(self,outf):
        """ output a frequency file suitable to be loaded in aCLIMAX"""

        buf='BEGIN ATOMS\n'
        for atom in self.atoms: buf+=atom.print2seq()+'\n'
        buf+='END ATOMS\n\n'
        
        buf+='BEGIN FREQUENCIES\n\n'
        for freq in self.freqs: buf+=freq.print2seq(self.atoms)+'\n'
        buf+='END FREQUENCIES\n'

        open(outf,'w').write(buf)
