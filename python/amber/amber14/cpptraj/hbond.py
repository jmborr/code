from pdb import set_trace as tr

types=('Acceptor', 'DonorH', 'Donor')

class Atom(object):
    '''Atom class'''
    def __init__(self,other=None,type=None):
        self.type=None
        self.resname=None
        self.resn=None
        self.atomname=None
        if type and isinstance(type,basestring) and type in types:
            self.type=type
        if other:
            if isinstance(other,basestring):
                '''assign by passing a string'''
                self.resname=other[0:3]
                self.resn=int( other.split('_')[1].split('@')[0] )
                self.atomname=other.split('_')[1].split('@')[1]

    def __repr__(self):
        return '{0}_{1}@{2}'.format(self.resname, self.resn, self.atomname)

    def __eq__(self, other):
        '''overload equality comparison'''
        return (isinstance(other, self.__class__) and self.__dict__==other.__dict__)
    
class Bond(object):
    '''Hydrogen bond class'''
    def __init__(self):
        self.acceptor=Atom(type='Acceptor')
        self.donorH=Atom(type='DonorH')
        self.donor=Atom(type='Donor')

    def __repr__(self):
        return '{0}-{1}-{2}'.format(self.acceptor, self.donor, self.donorH.atomname)

    def __eq__(self, other):
        '''overload equality comparison'''
        return (isinstance(other, self.__class__) and self.__dict__ == other.__dict__)

class AvgInfo(object):
    '''The average info stored in the output average file after running hbond'''
    def __init__(self, other=None):
        self.bond=Bond()
        self.frames=None
        self.frac=None
        self.avgDist=None
        self.avgAng=None
        if other:
            if isinstance(other,basestring):
                '''assign by passing a line from an "average file"'''
                items=other.split()
                self.bond.acceptor=Atom(items[0], type='Acceptor')
                self.bond.donorH=Atom(items[1], type='DonorH')
                self.bond.donor=Atom(items[2], type='Donor')
                self.frames=int(items[3])
                self.frac=float(items[4])
                self.avgDist=float(items[5])
                self.avgAng=float(items[6])
    
    def __repr__(self):
        acc=str(self.bond.acceptor)
        dh=str(self.bond.donorH)
        d=str(self.bond.donor)
        fm=self.frames
        fr=self.frac
        ad=self.avgDist
        aa=self.avgAng
        return '{0:13}{1:>14}{2:>14}{3:9d}       {4:6.4f}       {5:6.4f}     {6:8.4f}'.format(acc,dh,d,fm,fr,ad,aa)

class AvgInfoFile(object):
    '''list of AvgInfo objects plus header file'''
    def __init__(self,other=None):
        self.header=None
        self.records=[]
        if other:
            self.load(other)
    
    def load(self,other):
        local_other=other
        if isinstance(local_other,basestring):
            local_other = open(other,'r')
        if isinstance(local_other,file):
            for line in local_other.readlines():
                if line[0]=='#':
                    self.header=line.strip()
                else:
                    row=AvgInfo(line)
                    self.records.append(row)

    def __repr__(self):
        buf=self.header+'\n'
        buf += '\n'.join([str(row) for row in self.records])
        return buf


            
            
