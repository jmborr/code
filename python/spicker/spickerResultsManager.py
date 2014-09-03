#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp
from random import randint

pre0=re.compile('rank size percent density')
pre1=re.compile('\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)')
pre2=re.compile('Target:\s*(\w+)')

def TMscore(f1,f2):
    return os.popen('TMscore '+f1+' '+f2+' | grep "TM-score    ="').readline().split()[2]

def TMverbose(f1,f2,d0=''):
    return chomp(os.popen('tm_score.py -a '+f1+' -b '+f2+' -d y -e '+`d0`).readline())

def TMonlySecStr(f1,f2,seqd):
    return chomp(os.popen('tm_score.py -a '+f1+' -b '+f2+' -f '+seqd).readline())

class spickOut:
    
    def __init__(self,dir='.',target='',nat=''):        
        """If no default arguments are passed, it's assumed the
        following files are located in dir:
        rmsinp, seq.dat, combo*pdb closc*pdb
        """

        #check all files exist
        self.readError=''
        if not os.path.exists(dir+'/centroids.pdb'): self.readError='no centroids.pdb in '+dir
        if not os.path.exists(dir+'/closests.pdb') : self.readError='no closests.pdb in '+dir
        if not os.path.exists(dir+'/summary.txt')  : self.readError='no summary.txt in '+dir
        if nat and not os.path.exists(dir+'/'+nat) : self.readError='no '+nat+' in '+dir
        if self.readError: return None

        self.dir=dir
        self.target=''
        self.nat=None  
        self.ids=[]
        self.combo={}
        self.closc={}
        self.nc=0      #number of clusters (actually, number of combo files)
        self.size={}
        self.dens={}
        self.biggest=0 #cluster id of most populated cluster
        self.densest=0 #cluster id of most dense cluster
        self.densRanks={}
        self.rmsds={}
        self.zrrmsds={}
        self.tms={}
        self.Rg={}
        
        #get target name
        if not target:
            if os.path.exists(dir+'/rmsinp'):
                self.target=chomp(os.popen('tail -1 '+dir+'/rmsinp').readline())
            else:
                p=pre2.match(open(dir+'/summary.txt','r').readline())
                if p: self.target=p.group(1)
                else: self.target=''
        else: self.target=target

        #read centroids.pdb, populate also self.ids
        self.nc=1
        self.ids.append(self.nc)
        p=open(dir+'/centroids.pdb','r')
        self.combo[self.nc]=p.tell() #stores position where centroid begins
        line=p.readline()
        while line:
            if line[0:3]=='TER':
                pos=p.tell()                
                line=p.readline()
                if line:
                    self.nc+=1
                    self.ids.append(self.nc)
                    self.combo[self.nc]=pos #stores position where centroid begins
            line=p.readline()
        p.close()

        #read closc files
        nc=1
        p=open(dir+'/closests.pdb','r')
        self.closc[nc]=p.tell() #stores position where centroid begins
        nc+=1
        line=p.readline()
        while line:
            if line[0:3]=='TER':
                pos=p.tell()                
                line=p.readline()
                if line:
                    nc+=1
                    self.closc[nc]=pos #stores position where centroid begins
            line=p.readline()
        p.close()

        #read native file
        if nat:
            self.importNat(nat)

        #read summary.txt file
        p=open(dir+'/summary.txt','r')
        while not pre0.match(p.readline()): pass
        for i in range(0,self.nc):
            id,size,percent,dens=pre1.match( p.readline() ).groups()
            id=int(id) ; self.size[id]=int(size) ; self.dens[id]=float(dens)
        p.close()

        
        #biggest and densest centroids
        self.biggest=self.ids[0] #assumed first id corresponds to most populated cluster
        self.densest=self.ids[0]
        for id in self.ids:
            if self.dens[id] > self.dens[self.densest] : self.densest=id

    #get id of densest cluster on the 'top' best clusters
    def getDensestID(self,top=''):
        if not top: return self.densest
        idd=self.ids[0]
        for id in self.ids:
            if id<=top: #required id of cluster above (or equal) to top
                if self.dens[id] > self.dens[idd] : idd=id
        return idd
    
    #output a centroid. If passed range, only a segment, not all
    #coordinated are passed. Range must be a two-number list    
    def printCentroid(self,file,id='',segment='', set='centroids.pdb'):
        if not id: id=self.densest
        pf=open(file,'w')
        p=open(self.dir+'/'+set,'r')
        p.seek(self.combo[id])
        nread=0
        if segment:
            for i in range(1,segment[0]):
                p.readline() #skip N-terminal residues outside the range
                nread+=1
        while 1:
            line=p.readline()
            nread+=1
            if segment:
                if nread > segment[1]:
                    pf.write('TER\n')
                    break #skip C-terminal residues outside the range
            pf.write(line)
            if line[0:3]=='TER': break        
        p.close()
        pf.close()

    def getRg(self,id=''):
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        command='Rg.x '+file
        Rg=chomp(os.popen(command).readline())
        os.system('/bin/rm '+file)
        return float(Rg)

    def initRg(self):
        for id in self.ids: self.Rg[id]=self.getRg(id=id)
        
    def TMtoNat(self,id=''):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        tm= TMscore(file,self.nat)
        os.system('/bin/rm '+file)
        return float(tm)

    def TMidId(self,id1,id2):
        file1='junk'+`randint(0,99999)`  ;  self.printCentroid(file1,id=id1)
        file2='junk'+`randint(0,99999)`  ;  self.printCentroid(file2,id=id2)
        tm= TMscore(file1,file2)
        os.system('/bin/rm '+file1+' '+file2)
        return float(tm)

    def TMverboseIdId(self,id1,id2,d0=''):
        file1='junk'+`randint(0,99999)`  ;  self.printCentroid(file1,id=id1)
        file2='junk'+`randint(0,99999)`  ;  self.printCentroid(file2,id=id2)
        tmv= TMverbose(file1,file2,d0=d0)
        os.system('/bin/rm '+file1+' '+file2)
        return tmv
        
    #PDB WITH superposition of centroid with passed id to native.
    def alignToNatByTM(self,outfile,id='',d0=0.0):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        os.system('tm_score.py -a '+file+' -b '+self.nat+' -c '+outfile+' -e '+`d0`)
        os.system('/bin/rm '+file)

    #print the first line of REMARKS section of the superposition by TM. Somethink like:
    #REMARK lseq=  114 rmsd=  5.59 tm= 0.530 naliged=  92 cov= 0.81 rmsd_cov=  3.72 lext= 106 aligdens= 0.86
    def outFirstLineFromAlignToNatByTM(self,id='',d0=0.0):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        #first line of standard output contain the precious info
        rk=chomp(os.popen('tm_score.py -a '+file+' -b '+self.nat+' -d yes -e '+`d0`).readline())
        os.system('/bin/rm '+file)
        return rk
        
    def initTMtoNat(self):
        if not self.nat:
            sys.stderr.write('ERROR: from initTMtoNat: first use  importNat(nat=nativefile)\n')
            return
        for id in self.ids: self.tms[id]=self.TMtoNat(id=id)

    #rank a list of id's by decreasing TM
    def rankIDsByTMtoNat(self,list=[]):
        dict={}
        if not list: list=self.ids[:]
        if self.tms:
            for id in list: dict[id]=self.tms[id]
        else:
            for id in list: dict[id]=self.TMtoNat(id)            
        items=[(v, k) for k, v in dict.items()]
        items.sort()
        items.reverse()
        return [k for v, k in items]

    def TMonlySecStrToNat(self,seqd,id=''):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        tm= TMonlySecStr(file,self.nat,seqd)
        os.system('/bin/rm '+file)
        return float(tm)

    def rankIDsByTMonlySecStrToNat(self,seqd,list=[]):
        dict={}
        if not list: list=self.ids[:]
        for id in list: dict[id]=self.TMonlySecStrToNat(seqd,id)
        items=[(v, k) for k, v in dict.items()]
        items.sort()
        items.reverse()
        return [k for v, k in items]
        
    def ZrRMSDtoNat(self,id=''):
        if not self.nat:
            sys.stderr.write('ERROR: from ZrRMSDtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        zrrmsd=chomp(os.popen('zrrmsd.x '+file+' '+self.nat).readline())
        os.system('/bin/rm '+file)
        return float(zrrmsd)

    def initZrRMSDtoNat(self):
        if not self.nat:
            sys.stderr.write('ERROR: from initZrRMSDtoNat: first use  importNat(nat=nativefile)\n')
            return
        for id in self.ids: self.zrrmsds[id]=self.ZrRMSDtoNat(id=id)

    #rank a list of id's by incresing Z-rRMSD
    def rankIDsByZrRMSDtoNat(self,list=[]):
        dict={}
        if not list: list=self.ids[:]
        if self.zrrmsds:
            for id in list: dict[id]=self.zrrmsds[id]
        else:
            for id in list: dict[id]=self.zrrmsds[id]
        items=[(v, k) for k, v in dict.items()]
        items.sort()
        return [k for v, k in items]

    def rmsdToNat(self,id=''):
        if not self.nat:
            sys.stderr.write('ERROR: from rmsdToNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        file='junk'+`randint(0,99999)`
        self.printCentroid(file,id=id) #print centroid to temporary file
        rmsd=chomp(os.popen('rmsd.x '+file+' '+self.nat).readline())
        os.system('/bin/rm '+file)
        return float(rmsd)

    def rmsdIdId(self,id1,id2):
        file1='junk'+`randint(0,99999)`  ;  self.printCentroid(file1,id=id1)
        file2='junk'+`randint(0,99999)`  ;  self.printCentroid(file2,id=id2)
        rmsd=chomp(os.popen('rmsd.x '+file1+' '+file2).readline())
        os.system('/bin/rm '+file1+' '+file2)
        return float(rmsd)
    
    def initRmsdsToNat(self):
        if not self.nat:
            sys.stderr.write('ERROR: from initRmsdsToNat: first use  importNat(nat=nativefile)\n')
            return
        for id in self.ids: self.rmsds[id]=self.rmsdToNat(id=id)

    #rank a list of id's by incresing rmsd
    def rankIDsByRmsdtoNat(self,list=[]):
        dict={}
        if not list: list=self.ids[:]
        if self.rmsds:
            for id in list: dict[id]=self.rmsds[id]
        else:
            for id in list: dict[id]=self.rmsdToNat(id)
        items=[(v, k) for k, v in dict.items()]
        items.sort()
        return [k for v, k in items]

    def rankIDsByDens(self):
        items=[(dens, id) for id, dens in self.dens.items()] #list of (value,key) pairs
        #sort by first element of (value,key), that is, by value. reverse() makes decreasing dens.
        items.sort()
        items.reverse()
        rank=0
        for (dens,id) in items:
            rank+=1
            self.densRanks[id]=rank
        return [id for dens, id in items] #report key, which is cluster id

        
    def importNat(self,nat='./CA'):
        if nat[0]!='/': #if not passed absolute path, assume nat in self.dir
            self.nat=self.dir+'/'+chomp(os.popen('basename '+nat).readline())
        
if __name__=='__main__':
    ih=inpHand('Usage: resultsManager.py [options]\nwill output TMscore of densest cluster to native',
               ' -a _RA_dir directory where spicker results are stored',
               ' -b _RA_nat native pdb file')
    ih.parse(locals(),sys.argv)

    p=spickOut(dir=dir)
    p.importNat(nat)
    print p.TMtoNat()
