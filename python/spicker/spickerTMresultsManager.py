#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp
from random import randint

pre1=re.compile('\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)')

def TMalign(f1,f2):
    return chomp(os.popen('TMalign '+f1+' '+f2+' |grep "TM-score=" |cut -d\'=\' -f4 |cut -d\',\' -f1').readline())

def TMscore(f1,f2):
    return os.popen('TMscore '+f1+' '+f2+' | grep "TM-score    ="').readline().split()[2]

def TMverbose(f1,f2,d0=''):
    return chomp(os.popen('tm_score.py -a '+f1+' -b '+f2+' -d y -e '+`d0`).readline())

class spickOut:
    
    def __init__(self,dir='.',nat='',rebuilt=''):
        """Default arguments:
        dir: directory where files are located
        nat: pdb native file
        carange: subset of residues in the native file used to compare to generated models.
                 For instance:  34-99 means count only residues 34 to 99 (both included)
        rebuilt: read combc*.pdb.rebuilt and closc*.pdb rebuilt instead
        If no default arguments are passed, it's assumed the
        following files are located in dir:
        combc*.pdb closc*.pdb, rst.dat
        """

        #check files exist
        self.readError=''
        if not os.path.exists(dir+'/rst.dat'):
            self.readError='no rst.dat in '+dir
        if not rebuilt:
            if not os.path.exists(dir+'/combc1.pdb'): self.readError='no combc1.pdb in '+dir
            if not os.path.exists(dir+'/closc1.pdb'): self.readError='no closc1.pdb in '+dir
        else:
            if os.path.exists(dir+'/combc1.pdb.rebuilt'):
                self.readError='no combc1.pdb.rebuilt in '+dir
            if os.path.exists(dir+'/closc1.pdb.rebuilt') :
                self.readError='no closc1.pdb.rebuilt in '+dir
        if nat:
            if nat[0]!='/': #if not passed absolute path, assume nat in self.dir
                nat=dir+'/'+chomp(os.popen('basename '+nat).readline())
            if not os.path.exists(nat) : self.readError='no '+nat+' in '+dir
        if self.readError: return None
        if self.readError: return None

        #intialize several variables
        self.dir=dir
        self.target='' #target name
        self.length='' #number of residues
        self.nat=None  
        self.ids=[]
        self.combc={}
        self.closc={}
        self.nc=0      #number of clusters (actually, number of combc files)
        self.size={}
        self.dens={}
        self.biggest=0 #cluster id of most populated cluster
        self.densest=0 #cluster id of most dense cluster
        self.densRanks={}
        self.rmsds={}
        self.zrrmsds={}
        self.tms={}
        self.tmAls={}   #store TMalign, NOT TMscore
        self.Rg={}
        self.rebuilt='' #appendix name for rebuilt models
        self.nterCA=''
        self.cterCA=''
        
        #appendix if using rebuilt models
        if rebuilt: self.rebuilt='.rebuilt'
        
        #read number of residues from rst.dat (remove third and fourth lines in the future)
        line=os.popen('grep -e "Model length: " '+dir+'/rst.dat').readline();
        if line: self.length=int(re.compile(':\s+(\w+)').search(line).group(1))            
        line=os.popen('grep -e "Modeling Length: " '+dir+'/rst.dat').readline();
        if line: self.length=int(re.compile(':\s+(\w+)').search(line).group(1))

        #read target name from rst.dat
        line=os.popen('grep -e "Target name: " '+dir+'/rst.dat').readline();
        self.target=re.compile(':\s+(\w+)').search(line).group(1)

        #read number of output clusters
        line=os.popen('grep -e "Number of clusters: " '+dir+'/rst.dat').readline();
        self.nc=int(re.compile(':\s+(\w+)').search(line).group(1))

        #read native range used, if available
        line=os.popen('grep -e "Native range used: " '+dir+'/rst.dat').readline();
        if line:
            [self.nterCA,self.cterCA]=re.compile(':\s+(\d+)\s-\s+(\d+)').search(line).groups()
            self.nterCA=int(self.nterCA)
            self.cterCA=int(self.cterCA)

        #read native file
        if nat:
            self.importNat(nat)

        #fill densities array
        p=open(dir+'/rst.dat','r')
        line=p.readline()
        while line and not re.compile(' B---').match(line): line=p.readline()
        if line:
            for i in range(0,self.nc):
                id,percent,dens,avtm=pre1.match( p.readline() ).groups()
                id=int(id)
                self.size[id]=int(percent)
                self.dens[id]=float(dens)  
        p.close()

        #populate ids, combc and closc arrays
        for i in range(1,self.nc+1):
            self.ids.append(i)
            self.combc[i]=dir+'/combc'+`i`+'.pdb'+rebuilt
            self.closc[i]=dir+'/closc'+`i`+'.pdb'+rebuilt        

        #find biggest and densest centroids, out of all output clusters
        self.biggest=self.ids[0] #assumed first id corresponds to most populated cluste
        self.densest=self.ids[0]
        for id in self.ids:
            if self.dens[id] > self.dens[self.densest] : self.densest=id

        #read RMSD comparison to native from rst.dat, if available
        p=open(dir+'/rst.dat','r')
        line=p.readline()
        while line and not re.compile(' A---').match(line):  line=p.readline()
        if line:
            for i in range(0,self.nc):
                id=self.ids[i]
                self.rmsds[id]=float(p.readline().split()[1]) #second item
        p.close()
        
        #read TM comparison to native from rst.dat, if available
        p=open(dir+'/rst.dat','r')
        line=p.readline()
        while line and not re.compile(' AB---').match(line):  line=p.readline()
        if line:
            for i in range(0,self.nc):
                id=self.ids[i]
                self.tms[id]=float(p.readline().split()[1]) #second item
        p.close()
        
    #get id of densest cluster on the 'top' biggest clusters
    def getDensestID(self,top=''):
        if not top: return self.densest
        idd=self.ids[0]
        for id in self.ids:
            if id<=top: #required id of cluster above (or equal) to top
                if self.dens[id] > self.dens[idd] : idd=id
        return idd

    def getRg(self,id=''):
        if not id: id=self.densest
        command='Rg.x -a '+self.combc[id]
        Rg=os.popen(command).readline().split()[0] #only first number
        return float(Rg)

    def initRg(self):
        for id in self.ids: self.Rg[id]=self.getRg(id=id)
        
    def TMtoNat(self,id=''):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        tm= TMscore(self.combc[id],self.nat)
        return float(tm)

    def TMidId(self,id1,id2):
        tm= TMscore(self.combc[id1],self.combc[id2])
        return float(tm)

    def TMverboseIdId(self,id1,id2,d0=''):
        tmv= TMverbose(self.combc[id1],self.combc[id2],d0=d0)
        return tmv
        
    #PDB WITH superposition of centroid with passed id to native.
    def alignToNatByTM(self,outfile,id='',d0=0.0):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        os.system('tm_score.py -a '+self.combc[id]+' -b '+self.nat+' -c '+outfile+' -e '+`d0`)

    #print the first line of REMARKS section of the superposition by TM. Somethink like:
    #REMARK lseq=  114 rmsd=  5.59 tm= 0.530 naliged=  92 cov= 0.81 rmsd_cov=  3.72 lext= 106 aligdens= 0.86
    def outFirstLineFromAlignToNatByTM(self,id='',d0=0.0):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        #first line of standard output contain the precious info
        rk=chomp(os.popen('tm_score.py -a '+self.combc[id]+' -b '+self.nat+' -d yes -e '+`d0`).readline())
        return rk
        
    def initTMtoNat(self):
        if not self.nat:
            sys.stderr.write('ERROR: from initTMtoNat: first use  importNat(nat=nativefile)\n')
            return
        for id in self.ids: self.tms[id]=self.TMtoNat(id=id)

    #rank a list of id's by decreasing TM. Optional list parameter is a list of id's
    def rankIDsByTMtoNat(self,list=[]):
        dict={}
        if not list: list=self.ids[:]
        if self.tms: #use self.tms if the dictionary is filled
            for id in list: dict[id]=self.tms[id]
        else:
            for id in list: dict[id]=self.TMtoNat(id)
        items=[(v, k) for k, v in dict.items()]
        items.sort()
        items.reverse()
        return [k for v, k in items]

    def ZrRMSDtoNat(self,id=''):
        if not self.nat:
            sys.stderr.write('ERROR: from ZrRMSDtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        zrrmsd=chomp(os.popen('zrrmsd.x '+self.combc[id]+' '+self.nat).readline())
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
        rmsd=chomp(os.popen('rmsd.x '+self.combc[id]+' '+self.nat).readline())
        return float(rmsd)

    def rmsdIdId(self,id1,id2):
        rmsd=chomp(os.popen('rmsd.x '+self.combc[id1]+' '+self.combc[id2]).readline())
        return float(rmsd)
    
    def initRmsdsToNat(self):
        if not self.nat:
            sys.stderr.write('ERROR: from initRmsdsToNat: first use  importNat(nat=nativefile)\n')
            return
        for id in self.ids: self.rmsds[id]=self.rmsdToNat(id=id)

    #rank a list of id's by increasing rmsd
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
        else: self.nat=nat


    def tmAltoStr(self,structure,id=''):
        if not id: id=self.densest
        return float(TMalign(structure,self.combc[id]))  #d0 determined by combc
    
    def filltmAls(self,structure):
        for id in self.ids: self.tmAls[id]=self.tmAltoStr(structure,id)
  
    #rank a list of id's by decreasing TM. Optional list parameter is a list of id's
    def rankIDsBytmAlstoStrt(self,list=[],str=''):
        dict={}
        if not list: list=self.ids[:]
        if self.tmAls: #use self.tms if the dictionary is filled
            for id in list: dict[id]=self.tmAls[id]
        else:
            if not str: print 'ERROR: structure needed by rankIDsBytmAlstoStrt'
            for id in list: dict[id]=self.tmAltoStr(str,id)
        items=[(v, k) for k, v in dict.items()]
        items.sort()
        items.reverse()
        return [k for v, k in items]

if __name__=='__main__':
    ih=inpHand('Usage: resultsManager.py [options]\nwill output TMscore of densest cluster to native',
               ' -a _RA_dir directory where spicker results are stored',
               ' -b _RA_nat native pdb file')
    ih.parse(locals(),sys.argv)

    p=spickOut(dir=dir)
