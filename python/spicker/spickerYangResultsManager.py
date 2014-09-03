#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,junkName,Bye
from random import randint
from TMalign.findBestPDBhitByTMalign import findBestPDBhitByTMalign
from math import log #natural logarithm

difficulty_types={'easy':1,'medm':2,'hard':3}
pre1=re.compile('\s*(\d+)\s+(\d+)\s+(\d+.\d+)\s+(\d+.\d*)')
pre2=re.compile('\s*(\d+)\s+(\d+)\s+(\d+.\d+)\s+(\d+.\d+)')

def RMSD(f1,f2):
    return chomp(os.popen('rmsd '+f1+' '+f2+' | grep "RMSD="').readline()).split()[1]

def TMscore(f1,f2):
    #print f1,f2
    return os.popen('TMscore '+f1+' '+f2+' | grep "TM-score    ="').readline().split()[2]

def TMverbose(f1,f2,d0=''):
    return chomp(os.popen('tm_score.py -a '+f1+' -b '+f2+' -d y -e '+`d0`).readline())

#read a template from file pointer object "pt"
def readTemplate(pt):
    info={} #dictionary with all interesting quantities
    templ=[]
    line=pt.readline()
    if line:
        [info['libPDB'],L,Z,i]=line.split() #template-from-library
        info['L']=int(L)
        info['Z']=float(Z)
    else: return () #we reached end of file
    while(1):
        line=pt.readline() #;print line
        if not line: sys.exit(1) #break
        elif line.find('ATOM')>=0: templ.append(line)
        elif line.find('TER')>=0:
            seqid=line.split()[1]
            info['seqid']=float(seqid) #sequence identity to template-from-library
            break
    info['templ']=templ
    return info


def generateJunkFileName():
    junk='junk'+`randint(0,9999)`
    while os.path.exists(junk): junk='junk'+`randint(0,9999)`
    return junk

#calculate ZscoreOfCluster as average of all Z's with TM scores in the
#range [tmM,tmm] where tmM is bestTM and tmm=bestTM-0.1.  if tmm<0.45,
#then do the average of all Z's with TM scores in the [0.45,0.17]
#range, and use this average as another Z to be added to the average
#of all Z's with TM scores in the [tmM,0.45] range
def averageZs(zs):
    tmco=0.45 ; tmco2=0.17
     #sort zs by decreasing TM scores
    tms=zs.keys() ; tms.sort() ; tms.reverse() ;maxn=len(tms)-1
    tmM=tms[0] #;print 'tmM='+`tmM`

    #no relation between tasser model and templates
    if tmM<=tmco2: return zs[tmM]

    #very weak relation between tasser model and templates
    if tmM<tmco:
        sumZ=0 ; n=0
        while n<maxn:
            tm=tms[n]
            if tm<=tmco2: break
            sumZ+=zs[tm] ; n+=1
        return sumZ/n

    #significant relation between tasser model and some templates
    sumZ=0 ; n=0 ; sumZ2=0 ; m=0
    tmm=tmM-0.1 ; tm=tmM
    while tm>=tmco and tm>tmm:
        sumZ+=zs[tm] ; n+=1 #;print n
        if n>=maxn: break
        tm=tms[n]
    #do average of all templates in the (tmco,tmco2) window
    if tmm<tmco and n<maxn:
        m=n ; tm=tms[m]
        while tm>tmco2:
            sumZ2+=zs[tm] ; m+=1
            if m==maxn: break
            tm=tms[m]
    if m>n: return ( sumZ + sumZ2/(m-n) )/(n+1)
    return sumZ/n
    
class spickOut:    
    def __init__(self,dir='.',nat='',rebuilt=''):
        """Default arguments:
        dir: directory where files are located
        nat: pdb native file
        rebuilt: read combo*.pdb+rebuilt and closc*.pdb+rebuilt instead (rebuilt is a suffix)
        If no default arguments are passed, it's assumed the
        following files are located in dir:
        combo*.pdb closc*.pdb, rst.dat
        """

        #check files exist
        self.readError=''
        if not os.path.exists(dir+'/rst.dat'):
            self.readError='no rst.dat in '+dir
        if not os.path.exists(dir+'/combo1.pdb'+rebuilt):
            self.readError='no combo1.pdb'+rebuilt+' in '+dir
        if not os.path.exists(dir+'/closc1.pdb'+rebuilt):
            self.readError='no closc1.pdb'+rebuilt+' in '+dir
        if nat:
            if nat[0]!='/': #if not passed absolute path, assume nat in self.dir
                nat=dir+'/'+chomp(os.popen('basename '+nat).readline())
                if not os.path.exists(nat) : self.readError='no '+nat+' in '+dir
            else:
                if not os.path.exists(nat) : self.readError='no '+nat+' found'
        if self.readError: return None

        #intialize several variables
        self.nuse=0    #number of snapshots used in doing the spicker clustering
        self.dir=dir
        self.target='' #target name
        self.length='' #number of residues
        self.nat=None  
        self.ids=[]
        self.combo={}
        self.closc={}
        self.nc=0      #number of clusters (actually, number of combo files)
        self.size={}
        self.numberdens={} #number density (cluster-size/self.nuse)
        self.dens={}       #number_dessity/average_rmsd_to_centroid
        self.Rc_in={}  #average rmsd of cluster to centroid
        self.biggest=0 #cluster id of most populated cluster
        self.densest=0 #cluster id of most dense cluster
        self.densRanks={}
        self.rmsds={}
        self.zrrmsds={}
        self.tms={}
        self.Rg={}
        self.rebuilt='' #appendix name for rebuilt models
        self.nterCA=''
        self.cterCA=''
        self.templates=[]          #list storing templates from prospector
        self.Zscores=[]            #Zscore of each template from prospector
        self.templateOfCluster={}  #associated template of each cluster
        self.coverageSeqOfCombo={} #coverage sequence of each combo by its closest template
        self.ZscoreOfCluster={}    #associated Zscore of each cluster (it's an average of Z's)
        self.Cscores={}            #Cscore of each cluster
        self.Cest=0                #cluster id of cluster with higher C-score
        self.CscoresRanks={}       #rank id's by decreasing Cscore
        self.YingScoresRanks={}    #rank id's by decreasing Ying's score
        #appendix if using rebuilt models
        if rebuilt: self.rebuilt='.rebuilt'
        
        #read number of residues from rst.dat (remove third and fourth lines in the future)
        line=os.popen('grep -e "Modeling Length: " '+dir+'/rst.dat').readline();
        if line: self.length=int(re.compile(':\s+(\w+)').search(line).group(1))

        #read target name from rst.dat
        line=os.popen('grep -e "Target: " '+dir+'/rst.dat').readline();
        self.target=re.compile(':\s+(\w+)').search(line).group(1)

        #read number of output clusters
        line=os.popen('grep -e "Number of clusters: " '+dir+'/rst.dat').readline();
        self.nc=int(re.compile(':\s+(\w+)').search(line).group(1))

        #read number of used snapshots
        line=os.popen('grep -e "Number of structure in use=" '+dir+'/rst.dat').readline();
        self.nuse=int(re.compile('=\s+(\w+)').search(line).group(1))

        #read native range used, if available
        line=os.popen('grep -e "Native range used: " '+dir+'/rst.dat').readline();
        if line:
            [self.nterCA,self.cterCA]=re.compile(':\s+(\d+)\s-\s+(\d+)').search(line).groups()
            self.nterCA=int(self.nterCA)
            self.cterCA=int(self.cterCA)

        #read native file, if passed as argument
        if nat: self.nat=nat

        #fill several dictionaries
        p=open(dir+'/rst.dat','r')
        line=p.readline()
        while line and not re.compile(' B---').match(line): line=p.readline()
        if line:
            for i in range(0,self.nc):
                id,size,R_cut,dens=pre1.match( p.readline() ).groups()
                id=int(id)
                size=int(size)
                self.size[id]=size
                self.numberdens[id]=(1.0*size)/self.nuse            
        while line and not re.compile(' C---').match(line): line=p.readline()
        if line:
            for i in range(0,self.nc):
                id,size,R_in,Rc_in=pre2.match( p.readline() ).groups()
                id=int(id)
                Rc_in=float(Rc_in)
                if Rc_in<0.5: Rc_in=0.5 #fix for very small clusters
                self.Rc_in[id]=Rc_in
                self.dens[id]=self.numberdens[id]/Rc_in
        p.close()


        #populate ids, combo and closc arrays
        for i in range(1,self.nc+1):
            self.ids.append(i)
            self.combo[i]=dir+'/combo'+`i`+'.pdb'+rebuilt
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
        command='Rg.x -a '+self.combo[id]
        Rg=os.popen(command).readline().split()[0] #only first number
        return float(Rg)

    def initRg(self):
        for id in self.ids: self.Rg[id]=self.getRg(id=id)
        
    def TMtoNat(self,id='',str='combo'):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        if str=='combo':
            tm= TMscore(self.combo[id],self.nat)
        elif str=='closc':
            tm= TMscore(self.closc[id],self.nat)
        return float(tm)

    def TMidId(self,id1,id2,str='combo'):
        if str=='combo':
            tm= TMscore(self.combo[id1],self.combo[id2])
        else:
            tm= TMscore(self.closc[id1],self.closc[id2])
        return float(tm)

    def avTMid(self,id='',list=[]):
        avtm=0.0 ; n=0
        if not id: id=self.densest
        if not list: list=self.ids[:]
        for id2 in list:
            if id != id2:
                avtm+=self.TMidId(id,id2)
                n+=1
        return avtm/n

    def TMverboseIdId(self,id1,id2,d0=''):
        tmv= TMverbose(self.combo[id1],self.combo[id2],d0=d0)
        return tmv
        
    #PDB WITH superposition of centroid with passed id to native.
    def alignToNatByTM(self,outfile,id='',d0=0.0):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        os.system('tm_score.py -a '+self.combo[id]+' -b '+self.nat+' -c '+outfile+' -e '+`d0`)

    #print the first line of REMARKS section of the superposition by TM. Somethink like:
    #REMARK lseq=  114 rmsd=  5.59 tm= 0.530 naliged=  92 cov= 0.81 rmsd_cov=  3.72 lext= 106 aligdens= 0.86
    def outFirstLineFromAlignToNatByTM(self,id='',d0=0.0):
        if not self.nat:
            sys.stderr.write('ERROR: from TMtoNat: first use  importNat(nat=nativefile)\n')
            return
        if not id: id=self.densest
        #first line of standard output contain the precious info
        rk=chomp(os.popen('tm_score.py -a '+self.combo[id]+' -b '+self.nat+' -d yes -e '+`d0`).readline())
        return rk
        
    def initTMtoNat(self,str='combo'):
        if not self.nat:
            sys.stderr.write('ERROR: from initTMtoNat: first use  importNat(nat=nativefile)\n')
            return
        for id in self.ids: self.tms[id]=self.TMtoNat(id=id,str=str)

    def mutual_avTM(self,str='combo',list=[]):
        if not self.tms: self.initTMtoNat(str)
        if not list: list=self.ids[:]
        N=len(list)
        avTM=0.0
        for i in range(0,N-1):
            for j in range(i+1,N):
                avTM+=self.TMidId(list[i],list[j])
        return 2*avTM/(N*(N-1))
                    
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
        zrrmsd=chomp(os.popen('zrrmsd.x '+self.combo[id]+' '+self.nat).readline())
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
        return float(RMSD(self.combo[id],self.nat))

    def rmsdIdId(self,id1,id2):
        return float(RMSD(self.combo[id1],self.combo[id2]))
    
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

    #rank a list of id's by decreasing density
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

    #search against the filtered PDB database, and return the following line
    #modelHeader, pdbHeader, tmHit, Z_TM, gd, 1-gd, tmHit*(1-gd), (tmnat, tmHit-tmNat)
    def findBestPDBbyTMalign(self,id='',modelHeader='',pdbent='',pdblistf='',nat=''):
        if not id: id=self.densest
        return findBestPDBhitByTMalign(self.combo[id],modelHeader=modelHeader,pdbent=pdbent,pdblistf=pdblistf,nat=nat)

    def loadTemplates(self,templatefile):
        try: pt=open(templatefile,'r')
        except ValueError: print 'ERROR: template file',templatefile,'not found'
        n=0
        x=readTemplate(pt) #;print x; sys.exit(1)
        while x:
            self.Zscores.append(x['Z']) #;print self.Zscores[n]  #Z-score of the template
            self.templates.append(x) #all info on the template
            n+=1
            x=readTemplate(pt) #;print 'n=',n
        pt.close() #;print self.Zscores;print self.templates;sys.exit(1)
        return n #number of templates read

    def calcZscoreOfCluster(self,templatefile,id=''):

        """For every centroid, obtain PROSPECTOR templates that are
        structurally similar, then average their Zscores and the
        resulting average is the Zscore of the cluster associated to
        the centroid"""
        
        if not id: id=self.densest
        if not self.templates: self.loadTemplates(templatefile)
        zs={} #store pairs of (TM,Z)
        junk=generateJunkFileName() #will temporarily store each of the templates
        besttm=0.0 #;print id
        for i in range(0,len(self.templates)):
            open(junk,'w').writelines(self.templates[i]['templ']) #dump PDB lines to file
            tm=float(TMscore(junk,self.combo[id]))
            zs[tm]=self.Zscores[i]
            if tm>besttm:
                besttm=tm
                self.templateOfCluster[id]=self.templates[i]
        self.ZscoreOfCluster[id]=averageZs(zs)
        os.system('/bin/rm '+junk) #don't forget to clean up your mess!
        return self.ZscoreOfCluster[id]

    def initZscoreOfCluster(self,templatef):
        for id in self.combo.keys(): self.calcZscoreOfCluster(templatef,id=id)


    def calcCoverageSeqOfCombo(self):

        """for each combo file, create a sequence of 1's and 0's. '1'
        if the residue is covered by the PROSPECTOR template most
        similar to the combo, '0' otherwise"""

        if 'coverageSeqOfCombo' not in self.__dict__.keys(): self.coverageSeqOfCombo={}
            
        for id in self.combo.keys():
            template=self.templateOfCluster[id]['templ'] #template file as list
            seq='0'*self.length #init sequence as not covered by the template
            for xyz in template:
                if ' CA ' in xyz: #check a CA atom
                    i=int(xyz.split()[4])-1 #residue position, notice the shift
                    #mark position as covered
                    if i==self.length-1: #last index
                        seq=seq[:i]+'1'
                    else:
                        seq=seq[0:i]+'1'+seq[i+1:]
            self.coverageSeqOfCombo[id]=seq
            
    def initCscores(self,templatef):
        if not self.ZscoreOfCluster: self.initZscoreOfCluster(templatef)
        for (id,Zscore) in self.ZscoreOfCluster.items():
            self.Cscores[id]=log(self.dens[id]*self.ZscoreOfCluster[id])


    def rankIDsByCscores(self,templatef=''):

        """rank a list of id's by decreasing Cscore"""
        
        if not self.Cscores: self.initCscores(templatef)
        items=[(score, id) for id, score in self.Cscores.items()] #list of (value,key) pairs
        #sort by first elem. of (value,key),that is,by value. reverse() makes decreasing score
        items.sort()
        items.reverse()
        self.Cest=items[0][1]  #ID of cluster with highest C-score
        rank=0
        for (score,id) in items:
            rank+=1
            self.CscoresRanks[id]=rank
        return [id for score, id in items] #report key, which is cluster id


    def rankIDsByYingPredictedTM(self,templatef='',type='easy'):
        """Est_TMscore = 0.1060727 * fraction - 0.0480012* Average
        -0.0635406*type + 0.0085284 * Z.score + 0.0024590 * len +
        0.4592112.
        NOTE: the term 0.0024590 * len is not usefult towards ranking
        because all of the clusters have same length
        """
        diff=difficulty_types[type]
        if not self.ZscoreOfCluster: self.initZscoreOfCluster(templatef)
        #print 'Zscore=',self.ZscoreOfCluster
        #predict the TM-score of each model
        tms={}
        for (id,Zscore) in self.ZscoreOfCluster.items():
            #print id,Zscore,self.numberdens[id],self.Rc_in[id]
            tms[id]=0.1060727*self.numberdens[id] - 0.0480012*self.Rc_in[id] - 0.0635406*diff + 0.0085284*Zscore + 0.0024590*self.length + 0.4592112
        #sort scores from bigger to smaller
        items=[(tm, id) for id, tm in tms.items()] #list of (value,key) pairs
        items.sort()
        items.reverse() #;print 'tms=',items
        rank=0
        #init YingScoresRanks dictionary
        for (score,id) in items:
            rank+=1
            self.YingScoresRanks[id]=rank
        return [id for tm, id in items] #report key, which is cluster id
        
    
    def rankIDsByYingPredictedTMforAbInitio(self,templatef=''):
        """-0.67766 + 0.65850*TMotherModels - 0.05218*clusterRMSD +
        0.15072* fraction + 0.23459*log(len )"""
        pass
           
if __name__=='__main__':
    ih=inpHand('Usage: resultsManager.py [options]\nwill output TMscore of densest cluster to native',
               ' -a _RA_dir directory where spicker results are stored',
               ' -b _RA_nat native pdb file')
    ih.parse(locals(),sys.argv)

    p=spickOut(dir=dir)
