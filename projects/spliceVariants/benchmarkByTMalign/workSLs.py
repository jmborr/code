#!/usr/bin/python

import os,sys,re,pdb
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand

"""Manage the benchmarkByTMalign results"""


#global variables
joink=os.path.join
isthere=os.path.exists
pcd=joink(codedir,'projects','spliceVariants','benchmarkByTMalign')
scd=joink(scratchdir,'spliceVariants','benchmarkByTMalign','out')
pdbdatd=joink(scratchdir,'adrian/dat/pdb_jul07')
collf=joink(scd,'varsets.can')

#files for targets (prefix "t")
tf=joink(scd,'../filtered.list')
trepod=joink(scd,'trepod')
tcollf=joink(trepod,'collection.can')
#files for pairs of short_long targets (prefix "sl")
slf=joink(scd,'../short_long.list')
slrepod=joink(scd,'slrepod')
slcollf=joink(slrepod,'collection.can')


#################################################################

class classT(CAN):

    
    pass

#################################################################

class classTl(CANl):

    """inherited class from CANl for targets, so that we can add
    methods to fuzz around the attributes of our particular CANl
    object implementation"""

    pass

objTl=classTl(tcollf,tf,trepod) #holds info for targets

#################################################################

class classSL(CAN):


    def addJobname(self,jobname0):

        """add an identifier for simulations"""

        for typesl in ('short','long'):
            self.__dict__[typesl]['jobname']=jobname0+typesl[0]
        return True
    
            
    def addSimOut(self,target=''):

        """add tasser, spk dictionaries"""
        
        if target=='non interacting loop':
            for typesl in ('short','long'):
                pt0=self.__dict__[typesl]
                if 'nil' not in pt0.keys():
                    simdir=joink(pt0['simdir'],'nil')
                    pt0['nil']={'simdir':simdir,}
                pt=pt0['nil']
                simdir=pt['simdir']
                pt['tasser']={'outfile':''}
                outfile=joink(simdir,pt0['id']+'.out.tasser.tar.bz2')
                if isthere(outfile): pt['tasser']['outfile']=outfile
                pt['spk']={'outfile':''}
                outfile=joink(simdir,pt0['id']+'.out.spicker.tar.bz2')
                if isthere(outfile): pt['spk']['outfile']=outfile
        return True
        
        
    def AddSpkobj(self,target=''):

        """add spicker object attribute for each variant"""
        
        from utilities.small_utilities import unTARme
        from spicker.spickerYangResultsManager import spickOut
        for typesl in ('short','long'):
            pt1=self.__dict__[typesl]
            if target=='non interacting loop': pt1=pt1['nil']
            pt2=pt1['spk']
            pt2['spkobj']=None #assume no spicker results
            if not pt2['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            wd=unTARme(pt2['outfile']) #untar output file to simdir directory
            pt2['spkobj']=spickOut(dir=wd)
            os.system('/bin/rm -r '+wd)
        return False


    def addLoopdat(self,target=''):

        """add loop properties from loop.dat file"""

        from utilities.small_utilities import unTARme
        
        for type in ('short','long'):

            ct=self.__dict__[type]
#            if target=='non interacting loop':ct=ct['nil']            
            if not ct['spk']['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+type+'\n')
                continue
            wd=unTARme(ct['spk']['outfile']) #untar to some working directory
            if 'loopdat' not in ct['spk'].keys(): ct['spk']['loopdat']={}
            ct=ct['spk']['loopdat']
            pt=open(joink(wd,'loop.dat'),'r')
            line=''
            while '#lid begin end' not in line: line=pt.readline()
            line=pt.readline()
            while '#comparison to native' not in line:
                items=line.split()
                if len(items)<5: #loop is short, thus we did not compute info
                    ct[items[0]]={'id':items[0],'begin':int(items[1]),'end':int(items[2]),
                                  'L':int(items[3]),'istooshort':True }
                else:
                    ct[items[0]]={'id':items[0],'begin':int(items[1]), 'end':int(items[2]),
                                  'L':int(items[3]),'exp':items[4],'h':items[5],'e':items[6],
                                  'c':items[7],'rl':items[8],'lst':items[9],'rg':items[10],
                                  'rg2':items[11],'istooshort':False }
                line=pt.readline()
            while '#' in line: line=pt.readline()
            while line:
                items=line.split()
                if len(items)>4:ct[items[0]]['tonat']={'rl':items[4],'rg':items[5],'rg2':items[6]}
                line=pt.readline()
            os.system('/bin/rm -r '+wd)
        return


    def submitSpklg(self,joblist='',switch='',ibatch=0):

        """Submit spicker job only"""
        
        from jobs.job import job,genJobList

        wt=0.249
        if switch=='normal run':
            for typesl in ('short','long'):
                p0=self.__dict__[typesl]
                outdir=p0['simdir']
                jobname=p0['jobname']+'.spklg'
                inputs=[joink(p0['simdir'],'chain.dat'),p0['CA']]  #files to copy
                for f in ('exp.dat','rmsinp','seq.dat'): inputs.append(joink(scd,p0['id'],f))
                incmd='tar jxf '+p0['tasser']['outfile']+' -C . && '+\
                       '/bin/mv '+os.path.basename(p0['CA'])+' CA ' #input commands
                outf=p0['id']+'.out.spicker.tar.bz2'
                outcmd='tar jcf '+outf+' closc* combo* loop* rst* && '+\
                        '/bin/mv '+outf+' '+outdir #output command
                Job=job(exe='spk_rp',exed=pcd,inputs=inputs,incmd=incmd,outcmd=outcmd)
                Job.qsub(jobname,outdir,wallt=wt,mem_limit='890',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        elif switch=='nil':
            for typesl in ('short','long'):
                p0=self.__dict__[typesl]
                p1=p0['nil']
                outdir=p1['simdir']
                jobname=p0['jobname']+'.nil'
                tassout=p0['id']+'.out.tasser.tar.bz2' #basename of tasser output
                inputs=[joink(p0['simdir'],'chain.dat'),p0['CA'],\
                        joink(p1['simdir'],tassout)]  #files to copy
                for f in ('exp.dat','rmsinp','seq.dat'): inputs.append(joink(scd,p0['id'],f))
                incmd='tar jxf '+tassout+' && '+\
                       '/bin/mv '+os.path.basename(p0['CA'])+' CA ' #input commands
                outf=p0['id']+'.out.spicker.tar.bz2'
                outcmd='tar jcf '+outf+' closc* combo* loop* rst* && '+\
                        '/bin/mv '+outf+' '+outdir #output command
                Job=job(exe='spk_rp',exed=pcd,inputs=inputs,incmd=incmd,outcmd=outcmd)
                Job.qsub(jobname,outdir,wallt=wt,mem_limit='890',joblist=joblist)
                ibatch=Job.dormant(ibatch)
            
        return ibatch,False

    
    def addSpklg(self):

        """add spicker clustering of loop conformations using local or global rmsd"""

        from utilities.small_utilities import unTARme
        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]['spk'] #; Bye(ct.keys())
            if not  ct['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            ct['spklocal']={} ; ct['spkglobal']={}
            dstg=ct['spkobj'].densest #ID of densest global cluster
            wd=unTARme(ct['outfile'])
            for n in ct['loopdat'].keys(): #cycle through all loops
                if ct['loopdat'][n]['istooshort']: continue #loop is too short
                for rmsdtype in ('local','global'):
                    ct['spk'+rmsdtype][n]={}
                    st=ct['spk'+rmsdtype][n] #reference
                    rstf='rst_%s_%02d.dat'%(rmsdtype,int(n)) #;Bye(rstf)
                    pt=open(joink(wd,rstf),'r') #; print joink(wd,rstf)
                    line=pt.readline()
                    while 'Native Length:' not in line: line=pt.readline()
                    st['L']=int(line.split()[-1]) #loop size
                    st['nc']=int(pt.readline().split()[-1]) #number of clusters
                    st['availsnap']=int(pt.readline().split('=')[-1]) #total #avail. structures
                    st['usedsnap']=int(pt.readline().split('=')[-1]) #used snapshots
                    #init combo2nat and closc2nat
                    st['combo2nat']={} ; st['closc2nat']={}
                    while 'A--------------' not in line: line=pt. readline()               
                    while True: #cycle through comparisons to native
                        line=pt.readline()
                        items=line.split() #; print 'items=',items
                        if not items: break
                        st['combo2nat'][items[0]]=items[1]
                        st['closc2nat'][items[0]]=items[4]
                    #init size, R_cut, density
                    st['size']={} ; st['R_cut']={} ; st['density']={}
                    while 'B--------------' not in line: line=pt.readline()
                    while True:
                        line=pt.readline()
                        items=line.split() #; print 'items=',items
                        if not items: break
                        if 'include used structure' in line: break
                        st['size'][items[0]]=int(items[1])
                        st['R_cut'][items[0]]=items[2]
                        st['density'][items[0]]=items[3]
                    #init Rc_in
                    st['Rc_in']={} ; st['N_ex']={} ; st['Rc_ex']={}
                    while 'C--------------' not in line: line=pt.readline()
                    while True:
                        items=pt.readline().split() #; print 'items=',items
                        if not items: break
                        st['Rc_in'][items[0]]=items[3]
                        st['N_ex'][items[0]]=int(items[4])
                        st['Rc_ex'][items[0]]=items[6]
                    #find clusters most similar to loop of combo structe of densest global cluster
                    while 'horizontal-axis:' not in line: line=pt.readline()
                    for i in range( int(dstg)): line=pt.readline()#jump to row of densest glob clus
                    rmsds=line.split()[1:] #rmsd vals betw combos and combo of densest glob. clus
                    tmp={}
                    for i in range(len(rmsds)): tmp[`i+1`]=float(rmsds[i])
                    #order tmp from lowest to higest rmsd value
                    items=[(v, k) for k, v in tmp.items()]
                    items.sort() #; Bye(items)
                    #find similar combos, that is, clust id's having combos with rmsd to combo of
                    #densest glob. clust below 1.25 times the lowest rmsd value found
                    rmsdcutoff=items[0][0]*1.25
                    st['similars']=[]
                    for item in items:
                        rmsd=item[0]
                        if rmsd<rmsdcutoff: st['similars'].append( item[1] )
                        else: break
                    #find average size of all similar combos
                    st['avsize']=0
                    for id in st['similars']: st['avsize']+=st['N_ex'][id]
                    st['avsize']=float(st['avsize'])/st['usedsnap']
                    st['av_Rc_in']=0.0
                    overlappingsize=0 #different than st['usedsnap']
                    for id in st['similars']:
                        st['av_Rc_in']+=float(st['Rc_in'][id])*st['size'][id]
                        overlappingsize+=st['size'][id]
                    st['av_Rc_in']=st['av_Rc_in']/overlappingsize
            os.system('/bin/rm -r '+wd)


    def addLength(self):

        """add sequence lengths for short and long sequences"""

        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            rmsinpf=joink(scd,ct['id'],'rmsinp')
            if isthere(rmsinpf):
                ct['L']=int(open(rmsinpf,'r').readline().split()[1])
            else:
                ct['L']=''
                sys.stderr.write('no '+rmsinpf+' found\n')
        
    
    def addProxTerm(self):

        """add variable indicating proximity of each unaligned region
        to either termini

        addProxTerm()
        Example:                  ******** (unaligned region)
                <----------- b ---------->
                *******************************************************
                <----------------------- L --------------------------->
        The unaligned region is closer to N-terminus, thus we measure
        from here. We report the ratio prox=2*b/L, 0 < prox <= 1"""

        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]['spk'] #; Bye(ct.keys())
            if not  ct['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            for n in ct['loopdat'].keys(): #cycle through all loops
                cn=ct['loopdat'][n]
                cn['prox']=''
                if cn['istooshort']: continue #loop is too short
                middle=float(cn['begin']+cn['end'])/2 #middle point of loop
                #If loop closer to C-terminus than to N-terminus, count distance from C-term
                if middle>self.__dict__[typesl]['L']/2: middle=self.__dict__[typesl]['L']-middle
                cn['prox']=middle/self.__dict__[typesl]['L']/2
                if cn['prox']>1:
                    print self.id,cn['prox'],cn['begin'],cn['end'],self.__dict__[typesl]['L']
        

    def readDispro(self):

        """add disorder predictions from dispro

        dispro has three different types of disorder: 'COILS', 'REM465', and 'HOTLOOPS'.
        For each type, we have a line of 0's and 1's with same length as sequence length.
        One '0' represents ordered prediction, and 1 for disordered."""

        for typesl in ('short','long'): #cycle through short and long sequences
            pattern=re.compile('[A-Z]')
            ct=self.__dict__[typesl]
            ct['disorder']={'dispro':{'outfile':'','COILS':'','REM465':'','HOTLOOPS':'','cummulative':''}}
            disfile=joink(scd,ct['id'],ct['id']+'.dispro') #;print disfile
            if not isthere(disfile):
                sys.stderr.write('file '+disfile+' does not exists\n')
                continue
            cd=ct['disorder']['dispro']
            cd['outfile']=disfile
            pt=open(disfile,'r')
            if not pt.readline(): continue #file exists but is empty
            for x in ('COILS','REM465','HOTLOOPS'):
                while x not in pt.readline(): pass
                for y in pt.readline().strip():
                    if pattern.match(y): cd[x]+='1'
                    else: cd[x]+='0'
            #cummulative is adding the '1' for the three different methods
            lista=[ int(x) for x in list(cd['COILS'])]
            for x in ('REM465','HOTLOOPS'):
                for i in range(len(cd[x])):
                    if cd[x][i]=='1': lista[i]+=1
            cd['cummulative']=''.join([`x` for x in lista])


    def avDisorder(self):

        """compute average intrinsic disorder for every loop

        For every intrinsic disorder predictor, average the prediction over all residues
        of the loop
        """

        for typesl in ('short','long'): #cycle through short and long sequences
            if not self.__dict__[typesl]['disorder']['dispro']['outfile']: continue
            if not  self.__dict__[typesl]['spk']['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            cl=self.__dict__[typesl]['spk']['loopdat']
            for n in cl.keys():
                cn=cl[n]
                cn['disorder']={'dispro':{'average':{'COILS':0.0,'REM465':0.0,'HOTLOOPS':0.0,'one':0.0,'semi-consensus':0.0}}}
                begin=cn['begin']-1
                end=cn['end']-1
                cav=cn['disorder']['dispro']['average']
                for x in ('COILS','REM465','HOTLOOPS'):
                    prediction=self.__dict__[typesl]['disorder']['dispro'][x]
                    for y in prediction[begin:end+1]:
                        if y=='1': cav[x]+=1
                prediction=self.__dict__[typesl]['disorder']['dispro']['cummulative']
                for y in prediction[begin:end+1]:
                    if int(y)>0: cav['one']+=1
                    if int(y)>1: cav['semi-consensus']+=1
                for x in ('COILS','REM465','HOTLOOPS','one','semi-consensus'):
                    cav[x]=cav[x]/(end+1-begin)
        #return False !prevents writing of pickle file
                
    def correctHEC(self):

        """correct erroneous secondary structure content values imported from loop.dat files"""

        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            pdbid=ct['id']
            ss=['1',] #will contain secondary structure. Shift by one
            for line in open(joink(scd,pdbid,'seq.dat'),'r').readlines():
                ss.append( line.split()[2] )
            if not  ct['spk']['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            cl=ct['spk']['loopdat']
            for n in cl.keys():
                cn=cl[n]
                begin=cn['begin']
                end=cn['end']
                h=0.0 ; e=0.0
                for x in ss[begin:end+1]:
                    if x=='2': h+=1
                    elif x=='4': e+=1
                h=h/(end+1-begin) ; e=e/(end+1-begin) ; c=1-h-e
                cn['h']=h ; cn['e']=e ; cn['c']=c


    def addHigestDensityModel(self):
        
        """add attribute containing info on TASSER densest model"""

        from utilities.small_utilities import unTARme
        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            cs=ct['spk']
            cs['Dmodel']={}
            if not cs['outfile']: continue #no spicker job
            cd=cs['Dmodel']
            cd['id']=cs['spkobj'].densest #id of densest cluster
            cd['combof']='combo'+`cd['id']`+'.pdb'
            if not isthere(joink(ct['simdir'],cd['combof'])):
                unTARme( cs['outfile'],wd=ct['simdir'], fileL=[cd['combof'],] ) #extract combo
            cd['combof']=joink(ct['simdir'],cd['combof']) #full path
    

    def amberminimize(self,jobname='',joblist=''):

        """run pulchra+dukka+amber minimization on combo with highest C-score

        We do not store new info on the pickle file"""

        from jobs.job import job,genJobList
        from utilities.small_utilities import junkName
        
        if not jobname: jobname='min'+junkName()
        if not joblist: joblist=genJobList()

        from utilities.codedir import nodecodedir
        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            cs=ct['spk']
            if not cs['outfile']: continue #no spicker job
            cd=cs['Dmodel']
            inputfile=cd['combof']
            outputfile='combo%02d.pdb.minimized'%(cd['id'])
            cd['combof.min']=joink(ct['simdir'],outputfile)
            jobnameII=jobname+typesl[0] #specify short and long
            Job=job(name='ca2at',exe='ca2at.py',
                    exed=joink(nodecodedir,'python/tertiary'),
                    args='-a '+inputfile+' -b '+outputfile,
                    outcmd='/bin/mv '+outputfile+' '+ct['simdir'])
            Job.qsub(jobnameII,ct['simdir'],wallt='0.249',mem_limit='256',ddisk='300',
                     joblist=joblist)
    

    def surface_area(self):

        """Calculate exposed surface area for every unaligned region.

        Using amber-minimized densest centroid (see amberminimize
        method) and 'surface' executable from jackal package
        NOTE: surfaces, expressed as a percentage of the standard amino
        acid surface (GLY-XXX-GLY), can be bigger than one        
        """

        from tertiary.surface import area3 #standard surface areas of residues
        from utilities.small_utilities import junkName
        surfacex='/library/jose/code/bin/jackal/bin/surface'
        junkf=junkName()
        
        for typesl in ('short','long'): #cycle through short and long sequences
            ptv=self.__dict__[typesl]
            pts=ptv['spk']
            if not pts['outfile']: continue
            ptc=pts['Dmodel'] #;print ptc
            ptc['area']={}
            pta=ptc['area']
            pta['value']=[-1.0,]*ptv['L'] #init surface area to -1
            pta['perctg']=[-1.0,]*ptv['L']
            modelf=ptc['combof.min'] #minimized densest centroid
            if not isthere(modelf): continue
            cmd=surfacex+' -prm 2  -out '+junkf+' '+modelf #;Bye(cmd)
            os.system(cmd)
            lines=open(junkf,'r').readlines()
            nline=0 ;  line=lines[nline]
            while '!RES' not in line:
                nline+=1 ; line=lines[nline]
            nline+=1 ; line=lines[nline] #;Bye(line)
            while line:
                items=line.split()
                if len(items)!=3:break
                res,n,s=items
                n=int(n)-1
                s=float(s) #;print res,n,s;Bye('')
                pta['value'][n]=s
                pta['perctg'][n]=s/area3[res] #as a percentage of total residue area
                if pta['perctg'][n]>1:
                    print s,area3[res],res
                nline+=1 ; line=lines[nline]
            while line:
                if 'total area' in line:
                    pta['total']=float(line.split(':')[1])
                    break
                nline+=1 ; line=lines[nline]
        os.system('/bin/rm '+junkf)


    def add_sequence(self):

        """insert one-letter sequence for short, long, and each unaligned region"""

        from seq.fastaManager import importFastaEntry

        pdbds=('/gpfs1/scratch/adrian/dat/pdb_jul07/input',
               '/gpfs1/scratch/adrian/dat/pdb_mar07/input',
               '/local/library/pdb_jul05/input/',
               )
        for typesl in ('short','long'):  #cycle through short and long sequences
            ptv=self.__dict__[typesl]    #;print ptv['id']
            ptv['sequence']=''
            for pdbd in pdbds:
                fastaf=joink(pdbd,ptv['id']) #fasta file
                if isthere(fastaf): #we found the file
                    header,ptv['sequence']=importFastaEntry(fastaf)
                    break
            if not ptv['sequence']:
                sys.stderr.write('ERROR (add_sequence) could not import sequence\n')
                continue
            pts=ptv['spk']
            if not pts['outfile']: continue
            ptl=pts['loopdat'] #;print ptl.keys()
            for n in ptl.keys():
                ptn=ptl[n] #;print ptn['L']
                begin=ptn['begin']-1
                end=ptn['end']
                ptn['sequence']=ptv['sequence'][begin:end]

        
    def surface_area_loops(self):

        """Calculate surface area of unaligned region"""

        from tertiary.surface import area1
        
        for typesl in ('short','long'): #cycle through short and long sequences
            ptv=self.__dict__[typesl]
            pts=ptv['spk']
            if not pts['outfile']: continue
            ptl=pts['loopdat']
            pta=pts['Dmodel']['area']
            if not pta: continue
            for n in ptl.keys(): #cycle through all loops
                ptn=ptl[n]
                if ptn['istooshort']: continue #loop is too short
                begin=ptn['begin']-1
                end=ptn['end']
                ptn['area']={}
                ptna=ptn['area']
                ptna['value']=pta['value'][begin:end][:]
                ptna['perctg']=pta['perctg'][begin:end][:]
                ptna['total']=0.0                
                for value in ptna['value']: ptna['total']+=value
                total=0.0
                for X in ptn['sequence']: total+=area1[X]
                ptna['total_perctg']=ptna['total']/total

                
    def printLoopSeqs(self):

        """print sequence of unaligned region

        Each line contains:
        id pdbid loopID L sequence"""

        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            cs=ct['spk'] #; Bye(ct.keys())
            if not  cs['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            cl=cs['loopdat']
            for n in cl.keys(): #cycle through all loops
                cn=cl[n]
                if cn['istooshort']: continue #do not output short loops
                sequence=cn['sequence']
                print self.id,ct['id'],'%2s %3d'%(n,cn['L']),sequence


    def loopClosureEntropy(self):

        """calculate the average closure entropy per residue for every loop

        Because we only have one template in the benchmark, there are
        no constraints in comb.dat and combCA.dat between two residues
        in the same unaligned region. Consequently, no loop closure
        entropy can be predicted."""

        from math import log
        
        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            combf=joink(ct['simdir'],'comb.dat')
            combCAf=joink(ct['simdir'],'combCA.dat')
            cs=ct['spk'] #; Bye(ct.keys())
            if not  cs['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            for n in cs['loopdat'].keys(): #cycle through all loops
                cn=cs['loopdat'][n]
                lce=0.0 #initialization
                if cn['istooshort']:
                    cn['loop_closure_entropy']=lce                    
                    continue #do not output short loops
                begin=cn['begin']
                end=cn['end']
                for f in (combf,combCAf):
                    if not isthere(f): continue
                    ptin=open(f,'r')
                    l=ptin.readline() #first line is just number of contacts
                    while True:
                        l=ptin.readline()
                        if not l: break #end-of-file reached
                        i,j,p=l.split()
                        if i>=begin and i<=end and j>=begin and j<=end: #contact within loop
                            lce+=float(p)*log(abs(i-j))
                    ptin.close()
                lce=lce/(end+1-begin) #divide by number of residues
                cn['loop_closure_entropy']=lce ;print lce
        return False

        
    def printLoopProps(self):

        """print properties for each unalig reg maybe relevant to rmsd2nat prediction

        printLoopProp()
        Each line contains:
        id pdbid loopid L prox exp h e c loopdat.n.rl loopdat.n.rg loopdat.n.rg2 loopdat.n.tonat.rl loopdat.n.tonat.rg loopdat.n.tonat.rg2 spklocal.n.avsize spklocal.n.av_Rc_in spkglobal.n.avsize spkglobal.n.av_Rc_in dispro_av_one dispro_av_semi-consensus
        Look to 'datatree' function for explanation of these quantities
        """

        for typesl in ('short','long'): #cycle through short and long sequences
            ct=self.__dict__[typesl]
            cs=ct['spk'] #; Bye(ct.keys())
            if not  cs['outfile']:
                sys.stderr.write('no out.spicker file for '+self.id+' _'+typesl+'\n')
                continue
            for n in cs['loopdat'].keys(): #cycle through all loops
                cn=cs['loopdat'][n]
                if cn['istooshort']: continue #loop is too short
                cnl=cs['spklocal'][n]
                cng=cs['spkglobal'][n]
                cnd=cn['disorder']['dispro']['average']
                cna=cn['area']
                print self.id,ct['id'],'%2s %3d %4.2lf %s %4.2lf %4.2lf %4.2lf %5s %5s %5s %5s %5s %5s %4.2lf %5.2lf %4.2lf %5.2lf %4.2lf %4.2lf %4.2lf'%(n,cn['L'],cn['prox'],cn['exp'],cn['h'],cn['e'],cn['c'],cn['rl'],cn['rg'],cn['rg2'],cn['tonat']['rl'],cn['tonat']['rg'],cn['tonat']['rg2'], cnl['avsize'],cnl['av_Rc_in'],cng['avsize'],cng['av_Rc_in'],cnd['one'],cnd['semi-consensus'],cna['total_perctg'])


    def addCA(self):

        """add path to native file"""

        from inputArgs.inputArgs import deglobb
        
        benchd='/gpfs1/scratch/adrian/dat/pdb_mar07'
        CA0=joink(benchd,'CA','xxxxx.pdb')    #globbed PDB file
        for typesl in ('short','long'):
            p0=self.__dict__[typesl]
            CA=deglobb(CA0,p0['id'])
            if not isthere(CA): CA=joink(benchd,'aux',header+'.pdb')
            if not isthere(CA): sys.stderr.write('ERROR no native file for '+self.id+'/'+p0['id'])
            p0['CA']=CA
        return True
    
        
    def tassSpicSim (self,switch='',joblist=''):

        """do TASSER and SPICKER simulations
        tassSpicSim(id='',switch='')
        switch=='non interacting loops'. Perform TASSER simulations with non interacting
        loops (cas_noninteracting_loops.f) followed by extended SPICKER (spk_rp)"""
        
        from tasser.tasser_utilities import calc_wall_time
        from inputArgs.inputArgs import deglobb
        from jobs.job import job,pastry,genJobList
        from seq.fastaManager import importFastaEntry
        
        if not joblist: joblist=genJobList()
        benchd='/gpfs1/scratch/adrian/dat/pdb_mar07'
        fastaf0=joink(benchd,'input','xxxxx') #globbed fasta file
        CA0=joink(benchd,'CA','xxxxx.pdb')    #globbed PDB file
        ibatch=0
        if switch=='non interacting loops':
            cas=joink(pcd,'cas_noninteracting_loops') #TASSER version for non interacting loops
            for typesl in ('short','long'):
                p0=self.__dict__[typesl]
                outdir=joink(p0['simdir'],'nil')
                p0['nil']={'simdir':outdir,} #add to the CAN object
                spkout=joink(p0['nil']['simdir'],p0['id']+'.out.spicker.tar.bz2')
                if isthere(spkout):
                    sys.stdout.write(self.id+'/'+p0['id']+' already completed\n')
                    continue #we already did the simulation
                jobname='nil'+p0['jobname']
                header=p0['id']
                fastaf=deglobb(fastaf0,header)
                CA=deglobb(CA0,header)#;print CA
                if not os.path.exists(CA):
                    CA=joink(benchd,'aux',header+'.pdb') #; print CA
                    fastaf=joink(benchd,'aux',header)
                header,seq=importFastaEntry(open(fastaf,'r'))
                inpdirA=p0['simdir']
                inpdirB=joink(scd,header) #contains exp.dat seq.dat rmsinp pair1.dat pair3.dat
                wt=calc_wall_time(len(seq))
                Job=job(exe='tas_spk.py',exed=pcd,args='-a '+inpdirB+' -b '+inpdirA+' -c '+outdir+' -d '+CA+' -e '+header+' -f '+cas)
                Job.qsub(jobname,outdir,wallt=wt,mem_limit='890',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        return True

    
#################################################################

class classSLl(CANl):

    """inherited class from CANl for short_long pairs, so that we can
    add methods to fuzz around the attributes of our particular CANl
    object implementation"""

    pass

objSLl=classSLl(slcollf,slf,slrepod) #holds info for target pairs


#################################################################

def datatree():

    """return the datatree of one short_long pair (a classSL(CAN) object)
        
    short_long (a classSL(CAN) object)
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory containing files
    |_'short','long': dictionary for info on 'short' and 'long' sequences
    | |_'id': PDB ID of sequence
    | |_'jobname': unique identifier for jobsubmission
    | |_'sequence': amino acid sequence
    | |_'simdir': directory holding files
    | |_'L': sequence length
    | |_'disorder': intrinsic-disorder related
    | | |_'dispro': predictions with dispro. '0' for order, '1' for disorder
    | |   |_'outfile': valid output file containing predictions from dispro
    | |   |_'COILS':  disordered coils
    | |   |_'REM465': predicted amino acids that would not show in an X-ray structure
    | |   |_'HOTLOOPS': disordered loops
    | |   |_'cummulative': summ of each entry of the previous three
    | |_'tasser': dictionary holding info on TASSER run
    | | |_'outfile':TASSER output file (*.tar.bz2)
    | |_'spk': dictionary holding info on SPICKER run
    | |  |_'outfile': spicker output file
    | |  |_'spkobj': spicker object
    | |  |_'Dmodel': dictionary containing info of densest model
    | |  |  |_'id': cluster ID
    | |  |  |_'combof': combo file (full path)
    | |  |  |_'combof.min': combo file rebuilt by dukka+pulchra+amberminimiz (full path)
    | |  |    |_'area': dictionary containing info on sovent accesible area
    | |  |      |_'total': total surface area
    | |  |      |_'value': list with absolute values of solvent accesible area per residue
    | |  |      |          ( residue one is first item in list (index list zero) )
    | |  |      |_'perctg': same as 'value', but as a percentage of total residue area
    | |  |_'loopdat': dictionary containing info on loop.dat file
    | |  |  |_n: loop id
    | |  |    |_'begin': beginning of loop. Index begin at 1, not 0
    | |  |    |_'end': end of loop
    | |  |    |_'L': loop length
    | |  |    |_'istooshort': logical variable marks small loop
    | |  |    |_'sequence': amino acid sequence of the loop
    | |  |    |_'exp': predicted exposed surface
    | |  |    |_'h': percentage predicted helix
    | |  |    |_'e': percentage predicted strand
    | |  |    |_'c': percentage predicted coil
    | |  |    |_'loop_closure_entropy': similar to contact order, (comb.dat and combCA.dat)
    | |  |    |_'disorder': intrinsic disorder related
    | |  |    | |_'dispro':  dispro disorder predictor tool
    | |  |    |    |_'average': average disorder content
    | |  |    |      |_'COILS':  disordered coils
    | |  |    |      |_'REM465': predicted residues that would not show in an X-ray structure
    | |  |    |      |_'HOTLOOPS': disordered loops
    | |  |    |      |_'one': average for 'COILS'||'REM465'||'HOTLOOPS'
    | |  |    |      |_'semi-consensus': average requiring more than one disorder type equal
    | |  |    |_'rl': local rmsd for densest cluster
    | |  |    |_'lst': length of stems used for global rmsd
    | |  |    |_'rg': global rmsd (align both stems and loop)
    | |  |    |_'rg2': global rmsd (align only stems)
    | |  |    |_'tonat': dictionary holding comparisons to native
    | |  |    | |_'rl': local rmsd
    | |  |    | |_'rg': global rmsd (align both stems and loop)
    | |  |    | |_'rg2': second measure of global rmsd (align only stems)
    | |  |    |_'area': dictionary containing info on sovent accesible area
    | |  |      |_'value': list with absolute values of solvent accesible area per residue
    | |  |      |          ( residue one is first item in list (index list zero) )
    | |  |      |_'perctg': same as 'value', but as a percentage of total residue area
    | |  |      |_'total': total solvent accesible area
    | |  |      |_'total_perctg': same as total, but as a percentage of total residue area
    | |  |_'spklocal','spkglobal': clustering of loop conformations using local or global rmsd
    | |    |_n: loop id
    | |      |_'L': loop length
    | |      |_'nc': number of clusters
    | |      |_'availsnap': number of structures in the trajectory files
    | |      |_'usedsnap': number of used structures in the clustering, corresponding to the
    | |      |             population of the densest global cluster
    | |      |_'combo2nat': dictionary with rmsd between native and combo
    | |      |_'closc2nat': dictionary with rmsd between native and closc
    | |      |_'size': dictionary with cluster sizes
    | |      |_'R_cut': dictionary with rmsd cutoff for clustering
    | |      |_'density': dictionary with cluster densities
    | |      |_'Rc_in': dictionary with cluster size measure as rmsd value
    | |      |_'similars': list of clusters id's with combos similar to loop piece of combo of
    | |      |             densest cluster
    | |      |_'avsize': average size of similar clusters
    | |      |_'av_Rc_in': weighted average of Rc_in of similar clusters
    | |_'nil': info on simulation with non-interacting loops
    |   |_'simdir': directory to store results
    |   |_'tasser':
    |   | |_'outfile':
    |   |_'spk':
    |     |_'outfile':
    |     |_'spkobj':
    |_'templ': dictionary for info on common template
      |_id: PDB ID of selected template
    """

    print globals()['datatree'].__doc__

def initOne(id,tmpls={}):

    if not tmpls:
        file=joink(scd,'../findCommonTemplate.dat')
        for line in open(file,'r').readlines():
            short,long,templ=line.split()
            tmpls[short+'_'+long]=templ
    
    c=classSL(id=id)
    c.genDumpf(objSLl.repod)
    c.simdir=joink(scd,id)
    c.templ={'id':tmpls[id],}

    for type in ('short','long'):
        pdbid=id[:5]
        if type=='long': pdbid=id[6:]
        c.__dict__[type]={'id':pdbid,'tasser':{},'spk':{},'simdir':joink(scd,id,pdbid)}
        ct=c.__dict__[type] #a handy reference
        ct['tasser']['outfile']=''
        outfile=joink(ct['simdir'], ct['id']+'.out.tasser.tar.bz2')
        if isthere(outfile): ct['tasser']['outfile']=outfile
        ct['spk']['outfile']=''
        outfile=joink(ct['simdir'], ct['id']+'.out.spicker.tar.bz2')
        if isthere(outfile): ct['spk']['outfile']=outfile
        
    c.pickleDump()
    return c

    
def initColl():

    """initialize the dump files for (short,long) pairs"""

    #first find template for each short_long
    tmpls={}
    file=joink(scd,'../findCommonTemplate.dat')
    for line in open(file,'r').readlines():
        short,long,templ=line.split()
        tmpls[short+'_'+long]=templ
    #now initialize short_long objects
    remaining=len(objSLl.canlist)
    for id in objSLl.canlist:
        initOne(id,tmpls=tmpls)
        print remaining,id
        remaining-=1

    #objSLl.updColletion()


def doitall(id):

    """a recount of everything we had done on each varset"""

    methods=('AddSpkobj','addLoopdat','addSpklg','addLength','addProxTerm','readDispro',
             'correctHEC','addHigestDensityModel','amberminimize','surface_area',
             'add_sequence','surface_area_loops','loopClosureEntropy','addSimOut')
    
    objSL=initOne(id)
    for method in methods:
        print method
        getattr(objSL,method)()
    objSL.pickleDump()


def AddSpkobj(id='',target=''):

    """add spicker object attribute for each variant"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.AddSpkobj(target=target)!=False: objSL.pickleDump()
    else:
        remaining=len(objSLl.canlist)
        for objSL in objSLl.iterate():
            print remaining,objSL.id
            if objSL.AddSpkobj(target=target)!=False: objSL.pickleDump()
            remaining=remaining-1
 

def addLoopdat(id=''):

    """add loop properties from loop.dat file"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addLoopdat()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('addLoopdat',silent=True)
    

def addSpklg(id=''):

    """add spicker clustering of loop conformations using local or global rmsd"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addSpklg()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('addSpklg',silent=True)


def addLength(id=''):

    """add sequence lengths for short and long sequences"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addLength()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('addLength',silent=True)


def addProxTerm(id=''):

    """add variable indicating proximity of each unaligned region
    to either termini

    addProxTerm(id='')
    Example:                  ******** (unaligned region)
            <----------- b ---------->
            *******************************************************
            <----------------------- L --------------------------->
    The unaligned region is closer to N-terminus, thus we measure
    from here. We report the ratio prox=2*b/L, 0 < prox <= 1"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addProxTerm()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('addProxTerm',silent=False)


def readDispro(id=''):

    """add disorder predictions from dispro

    dispro has three different types of disorder: 'COILS', 'REM465', and 'HOTLOOPS'.
    For each type, we have a line of 0's and 1's with same length as sequence length.
    One '0' represents ordered prediction, and 1 for disordered."""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.readDispro()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('readDispro',silent=True)


def avDisorder(id=''):
    
    """compute average intrinsic disorder for every loop

    For every intrinsic disorder predictor, average the prediction over all residues
    of the loop
    """

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.avDisorder()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('avDisorder',silent=True)


def correctHEC(id=''):

    """correct erroneous secondary structure content values imported from loop.dat files"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.correctHEC()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('correctHEC',silent=True)


def addHigestDensityModel(id=''):
        
    """add attribute containing info on TASSER densest model"""
    
    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addHigestDensityModel()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('addHigestDensityModel',silent=False)
        

def amberminimize(id='',jobname=''):

    """run pulchra+dukka+amber minimization on combo with highest C-score"""

    from jobs.job import genJobList
    from utilities.small_utilities import junkName
    joblist=genJobList()
    if id:
        if not jobname: jobname='min'+junkName()
        objSL=objSLl.loadCAN(id)
        if objSL.amberminimize(jobname=jobname,joblist=joblist)!=False: objSL.pickleDump()
    else:
        ijob=0
        for objSL in objSLl.iterate():
            if objSL.amberminimize(jobname='min'+'%05d'%(ijob,),joblist=joblist)!=False:
                objSL.pickleDump()
            ijob+=1


def surface_area(id=''):

    """Calculate exposed surface area for every unaligned region.

    Using amber-minimized densest centroid (see amberminimize
    method) and 'surface' executable from jackal package"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.surface_area()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('surface_area',silent=True)


def add_sequence(id=''):

    """insert one-letter sequence for short, long, and each unaligned region"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.add_sequence()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('add_sequence',silent=False)

       
def surface_area_loops(id=''):

    """Calculate surface area of unaligned region"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.surface_area_loops()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('surface_area_loops',silent=False)


def printLoopSeqs(id=''):

    """print sequence of unaligned region
    
    Each line contains:
    id pdbid loopID L sequence"""
    print '#    id pdbid loopID L sequence'
    
    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.printLoopSeqs()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('printLoopSeqs',silent=True)


def loopClosureEntropy(id=''):

    """calculate the average closure entropy per residue for every loop

    Because we only have one template in the benchmark, there are no
    constraints in comb.dat and combCA.dat between two residues in the
    same unaligned region. Consequently, no loop closure entropy can
    be predicted."""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.loopClosureEntropy()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('loopClosureEntropy',silent=True)


def printLoopProps(id=''):

    """print properties for each unalig reg maybe relevant to rmsd2nat prediction

    printLoopProp()
    Each line contains:
    id pdbid L prox exp h e c loopdat.n.rl loopdat.n.rg loopdat.n.rg2 loopdat.n.tonat.rl loopdat.n.tonat.rg loopdat.n.tonat.rg2 spklocal.n.avsize spklocal.n.av_Rc_in spkglobal.n.avsize spkglobal.n.av_Rc_in dispro_av_one dispro_av_semi-consensus
    Look to 'datatree' function for explanation of these quantities"""
    print '############################################################################'
    print '#prox: proximity of loop to either termini'
    print '#exp: predicted fraction of surface that is exposed to solvent'
    print '#h e c: secondary structure content for helix (h), strand (e) and coil (c)'
    print '#loopdat.n.rl: average local rmsd of global cluster member to cluster center'
    print '#loopdat.n.rg: first measure of averag global rmsd of global cluster member to cluster center'
    print '#loopdat.n.rg2: second meas of averag global rmsd of global cluster member to cluster center'
    print '#loopdat.n.tonat.rl: local rmsd of loop model to native, for most dense global cluster'
    print '#loopdat.n.tonat.rg: first meas of global rmsd of loop to native, for most dense global clust'
    print '#loopdat.n.tonat.rg2:second  meas of global rmsd of loop to nat, for most dense global clust'
    print '#spklocal.n.avsize: number density of loop confs similar to loop conf in most dense glob cl\n#                   after spicker clustering using local rmsd'
    print '#spklocal.n.av_Rc_in: after spicker clustering using local rmsd, average Rc_in of clusters\n#                     having centroids similar to loop conf in most dense global cluster'
    print '#spkglobal.n.avsize: number density of loop confs similar to loop conf in most dense glob cl\n#                   after spicker clustering using global rmsd'
    print '#spkglobal.n.av_Rc_in: after spicker clustering using global rmsd, average Rc_in of clusters\n#                     having centroids similar to loop conf in most dense global cluster'
    print '#dispro_av_one: intrinsic disorder prediction of dispro program averaged over the loop. Dispro has three\n#               predictors and we assign a residue as disordered if at least one predictor assigns it'
    print '#dispro_av_semi-consensus: we assign a residue as disordered if more than one predictor assigns it'
    print '#sas: solvent accessible surface'
    print '############################################################################'
    print '#(1)id (2)pdbid (3)loopID (4)L (5)prox (6)exp (7)h (8)e (9)c (10)loopdat.n.rl'
    print '#(11)loopdat.n.rg (12)loopdat.n.rg2 (13)loopdat.n.tonat.rl (14)loopdat.n.tonat.rg'
    print '#(15)loopdat.n.tonat.rg2 (16)spklocal.n.avsize (17)spklocal.n.av_Rc_in'
    print '#(18)spkglobal.n.avsize (19)spkglobal.n.av_Rc_in (20)dispro_av_one'
    print '#(21)dispro_av_semi-consensus (22) sas'
    print '#    1       2     3   4   5     6     7    8   9     10    11    12    13    14    15    16   17   18    19   20   21  22'
    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.printLoopProps()!=False: objSL.pickleDump()
    else:
        objSLl.cycleOverList('printLoopProps',silent=True)
    

def tassSpicSim(id='',switch=''):

    """do TASSER and SPICKER simulations

       tassSpicSim(id='',switch='')

       switch=='non interacting loops'. Perform TASSER simulations with non interacting
    loops (cas_noninteracting_loops.f) followed by extended SPICKER (spk_rp)"""
   
    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.tassSpicSim(switch=switch)!=False: objSL.pickleDump()
    else:
        from jobs.job import genJobList
        joblist=genJobList()
        remaining=len(objSLl.canlist)
        for objSL in objSLl.iterate():
            print remaining,objSL.id
            if objSL.tassSpicSim(switch=switch,joblist=joblist)!=False: objSL.pickleDump()
            remaining=remaining-1
             

def view(id):
    
    """view object properties

    view(id)"""

    objSLl.loadCAN(id).view()
        

def addSimOut(id='',target=''):
    
    """add tasser, spk dictionaries"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addSimOut(target=target)!=False: objSL.pickleDump()
    else:
        remaining=len(objSLl.canlist)
        for objSL in objSLl.iterate():
            print remaining,objSL.id
            if objSL.addSimOut(target=target)!=False: objSL.pickleDump()
            remaining=remaining-1
            

def addJobname():

        """add an identifier for simulations"""
        
        remaining=len(objSLl.canlist)
        for objSL in objSLl.iterate():
            print remaining,objSL.id
            if objSL.addJobname('j%05d'%(remaining,))!=False: objSL.pickleDump()
            remaining=remaining-1
            

def addCA(id=''):

    """add path to native file"""

    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.addCA()!=False: objSL.pickleDump()
    else:
        remaining=len(objSLl.canlist)
        for objSL in objSLl.iterate():
            print remaining,objSL.id
            if objSL.addCA()!=False: objSL.pickleDump()
            remaining=remaining-1


def submitSpklg(id='',switch='normal run'):

    """Submit spicker job only"""

    from jobs.job import job,genJobList
    
    if id:
        objSL=objSLl.loadCAN(id)
        if objSL.submitSpklg(switch=switch)!=False: objSL.pickleDump()
    else:
        from jobs.job import genJobList
        joblist=genJobList() ; ibatch=0 ; remaining=len(objSLl.canlist)
        for objSL in objSLl.iterate():
            print remaining,objSL.id
            ibatch,record=objSL.submitSpklg(joblist=joblist,switch=switch,ibatch=ibatch)
            if record!=False: objSL.pickleDump()
            remaining=remaining-1
 

def help(*kargs):

    """list available services"""
    
    services=['initOne','datatree','initColl','AddSpkobj','addLoopdat',
              'addSpklg','addProxTerm','addLength','readDispro','correctHEC',
              'avDisorder','addHigestDensityModel','amberminimize','surface_area',
              'add_sequence','surface_area_loops','loopClosureEntropy',
              'tassSpicSim','view','viewall','addSimOut','addJobname','addCA',
              'submitSpklg']
    os.system('/usr/bin/clear')
    if not kargs:
        print 'Services='
        for service in services:
            print '  '+service+' : '+globals()[service].__doc__.split('\n')[0]
        print 'Type "workSLs.py -a help -b service" for more detail on "service"'
    else:
        sv=kargs[0]
        print globals()[sv].__doc__        


if __name__=='__main__':
    inpHand('Usage: workSLs.py [options]',
            ' -a _R_service required service. Type "help" for list of services.',              
            ' -b __servargs service required arguments in python syntax, and separate them with ",". Escape symbols \' and " (def:None)',
            ' -c __servoptargs service optional arguments (Example: "key=val" for one argument and "key1=val1, key2=val2,..." for two or more arguments (def:None)',
            ).parse(locals(),sys.argv)

    #pdb.set_trace() #for debugging purposes only
    #parse required arguments [arg1,arg2,...]
    reqargs=[]
    if servargs:
        reqargs=servargs.split(',')
    #parse optional arguments {key1:val1, key2:val2,...}
    optsargs={}
    if servoptargs:
        for pair in servoptargs.split(','):
            key,val=pair.split('=')
            optsargs[key]=val
    #execute the service
    locals()[service](*reqargs,**optsargs)

