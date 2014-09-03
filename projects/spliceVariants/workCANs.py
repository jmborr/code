#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,nodecodedir
from utilities.small_utilities import chomp,Bye
from utilities.readingWritingFiles import arethere
from inputArgs.inputArgs import inpHand
from jobs.job import job

"""Manage the splice variants"""

#global variables
joink=os.path.join
isthere=os.path.exists
basename=os.path.basename
pcd=joink(codedir,'projects','spliceVariants')
scd=joink(scratchdir,'spliceVariants','out2')
collf=joink(scd,'varsets.can')
#canlistf=joink(scd,'../preparing_database/refs300_nonduplicated.list')
#Next canlist does not include the 45 varsets with discrepancies
#between sequences derived from the Swiss-Prot entries and what I had
#derived from uniprot_sprot.fasta and uniprot_sprot_varsplic.fasta
#files
canlistf=joink(scd,'../preparing_database/refs300_nonduplicated.reduced.list')
repod=joink(scd,'repod')
alignx=codedir+'/f77/seq/align'


class PEPSIl(CANl):

    """inherited class from CANl, so that we can add methods to fuzz
    around the attributes of our particular CANl object
    implementation"""

    pass

cl=PEPSIl(collf,canlistf,repod) #list of PEPSI objects
canlist=chomp(open(cl.canlistf,'r').readlines())
    
def readSwissProtXML(acc):

    """retrieve XML entry for a Swiss Prot accesion number

    readSwissProtXML(acc) returns the Swiss Prot entry for 'acc'
    accession number as a single string"""

    db=joink(scd,'../preparing_database/swissProtDB/uniprot_sprot.xml') #database
    idx=joink(scd,'../preparing_database/swissProtDB/uniprot_sprot.xml.idx') #index file

    for l in open(idx,'r').readlines(): #cycle over the index file
        if acc in l: #found accession number
            byte=int(l.split(':')[0]) #byte position in database
            break
    pin=open(db,'r') #open database
    pin.seek(byte)   #goto byte position
    entry=''
    l=pin.readline()
    while True:
        entry+=l
        if  '</entry>' in l: break
        l=pin.readline()
    return entry

    
def readVarsetsDat():

    """read varsets.dat file into dictionary object"""
    
    tmpl={}
    tmpl2={}
    for line in chomp(open( scd+'/../preparing_database/varsets.dat' ,'r').readlines()):
        if line[0]=='#':continue
        list=line.split()
        header=list[0]
        tmpl[header]=list[1:]
        tmpl2[header]=[ header+var for var in list[1:] ]
    return tmpl,tmpl2

 
def feature2vars(varseqf):

    """ extract variant ID's in the VAR_SEQ feature list

    return the variant ID's as two-digits ID's in a list"""
    
    pattern=re.compile('isoform (\d+)')
    description=varseqf[3]
    return [ '%02d'%(var) for var in pattern.findall(description) ]


#%%%%%%%%%%%%%%%%%% BEGIN CLASS %%%%%%%%%%%%%%%%%%%%%%%%


class PEPSI(CAN):

    """inherited class from CAN, so that we can add methods to fuzz
    around the attributes of our particular CAN object
    implementation"""


    def removeAtribute(self):

        """remove some attribute that we mistakenly create"""
        for var in self.varid2:
            del self.__dict__['vars'][var]['addSpk']


    def addSeqAndL(self):

        """add sequence as single string, as well as sequence length,
        for each variant"""

        from seq.fastaManager import importFastaEntry
        for var in self.varid8: #cycle through every variant
            h,seq=importFastaEntry( open(scd+'/../input/'+var,'r') )
            var2=var[6:]
            self.vars[var2]['seq']=seq    #sequence
            self.vars[var2]['L']=len(seq) #sequence length

    def addVar8(self):

        """add var8 field to each variant"""

        for var in self.varid2:
            self.vars[var]['var8']=self.id+var
            

    def addIsoname(self):

        """find the name given to each splice variant in file uniprot_sprot_varsplic.fasta"""
        
        db=joink(scd,'../preparing_database/swissProtDB/uniprot_sprot_varsplic.fasta')
        idx=joink(scd,'../preparing_database/swissProtDB/uniprot_sprot_varsplic.fasta.idx')
        
        for var in self.varid2:
            varref=self.vars[var] #handy reference
            if var=='00': #there is no name in uniprot_sprot_varsplic.fasta for the reference
                varref['isoname']='0' #assign by hand
                continue
            varref['isoname']='' #assume no name
            #search for byte offset in index file
            varraw=self.id+'-'+`int(var)`
            for l in open(idx,'r').readlines():
                if varraw in l:
                    byte=int(l.split(':')[0])
                    break
            #search for name in uniprot_sprot_varsplic.fasta file
            pin=open(db,'r')
            pin.seek(byte) #jump to appropriate byte offset
            match=re.compile('Splice isoform\s(.+)\sof\s'+self.id).search(pin.readline())
            if match:
                varref['isoname']=match.group(1)
                print '    '+varref['isoname']
            else:
                sys.stderr.write('no match found for '+varref['var8']+'\n')
            pin.close()

            
    def addSwissp(self):

        """Retrieve Swiss Prot entry (XML format) into a single string"""

        self.swissp=readSwissProtXML(self.id)


    def getSwissp(self):

        """return swissp field"""

        return self.swissp
        

    def addAlignmentToRef_VAR_SEQ(self):
        
        """seq. alig. object. Alignment done with VAR_SEQ from Swiss-Prot

        add sequence alignment of each variant to its reference as an
        alignObj sequence alignment object. The alignment is done from
        the VAR_SEQ annotations of the Swiss-Prot entry for the
        varset."""

        from xml.dom.minidom import parseString,Node
        from seq.alignYangManageResults import gappedSeq,alignment

        if not self.checkRef(): return #discrepancy
        
        #print self.swissp
        #list of splice variant features
        spfs=[]
        doc=parseString(self.swissp) #read swissp string
        for node in doc.getElementsByTagName("feature"):
            if node.getAttribute("type")=="splice variant": spfs.append(node)
                
        #init dictionary that assigns features to variants. Every
        #variant has a list of features that affect it
        vars2f={}
        for var in self.varid2: #cycle through list of variants
            varref=self.vars[var] #a handy reference            
            vars2f[var]=[]     #init list of associated features            
            for spf in spfs:   #cycle through splice variant features
                insert=re.escape(varref['isoname']) #escape non alphanumeric characters
                pattern=re.compile('isoform\s'+insert+'[\s,\)]')
                if pattern.search(spf.getAttribute("description")):
                    vars2f[var].append(spf)

        #implement the splicing events for each variant
        refseq=self.vars['00'] #reference sequence
        #print refseq['seq']
        for var in self.varid2: #cycle through list of variants
            #print var, self.vars[var]['isoname']
            #init alignment. Both reference and variant same initially
            alg=alignment(gappedSeq(refseq['seq']),gappedSeq(refseq['seq']))
            for spf in vars2f[var]: #cycle through its list of splice variant features
                #retrieve beginning and ending positions of the splicing event
                location=spf.getElementsByTagName("location")[0]
                #print 'location=\n',location.__dict__
                pos=location.getElementsByTagName("position")
                begin=-1
                end=-1
                if pos:
                    begin=int(pos[0].getAttribute("position"))
                    end=begin
                else:
                    begin=int(location.getElementsByTagName("begin")[0].getAttribute("position"))
                    end=int(location.getElementsByTagName("end")[0].getAttribute("position"))
                #print 'begin=',begin,' end=',end
                #retrieve the variant sequence, if any
                #print 'spf=',spf
                #print spf.__dict__
                variation=spf.getElementsByTagName("variation")
                #print 'VARIATION=',variation
                if variation: #variant sequence is a mutation
                    mutaseq=variation[0].__dict__['childNodes'][0].data
                else:         #the original sequence is missing altogether
                    mutaseq=''
                #print begin,end
                #print 'mutaseq='+mutaseq                
                #add splicing event to current alignment. Splice first
                #sequence. Notice shift in the indexes
                alg=alg.addSplicingEvent(begin-1,end-1,variation=mutaseq,i=1)
                #print alg
            self.vars[var]['alg2refSP']=alg
            Ltrue=alg.gs[2].sL
            Lfalse=self.vars[var]['L']
            #if Ltrue!=Lfalse:
            if alg.gs[2].s != self.vars[var]['seq']:
                sys.stderr.write(self.vars[var]['var8']+' '+self.vars[var]['isoname']+' variant correct:'+`Ltrue`+' wrong:'+`Lfalse`+'\n')
                #os.system('align '+alg.gs[2].s+' '+self.vars[var]['seq']+' 3')
                

    def checkRef(self):

        """check reference sequence is that one given by the Swiss-Prot entry"""

        from xml.dom.minidom import parseString,Node
        from seq.alignYangManageResults import gappedSeq,alignment
        
        doc=parseString(self.swissp).getElementsByTagName("entry")[0]
        for node in doc.childNodes:
            if node.localName=='sequence':
                seq_sp=node.childNodes[0].data.replace('\n','') #sequence from xml file
                break

        if seq_sp != self.vars['00']['seq']:
            Ltrue=len(seq_sp)
            Lfalse=len(self.vars['00']['seq'])
            #print self.swissp
            sys.stderr.write(self.id+" (correct:"+`Ltrue`+" wrong:"+`Lfalse`+")\n")
            return False
            #if Ltrue==Lfalse:
            #    print alignment(gappedSeq(self.vars['00']['seq']),gappedSeq(seq_sp))
        return True

        
    def addAlignmentToRef(self):

        """add sequence alignment of each variant to its reference as
        an alignObj sequence alignment object"""

        from seq.alignYangManageResults import alignOut
        seq0=self.vars['00']['seq'] #reference sequence
        for var in self.varid2: #cycle through every variant
            seq1=self.vars[var]['seq']
            cmd=alignx+' '+seq0+' '+seq1+' 3' #alignment command
            self.vars[var]['alg2ref']=alignOut( os.popen(cmd).readlines() ).alg #alignment object
            self.pickleDump() #update dump file


    def AddSeqIDtoRef(self):

        """add sequence identity between reference and each var"""

        for var in self.varid2: #cycle through every variant
            aligment=self.vars[var]['alg2ref']
            self.vars[var]['seqid']=aligment.seqIdentity()


    def reportSeqID(self):

        """print sequence identity between reference and each var"""

        for var in self.varid2: #cycle through every variant
            print '%s %4.2lf'%(self.vars[var]['var8'],self.vars[var]['seqid'])
            

    def addProspector(self):

        """add prospector dictionary to each variant. It will hold
        info related to PROSPECTOR run"""

        for var in self.varid2: #cycle through every variant
            self.vars[var]['prosp']={}


    def addProspOutFile(self):

        """add the path to the prospector output file, if extant"""

        for var in self.varid2: #cycle through every variant
            f=joink(scd,self.id,var,self.id+var+'.out.prospector.tar.bz2')
            self.vars[var]['prosp']['outfile']='' #init as not extant
            if isthere(f): self.vars[var]['prosp']['outfile']=f
       
        
    def addEasyMedmHard(self):
        
        """add PROSPECTOR easy, medm, hard classification

        add an attribute to classify PROSPECTOR difficulty for each
        variant (easy, medm, hard)"""
        from utilities.small_utilities import unTARme
        for var in self.varid2: #cycle through every variant
            pd=self.vars[var]['prosp'] #handy reference
            if 'emh' in pd.keys() and pd['emh']: continue #we already made the assignment
            pd['emh']='' #no prospector output, thus no type
            if pd['outfile']: #there's a prospector output file
                wd=unTARme(pd['outfile']) #untar to some working directory
                sumd=joink(wd,'pdbbsummary') #summary directory
                pd['emh']='hard' #assume a hard case
                for f in ('listbest','list_verygood','list_fair'):
                    x=open(joink(sumd,f),'r').readline() #read first line
                    if '1' in x: #we found assignment
                        pd['emh']='easy' #assume easy (listbest and list_verygood cases)
                        if f=='list_fair': pd['emh']='medm'
                        break
                os.system('/bin/rm -r '+wd) #remove temporary working directory

    def addTasser(self):

        """add tasser dictionary to each variant.

        It will hold info related to TASSER run"""

        for var in self.varid2: #cycle through every variant
            self.vars[var]['tasser']={}

    def addTasserInFile(self):

        """add the path to the tasser input file, if extant"""

        for var in self.varid2: #cycle through every variant
            f=joink(scd,self.id,var,self.id+var+'.in.tasser.tar.bz2')
            self.vars[var]['tasser']['infile']='' #init as not extant
            if isthere(f): self.vars[var]['tasser']['infile']=f


    def tassSpicSim(self,switch='',joblist='',ibatch=0):

        """submit TASSER & spicker simulations"""

        from utilities.small_utilities import unTARme,junkName
        from tasser.chainDatReader import chainDatReader
        from prospector.prospManager import templ_file
        from jobs.job import job,pastry,genJobList
        from tasser.tasser_utilities import calc_wall_time
        
        if not joblist: joblist=genJobList()
        if switch=='non interacting loops':
            for var in self.varid2: #cycle through every variant
                try:
                    pt=self.vars[var]
                    #we will identify which template in chain.dat is most similar to the
                    #PROSPECTOR template identified in the Cmodel
                    if not pt['tasser']['infile'] or not pt['prosp']['outfile']: continue
                    wd=unTARme(pt['tasser']['infile'],fileL=['chain.dat',]) #untar to some working directory
                    chaindat=chainDatReader(chainf=wd+'/chain.dat',fastaf=pt['fastafile'])
                    unTARme(pt['prosp']['outfile'],wd=wd,fileL=['pdbbpdborienrev',])
                    templobj=templ_file( joink(wd,'pdbbpdborienrev/'+pt['var8']+'rap3orienrev5s.pdb') )
                    templobj.loadTemplates()
                    pt['nil']={'simdir':joink(pt['simdir'],'nil'),}
                    os.system('/bin/mkdir -p '+pt['nil']['simdir'])
                    for header in templobj.tpls.keys():
                        if header==pt['spk']['Cmodel']['pdbtempl']:
                            ptempl=junkName() #file to store the prospector template
                            templobj.templ2file(header=header,fname=ptempl)
                            ttempl=joink(pt['nil']['simdir'],'chain.%s.dat'%(header,))
                            tmmax=0.0
                            while chaindat.readTempl():
                                tm=chaindat.TMscore(struct=ptempl) #TM-score between TASSER and PROSP. templ
                                if tm>tmmax:
                                    tmmax=tm  ;  chaindat.templToChaindat(fname=ttempl)
                    pt['nil']['tasser']={'tm_prospTempl2tassTempl':tmmax,'chaindat':ttempl}
                    os.system('/bin/rm -r %s %s'%(ptempl,wd))
                    #gather necessary TASSER input in a compressed file
                    wd=unTARme(pt['tasser']['infile'])
                    os.system('cd %s && /bin/cp %s ./chain.dat && tar jcf %s.in.tasser.tar.bz2 * && /bin/mv %s.in.tasser.tar.bz2 %s'%(wd,pt['nil']['tasser']['chaindat'],pt['var8'],pt['var8'],pt['nil']['simdir']))
                    pt['nil']['tasser']['infile']=joink(pt['nil']['simdir'],pt['var8']+'.in.tasser.tar.bz2')
                    os.system('/bin/rm -r '+wd)
                    #submit
                    cas=joink(pcd,'benchmarkByTMalign/cas_noninteracting_loops')
                    spicker=joink(pcd,'benchmarkByTMalign/spk_rp') #spicker evaluation of unaligned regions
                    spkout=joink(pt['nil']['simdir'],'out.spicker.tar.bz2')
                    jobname='nil'+pt['var8'] ; wt=calc_wall_time(len(pt['seq']))
                    outdir=pt['nil']['simdir']
                    Job=job(exe='tas_spk2.py',exed=pcd,args='-a %s -b %s -f %s -g %s -i %s'%\
                            (pt['fastafile'],pt['nil']['tasser']['infile'],outdir,cas,spicker))
                    Job.qsub(jobname,outdir,wallt=wt,mem_limit='890',joblist=joblist)
                    ibatch=Job.dormant(ibatch)
                except:
                    pass
            return True,ibatch
            
            
    def addTasserEasyMedmHard(self):
    
        """add TASSER easy, medm, hard classification
        
        add an attribute to classify PROSPECTOR difficulty for each
        variant (easy, medm, hard)"""

        from utilities.small_utilities import unTARme
        for var in self.varid2: #cycle through every variant
            pt=self.vars[var]['tasser'] #handy reference
            if 'emh' in pt.keys() and pt['emh']: continue #we already made the assignment
            pt['emh']=''
            if pt['infile']: #there's a tasser input file
                wd=unTARme(pt['infile'],fileL=['chain.dat',]) #untar to some working directory
                pt['emh']=open(joink(wd,'chain.dat'),'r').readline().strip().split()[1]
                os.system('/bin/rm -r '+wd) #remove temporary working directory
                
        
    def addTasserOutFile(self):

        """add the path to the tasser output file, if extant"""

        for var in self.varid2: #cycle through every variant
            f=joink(scd,self.id,var,self.id+var+'.out.tasser.tar.bz2')
            self.vars[var]['tasser']['outfile']='' #init as not extant
            if isthere(f): self.vars[var]['tasser']['outfile']=f


    def addZ1(self):

        """add Z score of first PROSPECTOR template

        This field will help us later to select variants according to
        a certain Z value"""

        from utilities.small_utilities import unTARme
        for var in self.varid2: #cycle through every variant
            var8=self.id+var
            pd=self.vars[var]['prosp'] #handy reference
            if 'Z1' in pd.keys() and pd['Z1']: continue #we already made the assignment
            pd['Z1']=None #no prospector output, thus no value
            if pd['outfile']: #there's a prospector output file
                wd=unTARme(pd['outfile'],fileL=['pdbbpdborienrev',])
                pdbfile=joink(wd,'pdbbpdborienrev',var8+'rap3orienrev5s.pdb')
                pd['Z1']=float(open(pdbfile,'r').readline().strip().split()[2])
                os.system('/bin/rm -r '+wd) #remove temporary working directory

            
    def addSpk(self,switch='',ibatch=0,joblist=''):

        """submit SPICKER simulations"""

        if not joblist:
            from jobs.job import genJobList
            joblist=genJobList()
            
        if switch=='loop info': #submit spicker for loop info                
            for var in self.varid2: #cycle through every variant
                p=self.vars[var]
                p['spk']['loopdat']={}
                p['spk']['outfileloopdat']=''
                if not arethere([p['prosp']['outfile'],p['tasser']['infile'],
                                 p['spk']['infile']]):
                    sys.stderr.write('ERROR some input file not available for '+p['var8']+'\n')
                    continue
                templf='pdbbpdborienrev/'+p['var8']+'rap3orienrev5s.pdb' #PROSPECTOR templf file
                pdbid=p['spk']['Cmodel']['pdbtempl'] #Cmodel most similar to PROSPECTOR template
                comboid=p['spk']['Cmodel']['id']     #cluster id of Cmodel
                outfile=p['var8']+'.out.spicker.tar.bz2'
                outdir=joink(p['simdir'],'loopinfo')
                jobname=p['var8']+'.spk'
                Job=job(exe='spicker',
                        exed=joink(codedir,'f77/spicker/spicker0.0/report_loops/from_prospector'),
                        incmd='tar jxf '+p['prosp']['outfile']+' -C . '+templf+\
                        ' && tar jxf '+p['tasser']['infile']+\
                        ' -C . exp.dat seq.dat rmsinp &&  tar jxf '+p['spk']['infile']+' -C .',
                        args='-a '+templf+' -b '+pdbid+' -c '+`comboid`,
                        outcmd=' tar jcf '+outfile+' combo* closc* loop* rst.dat '+\
                        ' && /bin/mv '+outfile+' '+outdir)
                Job.qsub(jobname,outdir,wallt='0.249',mem_limit='890',ddisk='200',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        elif switch=='non interacting loops':
            for var in self.varid2:
                pt=self.vars[var]  ;  ptnil=pt['nil']  ;  outdir=ptnil['simdir']  ; jobname='spk'+pt['var8']
                inspicker=joink(outdir,pt['var8']+'.in.tasser.tar.bz2')
                if not os.path.exists(inspicker): continue #no input
                outspicker=joink(outdir,pt['var8']+'.out.tasser.tar.bz2')
                if os.path.exists(outspicker): continue #output already done                
                Job=job(exe='spk_rp',exed=joink(pcd,'benchmarkByTMalign'),args='',
                        outcmd='tar jcf %s combo* closc* str.txt rst*.dat && /bin/mv %s %s'%\
                        (outspicker,outspicker,outdir)
                        )
                Job.qsub(jobname,outdir,wallt=0.2499,mem_limit='890',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        return ibatch,True
    
 
    def addSpkInFile(self):
        
        """add the path to the spicker input file, if extant"""

        for var in self.varid2: #cycle through every variant
            f=joink(scd,self.id,var,self.id+var+'.in.spicker.tar.bz2')
            self.vars[var]['spk']['infile']='' #init as not extant
            if isthere(f): self.vars[var]['spk']['infile']=f
        return True
    
                    
    def addSpkOutFile(self):
        
        """add the path to the spicker output file, if extant"""

        for var in self.varid2: #cycle through every variant
            f=joink(scd,self.id,var,self.id+var+'.out.spicker.tar.bz2')
            self.vars[var]['spk']['outfile']='' #init as not extant
            if isthere(f): self.vars[var]['spk']['outfile']=f
       
                    
    def addConReg(self):

        """add constitutive region (empty) dictionary"""

        for var in self.varid2: #cycle through every variant
            self.vars[var]['conreg']={}
 
            
    def addVarReg(self):

        """add variant region (empty) dictionary"""

        for var in self.varid2: #cycle through every variant
            self.vars[var]['varreg']={}


    def addSimdir(self):

        """add home directory for the variant where we keep simulation results"""

        for var in self.varid2: #cycle through every variant
            self.vars[var]['simdir']=joink(scd,self.id,var)

        
    def AddSpkobj(self):

        """add spicker object attribute for each variant"""
        
        from utilities.small_utilities import unTARme
        from spicker.spickerYangResultsManager import spickOut
        for var in self.varid2: #cycle through every variant
            pt1=self.vars[var]
            pt2=pt1['spk']
            pt2['spkobj']=None #assume no spicker results
            if pt2['outfile']:
                unTARme(pt2['outfile'],wd=pt1['simdir']) #untar output file to simdir directory
                pt2['spkobj']=spickOut(dir=pt1['simdir'])
            
        
    def addConregTmMatrix(self):

        """TM-score sup. between var and ref TASSER models for constitutive region

        TM matrix (a dict of dicts) between variant and reference TASSER models.
        Thus, tm[1][3] is detailed TM-score of the constitutive
        region between variant model 1 and reference model 3"""

        from utilities.small_utilities import junkName

        algfile=junkName() #store alignment between variant and reference        
        spkref=self.vars['00']['spk']['spkobj'] #spicker object for reference
        if not spkref: sys.stdout.write('no spicker object for reference '+self.id+'\n')
        for var in self.varid2:
            ptc=self.vars[var]['conreg']
            ptc['tmsup']={}        #assume no matrix
            if var=='00': continue #unnecessary to superimpose the reference to itself
            spkvar=self.vars[var]['spk']['spkobj'] #spicker object for variant
            if not spkvar: sys.stdout.write('no spicker object for variant '+var+'\n')
            #alignment according to Swiss-Prot entry. Upper sequence is reference in the alignment
            alg=self.vars[var]['alg2refSP']
            if not alg: sys.stdout.write('no alg2refSP alignment\n')
            if spkref and spkvar and alg:
                open(algfile,'w').write( alg.info() ) #write alignment to file
                for i in range(1,spkvar.nc+1): #cycle over variant combo files
                    ptc['tmsup'][i]={}
                    cvar=spkvar.combo[i] #combo file for cluster "i" of variant
                    for j in range(1,spkref.nc+1):  #cycle over reference combo files
                        ptx=ptc['tmsup'][i][j]={}
                        cref=spkref.combo[j] #combo file for cluster "j" of reference
                        #print cref+'\n'+cvar
                        cmd='tm_score.py -a '+cref+' -b '+cvar+' -d yes -j identical -i '+algfile
                        lines=os.popen(cmd).readlines() #read TM-score output
                        #print lines                       
                        for item in ('rmsd','tm','cov','rmsd_cov'):
                            match=re.compile(item+'=\s*(\d+\.\d+)').search(lines[0])
                            if match:
                                value=match.group(1)
                                #print value
                                ptx[item]=value
                        #next read alignment from TM-score superposition. Alignment of the form
                        #MIYVstgawgwagSCFWQWSFcppac where uppercase denote aligned residues. Here
                        #we represent the constitutive region as a single continuous chunk
                        tmalignment=''.join( [ l[7:-1] for l in lines[2:] ] )
                        ptx['alignment']=''
                        #print self.id,'var=',var,'i=',i,'j=',j,'tm=',ptx['tm']
                        ptx['alignment']=tmalignment
        os.system('/bin/rm '+algfile+' &>/dev/null') #clean up alignment file


    def correctSpkobj(self):

        """Correct the typographic error spkojb --> spkobj"""

        for var in self.varid2:
            spkref=self.vars[var]['spk']
            spkref['spkobj']=spkref['spkojb'] #create new correct reference
            del spkref['spkojb'] #delete old incorrect reference
        

    def addTmsupdd(self,cc=False):

        """TM-score superposition between densest clusters for reference and variant

        This is a reference to the field c.vars[var]['conreg'][i][j]['tmsup'] where
        'i' is ID of densest cluster for variant and 'j' is ID of densest cluster
        for reference"""


        if not self.vars['00']['spk']['outfile']: return
        spkref=self.vars['00']['spk']['spkobj']
        dref=''
        if spkref:
            if cc: dref=spkref.rankIDsByCscores()[0] #ID of cluster with higest C-score for ref
            else:  dref=spkref.densest #ID of densest cluster for reference
        for var in self.varid2: #cycle through every variant
            ptc=self.vars[var]['conreg']
            if  cc: ptc['tmsupcc']={}
            else:   ptc['tmsupdd']={}
            if not self.vars[var]['spk']['outfile']: continue
            spkvar=self.vars[var]['spk']['spkobj']
            if spkref and spkvar and var != '00':                
                if cc:
                    dvar=spkvar.rankIDsByCscores()[0]
                    ptc['tmsupcc']=ptc['tmsup'][dvar][dref]
                else:
                    dvar=spkvar.densest ##ID of densest cluster for variant
                    ptc['tmsupdd']=ptc['tmsup'][dvar][dref] #ref to TM-score superposition info


    def addTmsupcc(self):

        """TM-score superposition between densest clusters for reference and variant

        same as addTmsupdd(), but use clusters with highest C-score
        instead of densest clusters"""

        self.addTmsupdd(cc=True)
        
        
    def addProspTemplf(self):

        """For every variant, add filename for prospector template file

        addProspTemplf()

        Will create self.vars[var]['prosp']['templf'] containing absolute file
        name for PROSPECTOR template file"""

        from utilities.small_utilities import unTARme

        for var in self.varid2: #cycle through every variant
            pd=self.vars[var]['prosp']
            pd['templf']=''
            if not pd['outfile']: continue #jump to next variant (no PROSPECTOR output)
            wd=unTARme(pd['outfile'],fileL=['pdbbpdborienrev',]) #untar to some working directory
            var8=self.vars[var]['var8']
            pdbfile=joink(wd,'pdbbpdborienrev',var8+'rap3orienrev5s.pdb')
            templatef=joink(self.vars[var]['simdir'],var8+'rap3orienrev5s.pdb')
            os.system('/bin/mv '+pdbfile+' '+templatef) #move template file to proper place
            pd['templf']=templatef
            os.system('/bin/rm -r '+wd) #remove temporary working directory


    def addSpkProps(self):

        """for every variant, add several properties to the spicker object

        addSpkProps()
        
        We add average Zscores, Cscores, ZscoreOfCluster, CscoresRanks,
        coverageOfCombo, and rankIDsByYingPredictedTM fields to
        self.vars[var]['spk']['spkobj']. We purposedly delete the
        coordinates from the 'templates' field because it would
        contain all coordinates from the templates making the dump
        file too big"""

        for var in self.varid2: #cycle through every variant
            templatef=self.vars[var]['prosp']['templf']
            if not templatef: continue #jump to next variant (no template file)
            spkobj=self.vars[var]['spk']['spkobj'] #a reference
            if not spkobj: continue #jump to next variant (no spicker object)
            if hasattr(spkobj,'coverageSeqOfCombo'): continue #we already processed it
            #initialize fields ZscoreOfCluster, Cscores, and CscoresRanks
            spkobj.rankIDsByCscores(templatef=templatef)
            type=self.vars[var]['tasser']['emh'] #easy,medm,hard evaluated by TASSER
            if not type: continue #jump to next variant 
            spkobj.rankIDsByYingPredictedTM(templatef=templatef,type=type)
            spkobj.calcCoverageSeqOfCombo()
            #remove template coordinates to save space
            for template in spkobj.templates: template['templ']=[]
           

    def addIsconst(self):

        """for every var, add a seq. of 1's and 0's for const. or var residues

        addIisconst()

        add field self.vars[var]['isconst']
        For instance, 0000011111111111000001111 indicates a variant sequence of
        length 25 where the first five residues are variant, the next 11
        residues are constitutive, and so on."""

        for var in self.varid2: #cycle through every variant
            pvar=self.vars[var]
            pvar['isconst']=''
            #alignment between reference and variant according to Swiss-Prot entry file
            alg=pvar['alg2refSP']
            if not alg: continue           
            pvar['isconst']=alg.isAlignedLines(algtype='identical')[2]

       
    def addConstCovTypes(self):

        """for each var, find type of coverage for residues in the constitutive region

        addConstCovTypes()
        
        Given a residue in the variant sequence, we give it a type:
          0: residue is variant
          1: res is constitutive, not covered by a PROSPECTOR templ in both ref and var
          2: res is constitutive, covered by a PROSPECTOR templ only in ref
          3: res is constitutive, covered by a PROSPECTOR templ only in var
          4: res is constitutive, covered by a PROSPECTOR templ in both ref and var

        The point is to use this info in order to answer the following question: What
        are the reasons for two residues (one in the ref, the other in the var) being
        unaligned in the constitutive region? Is it because:
        (i)  corresponding PROSPECTOR templates are different?
        (ii) one or both residues are not covered by a PROSPECTOR template?"""

        spkref=self.vars['00']['spk']['spkobj']
        for var in self.varid2: #cycle through every variant
            pv=self.vars[var] ;print var
            #retrieve constitutive alignment as a list of pairs of ungapped indexes
            alg2refSP=pv['alg2refSP']
            if not alg2refSP:
                sys.stderr.write('ERROR(addConstCovTypes) no alg2refSP\n')
                continue #jump to next variant
            pairs=alg2refSP.listIndexPairsIdentical() #list of pairs
            pv['covtypes']={}
            if var=='00':continue #unnecessary to compare to itself
            isconst=pv['isconst']
            #print 'isconst=\n',isconst
            if not isconst:
                sys.stderr.write('no pv[\'isconst\'] found\n')
                continue
            spkvar=pv['spk']['spkobj']            
            if not spkref or not spkvar:
                sys.stderr.write('no [\'spk\'][\'spkobj\'] for variant or reference\n')
                continue
            for idv in spkvar.combo.keys():
                pv['covtypes'][idv]={}
                #print 'covv=\n',spkvar.coverageSeqOfCombo[idv]
                covv=[ int(x) for x in spkvar.coverageSeqOfCombo[idv] ]
                for idr in spkref.combo.keys():
                    conl=[ int(x) for x in isconst ]
                    covr=[ int(x) for x in spkref.coverageSeqOfCombo[idr] ]
                    #print 'covr=\n',spkref.coverageSeqOfCombo[idr]
                    #constitutive residue 'j' in ref is aligned to const. res. 'i' in var
                    for (j,i) in pairs: conl[i]+=conl[i]*(covr[j]+2*covv[i])
                    pv['covtypes'][idv][idr]=''.join([`x` for x in conl])
                    #print 'pv[covtypes][',idv,'][',idr,']',pc[idv][idr]
        return True
    

    def addHigestCscoreModel(self):

        """add attribute containing info on TASSER model with highest C-score"""

        for var in self.varid2:
            p0=self.vars[var]
            p=p0['spk']
            if not p['outfile']: continue
            p['Cmodel']={}
            spkobj=p['spkobj']
            if not hasattr(spkobj,'Cest'): spkobj.__dict__['Cest']=spkobj.rankIDsByCscores()[0]
            p['Cmodel']['id']=spkobj.Cest #id of cluster with highest C-score
            p['Cmodel']['combof']=spkobj.combo[spkobj.Cest] #combo-file
            p['Cmodel']['combof.reb']=p0['simdir']+'/combo%02d.pdb.rebuilt'%(p['Cmodel']['id'])
            p['Cmodel']['closcf']=spkobj.closc[spkobj.Cest] #closc-file
            p['Cmodel']['C']=spkobj.Cscores[spkobj.Cest]    #C-score
            p['Cmodel']['Z']=spkobj.ZscoreOfCluster[spkobj.Cest]    #average Z-score

    
    def addTmsupGapAlg(self):

        """add ...['tmsup'][i][j]['galignment'] field

        The 'galignment' field is same as the 'alignment' field, but we
        include '-' symbols for each variant residue, so that the total
        length of the alignment is same as the variant length.

        The 'alg012' field is same as the 'galignment' field with the
        substitutions: '-'-->'0', [a-z]-->'1', [A-Z]-->'2'

        Also add fields 'i' and 'j' to ...['tmsup'][i][j], indicating
        the respective model ID's"""

        patt=re.compile('[a-z]')
        spkobjref=self.vars['00']['spk']['spkobj']
        if not spkobjref: return
        
        for var in self.varid2:
            if var=='00': continue
            pv=self.vars[var]
            #add fields 'i' and 'j'
            try:
                spkobjvar=pv['spk']['spkobj']  #;print spkobjvar
                if not spkobjvar: continue
            except: continue
            pt=pv['conreg']['tmsup']
            for i in spkobjvar.ids: #cycle through variant models
                for j in spkobjref.ids: #cycle through reference models
                    pt[i][j]['i']=i
                    pt[i][j]['j']=j
            #add fields 'galignment' and 'alg012'
            try:
                isconst=pv['isconst'] #sequence of 1's (constitutive residue) and 0's (var res)
            except: continue #jump to next variant
            for i in spkobjvar.ids: #cycle through variant models
                for j in spkobjref.ids: #cycle through reference models
                    alignment=pt[i][j]['alignment']
                    galignment=''
                    alg012=''
                    n=0
                    for x in isconst:
                        if x=='0':
                            galignment+='-'
                            alg012+='0'
                        else:                            
                            galignment+=alignment[n]
                            if patt.match(alignment[n]): alg012+='1'
                            else: alg012+='2'
                            n+=1
                    pt[i][j]['galignment']=galignment
                    pt[i][j]['alg012']=alg012
                    if n!=len(alignment): print n,len(alignment)


    def addUnalgTMpercents(self,modelsel='C'):

        """report info on the TM-score-unaligned residues of the constitutive region

        After doing TM-score on the constitutive region between a model
        for the reference and a model for the variant, we find some
        residues are unaligned. Possible reasons are:

        (1) PROSPECTOR templates for ref and var are different for these regions
        
        (2) These residues are not covered by a PROSPECTOR template in
        either the ref, var, or both sequences. One can argue that
        (obviating the predicted contacts), TASSER operated in
        abinitio mode for these residues and provided different
        results

        report the length of the unaligned region, as well as the
        percentage of unaligned residues covered by a PROSPECTOR
        template in the reference, variant, or both sequences

        modelsel: select which models in ref and var to use for
        alignment. 'C' stands for models with highest C-score"""

        for var in self.varid2:
            if var=='00': continue
            try:
                if modelsel=='C':
                    p=self.vars[var]['conreg']['tmsupcc']
                    alignment=p['alignment']
                    alg012=p['alg012'] #TM-score alignment
                    covtypes=self.vars[var]['covtypes'][ p['i'] ][ p['j'] ]#PROSP. coverages
            except: continue
            n=0
            props={'1':0.,'2':0.,'3':0.,'4':0.} #see explanation of 'covtypes' field
            for x in alg012: #cycle over alignment
                if x=='1': #'1' unaligned residue after TM-score superposition
                    props[ covtypes[n] ]+=1
                n+=1
            sum=0
            for key in props.keys(): sum+=props[key]
            if sum:
                for key in props.keys(): props[key]=props[key]/sum
            self.vars[var]['conreg']['unalgTMpercents']=props
        return True


    def superimpose(self,var,file):

        """superimpose variant and reference constitutive portion

        superimpose(self,var,file)

        var: pass the two-digit variant. Otherwise do superposition for all variants.
        file: output file for superposition
        """

        from utilities.small_utilities import junkName

        print var,file
        algfile=junkName() #store alignment between variant and reference
        try:
            alg=self.vars[var]['alg2refSP']
            if not alg: return False
            cref=self.vars['00']['spk']['Cmodel']['combof']
            cvar=self.vars[var]['spk']['Cmodel']['combof']
        except: return False
        open(algfile,'w').write( alg.info() )
        cmd='tm_score.py -a '+cref+' -b '+cvar+' -c '+file+' -d yes -j identical -i '+algfile
        lines=os.popen(cmd).readlines()
        os.system('/bin/rm '+algfile+' &>/dev/null') #clean junk files
        return lines


    def amberminimize(self,joblist):

        """run pulchra+dukka+amber minimization on combo with highest C-score

        We do not store new info on the pickle file"""

        from jobs.job import job
        from utilities.codedir import nodecodedir
        for var in self.varid2:
            p=self.vars[var]
            q=p['spk']
            if not q['outfile']: continue
            inputfile=q['Cmodel']['combof']            
            outputfile='combo%02d.pdb.minimized'%(q['Cmodel']['id'])
            q['Cmodel']['combof.min']=joink(p['simdir'],outputfile) #;Bye(q['Cmodel'])
            Job=job(name='ca2at',exe='ca2at.py',
                    exed=joink(nodecodedir,'python/tertiary'),
                    args='-a '+inputfile+' -b '+outputfile,
                    outcmd='/bin/mv '+outputfile+' '+p['simdir'])
            Job.qsub('min'+p['var8'],p['simdir'],wallt='0.249',mem_limit='256',ddisk='300',
                     joblist=joblist)


    def surface_area(self):

        """Calculate exposed surface area for every unaligned region.

        Using amber-minimized densest centroid (see amberminimize
        method) and 'surface' executable from jackal package"""

        from tertiary.surface import area3 #standard surface areas of residues
        from utilities.small_utilities import junkName
        surfacex='/library/jose/code/bin/jackal/bin/surface'
        junkf=junkName()
        for var in self.varid2:
            ptv=self.vars[var]
            pts=ptv['spk']
            if not pts['outfile']: continue
            ptc=pts['Cmodel'] #;print ptc
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
                nline+=1 ; line=lines[nline]
            while line:
                if 'total area' in line:
                    pta['total']=float(line.split(':')[1])
                    break
                nline+=1 ; line=lines[nline]
        os.system('/bin/rm '+junkf)           


    def addConstL(self):

        """calculate length of constitutive region"""

        for var in self.varid2:
            ptv=self.vars[var] #;print var
            L=0
            for x in ptv['isconst']:
                if x=='1': L+=1
            ptv['conreg']['L']=L #;print L
        return True
    

    def add_seq_chunks(self):

        """one-letter sequence for constitutive and variant regions"""

        for var in self.varid2:
            ptv=self.vars[var] #;print var
            seq=ptv['seq']
            isconst=ptv['isconst']
            seqC='' #init sequence for constitutive region
            seqV='' #init sequence for variant region
            for i in range(ptv['L']):
                if isconst[i]=='1': seqC+=seq[i]
                else: seqV+=seq[i]
            ptv['conreg']['seq']=seqC #;print seqC
            ptv['varreg']['seq']=seqV #;print seqV            
        return False


    def genTMsupCmodelCmodel(self):

        """isoform-reference TM-score superposition output.

        For each isoform, create superposition PDB to its reference, where
        the superposition is determined by the TM-score alignment of the
        constitutive region"""
    
        from utilities.small_utilities import junkName

        algfile=junkName() #store alignment between variant and reference        
        spkref=self.vars['00']['spk']['spkobj'] #spicker object for reference
        if not spkref:
            sys.stdout.write('no spicker object for reference '+self.id+'\n')
            return False
        cref=self.vars['00']['spk']['Cmodel']['combof'] #combo model with highest C-score
        for var in self.varid2:
            ptc=self.vars[var]['conreg']
            if var=='00': continue #unnecessary to superimpose the reference to itself
            spkvar=self.vars[var]['spk']['spkobj'] #spicker object for variant
            if not spkvar:
                sys.stdout.write('no spicker object for variant '+var+'\n')
                continue
            #alignment according to Swiss-Prot entry. Upper sequence is reference in the alignment
            alg=self.vars[var]['alg2refSP']
            if not alg: sys.stdout.write('no alg2refSP alignment\n')
            if spkref and spkvar and alg:
                open(algfile,'w').write( alg.info() ) #write alignment to file
                cvar=self.vars[var]['spk']['Cmodel']['combof'] #combo model with highest C-score
                #supf=cvar+'.sup2ref' #superposition to reference sequence
                supf=cvar+'.sup2refII'
                os.system('tm_score.py -a '+cref+' -b '+cvar+' -j identical'+\
                          ' -i '+algfile+' -d yes -c '+supf)#+ '>/dev/null')
                self.vars[var]['spk']['Cmodel']['sup2ref']=supf
        os.system('/bin/rm '+algfile)
        return False
    
            
    def addCovtypeCmodel(self):

        """add covtype entry in Cmodel dictionary"""

        if not self.vars['00']['spk']['outfile']:
            sys.stderr.write('ERROR(addCovtypeCmodel) no spicker output for ref in '+self.id+'\n')
            return False
        idref=self.vars['00']['spk']['Cmodel']['id']
        for var in self.varid2:
            if var=='00': continue #has no sense for reference sequence
            if not self.vars[var]['spk']['outfile']: continue
            pta=self.vars[var]['spk']['Cmodel']
            idvar=pta['id']
            ptb=self.vars[var]['covtypes'] #;Bye(self.vars[var]['covtypes'])
            pta['covtypes']=ptb[idvar][idref]
        return True


    def addUnalgTMcov(self):

        """average template coverage for unaligned constitutive region

        Given a TM-score alignment between the constitutive region of
        one isoform and its reference, we look at the portion of the
        constitutive region that is unaligned by the TM-score
        procedure. This region has a different structure in isoform
        than in the reference. It could be an artifact of the
        models. Models are Cmodel of the reference and the
        isoform. For this unaligned region, we calculate the average
        of ['Cmodel']['covtypes'], which is an indication of how much
        of this region is covered by a PROSPECTOR template in the
        reference and/or isoform. Coverage ranges from 0 (both
        reference and isoform have no PROSPECTOR template) to 2(both
        have a PROSPECTOR template)"""

        if not self.vars['00']['spk']['outfile']:
            sys.stdout.write('no spicker object for reference '+self.id+'\n')
            return False
        for var in self.varid2:
            if var=='00': continue #has no sense for reference sequence
            pta=self.vars[var]
            ptb=pta['conreg']['tmsupcc']
            if not ptb: continue
            cvar=self.vars[var]['spk']['Cmodel']['covtypes']
            av=0.0
            n=0
            for i in range(pta['L']):
                if ptb['alg012'][i]=='1': #unaligned constitutive region
                    y=int(cvar[i])
                    if y==2 or y==3: av+=1
                    if y==4: av+=2
                    n+=1
            pta['conreg']['unalgTMcov']=''
            if n!=0: #there exists an unaligned region in the TM-score alignment
                pta['conreg']['unalgTMcov']=float(av)/n #;print av
        return True
 

    def addTempl2Cmodel(self):

        """add associated PDB template to Cmodel"""

        for var in self.varid2:
            p0=self.vars[var]['spk']
            if not p0['outfile']: continue
            spkobj=p0['spkobj']
            if not hasattr(spkobj,'Cest'): spkobj.__dict__['Cest']=spkobj.rankIDsByCscores()[0]
            p0['Cmodel']['pdbtempl']=spkobj.templateOfCluster[spkobj.Cest]['libPDB']
        return True


    def addFasta(self):

        """add fasta file property"""

        for var in self.varid2:
            p0=self.vars[var]
            if not 'fastafile' in p0.keys():
                p0['fastafile']=joink(scratchdir,'spliceVariants/input',p0['var8'])
        return True
    
            
    def submitBlast(self,switch='',joblist='',ibatch=0,var='',nalg=9999):

        """submit/check blast simulations"""

        if switch=='check': #add outfile
            vars=self.varid2
            if var: vars=[var,]
            for var in vars:
                p0=self.vars[var]
                if 'blast' not in p0.keys(): p0['blast']={'outfile':'',}
                outfile=joink(p0['simdir'],p0['var8']+'.out.blast.tar.bz2')
                if isthere(outfile): p0['blast']['outfile']=outfile
                else: sys.stderr.write('ERROR no blast output for %s\n'%(self.vars[var]['var8'],))
        elif switch=='submit': #submit blast job
            vars=self.varid2
            if var: vars=[var,]
            for var in vars:
                p0=self.vars[var]
                if 'blast' not in p0.keys(): p0['blast']={'outfile':'',}
                outf=p0['var8']+'.out.blast.tar.bz2'
                #if not p0['blast']['outfile']:
                Job=job(name='blast',
                        exe='blastManager.py',
                        exed=joink(codedir,'python/seq'),
                        inputs=[p0['fastafile'],],
                        outcmd='tar jcf %s %s.blast %s.chk %s.fasta && /bin/mv %s %s'%\
                        (outf,p0['var8'],p0['var8'],p0['var8'],outf,p0['simdir']),
                        args=' -a blastRun -b %s -c "opt=-b %s"'%\
                        ( basename(p0['fastafile']), repr(nalg) )
                        )
                Job.qsub(p0['var8']+'.bl',p0['simdir'],wallt=0.249,mem_limit='1500',
                         ddisk='200',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        return False,ibatch
                

    def submitProspector(self,switch='',joblist='',ibatch=0,var=''):

        """submit/check prospector simulations"""

        if switch=='check input': #add outfile
            vars=self.varid2
            if var: vars=[var,]
            for var in vars:
                p0=self.vars[var]
                if 'prosp' not in p0.keys(): p0['prosp']={'infile':'','outfile':''}
                infile=joink(p0['simdir'],p0['var8']+'.in.prospector.tar.bz2')
                if isthere(infile): p0['prosp']['infile']=infile
                else: sys.stderr.write('ERROR no prospector input for %s\n'%(self.vars[var]['var8']))
        if switch=='check output': #add outfile
            vars=self.varid2
            if var: vars=[var,]
            for var in vars:
                p0=self.vars[var]
                if 'prosp' not in p0.keys(): p0['prosp']={'infile':'','outfile':''}
                outfile=joink(p0['simdir'],p0['var8']+'.out.prospector.tar.bz2')
                if isthere(outfile): p0['prosp']['outfile']=outfile
                else: sys.stderr.write('ERROR no prospector output for %s\n'%(self.vars[var]['var8']))
        elif switch=='submit input': #submit production of input for prospector
            vars=self.varid2
            if var: vars=[var,]
            for var in vars:
                p0=self.vars[var]
                fastaf=p0['fastafile']
                blastout=basename(p0['blast']['outfile'])
                outdir=p0['simdir']
                inprosp=p0['var8']+'.in.prospector.tar.bz2'
                Job=job(name='inprospector',
                        exe='inprosp.py',
                        exed=joink(codedir,'python/combo_jobs/predictor'),
                        inputs=[p0['blast']['outfile'],],
                        incmd='tar jxf '+blastout,
                        args='-a '+fastaf+' -b '+blastout+' -c '+outdir,
                        )
                Job.qsub(p0['var8']+'.ip',p0['simdir'],wallt=0.249,mem_limit='255',
                             ddisk='200',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        elif switch=='submit':
            vars=self.varid2
            if var: vars=[var,]
            for var in vars:
                p0=self.vars[var]
                fastaf=p0['var8']+'.fasta' #hidden within blast output
                blastout=basename(p0['blast']['outfile'])
                inprosp=basename(p0['prosp']['inpfile'])
                outdir=p0['simdir']
                Job=job(name='prospector',
                        exe='prosp.py',
                        exed=joink(codedir,'python/combo_jobs/predictor'),
                        inputs=[p0['blast']['outfile'],p0['prosp']['inpfile']],
                        incmd='tar jxf '+blastout,
                        args='-a %s -b %s -c %s -d 1.0'%(fastaf,inprosp,outdir)
                        )
                Job.qsub(p0['var8']+'.p',p0['simdir'],wallt=0.249,mem_limit='499',
                             ddisk='500',joblist=joblist)
                ibatch=Job.dormant(ibatch)
            
        return False,ibatch
        
##################### END OF CLASS PEPSI ########################

            
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    varset ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'varid2': two-character id list of variants
    |_'varid8': full-name id list of variants
    |_'swissp': Swiss-Prot entry as single string
    |_'vars'  : dictionary of info for all isoforms. Keys are the two-char ids of variants
      |_varid2: dictionary to hold info on a particular variant.
        |_'var8': full-name of variant name
        |_'isoname': isoform name given in the headers of file uniprot_sprot_varsplic.fasta
        |_'seq': string holding single-letter amino-acid sequence
        |_'fastafile': absolute path to fasta file. Header is var8
        |_'L': sequence length
        |_'alg2ref': seq.alignYangManageResults.alignment object to the reference sequence
        |_'alg2refSP': as 'alg2ref', but alignment from VAR_SEQ annotation of SwissProt entry
        |_'isconst': sequence of 1's and 0's for either constitutive or variant residues
        |_'seqid':  sequence alignment of the 'alg2ref' object
        |_'simdir': home directory for the variant where we keep simulation results
        |_'covtypes': matrix of sequences of coverage types...['covtypes'][i][j] indicates
        | |           a sequence between TASSER model 'i' of variant and TASSER model 'j'
        | |           of reference. The sequence is as long as the variant sequence.
        | |           A position in the sequence is occupied by either of '0','1','2','3','4'
        | |           where:
        | |           '0': residue is variant
        | |           '1': residue is constitutive and is not covered by a PROSPECTOR
        | |                template in both variant and reference sequence
        | |           '2': res const., covered by a PROSP templ only for the ref. seq.
        | |           '3': res const., covered by a PROSP templ only for the var. seq.
        | |           '4': res const., covered by a PROSP templates in both var and ref seqs.
        | |_i: a TASSER model for the variant
        |   |_j: a TASSER model for the reference. We store the sequence of coverage types 
        |     
        |_'prosp': dictionary holding info related to PROSPECTOR run
        | |_'outfile': prospector output file (*.tar.bz2)
        | |_'emh': easy, medm, hard according to PROSPECTOR
        | |_'Z1': Z-score of first PROSPECTOR template
        |
        |_'tasser': dictionary holding info related to TASSER run
        | |_'infile': tasser input file (*.tar.bz2)
        | |_'outfile' tasser output file (*.tar.bz2)
        | |_'emh': easy, medm, hard according to TASSER
        |
        |_'spk': dictionary holding info related to SPICKER results
        | |_'outfileloopdat': spicker output file for spicker that does loop analysis
        | |_'infile': spicker input file (output from TASSER)
        | |_'outfile': spicker output file (*.tar.bz2)
        | |_'spkobj': spicker object
        | |_'loopdat':info on spicker on the loops        
        | |_'Cmodel': dictionary containing info of model with highest C-score
        |    |_'id': cluster ID
        |    |_'combof': combo file (full path)
        |    |_'closcf': closc file (full path)
        |    |_'combof.reb': combo file rebuilt by dukka+pulchra (full path)
        |    |_'combof.min': combo file rebuilt by dukka+pulchra+amberminimizat (full path)
        |    |_'sup2ref': superposition file to Cmodel of reference. Superposition given by
        |    |            TM-score alignment of the constitutive region. Only constitutive
        |    |            region is printed in the superposition file
        |    |_'C': C-score    
        |    |_'Z': average Z-score of the Cmodel 
        |    |_'pdbtempl': associated PDB file acting as template'
        |    |_'covtypes': sequence  of coverage types between Cmodel of variant and Cmodel
        |    |             of reference. The sequence is as long as the variant sequence.
        |    |             A position in the sequence is occupied by either of '0','1','2',
        |    |             '3','4'. See 'covtypes' under var for explanation of numbers.
        |    |_'area': dictionary containing info on surface area
        |      |_'total': total surface area
        |      |_'value': list with absolute values of surface area per residue
        |      |_'perctg': as a percentage of total residue area
        |_'nil': dictionary containing info of unaligned regions
        | |_'simdir': directory where we store results from simulation of non interacting unaligned regions
        | |_'tasser': dictionary containing info on the TASSER run
        | | |_'chaindat': chain.dat-like file used for TASSER with non interacting unaligned regions
        | | |_'tm_prospTempl2tassTempl': TM score between previous chaindat and PROSP template for Cmodel
        | | |_'infile': all necessary input for TASSER, except for common files
        | | |_'outfile':
        | |_'spicker': directory containing info on spicker results
        |   |_'outfile':
        |_'conreg': dictionary holding info on the constitutive region
        | |_'L': length of the region
        | |_'seq': one-letter amino acid sequence
        | |_'tmsup':tm detailed matrix (a dict of dicts) between reference and variant TASSER
        | | |       models. Thus, tmsup[1][3] is detailed TM of the constitutive
        | | |       region between variant model 1 and reference model 3
        | | |
        | | |_i: cluster i of variant
        | |   |
        | |   |_j: cluster j of reference.All info is stored as strings, not float.
        | |     |
        | |     |_'tm': tm-score between combo i of variant and combo j of reference
        | |     |_'rmsd':rmsd
        | |     |_'cov': coverage of the tm-score superposition
        | |     |_'rmsd_cov': rmsd of the covered region
        | |     |_'alignment': alignment of the form 'MIYVstgawgwagSCFWQWSFcppacvgc'.
        | |     |              We represent the constit region as a single continuous chunk.
        | |     |_'galignment': same as 'alignment', but include '-' for each variant residue
        | |     |_'alg012': same as 'galignment' with '-'-->'0',[a-z]-->'1',[A-Z]-->'2'
        | |_'tmsupdd': reference to some ['tmpsup'][i][j], with 'i' densest variant cluster
        | |            and 'j' is densest reference cluster.
        | |_'tmsupcc': same as 'tmsupdd' but using C-score instead of density
        | |_'unalgTMcov': a measure of confidence in the TM-score unaligned constit region.
        | |               Ranges from 0 (ref and isof have no PROSPECTOR template) up to
        | |               2 (both ref and isof have PROSPECTOR template)
        | |_'unalgTMpercents': dictionary holding percentages of TM-score-unaligned residues
        |      covered by a PROSPECTOR template in the reference, variant, or both sequences.
        |      We use Cmodels for TM-score alignment. Dictionary keys are '1','2','3','4'.
        |      Look in 'covtypes' for meaning of keys.
        |_'varreg' : dictionary holding info on the variant region
          |_'seq': one-letter amino acid sequence
    """
    print globals()['datatree'].__doc__


def doitall(id):

    """a recount of everything we had done on each varset

    Useful if we screw-up a particular dump file, because we can redo
    it all"""

    c=PEPSI(id=id)
    c.genDumpf(cl.repod)
    tmpl,tmpl2=readVarsetsDat()
    c.update('varid2', tmpl[id])
    c.update('varid8', tmpl2[id])
    c.update('vars',{})
    for varid in c.varid2: c.vars[varid]={} #prepare dictionaries to hold variant-specific info

    methods=('addSeqAndL','addSwissp','addVar8','addAlignmentToRef','AddSeqIDtoRef','reportSeqID','addProspector','addProspOutFile','addEasyMedmHard','addTasser','addTasserInFile','addTasserEasyMedmHard','addTasserOutFile','addZ1','addSpk','addSpkOutFile','addConReg','addVarReg','addSimdir','AddSpkobj','checkRef','addAlignmentToRef_VAR_SEQ','addIsoname','addConregTmMatrix','addTmsupdd','addTmsupcc','addProspTemplf','addSpkProps','addIsconst','addConstCovTypes','addHigestCscoreModel','addTmsupGapAlg','amberminimize','surface_area','add_seq_chunks','addConstL','genTMsupCmodelCmodel','addCovtypeCmodel','addUnalgTMcov','addUnalgTMpercents','addTempl2Cmodel')
    for method in methods:
        print method
        getattr(objSL,method)()
    #finally secure the beast!
    c.pickleDump()

    
def initCollection():

    """Initialize dump pickle files and collection file

    initCollection(id='')
    Initialize set of pickle files and collection file to store info
    on the splice variants. Remove previously stored pickle and
    collection file

    id: initialize only one dump file instead of all dump files and
    the collection file

    ***************************************************
    * WARNING!: this service will remove any existing *
    *           dump files and the collection file.   *
    ***************************************************"""

    print '**************************************************************'
    print '* WARNING!: this service will remove any existing dump files *'
    print '*         and the collection file. Do you want to proceed?   *'
    print '**************************************************************'
    sys.stdout.write('(y/n): ')
    reply=sys.stdin.readline().strip()
    if reply != 'y': return
    tmpl,tmp2=readVarsetsDat()

    remaining=len(canlist)
    for header in canlist:
        c=PEPSI(id=header)
        c.genDumpf(cl.repod)
        c.update('varid2', tmpl[header])
        c.update('varid8', tmpl2[header])
        c.update('vars',{})
        for varid in c.varid2: c.vars[varid]={} #prepare dictionaries to hold variant-specific info
        c.pickleDump()
        print remaining,header
        remaining-=1
    cl.updColletion()


def updateCollection():

    """update the collection file"""
    
    cl.updColletion()
    
    
def addSeqAndL():
    
    """add seq and seq length
    
    add sequence as single string, as well as sequence length, for
    each variant"""
    
    cl.cycleOverList('addSeqAndL')


def addVar8():

    """add var8 field to each variant"""

    cl.cycleOverList('addVar8')


def addIsoname():

     """find the name given to each splice variant in file uniprot_sprot_varsplic.fasta"""

     cl.cycleOverList('addIsoname')
     
    
def listIDs():

    """return a list with all varset id's"""

    list=[]
    for c in cl.iterate(): list.append(c.id)
    return list

        
def addSwissp():

    """connect to Local Swiss-Prot and retrieve entries"""
            
    cl.cycleOverList('addSwissp')


def printSwissp(id):

    """print swissp field of a given varset"""

    print cl.loadCAN(id).getSwissp()
    
    
def addAlignmentToRef_VAR_SEQ(id=''):

    """seq. alig. object. Alignment done with VAR_SEQ from Swiss-Prot

    add sequence alignment of each variant to its reference as an
    alignObj sequence alignment object. The alignment is done from
    the VAR_SEQ annotations of the Swiss-Prot entry for the
    varset."""

    if id:
        c=cl.loadCAN(id)
        if c.addAlignmentToRef_VAR_SEQ()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addAlignmentToRef_VAR_SEQ')


def addAlignmentToRef():
    
    """add sequence alignment of each variant to its reference

    add sequence alignment of each variant to its reference as an
    seq.alignYangManageResults.alignment object"""
    
    cl.cycleOverList('addAlignmentToRef')

    
def AddSeqIDtoRef():

    """check sequence identity of the constitutive parts"""

    cl.cycleOverList('AddSeqIDtoRef')


def reportSeqID():

    """print sequence identity between reference and each var"""
    for c in cl.iterate(): c.reportSeqID()


def addProspector():

    """add prospector dictionary to each variant

    It will hold info related to PROSPECTOR run"""

    cl.cycleOverList('addProspector')
    
    
def addProspOutFile():

    """add the path to the prospector output file, if extant"""

    cl.cycleOverList('addProspOutFile')


def reportMissingProspectorJobs(outf=joink(scd,'fail_prospector.list')):

    """output the list of variants with missing prospector jobs

    reportMissingProspectorJobs(outf=joink(scd,'fail_prospector.list')"""

    fail=[]
    success=[]
    for c in cl.iterate(): #cycle over all CAN objects
        for var in c.varid2:
            if not c.vars[var]['prosp']['outfile']: fail.append(c.id+var)
            else:  success.append(c.id+var)
    open(outf,'w').write( '\n'.join(fail) )
    open(joink(scd,'success_prospector.list'),'w').write( '\n'.join(success) )


def addEasyMedmHard():

    """add PROSPECTOR easy, medm, hard classification

        add an attribute to classify PROSPECTOR difficulty for each
        variant (easy, medm, hard)"""

    cl.cycleOverList('addEasyMedmHard')


def createListEasyCases( listf=joink(scd,'easy.list') ):

    """create a file with list of easy cases

    createListEasyCases() creates file
    $SC/spliceVariants/out2/easy.list listing all 8-character variants
    for which PROSPECTOR assigned easy case """

    easylist=[]
    for c in cl.iterate(): #cycle over all CAN objects
        for var in c.varid2:
            if c.vars[var]['prosp']['emh']=='easy':easylist.append(c.id+var)
    open(listf,'w').write( '\n'.join(easylist) )
            

def addTasser():

    """add tasser dictionary to each variant
    
    It will hold info related to TASSER run"""

    cl.cycleOverList('addTasser')


def addTasserInFile():

    """add the path to the tasser input file, if extant"""

    cl.cycleOverList('addTasserInFile')


def addTasserOutFile():

    """add the path to the tasser output file, if extant"""

    cl.cycleOverList('addTasserOutFile')


def  addTasserEasyMedmHard():
    
    """add TASSER easy, medm, hard classification
    
    add an attribute to classify PROSPECTOR difficulty for each
    variant (easy, medm, hard)"""

    cl.cycleOverList('addTasserEasyMedmHard')


def createListEasyCasesTasser( listf=joink(scd,'easy_tasser.list') ):

    """create a file with list of easy cases according to TASSER

    createListEasyCases(listf=joink(scd,'easy_tasser.list')) creates
    file $SC/spliceVariants/out2/easy_tasser.list listing all
    8-character variants for which TASSER assigned easy case """

    easylist=[]
    for c in cl.iterate(): #cycle over all CAN objects
        for var in c.varid2:
            if c.vars[var]['tasser']['emh']=='easy':easylist.append(c.id+var)
    open(listf,'w').write( '\n'.join(easylist) )


def addZ1():

    """add the path to the tasser output file, if extant"""

    cl.cycleOverList('addZ1')


def createListByZ1cutoff(Zco,listf):

    """list of variants by Z-score cutoff

    createListByZ1cutoff(Zco,listf)
    create a list of variants whose first template has Z-score above
    certain cutoff.
    Required parameters:
      Zco : Z-score cut-off
      listf : output file list """

    Zco=float(Zco)
    zlist=[]
    for c in cl.iterate(): #cycle over all CAN objects
        for var in c.varid2:
            Z=c.vars[var]['prosp']['Z1']
            if Z and Z>Zco: zlist.append(c.vars[var]['var8'])
    open(listf,'w').write( '\n'.join(zlist) )

    
def addSpk(id='',switch=''):

    """add spicker results (empty) dictionary"""

    if id:
        c=cl.loadCAN(id)
        if c.addSpk(switch=switch)!=False: c.pickleDump()
    else:
        from jobs.job import genJobList
        joblist=genJobList() ; ibatch=0 ; remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id
            ibatch,record=c.addSpk(switch=switch,joblist=joblist,ibatch=ibatch)
            if record!=False: c.pickleDump()
            remaining-=1

    
    
def addSpkOutFile():

    """add the path to the spicker output file, if extant"""

    cl.cycleOverList('addSpkOutFile')

    
def addSpkInFile():

    """add the path to the spicker input file, if extant"""

    cl.cycleOverList('addSpkInFile',silent=False)

    
def reportMissingSpickerJobs(outf=joink(scd,'fail_spicker.list')):

    """output the list of variants with missing spicker jobs

    reportMissingSpickerJobs(outf=joink(scd,'fail_spicker.list')"""

    fail=[]
    success=[]
    for c in cl.iterate(): #cycle over all CAN objects
        for var in c.varid2:
            if not c.vars[var]['spk']['outfile']: fail.append(c.id+var)
            else:  success.append(c.id+var)
    open(outf,'w').write( '\n'.join(fail) )
    open(joink(scd,'success_spicker.list'),'w').write( '\n'.join(success) )


def addConReg():

    """add constitutive region (empty) dictionary"""

    cl.cycleOverList('addConReg')

    
def addVarReg():

    """add variant region (empty) dictionary"""

    cl.cycleOverList('addVarReg')


def addSimdir():

    """add home directory for the variant where we keep simulation results"""

    cl.cycleOverList('addSimdir')


def AddSpkobj():

    """add spicker object attribute for each variant"""

    cl.cycleOverList('AddSpkobj')
    

def loadCAN(id):

    """load and print one varset"""

    c=cl.loadCAN(id)
    print c.__dict__
    
    
def removeAtribute():

    """remove some attribute that we mistakenly create"""
    
    cl.cycleOverList('removeAtribute')

def checkRef(id=''):

    """check reference sequence is that one given by the Swiss-Prot entry"""
    
    if id:
        c=cl.loadCAN(id)
        if c.checkRef()!=False: c.pickleDump()
    else:
        cl.cycleOverList('checkRef')

def addConregTmMatrix(id=''):

    """TM-score sup. between var and ref TASSER models for constitutive region"""

    if id:
        c=cl.loadCAN(id)
        if c.addConregTmMatrix()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addConregTmMatrix')


def addTmsupdd(id=''):

    """TM-score superposition between densest clusters for reference and variant

    addTmsupdd(id='')
    
    This is a reference to the field c.vars[var]['conreg'][i][j]['tmsup'] where
    'i' is ID of densest cluster for variant and 'j' is ID of densest cluster
    for reference.
    
    id: optional parameter to work on a particular varset. Otherwise will cycle
    through all varsets"""
    
    if id:
        c=cl.loadCAN(id)
        if c.addTmsupdd()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addTmsupdd')


def addTmsupcc(id=''):

    """TM-score superposition between densest clusters for reference and variant

    addTmsupdd(id='')
    
    This is a reference to the field c.vars[var]['conreg'][i][j]['tmsup'] where
    'i' is ID of densest cluster for variant and 'j' is ID of densest cluster
    for reference.
    
    id: optional parameter to work on a particular varset. Otherwise will cycle
    through all varsets"""
    
    if id:
        c=cl.loadCAN(id)
        if c.addTmsupcc()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addTmsupcc')


def correctSpkobj():

    """Correct the typographic error spkojb --> spkobj"""

    cl.cycleOverList('correctSpkobj')

    
def addProspTemplf(id=''):

    """For every variant, add filename for prospector template file

    addProspTemplf()

    Will create self.vars[var]['prosp']['templf'] containing absolute file
    name for PROSPECTOR template file

    id:  parameter to work on a particular varset. Default cycles
         through all varsets"""

    if id:
        c=cl.loadCAN(id)
        if c.addProspTemplf()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addProspTemplf')



def addSpkProps(id=''):
        
    """for every variant, add several properties to the spicker object
    
    addSpkProps()
       
    We add Cscores, ZscoreOfCluster, CscoresRanks, and rankIDsByYingPredictedTM
    fields to self.vars[var]['spk']['spkobj']

    id:  parameter to work on a particular varset. Default cycles
         through all varsets"""

    if id:
        c=cl.loadCAN(id)
        if c.addSpkProps()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addSpkProps')


def addIsconst(id=''):

    """for every var, add a seq. of 1's and 0's for const. or var residues

    addIisconst(id='')

    add field self.vars[var]['isconst']
    For instance, 0000011111111111000001111 indicates a variant sequence of
    length 25 where the first five residues are variant, the next 11
    residues are constitutive, and so on
    id:  parameter to work on a particular varset. Default cycles
         through all varsets"""

    if id:
        c=cl.loadCAN(id)
        if c.addIsconst()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addIsconst')


def addConstCovTypes(id=''):

    """for each var, find type of coverage for residues in the constitutive region

    addConstCovTypes()
        
    Given a residue in the variant sequence, we give it a type:
      0: residue is variant
      1: res is constitutive, not covered by a PROSPECTOR templ in both ref and var
      2: res is constitutive, covered by a PROSPECTOR templ only in ref
      3: res is constitutive, covered by a PROSPECTOR templ only in var
      4: res is constitutive, covered by a PROSPECTOR templ in both ref and var

    The point is to use this info in order to answer the following question: What
    are the reasons for two residues (one in the ref, the other in the var) being
    unaligned in the constitutive region? Is it because:
    (i)  corresponding PROSPECTOR templates are different?
    (ii) one or both residues are not covered by a PROSPECTOR template?"""

    if id:
        c=cl.loadCAN(id)
        if c.addConstCovTypes()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addConstCovTypes')


def addHigestCscoreModel(id=''):

    """add attribute containing info on TASSER model with highest C-score"""

    if id:
        c=cl.loadCAN(id)
        if c.addHigestCscoreModel()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addHigestCscoreModel')
    

def addTmsupGapAlg(id=''):

    """add ...['tmsup'][i][j]['galignment'] field

    The 'galignment' field is same as the 'alignment' field, but we
    include '-' symbols for each variant residue, so that the total
    length of the alignment is same as the variant length"""

    if id:
        c=cl.loadCAN(id)
        if c.addTmsupGapAlg()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addTmsupGapAlg')


def report(zrefco='',zvarco='',crefco='',cvarco='',typeref='',typevar='',maxtm='',minrmsd=''):

    """report info on reference variant pairs

    report(zrefco='',zvarco='',crefco='',cvarco='',typeref='',typevar='',maxtm='')

    zrefco: do not output references with average Z-score below this cutoff
    zvarco: do not output variants with average Z-score below this cutoff
    crefco: do not output references with C-score below this cutoff
    cvarco: do not output variants with C-score below this cutoff
    typeref: only output reference with TASSER type of this type (easy,medm,hard)
    typevar: only output variant with TASSER type of this type (easy,medm,hard)
    maxtm: do not output if TM-score above this cutoff
    """

    if zrefco: zrefco=float(zrefco)
    if zvarco: zvarco=float(zvarco)
    if crefco: crefco=float(crefco)
    if cvarco: cvarco=float(cvarco)
    if maxtm: maxtm=float(maxtm)
    if minrmsd: minrmsd=float(minrmsd)
    
    buf='#varset-id, var-id, average Z-score-reference, average Z-score-variant,'
    buf+=' C-score-reference, C-score-variant,'
    if typeref: buf+=' TASSER-type-reference,'
    if typevar: buf+=' TASSER-type-variant,'
    buf+=' TM-score-constitutive-region'
    print buf
    
    for c in cl.iterate(): #cycle over all CAN objects
        pv=c.vars['00']
        try:
            zref=pv['spk']['Cmodel']['Z']
            if zrefco and zref<zrefco: continue
            cref=pv['spk']['Cmodel']['C']
            if crefco and cref<crefco: continue
            idref=pv['spk']['Cmodel']['id']
            if typeref:
                typer=pv['tasser']['emh']
                if typeref!=typer: continue
        except: continue
        for var in c.varid2:
            if var=='00': continue
            pv=c.vars[var]
            try:
                tm=pv['conreg']['tmsupcc']['tm']
                rmsd=pv['conreg']['tmsupcc']['rmsd']
                if maxtm and float(tm)>maxtm: continue
                if minrmsd and float(rmsd)<minrmsd: continue
            except: continue
            try:
                zvar=pv['spk']['Cmodel']['Z']
                if zvarco and zvar<zvarco: continue
                cvar=pv['spk']['Cmodel']['C']
                if cvarco and cvar<cvarco: continue
                idvar=pv['spk']['Cmodel']['id']
                if typevar:
                    typev=pv['tasser']['emh']
                    if typevar!=typev: continue
            except: continue
            
            buf='%s %s %5.1lf %5.1lf %5.2lf %5.2lf'%(c.id,var,zref,zvar,cref,cvar)
            if typeref: buf+=' '+typer
            if typevar: buf+=' '+typev                
            buf+=' '+tm+' '+rmsd
            print buf
        

def superimpose(id,var,file):

        """superimpose variant and reference constitutive portion

        superimpose(id,file,var='')

        id: varset id
        var: pass the two-digit variant. Otherwise do superposition for all variants.
        file: output file for superposition"""

        print id,var,file
        cl.loadCAN(id).superimpose(var,file)
        

def amberminimize(id=''):

    """run pulchra+dukka+amber minimization on combo with highest C-score"""

    from jobs.job import genJobList
    joblist=genJobList()
    if id:
        c=cl.loadCAN(id)
        if c.amberminimize(joblist)!=False: c.pickleDump()
    else:
        for c in cl.iterate():
            if c.amberminimize(joblist)!=False: c.pickleDump() #cycle over all CAN objects
            

def surface_area(id=''):

    """Calculate exposed surface area for every unaligned region.

    Using amber-minimized densest centroid (see amberminimize
    method) and 'surface' executable from jackal package"""

    if id:
        c=cl.loadCAN(id)
        if c.surface_area()!=False: c.pickleDump()
    else:
        cl.cycleOverList('surface_area',silent=True)


def addConstL(id=''):

    """calculate length of constitutive region"""
    
    if id:
        c=cl.loadCAN(id)
        if c.addConstL()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addConstL',silent=True)


def add_seq_chunks(id=''):

    """insert one-letter sequence for constitutive and variant regions"""

    if id:
        c=cl.loadCAN(id)
        if c.add_seq_chunks()!=False: c.pickleDump()
    else:
        cl.cycleOverList('add_seq_chunks',silent=True)


def genTMsupCmodelCmodel(id=''):

    """isoform-reference TM-score superposition output.

    For each isoform, create superposition PDB to its reference, where
    the superposition is determined by the TM-score alignment of the
    constitutive region"""

    if id:
        c=cl.loadCAN(id)
        if c.genTMsupCmodelCmodel()!=False: c.pickleDump()
    else:
        cl.cycleOverList('genTMsupCmodelCmodel',silent=False)

    
def addCovtypeCmodel(id=''):

    """add covtype entry in Cmodel dictionary"""

    if id:
        c=cl.loadCAN(id)
        if c.addCovtypeCmodel()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addCovtypeCmodel',silent=True)


def addUnalgTMcov(id=''):

    """average template coverage for unaligned constitutive region

    Given a TM-score alignment between the constitutive region of one
    isoform and its reference, we look at the portion of the
    constitutive region that is unaligned by the TM-score
    procedure. This region has a different structure in isoform than
    in the reference. It could be an artifact of the models. Models
    are Cmodel of the reference and the isoform. For this unaligned
    region, we calculate the average of ['Cmodel']['covtypes'], which
    is an indication of how much of this region is covered by a
    PROSPECTOR template in the reference and/or isoform. Coverage
    ranges from 0 (both reference and isoform have no PROSPECTOR
    template) to 2 (both have a PROSPECTOR template)"""

    if id:
        c=cl.loadCAN(id)
        if c.addUnalgTMcov()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addUnalgTMcov',silent=True)


def addUnalgTMpercents(id=''):

    """report info on the TM-score-unaligned residues of the const. region

    After doing TM-score on the constitutive region between a model
    for the reference and a model for the variant, we find some
    residues are unaligned. Possible reasons are:

    (1) PROSPECTOR templates for ref and var are different for these regions
    
    (2) These residues are not covered by a PROSPECTOR template in
        either the ref, var, or both sequences. One can argue that
        (obviating the predicted contacts), TASSER operated in
        abinitio mode for these residues and provided different
        results
        
    report the length of the unaligned region, as well as the
    percentage of unaligned residues not covered by a PROSPECTOR
    template in the reference, variant, or both sequences"""

    if id:
        c=cl.loadCAN(id)
        if c.addUnalgTMpercents()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addUnalgTMpercents',silent=True)


def addTempl2Cmodel(id=''):

    """add associated PDB template to Cmodel"""

    if id:
        c=cl.loadCAN(id)
        if c.addTempl2Cmodel()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addTempl2Cmodel',silent=False)


def printRefVarTypeProspTMconst(id='',typeref='',typevar='',cref=None,cvar=None,tmco=''):

    """print ref var type-prosp TMconst-region

    printRefVarTypeProspTMconst(typeref='',typevar='',cref='',cvar='',tm=1.0)
    
    Ouput reference variant (easy,memd,hard)according to PROSPECTOR, and TM-score between
    constitutive region for corresponding densest clusters
    
    typeref:  will report cases with such type ('easy' or 'memd' or 'hard')
              for the reference according to PROSPECTOR labelling.
              Default is accept all cases
    typevar:  same as 'typeref' but for variant
    cref: only report cases with densest reference cluster having C-score above cref
    cvar: same as 'cref' but for the densest variant cluster
    tmco=1.0: will report cases with TM-score BELOW the passed cut-off"""

    #parse input
    if cref: cref=float(cref)
    if cvar: cvar=float(cvar)
    if not tmco: tmco=1.0
    else: tmco=float(tmco)

    print '#(1)id  (2)var  (3)Lref  (4)Lvar (5)PDB template for reference'
    
    print '#(6) PDB template for variant (7) average Z-score reference template'
    print '#(8) average Z-score of variant template (9) Z-score of ref templ as a template'
    print '#for variant sequence (10) Z-score of var templ as templ for ref seq.'
    print '#(11) type-ref  (12) C-score-ref  (13) cluster-ID-ref  (14) type-var'
    print '#(15) C-score-var  (16) cluster-ID-ref  (17) length-constitutive-region'
    print '#(18) TM-score-Cmodels (19) RMSD Cmodels after TM-score superposition'
    print '#(20 to 23) percentage of residues in the constitutive'
    print '#region that are unaligned after TM-score superposition of the const region'
    print '#and that are: (20)not covered by PROSPECTOR template, (21)only covered'
    print '#in the reference sequence,  (22)only covered in the isoform, (23)covered'
    print '#in both reference and isoform  (24)superposition PDB file showing'
    print '#superposition of the constitutive region'
    print '#(1)   (2)(3) (4)   (5)  (6)    (7)   (8)   (9)  (10) (11)   (12)(13)(14)   (15)(16)(17) (18)  (19)(20) (21) (22) (23)           (24)'
    idlist=cl.canlist
    if id: idlist=[id,]
    for id in idlist:
        c=cl.loadCAN(id)
        #check we have all info for reference
        pref=c.vars['00']
        #check prospector output and type
        if not pref['prosp']['outfile']: continue #move to next varset
        typeR=pref['prosp']['emh'] #PROSPECTOR difficulty type
        if typeref and typeR!=typeref: continue #move to next varset
        #check spicker output and C-score
        if not pref['spk']['outfile']: continue #move to next varset
        spkref=pref['spk']['spkobj']
        cR=spkref.Cscores[ spkref.densest ] #C-score of densest reference cluster
        idR=pref['spk']['Cmodel']['id']
        pdbR=pref['spk']['Cmodel']['pdbtempl'] #associated template
        Zref=pref['spk']['Cmodel']['Z']
        if isinstance(cref,float) and cref>cR:  continue #move to next varset
        #cycle through variants
        for var in c.varid2:
            if var=='00': continue #avoid comparing reference to itself
            pvar=c.vars[var]
            #check prospector output and type
            if not pvar['prosp']['outfile']: continue #move to next variant
            typeV=pvar['prosp']['emh'] #PROSPECTOR difficulty type
            if typevar and typeV!=typevar: continue #move to next variant
            #check spicker output and C-score
            if not pvar['spk']['outfile']: continue #move to next variant
            spkvar=pvar['spk']['spkobj']
            cV=spkvar.Cscores[ spkvar.densest ] #C-score of densest variant cluster
            if isinstance(cvar,float) and cvar>cV: continue #move to next varset
            #info on TM-score alignment for constitutive region between most dense clusters
            if not pvar['conreg']: continue #move to next variant
            tmsupcc=pvar['conreg']['tmsupcc']
            if not tmsupcc: continue
            tm=tmsupcc['tm'] #TM-score between the two most dense clusters
            rmsd=float(tmsupcc['rmsd'])
            tmf=float(tm)
            if tmf>tmco: continue #too similar, move to next variant           
            L=pvar['conreg']['L']
            idV=pvar['spk']['Cmodel']['id']
            pdbV=pvar['spk']['Cmodel']['pdbtempl'] #associated template
            Zvar=pvar['spk']['Cmodel']['Z']
            Zrefvar=-99.9 #average Z-score of reference template as  template for variant sequence
            for templ in spkvar.templates: #cycle over templates for variant
                if templ['libPDB']==pdbR: #the template is actually the selected templ for ref.
                    Zrefvar=templ['Z']
                    break
            Zvarref=-99.9
            for templ in spkref.templates:
                if templ['libPDB']==pdbV:
                    Zvarref=templ['Z']
                    break
            buf=c.id+' '+var
            buf+=' %3d %3d'%(pref['L'],pvar['L'])
            buf+=' %s %s'%(pdbR,pdbV)
            buf+=' %5.1lf %5.1lf %5.1lf %5.1lf'%(Zref,Zvar,Zrefvar,Zvarref)
            buf+=' %s %6.2lf %2d'%(typeR,cR,idR)
            buf+=' %s %6.2lf %2d'%(typeV,cV,idV)
            buf+=' %3d %s %4.1lf'%(L,tm,rmsd)
            pun=pvar['conreg']['unalgTMpercents']
            for x in ('1','2','3','4'): buf+=' %4.2lf'%(pun[x],)
            buf+=' %s'%'/'.join(pvar['spk']['Cmodel']['sup2ref'].split('/')[-3:])
            print buf


def print_number_fragments(id='',minsize=1):

    """for every isoform, print number of const and variant fragments

    print_number_fragments(id='',minsize=1)
    For instance, one variant insertion in the middle of the sequence
    would yield two constitutive fragments (the termini) and one
    variant region. For instance, one variant deletion would yield one
    constitutive fragment and zero variant fragments"""

    listC=[]
    if id: #only one element
        listC.append(cl.loadCAN(id))
    else:
        for c in cl.iterate():
            listC.append(c)
    
    for c in listC:
        for var in c.varid2:
            ptv=c.vars[var] #;print var
            ncon=0 #number of constitutive fragments
            nvar=0   #number of variant fragments
            isconst=ptv['isconst'] ;print isconst
            x=isconst[0]
            if x=='1': ncon+=1
            else: nvar+=1        
            for y in isconst[1:]:
                if y!=x: #we change identity of fragment
                    if y=='1': ncon+=1
                    else: nvar+=1
                    x=y #overwrite with new fragment identity
            print ptv['var8'],'%2d %2d'%(ncon,nvar)

    
def print_alg2refSP(id=''):

    """print alignment between reference and variants """

    listC=[]
    if id: #only one element
        listC.append(cl.loadCAN(id))
    else:
        for c in cl.iterate():
            listC.append(c)

    for c in listC:
        for var in c.varid2:
            ptv=c.vars[var] #;print var
            print ptv['alg2refSP']


def printSomeStuff(id):

    """small function to be edited for custom printing of some property"""

    c=cl.loadCAN(id)
    for var in c.varid2:
        ptv=c.vars[var] 
        if ptv['spk']['spkobj']:
            print '###################\n'+var
            print ptv['spk']['spkobj'].templates
    

def addFasta(id=''):

    """add fasta file property"""

    if id: #only one element
        c=cl.loadCAN(id)
        if c.addFasta()!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining, c.id
            if c.addFasta()!=False: c.pickleDump()
            remaining-=1
    return True


def view(id):
    
    """view object properties

    view(id)"""

    cl.loadCAN(id).view()
        

def submitBlast(id='',switch='',listvar8f='',nalg=9999):

    """submit blast simulations and check appropriate CAN attribute"""

    if id: #only one element
        c=cl.loadCAN(id)
        if c.submitBlast(switch=switch,nalg=nalg)!=False: c.pickleDump()
    else:
        from jobs.job import genJobList
        joblist=genJobList() ; remaining=len(cl.canlist) ; ibatch=0
        if listvar8f:
            listL=chomp(open(listvar8f,'r').readlines())
            remaining=len(listL)
            for h in listL:
                id=h[0:6] ; var=h[6:]
                c=cl.loadCAN(id)
                print remaining,id
                record,ibatch=c.submitBlast(switch=switch,ibatch=ibatch,var=var,nalg=nalg)
                if record!=False: c.pickleDump()
                remaining-=1
        else:
            for c in cl.iterate():
                print remaining, c.id
                record,ibatch=c.submitBlast(switch=switch,ibatch=ibatch,nalg=nalg)
                if record!=False: c.pickleDump()
                remaining-=1


def submitProspector(id='',switch='',listvar8f=''):

    """submit blast simulations and check appropriate CAN attribute"""

    if id: #only one element
        c=cl.loadCAN(id)
        if c.submitProspector(switch=switch)!=False: c.pickleDump()
    else:
        from jobs.job import genJobList
        joblist=genJobList() ; remaining=len(cl.canlist) ; ibatch=0
        if listvar8f:
            listL=chomp(open(listvar8f,'r').readlines())
            remaining=len(listL)
            for h in listL:
                id=h[0:6] ; var=h[6:]
                c=cl.loadCAN(id)
                print remaining,id
                record,ibatch=c.submitProspector(switch=switch,ibatch=ibatch,var=var)
                if record!=False: c.pickleDump()
                remaining-=1
        else:
            for c in cl.iterate():
                print remaining, c.id
                record,ibatch=c.submitProspector(switch=switch,ibatch=ibatch)
                if record!=False: c.pickleDump()
                remaining-=1


def tassSpicSim(id='',switch='',idlistf=''):
    """submit TASSER & spicker simulations"""
    from jobs.job import genJobList
    joblist=genJobList()  ;  remaining=len(cl.canlist)  ;  ibatch=0  ;
    if id: idlist=[id,]
    if idlistf: idlist=chomp(open(idlistf,'r').readlines())
    if not idlist:  idlist=cl.canlist
    remaining=len(idlist)
    for id in idlist:
        print remaining,id
        c=cl.loadCAN(id)
        record,ibatch=c.tassSpicSim(switch=switch,ibatch=ibatch,joblist=joblist)
        remaining-=1


def help(*kargs):

    """list available services"""
    
    services=['datatree','doitall','updateCollection','initCollection','addSeqAndL',
              'addAlignmentToRef','AddSeqIDtoRef','addProspector','addProspOutFile',
              'addEasyMedmHard','createListEasyCases','reportMissingProspectorJobs',
              'addTasser','addTasserInFile','addTasserEasyMedmHard',
              'createListEasyCasesTasser','addTasserOutFile','removeAtribute','addZ1',
              'createListByZ1cutoff','addSpk','addSpkOutFile','reportMissingSpickerJobs',
              'addVarReg','addConReg','addSimdir','AddSpkobj','addSwissp','listIDs',
              'addAlignmentToRef_VAR_SEQ','loadCAN','addIsoname','checkRef',
              'addConregTmMatrix','addTmsupdd','addProspTemplf','addSpkProps',
              'addIsconst','addConstCovTypes','addTmsupcc','addHigestCscoreModel',
              'report','addTmsupGapAlg','addUnalgTMpercents','superimpose','amberminimize',
              'surface_area','add_seq_chunks','addConstL','genTMsupCmodelCmodel',
              'addCovtypeCmodel','addUnalgTMcov','addTempl2Cmodel',
              'printSwissp','printRefVarTypeProspTMconst','print_number_fragments',
              'printSomeStuff','view','addSpkInFile','addFasta','tassSpicSim']
    os.system('/usr/bin/clear')
    if not kargs:
        print 'Services='
        for service in services:
            print '  '+service+' : '+globals()[service].__doc__.split('\n')[0]
        print 'Type "workCANs.py -a help -b service" for more detail on "service"'
    else:
        sv=kargs[0]
        print globals()[sv].__doc__        


if __name__=='__main__':
    inpHand('Usage: workCANs.py [options]',
            ' -a _R_service required service. Type "help" for list of services.',              
            ' -b __servargs service required arguments in python syntax, and separate them with ",". Escape symbols \' and " (def:None)',
            ' -c __servoptargs service optional arguments (Example: "key=val" for one argument and "key1=val1, key2=val2,..." for two or more arguments (def:None)',
            ).parse(locals(),sys.argv)

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

