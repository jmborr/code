#!/usr/bin/python

import os,re,sys,numpy
from utilities.codedir import codedir
from utilities.small_utilities import junkName,Bye

joink=os.path.join
basename=os.path.basename
dirname=os.path.dirname

class leap:

    """dealing with LEap"""
    
    def __init__(self):

        self.AAA=('ALA', 'ARG', 'ASN', 'ASP', 'ASH',\
                  'CYS', 'CYM', 'CYX', 'GLU', 'GLH',\
                  'GLN', 'GLY', 'HID', 'HIE', 'HIP',\
                  'ILE', 'LEU', 'LYS', 'MET', 'PHE',\
                  'PRO', 'SER', 'THR', 'TRP', 'TYR',\
                  'VAL')

        self.basenameAAA={
            'ALA':'ALA', 'ARG':'ARG', 'ASN':'ASN', 'ASP':'ASP',\
            'ASH':'ASP', 'CYS':'CYS', 'CYM':'CYS', 'CYX':'CYS',\
            'GLU':'GLU', 'GLH':'GLU', 'GLN':'GLN', 'GLY':'GLY',\
            'HID':'HIS', 'HIE':'HIS', 'HIP':'HIS', 'ILE':'ILE',\
            'LEU':'LEU', 'LYS':'LYS', 'MET':'MET', 'PHE':'PHE',\
            'PRO':'PRO', 'SER':'SER', 'THR':'THR', 'TRP':'TRP',\
            'TYR':'TYR', 'VAL':'VAL'}
        self.isstandarAAAname={
            'ALA':True, 'ARG':True, 'ASN':True, 'ASP':True, \
            'ASH':False,'CYS':True, 'CYM':False,'CYX':False,\
            'GLU':True ,'GLH':False,'GLN':True, 'GLY':True, \
            'HID':False,'HIE':False,'HIP':False,'ILE':True, \
            'LEU':True, 'LYS':True, 'MET':True, 'PHE':True, \
            'PRO':True, 'SER':True, 'THR':True, 'TRP':True, \
            'TYR':True, 'VAL':True}

class tleap(leap):

    """derived class from LEaP"""

    def __init__(self,cmd=[],seq=''):
        self.cmd=cmd #list of commands
        if seq: self.cmd+=seq.split('\n')
        
    def go(self,script=None):
        """execute the stored commands"""
        if not script: script=junkName()
        open(script,'w').write('\n'.join(self.cmd)+'\nquit')
        cmd='tleap -f  '+script
        os.system(cmd)

class top:
    """manage topology file"""

    def __init__(self,topf):
        self.top=topf

    def numberAtoms(self,mask='*'):
        """number of Atoms with a particular mask"""
        N=-1
        buf='strip '+mask+'\ngo'
        scriptf=junkName()
        open(scriptf,'w').write(buf)
        outf=junkName()
        cmd='ptraj '+self.top+' < '+scriptf+' > '+outf+' 2>/dev/null'
        os.system(cmd)
        buf=''.join(open(outf,'r').readlines())
        os.system('/bin/rm '+scriptf+' '+outf) #cleanup
        match=re.compile('Mask.*represents (\d+) atoms').search(buf)
        if match: N=int(match.group(1))
        return N
        
class crd:
    """manage output trajectory coordinate files"""
    def __init__(self,crdf,nat):
        """index the frames in the trajectory"""
        self.crdf=crdf
        self.nat=nat
        self.nlpf=nat/10 #number of lines per frame
        if nat%10: nlpf+=1
        self.nframe=0
        self.idx=[] #frame indexing of crd
        ptf=open(crdf,'r')
        while ptf:
            for iline in range(self.nlpf): ptf.readline()
            self.nframe+=1
            self.idx.append( ptf.tell() )
            self.frames=numpy.array([None]*(nframes))
        self.ptf=None
        self.currframe=0
        
    def opencrd(self,startframe=None):
        """open trajectory file"""
        self.ptf=open(self.crdf,'r')
        if not startframe: startframe=0
        self.currframe=startframe         #current frame index
        self.curr=self.idx[startframe]    #current file-frame index
        ptf.seek(self.curr)               #update to startframe
        
    def loadframe(self,iframe=None):
        """load frame to memory"""
        if iframe:
            self.currframe=iframe
            self.curr=self.idx[iframe]
            self.ptf.seek(self.curr)
        frame=[]
        for i in range(self.nlpf):
            frame+=self.ptf.readline().split()
        return numpy.array(frame).reshape(nat,3)
        
    def loadframes(self,startframe=None,endframe=None):
        """load to memory frames"""
        if not startframe: startframe=self.currframe
        if not endframe: endframe=self.nframe
        self.opencrd(startframe=startframe)
        iframe=startframe
        while iframe<=endframe:
            self.frames[iframe]=self.loadframe()
            iframe+=1
        
class sander:
    """dealing with sander input and output"""

    def __init__():
        self.exe=None  #path to executable
        self.inp=None  #path to input file


class ptraj:
    """dealing with ptraj"""

    def __init__(self,top,crd,vel=None):
        """crd can be one string or a list

        top   (str) topoly file
        cmd   (list)
        crds  (list)
        Ncs   (list, number of frames in every coordinate file)
        Nc    (total number of frames)
        vel
        outcrd output trajectory file
        """
        self.cmd=[]   #list of commands
        self.crds=[]  #list of coordinate files
        self.Ncs=[] 
        self.Nc=0   
        self.top=top
        if type(crd)==str:
            self.crds=[crd,]
        elif type(crd)==list:
            self.crds=crd
        self.crd=self.crds[0]
        if vel: self.vel=vel
        
    def iscommitted(self,command):
        """check if command keyword is already in self.cmd"""
        buf='\n'+'\n'.join(self.cmd)
        if '\n'+command in buf:
            sys.stderr.write(command+' already committed\n')
            True
        return False
        
    def trajin(self,start=None,end=None):
        #check first we did not do trajin first
        if self.iscommitted('trajin'): return {'self':self}
        tmp=[]
        for crd in self.crds:
            buf='trajin '+crd
            if start: buf+=' '+`start`
            if end: buf+=' '+`end`
            tmp.append(buf)
        self.cmd=tmp+self.cmd #trajin goes to the beginning
        return {'self':self}

    def trajout(self,outcrd):
        if self.iscommitted('trajout'): return {'self':self}
        self.cmd.append('trajout '+outcrd)
        return {'self':self}

    def strip(self,entities=[':WAT',]):
        """remove entities list"""
        for entity in entities: self.cmd.append('strip '+entity)
        return {'self':self,}

    def reference(self,reffile):
        """set reference structure"""
        self.cmd.append('reference '+reffile)
        return {'self':self}

    def rms(self,mode='first',outrms=''):
        self.cmd.append('rms '+mode+' out '+outrms)
        return {'self':self}

    def correl(self):
        """add correlation matrix"""
        self.cmd.append('matrix correl name correl out '+\
                        'correl.dat byres mass')
        return {'self':self}
    
    def quasiharmonic(self,nvecs=1000):
        """output quasiharmonic analysis"""
        list=['matrix mwcovar name mwcovar out mwcovar.dat byatom',
              'analyze matrix mwcovar out evecs.dat thermo vecs '+`nvecs`]
        self.cmd+=list
        return {'self':self}

    def go(self,script=None,outfile=None,outcrd=None):
        """execute the stored commands
        outfile: output from ptraj
        """
        rmscript=False
        rmoutfile=False
        if not script:
            script=junkName()
            rmscript=True
        if not outfile:
            outfile=junkName()
            rmoutfile=True
        self.trajin()
        if outcrd and not self.iscommitted('trajout'):
            self.outcrd=outcrd
            self.trajout(outcrd)
        open(script,'w').write('\n'.join(self.cmd))
        cmd='ptraj '+self.top+' < '+script+' &> '+outfile
        os.system(cmd)
        output='\n'.join(open(outfile,'r').readlines())
        if rmscript: os.system('/bin/rm '+script)
        if rmoutfile: os.system('/bin/rm '+outfile)
        return {'self':self,'output':output}

    def saveSlice(self,start,end,out):
        """save a slice of the first coordinate file"""
        traj=ptraj(self.top,self.crds[0]).trajin(start,end).trajout(out)
        
    def count_frames(self,crd=None):
        if not crd: crd=self.crds
        if type(crd).__name__!='list': crd=[crd,] #make it a list
        self.Ncs=[]
        for x in crd:
            output=''.join(ptraj(self.top,x).go()['output'])
            pattern=re.compile('Successfully read in (\d+) sets')
            match=pattern.search(output)
            self.Ncs.append( int(match.group(1)) )
        self.Nc=sum(self.Ncs)
        return {'self':self,'N':self.Nc}            
        
    def create_rst(self,rst,offset=0,vel=None):
        """create restart file from crd and vel files"""
        if not vel: vel=self.vel
        Nx=self.count_frames()['N']         #;print 'Nx='+`Nx`
        Nv=self.count_frames(crd=vel)['N']  #;print 'Nv='+`Nv`
        N=min(Nx,Nv)-offset                 #;print 'N ='+`N`
        if N<=0: return {'self':self,'N':0} #empty files or some error
        #do a restart of coordinates
        restart_crd=junkName()
        cmd='trajin '+self.crd+' '+`N`+' '+`N`+'\n'+\
             'trajout '+restart_crd+' restart'
        #print cmd
        script=junkName()
        open(script,'w').write(cmd)
        cmd='ptraj '+self.top+' < '+script+' &>/dev/null '+\
             ' && /bin/rm '+script
        os.system(cmd)
        restart_crd+='.'+`N`   #annoying twist by ptraj
        #do a restart of velocities
        restart_vel=junkName()
        cmd='trajin '+vel+' '+`N`+' '+`N`+'\n'+\
             'trajout '+restart_vel+' restart\n'
        #print cmd
        script=junkName()
        open(script,'w').write(cmd)
        cmd='ptraj '+self.top+' < '+script+' &>/dev/null '+\
             ' && /bin/rm '+script
        os.system(cmd)
        restart_vel+='.'+`N`
        #merge both restart files
        box=os.popen('tail -1 '+restart_crd).readline()  #;print 'box=',box
        buf=''
        pt=open(restart_crd,'r') ; line=pt.readline()
        while '   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000' not in line:
            buf+=line ; line=pt.readline()
        pt=open(restart_vel,'r')
        pt.readline() ; pt.readline() ; line=pt.readline() ;
        while '   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000' not in line:
            buf+=line ; line=pt.readline()
        buf+=box
        open(rst,'w').write(buf)
        os.system('/bin/rm '+restart_crd+' '+restart_vel)
        return {'self':self,'Nc':Nx,'Nv':Nv,'N':N}


class evecs:
    """modes as given by amber format, output of analyze matrix"""

    def __init__(self,file):
        self.file=file
        
    def create_movies(self,topfile,nmodes=100,outdir=None):
        """create PDB movies for each mode in evecs.dat

        nmodes=100 create movies only for the first 100 modes
        outdir=? directory to put the PDB movies"""
        currdir=os.getcwd()
        if not outdir: outdir=currdir
        workdir=joink(currdir,junkName())  ;  print workdir
        arvindfiles=joink(codedir,'software/amber9/arvind_movies/*')
        os.system('mkdir -p '+workdir)
        cmd= 'cd '+workdir+' && /bin/cp '+arvindfiles+' . && '+\
             'ln -s '+topfile+' topology ; '+\
             './vecs_to_binpos1.0 < '+self.file+' > binpos 2> junkX ; '+\
             'tail -1 junkX|tr -s " "|awk \'{print $6}\' > bytes'
        os.system(cmd)
        bytes=open(joink(workdir,'bytes'),'r').readline().strip()
        cmd= 'cd '+workdir+' && split -a 3 -b '+bytes+' binpos && '+\
             './doit '+`nmodes`+' '+outdir+' && '+\
             'cd '+currdir+' && '+\
             '/bin/rm -rf '+workdir
        os.system(cmd)
        
    
class mode:
    """mode, normal or quasi-harmonic"""
    pass

class normal_mode(mode):
    """normal mode, for harmonic systems"""
    pass
    
class quasiharmonic_mode(mode):
    """quasiharmonic mode"""
    pass

