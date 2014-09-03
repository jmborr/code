#!/usr/bin/python
import os,sys
from utilities.small_utilities import chomp
from seq.letters import three2one,one2three
from seq.alignYangManageResults import gappedSeq,alignment

class templ_file:
    """
    store info on memory of a PROSPECTOR output file containing the templates
    templ_file, class name
     |_:'parser', PROSPECTOR version that specifies the templates file format
     |_:'file', file name
     |_:'error', store error messages
     |_:'tpls', dictionary of templates
         |_header, template header
            |_'nagl', number of aligned residues
            |_'Z', Z-score
            |_'M'
            |_'xyz', list of coordinates
            |_'alg', seq.alignYangManageResults object
            |_'seqid', sequence identity over aligned region
    """
    
    def __init__(self,file,parser='genomesmay06'):
        self.error=''
        self.parser=parser
        self.file=''
        if not os.path.exists(file):
            self.error+='file '+file+' does not exists\n'
        else:
            self.file=file


    def loadTemplate(self,pin):
        t={}
        if self.parser=='genomesmay06': #prospector version genomesmay06
            while(1): #read header file
                line=pin.readline()   #;print line;sys.exit(1)
                if not line: return {} #reached end-of-file
                if line.find('ATOM')<0 and line.find('TER') and len(line)>1: #header line
                    t['header'],t['nalg'],t['Z'],t['M']=line.split()
                    t['nalg']=int(t['nalg']); t['Z']=float(t['Z'])
                    break
            t['xyz']=[]
            t['ATOM']=[]
            aligned_pairs=[]
            gp1=['X',]*10000 #initialize sequence 1 as a bunch of unknown amino acids
            gp2=['X',]*10000
            p1=0 ;m1=0 ;p2=0; m2=0; dx=[] ; nalg=0            
            while(1): #read coordinates and sequences
                line=pin.readline()
                if not line or line.find('ATOM')<0: break #we reached end-of-file
                t['ATOM'].append(line.strip())
                nalg+=1 #number of aligned residues
                m1=int(line[20:26]) ; m2=int(line[54:59])
                aligned_pairs.append( (m1,m2) )
                dx.append( m1-p1 - (m2-p2) ) #difference in number of 'X'
                p1=m1 ; p2=m2
                #note m1>0 and m2>0 allways, thus gp1[0] and gp2[0] are 'X' allways
                gp1[m1]=three2one[line[17:20]] #substitute unknown for amino acid
                gp2[m2]=three2one[line[60:63]]
                t['xyz'].append(line[29:54]) #coordinates
            #print nalg,t['nalg'] ;sys.exit(1)
            if nalg != t['nalg']: return {}  #some error reading coordinates
            del gp1[m1+1:] ; del gp2[m2+1:]  #clip tails of unknown amino acids
            #print dx ;print ''.join(gp1)+'\n'+''.join(gp2);sys.exit(1)
            #insert gaps to correctly align gp1 and gp2
            i1=0 ; i2=0 ; malg=0; prev=[0,0]
            hp1='' ;hp2=''
            for pair in aligned_pairs:
                nX=pair[0]-prev[0]-1
                hp1+='X'*nX
                hp2+='-'*nX
                nX=pair[1]-prev[1]-1
                hp1+='-'*nX
                hp2+='X'*nX
                hp1+=gp1[ pair[0] ]
                hp2+=gp2[ pair[1] ]
                prev=pair                
            #print ''.join(hp1)+'\n'+''.join(hp2);sys.exit(1)
            #store as alignment object
            t['alg']=alignment( gappedSeq( ''.join(hp1)), gappedSeq(''.join(hp2) ))
            #print t['alg'];sys.exit(1)
            if line.find('TER')>=0:
                junk,t['seqID']=line.split() ; t['seqID']=float(t['seqID'])
        return t

        
    def loadTemplates(self):
        if 'tpls' not in self.__dict__: self.tpls={}
        pin=open(self.file,'r') #open list of templates for reading
        while True:
            tmpl=self.loadTemplate(pin) #read another template
            if not tmpl: break            
            self.tpls[ tmpl['header'] ]=tmpl
     
        
    def templ2file(self,header='',fname=''):

        """write a template to file"""

        from copy import deepcopy

        if not header and not n:
            sys.stderr.write('ERROR templ2file: provide either header or n\n')
        template=None
        if header:
            for h in self.tpls.keys():
                if h==header:
                    template=self.tpls[h]
                    break
        if not fname:
            from utilities.small_utilities import junkName
            fname=junkName()
        open(fname,'w').write( '\n'.join(template['ATOM'])+'\nTER\n' )
        return fname

            
#######################################################

    
class prospOut:

    def __init__(self,tarf='',dir='.',parser='genomesmay06'):

        """

        tarf : prospector output is tarred instead of inside 'dir'
        """
        from utilities.small_utilities import unTARme
        self.parser=parser
        self.fs={} #list of output files from prospector
        self.tpls=[] #list of templates, every template is a dictionary of properties
        if tarf: self.dir=unTARme(tarf) #untar to temporary directory
        if dir: self.dir=dir
        self.filesystem(parser,dir)
        self.loadTemplates()
        if tarf: os.system('/bin/rm -r '+dir) #remove temporary directory

    def loadTemplate(self,pin):
        parser=self.parser
        t={}
        if parser=='genomesmay06': #prospector version genomesmay06
            while True: #read header file
                line=pin.readline()   #;print line;sys.exit(1)
                if not line: return {} #reached end-of-file
                if line.find('ATOM')<0 and line.find('TER') and len(line)>1: #header line
                    t['header'],t['nalg'],t['Z'],t['M']=line.split()
                    t['nalg']=int(t['nalg']); t['Z']=float(t['Z'])
                    break
            t['xyz']=[]
            gp1=['X',]*10000 #initialize sequence 1 as a bunch of unknown amino acids
            gp2=['X',]*10000
            p1=0 ;m1=0 ;p2=0; m2=0; dx=[] ; nalg=0            
            while True: #read coordinates and sequences
                line=pin.readline()
                if not line or line.find('ATOM')<0: break #we reached end-of-file
                nalg+=1 #number of aligned residues
                m1=int(line[20:26]) ; m2=int(line[54:59])
                dx.append( m1-p1 - (m2-p2) ) #difference in number of 'X'
                p1=m1 ; p2=m2
                #note m1>0 and m2>0 allways, thus gp1[0] and gp2[0] are 'X' allways
                gp1[m1]=three2one[line[17:20]] #substitute unknown for amino acid
                gp2[m2]=three2one[line[60:63]]
                t['xyz'].append(line[29:54]) #coordinates
            #print nalg,t['nalg'] ;sys.exit(1)
            if nalg != t['nalg']: return {}  #some error reading coordinates
            del gp1[m1+1:] ; del gp2[m2+1:]  #clip tails of unknown amino acids
            #print dx ;print ''.join(gp1)+'\n'+''.join(gp2);sys.exit(1)
            #insert gaps to correctly align gp1 and gp2
            i1=0 ; i2=0 ; malg=0
            while malg<nalg:
                while gp1[i1]=='X': i1+=1 #goto first known amino acid in sequence 1
                while gp2[i2]=='X': i2+=1
                #print malg,i1,i2,dx[malg]
                if dx[malg]!=0: #there are some gaps to insert
                    if dx[malg] > 0: #insert dx[i] gaps in the second sequence
                        for i in range(0,dx[malg]): gp2.insert(i2,'-')
                        i2+=dx[malg]
                    else:
                        for i in range(0,-dx[malg]): gp1.insert(i1,'-')
                        i1-=dx[malg]
                #print malg,i1,i2,dx[malg]
                #print ''.join(gp1)+'\n'+''.join(gp2);sys.exit(1)
                malg+=1 ; i1+=1 ; i2+=1
            del gp1[0] ; del gp2[1] #clip first 'X' (allways this issue with C-arrays numbering)
            #print ''.join(gp1)+'\n'+''.join(gp2);sys.exit(1)
            #store as alignment object
            t['alg']=alignment( gappedSeq( ''.join(gp1)), gappedSeq(''.join(gp2) ))
            #print t['alg'];sys.exit(1)
            if line.find('TER')>=0:
                junk,t['seqID']=line.split() ; t['seqID']=float(t['seqID'])
        return t
        
    def loadTemplates(self):
        parser=self.parser
        pin=open(self.fs['templf'],'r') #open list of templates for reading
        while(1):
            tmpl=self.loadTemplate(pin,parser) #read another template
            if not tmpl: break            
            self.tpls.append(tmpl)

    def templ2file(self,header='',n=0,fname=''):

        """write a template to file"""

        from copy import deepcopy

        if not header and not n:
            sys.stderr.write('ERROR templ2file: provide either header or n\n')
        template=None
        if header:
            for tpl in self.tpls:
                if header==template['header']:
                    template=tpl
                    break
        elif n: template=self.tpls[n]

        if not fname:
            from utilities.small_utilities import junkName
            fname=junkName()
        open(fname,'w').write( '\n'.join(template['ATOM'])+'\nTER\n' )
        return fname
            
        
    def filesystem(self,parser,dir=''):
        if not dir: dir=self.dir
        if parser=='genomesmay06': #prospector version genomesmay06
            cmd='basename `ls '+dir+'/pdbbpdborienrev/*rap3orienrev5s.pdb`'
            self.fs['templf']=dir+'/pdbbpdborienrev/'+chomp(os.popen(cmd).readline()) #templ file
            cmd='basename `ls '+dir+'/pdbbzpot4a/*.predictedrap3orienrev`'
            self.fs['contf']=dir+'/pdbbzpot4a//'+chomp(os.popen(cmd).readline()) #contact file
            
