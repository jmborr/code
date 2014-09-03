#!/usr/bin/python

import os,sys,re
from utilities.codedir import codedir,nrfiltdb,yanglib
from utilities.small_utilities import abort,Bye
from inputArgs.inputArgs import inpHand,deglobb,update_opts
from seq.fastaManager import retrieveHeader,fastaTMPfile
from seq.alignYangManageResults import gappedSeq,alignment
from inputArgs.inputArgs import inpHand


blastpgp=yanglib+'/blast/blastpgp'
#blastpgp=codedir+'/bin/blast/blastpgp'

opts={}
opts['blastRun']='-b 9999 -j 3 -h 0.001 -m 0 -d '+nrfiltdb+' -Q xxxxx.pssm -C xxxxx.chk -o xxxxx.blast'

res={}
res['Qname']=re.compile('Query=\s*(\S+)\n')
res['Ql']=re.compile('\((\d+)\sletters')
res['Sname']=re.compile('>(\S+)\n')
res['Sl']=re.compile('Length\s=\s*(\d+)')
res['Score']=re.compile('\sScore\s=\s*(\d+)')
res['Expect']=re.compile('Expect\s=\s*([\w,\.,-]+)')
res['Identities']=re.compile('Identities\s=\s*(\d+)')
res['Positives']=re.compile('Positives\s=\s*(\d+)')
res['Gaps']=re.compile('Gaps\s=\s*(\d+)')
res['length']=re.compile('Identities\s=\s*\d+/(\d+)')

#if the local alignment does not include all residues from the Query
#sequence, add those missing residues, and add corresponding gaps to
#the Sbjct sequence
def adjustToQuerySeq(querySeq,gs1,gs2,begin,end):
    '''querySeq: string containing whole Query sequence
    gs1:string containing local (gapped) alignment of Query sequence
    gs2:string containing local (gapped) alignment of Sbjct sequence
    begin, end: initial and ending indexes of querySeq string where alignment begins (C-array convention)
    '''
    Ql=len(querySeq) #;print querySeq
    if begin==0 and end==Ql-1: return gs1,gs2 #nothing to do
    if begin>0:
        gs1=querySeq[:begin]+gs1 #include leading residues not explicitly shown in the alignment
        gs2='-'*begin+gs2
    if end<Ql-1:
        gs1+=querySeq[end+1:]
        gs2+='-'*(Ql-1-end)
    return gs1,gs2

def readblastAlignment(pin):
    all=''
    
    #go to beginning of the alignment
    line=pin.readline()  ;  match=res['Score'].match(line)
    while not match and line:
        line=pin.readline()  ;  match=res['Score'].match(line)
        
    #add lines to string "all" until we reach two consecutive '\n'
    all=line
    stop=0
    while stop<2 and line:
        line=pin.readline()
        if line=='\n': stop+=1
        else: stop=0
        all+=line

    if not all: return all #we did not find any alignment
    position=pin.tell()
    
    #find the query sequence name and length, prepend appropriate lines to "all"
    pin.seek(0) #go to beginning of the file
    all2=''
    while all2.find('Query=')<0: all2=pin.readline()
    all2+=pin.readline() #record the "Query=..." line, as well as the next line

    #record the entry sequence name and length, by reading backwards. Prepend appropriate lines to "all"
    pin.seek(position)
    character=''
    while character!='>':
        pin.seek(-2,1)        #go backwards two bytes from current position
        character=pin.read(1) #read only one character (thus we'll advance by one position)
    pin.seek(-1,1) ; all2+=pin.readline()+pin.readline() #read name and length

    pin.seek(position)
    return all2+all

#A typical alignment entry will have the lines:
#
# Score =  404 bits (1037), Expect = e-111
# Identities = 228/228 (100%), Positives = 228/228 (100%), Gaps = 0/228 (0%)
#
#Query: 1   MTKKILIIEDEKNLARFVSLELQHEGYEVIVEVNGREGLETALEKEFDLILLDLMLPEMD 60
#           MTKKILIIEDEKNLARFVSLELQHEGYEVIVEVNGREGLETALEKEFDLILLDLMLPEMD
#Sbjct: 1   MTKKILIIEDEKNLARFVSLELQHEGYEVIVEVNGREGLETALEKEFDLILLDLMLPEMD 60
#
#Query: 61  GFEVTRRLQTEKTTYIMMMTARDSIMDVVAGLDRGADDYIVKPFAIEELLARIRATFRRQ 120
#           GFEVTRRLQTEKTTYIMMMTARDSIMDVVAGLDRGADDYIVKPFAIEELLARIRATFRRQ
#Sbjct: 61  GFEVTRRLQTEKTTYIMMMTARDSIMDVVAGLDRGADDYIVKPFAIEELLARIRATFRRQ 120
#
#
class blastAlignment:
        
    def __init__(self,pin,querySeq=''):
        """
        pin: file-handle
        """
        self.error=''
        #property list
        self.pl=['all','Qname','Ql','Sname','Sl','Score','Expect','Identities','Positives','Gaps','length','similarity','E','I','P','G']
        self.pd={}   #property dictionary
        self.gs={}   #the two gapped sequences
        self.Yalg='' #aligment object to manage alignments with the format from Yang's align.f
        
        for p in self.pl: self.pd[p]='' #initialize property dictionary

        self.pd['all']=readblastAlignment(pin) #dump alignment into a string

        all=self.pd['all']
        if not all:
            self.error='ERROR from blastAlignment::__init__ reading blast alignment'
            return

        for property in ['Qname','Sname']:
            self.pd[property]=''
            match=res[property].search(all)
            if match: self.pd[property]=match.group(1)
            
        for property in ['Ql','Sl','Score','Identities','Positives','Gaps','length']:
            self.pd[property]=0             #;print 'property=',property
            match=res[property].search(all)
            if match: self.pd[property]=int(match.group(1))

        self.pd['Expect']=res['Expect'].search(all).group(1)
        if self.pd['Expect'][0]=='e': self.pd['Expect']='1.0'+self.pd['Expect']
        self.pd['E']=float(self.pd['Expect'])

        for property in ['Identities','Positives','Gaps']:
            self.pd[ property[0] ]=(1.0*self.pd[property])/self.pd['length']

        #alignment object to manage output alignment from Yang's align.f object
        self.ix={'Query':1,'Sbjct':2} #translation of indexes in Yalg
        gs1=''.join( re.compile('Query:\s*\d+\s*([A-Z,-]+)').findall(all) )
        gs2=''.join( re.compile('Sbjct:\s*\d+\s*([A-Z,-]+)').findall(all) )
        #we may want to extend the alignment to include the full Query sequence
        if querySeq:
            #find beginnnig and ending amino acids indexes for ungapped
            #Query sequence where the alignment begins and ends
            indexes=re.compile('Query:\s*(\d+)\s*[A-Z,-]+\s*(\d+)').findall(all)
            begin=int(indexes[0][0])-1 ; end=int(indexes[-1][1])-1 #-1 because we follow C-array convention
            gs1,gs2=adjustToQuerySeq(querySeq,gs1,gs2,begin,end)
        self.Yalg=alignment(gappedSeq(gs1),gappedSeq(gs2))

        #handy REFERENCES to the gapped sequences objects
        for seq in ['Query','Sbjct']:
            self.gs[seq]=self.Yalg.gs[ self.ix[seq] ]


    def info(self):
        if not self.pd['all']:
            str+=' Score =%5d bits (000), Expect = %s\n'%(self.pd['Score'],self.pd['Expect'])
            str+=' Identities = %d/%d'%(self.pd['Identities'],self.pd['length'])
            str+=' (%d%),'%(int( (100.0*self.pd['Identities'])/self.pd['length'] ))
            str+=' Positives = %d/%d'%(self.pd['Positives'],self.pd['length'])
            str+='(%d%),'%(int( (100.0*self.pd['Positives'])/self.pd['length'] ))
            str+=' Gaps = %d/%d'%(self.pd['Gaps'],self.pd['length'])
            str+='(%d%),'%(int( (100.0*self.pd['Gaps'])/self.pd['length'] ))
            str+='\n'
            gs1=self.Yalg.gs[1].gs  ;  gs1=self.Yalg.gs[2].gs  ;  sim=self.pd['similarity']
            while(gs1):
                str+='\n'
                str+='Query:      '+gs1[0:61]+'    \n'  ;  gs1=gs1[61:]
                str+='            '+sim[0:61]+'\n'      ;  sim=gs1[61:]
                str+='Sbjct:      '+gs2[0:61]+'    \n'  ;  gs2=gs2[61:]
            self.pd['all']=str
        return self.pd['all']
    
    def __str__(self):
        return self.info()


################################################### end of class blastAlignment


def blastOutFile2Object(blastOutf,outf):
    """
    save the contents of a blast run, which are stored in a file with name blastOutf into
    object outf (file,file-handle,list)
    """
    if isinstance(outf,str):
        if outf!=blastOutf: os.system('/bin/cp '+blastOutf+' '+outf)
    elif isinstance(outf,list): outf+=open(blastOutf,'r').readlines()
    elif isinstance(outf,file): outf.writelines( open(blastOutf,'r').readlines() )
    
def blastRun(fastaf,outf=None,header='',opt=''):
    """
    do blast run
    
    input description:
    (1)fastaf: fasta entry (file,file-handle,list)
    (2)outf: store blast output (file,file-handle,list)
    (3)header: prefix to append to output files. If not passed, it will be the first word
               of the fasta entry header
    (4)opt: change default options to blast
    blast output:
    by default, files header.blast and header.chk are generated. if outf is passed, the contents of
    header.will also be stored in outf
    """
    
    if opt: opt=update_opts(opt,opts['blastRun'])
    else: opt=opts['blastRun']
    if not header:header=retrieveHeader(fastaf)
    opt=deglobb(opt,header)

    fastaFile=fastaTMPfile(fastaf,filename=header+'.fasta') #temp file if entry is not a file
    cmd=blastpgp+' -i '+fastaFile+' '+opt  #;print cmd;sys.exit(1)
    blastpgp_error=os.system(cmd) #run the blast
    if blastpgp_error!=0:
        raise OSError, 'ERROR blastpgp job failed'
    #save also contents of header+'.blast in object outf
    if outf: blastOutFile2Object(header+'.blast',outf)

    #return outf object and names of generated files
    return {'outf':outf, 'blast':header+'.blast', 'chk':header+'.chk', 'fasta':fastaFile, 'pssm':header+'.pssm','errors':{'blastpgp':blastpgp_error} }

#parse psi-blast output to gather all pairwise alignments.
def parse_blast(blastf,Eco,maxId,minID,red,querySeq='',m=0):
    '''Eco: maximum E-value
     maxID: maximum number of identical residues, divided by parent length
     minID: minimum  ---   --    ---       ---      ---   --  ---    ---
     red: redundancy, after we project all sequences onto the query
          sequence, no two sequences can have more than \"redundancy\"          
          sequence identity. This sequence identity is the number of
          identical residues divided by the average of the two
          (ungapped) sequences
     querySeq: we may want to extend the local alignments with the
          full Query sequence
     m: alignment view, same as flag -m from blastpgp
    '''
     
    al=[] #list of all projections of alignments onto query sequence

    #find sequence length
    line=os.popen('grep -A 1 -e "Query=" '+blastf+' | grep letters').readline()
    m=re.compile('(\d+)').search(line)
    if m: lch=int(m.group(1)) #length of query sequence
    else: abort('ERROR: parsing blast output file, could not find sequence length')

    #go to last round
    fpt=open(blastf,'r') ; line=fpt.readline()
    r=re.compile('\s*Results from round 3|CONVERGED\!')
    while not r.match(line):
        line=fpt.readline() #scan until match
        if not line:
            sys.stderr.write('ERROR seq.blastManager.parse_blast: blast file '+blastf+' does not contain line "Results from round 3" or line "CONVERGED"\n')
            return []
    #read all blast alignments
    nalg=0 #number of read alignments
    while(1):
        alg=blastAlignment(fpt,querySeq=querySeq) #;Bye(alg)#read one blast alignment
        if alg.error: break #reading error or just end-of-file reached
        nalg+=1
        seqId=(1.0*alg.pd['Identities'])/alg.pd['Ql'] #note we divide by ungapped Query sequence length
        if alg.pd['E']<Eco and  seqId<maxId and seqId>minID:
            #project Sbjct onto Query
            gs0=alg.Yalg.project( alg.ix['Sbjct'], alg.ix['Query'] ) #;print gs0.gs
            #check that the projection is not redundant to other previous projections
            nonred=1
            for gs in al:
                alg2=alignment(gs,gs0) #create alignment object with both gapped sequence objects
                #sequence identity is number of identical residues over average of two ungapped lengths
                seqId=(2.0*alg2.identicalL())/(gs.sL+gs0.sL)
                if seqId>red:
                    nonred=0
                    break
            if nonred: al.append(gs0)
    #transform all gapped sequence objects to strings
    for i in range(0,len(al)): al[i]=al[i].gs
    if not al:
        sys.stderr.write('blastManager.parse_blast found '+`nalg`+' local alignments but all were filtered out\n')
    return al


def help(*kargs):

    """list available services"""
    
    services=['blastRun',]
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

    inpHand('Usage: blastManager.py [options]',
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
    exitcode=0
    try:
        locals()[service](*reqargs,**optsargs)
    except Exception, e:
        sys.stderr.write(repr(e.args[0])+'\n')
        exitcode=1
    sys.exit(exitcode)

        
