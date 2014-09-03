#!/usr/bin/python

import os,re,sys
from inputArgs.inputArgs import inpHand,deglobb,update_opts
from utilities.small_utilities import chomp,junkName
from utilities.codedir import codedir,nrfiltdb

#some global variables
fasta=codedir+'/bin/fasta/fasta' #fasta executable
fastaout=codedir+'/f77/seq/fastaout.x' #Adrian's executable to parse fasta results
opts={}
opts['fastaRun']='-b 500 -E 0.001 -Q -m 10 -p -w 50' #default options
opts['fastaParse']='0.95 0.10 0.35 0.90' #for long and close alignments
db0=nrfiltdb #default nr-filter database

ph=re.compile('>\s*(\S+)\s*') #extract the header out of the fasta header line
paa=re.compile('[G,A,V,L,I,S,T,C,M,P,D,N,E,Q,K,R,H,F,Y,W,B,Z,X]') #check for amino acid letter
res={}
res['getHeader']=re.compile('>\s*(\w+)')

def fastaTMPfile(fastaf,filename=''):
    """
    create temporary file storing the fastaf entry. Name ir random unless \"filename\" is passed
    """
    if not filename: filename=junkName()
    if isinstance(fastaf,str):
        if filename!=fastaf: os.system('/bin/cp '+fastaf+' '+filename)
    else:
        if isinstance(fastaf,file):
            header,sequence=importFastaEntry(fastaf)
            writeFasta(header,sequence,filename=filename)
        elif isinstance(fastaf,list):
            writeFasta(fastaf[0],fastaf[1],filename=filename)
    return filename
    
#get header as the first word in the header line of the fasta entry
def retrieveHeader(fastaf):
    line=''
    if isinstance(fastaf,str): line=open(fastaf,'r').readline()
    elif isinstance(fastaf,file):
        pos=int(fastaf.tell()) #record current position
        line=fastaf.readline()
        while not res['getHeader'].match(line): line=fastaf.readline()
        fastaf.seek(pos)       #return to previous position
    elif isinstance(fastaf,list): line='>'+fastaf[0]
    return res['getHeader'].match(line).group(1)
    
#out can be output filename, a list where to store output, or a file handle
def fastaRun(fastaf,out,db=db0,opt=''):
    """
    fastaRun(fastaf,out,db=db0,opt='')
    fastaf: input fasta file name
    out: output, either a filename, a file handle, or a list    db: sequence database
    opt: options to pass to the fasta executable
    """
    if opt: opt=update_opts(opt,opts['fastaRun']) #change default options
    else: opt=opts['fastaRun']
    
    if isinstance(out,str):
        cmd=fasta+' '+fastaf+' '+db+' -O '+out+' '+opt+' >/dev/null' ; os.system(cmd)
    else:
        junk=junkName()
        cmd=fasta+' '+fastaf+' '+db+' -O '+junk+' '+opt ; os.system(cmd)
        if isinstance(out,list):
            out=open(junk).readlines()
        elif isinstance(out,file): out.writelines(open(junk).readlines())
        os.system('/bin/rm '+junk)

#parse fasta results. Generate header.close.aln, header.dist.aln, and fastaout.log files.
def fastaParse(fastaf,fastaOut,header='',opt=''):
    """
    fastaf: fasta file (or list)
    fastaOut: output file from fasta program (or list)
    header: prefix to output files
    opt: changes to default fastaParse options
    """
    if opt: opt=update_opts(opt,opts['fastaParse']) #change default options
    else: opt=opts['fastaParse']

    #determine the format of the fasta file
    junk1=''
    if isinstance(fastaf,str): junk1=fastaf
    elif isinstance(fastaf,list): #we passed a list, thus write to a file
        junk1=junkName()  ;  open(junk1,'w').writelines(fastaf)

    #determine the format of the output from running fasta on the fasta file
    junk2=''
    if isinstance(fastaOut,str): junk2=fastaOut
    elif isinstance(fastaOut,list):
        junk2=junkName()  ;  open(junk2,'w').writelines(fastaOut)    

    #determine header, the prefix for the output files, and run the parser
    if not header: header=retrieveHeader(fastaf)
    else:
        cmd=fastaout+' '+junk1+' '+junk2+' '+header+' '+opt  ;  os.system(cmd) #;print cmd
    
    #clean-up
    if isinstance(fastaf,list):   os.system('/bin/rm '+junk1)
    if isinstance(fastaOut,list): os.system('/bin/rm '+junk2)
    
#return sequence and header of current fasta entry where
def importFastaEntry(pinX,paa=paa):
    seq='' ;  header=''
    pin=pinX
    if isinstance(pinX,str):
        try:
            pin=open(pinX,'r')
        except:
            sys.stderr.write('ERROR (importFastaEntry) could not open '+pinX+'\n')
            return [header,seq]
    #check end-of-file nad beginning of fasta entry
    while not header:
        header=pin.readline()
        if not header: return ['',''] #check end-of-file
        header=ph.match(header).group(1) #check beginning of fasta entry
        #print 'header=',header
    #read fasta sequence
    n=pin.tell()  ;  line=chomp( pin.readline() )
    while paa.match(line): #correspond to an amino acid letter
        seq+=line
        n=pin.tell()
        line=chomp( pin.readline() )

    if line: pin.seek(n) #go back one line, since we over-read in the previous loop
    return [header,seq]

def writeFastaToString(header,seq,n=60):
    buf='>'+header
    tmp=seq
    while tmp:
        buf+='\n'+tmp[0:n]
        tmp=tmp[n:]
    return buf
    
#write fasta file to current working directory
def writeFasta(header,seq,filename='',n=60):
    tmp=seq
    if not filename: filename=header
    if filename:
        if isinstance(filename,file): pout=filename
        else: pout=open(filename,'w')
    pout.write('>'+header+'\n')
    while tmp:
        pout.write(tmp[0:n]+'\n')
        tmp=tmp[n:]
    pout.close()

def alnGen(inpf,outd,listH=[],optsfasta='',optsparse=''):
    """
    input description:
    (1) inpf (globbed) file containing one or more fasta entries
    (2) outd (globbed) directory where to put output files
    (3) listH list of headers, to be used if we pass a globbed inpf
    (4) optsfasta override default fastaRun options
    (5) optsparse override default fastaParse options
    create following output:
    (1) outd/fastaout/header.out file containing output from fasta
    (2) outd/3590/header.aln close multiple sequence alignment
    (3) outd/e10/header.aln  distant multiple sequence aligment
    (4) outd/fastaout/fastaout.log log file
    """
    junk=junkName() #;print 'junk',junk #temporary file to write fasta entry to hard drive
    if not listH: listH=['xxxxx',]
    for h in listH:
        inpf2=deglobb(inpf,h) ; outd2=deglobb(outd,h) ;print inpf2 #deglobbing
        pin=open(inpf2,'r') #file with one or more fasta files
        while(1): #process all fasta entries in inpf2
            header,seq=importFastaEntry(pin)
            if not header: break
            header=retrieveHeader(inpf2) #retain only first "word" as header
            writeFasta(header,seq,filename=junk) #write to temporary file
            out=header+'.out' #out file to store fasta results
            fastaRun(junk,out,opt=optsfasta) #passing an empty "optsfasta" option is OK
            fastaParse(junk,out,header=header,opt=optsparse) #parse fasta results
            closeD=outd2+'/3590' ; distD=outd2+'/e10' ; fastaOutD=outd2+'/fastaout' #subdirs
            os.system('/bin/mkdir -p '+closeD+' '+distD+' '+fastaOutD)            
            os.system('/bin/mv '+header+'.close.aln '+closeD+'/'+header+'.aln')
            os.system('/bin/mv '+header+'.dist.aln ' +distD +'/'+header+'.aln')
            os.system('/bin/mv '+header+'.out fastaout.log '+fastaOutD)
        pin.close()
    os.system('/bin/rm '+junk)
    
    
if __name__=='__main__':
    ih=inpHand('Usage: fastaManager.py [options]\ncreate close and distant multiple sequen alignments\nfor a set of fasta entries. Entries can be all in one file, or in many files',
               ' -a _AR_inpf (globbed) file containing one or more fasta entries. Header will be generated from the first word in the fasta header line',
               ' -b _A_outd (globbed) root directory where I will create subdirectories "3590", "e10", and "fastaout" directories (def: current dir)',
               ' -c _A_listf list of headers, used only to deglob inpf and outf',
               ' -d __optsfasta pass your own fasta options',
               ' -e __optsparse pass your own option to parse fasta results')
    ih.parse(locals(),sys.argv)

    if not outd: outd=os.getcwd()
    listH=[]
    if listf: listH=chomp(open(listf,'r').readlines())

    alnGen(inpf,outd,listH,optsfasta,optsparse)
    sys.exit(0)
