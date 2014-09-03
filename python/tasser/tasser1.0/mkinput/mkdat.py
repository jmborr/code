#!/usr/bin/python
import sys,os,re,copy
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.codedir import codedir,libseq,yangpdb
from utilities.small_utilities import junkName,Bye
from jobs.job import pastry
from seq.fastaManager import importFastaEntry,writeFasta
from seq.letters import valid_A
from jobs.job import pastry
from seq.blastManager import blastRun,blastAlignment

zal=codedir+'/perl/tasser/tasser1.0/mkinput/bin/zal'
dat=codedir+'/perl/tasser/tasser1.0/mkinput/bin/dat'

def formatTempl(before,after):
    p0=re.compile('TER')
    p1=re.compile('ATOM')
    inf=open(before,'r')
    outf=open(after,'w')
    n=os.popen('grep TER '+before+'|wc|tr -s " "|cut -d" " -f2').readline()
    outf.write(n)
    c=1 #template number
    for line in inf.readlines():
        if p0.search(line):
            outf.write('TER\n')
            c+=1 #another template read
        else:
            if p1.search(line):
                aa=line[17:20]
                n1=int(line[22:26]) ; n2=int(line[55:59])
                rx=float(line[30:38]) ; ry=float(line[38:46]) ; rz=float(line[46:54])
                outf.write('ATOM %6d  CA  %s%6d    %8.3f%8.3f%8.3f%5d\n'%(n1,aa,n1,rx,ry,rz,n2))
            else:
                it=line.split()
                it[1]=int(it[1])  ;  it[2]=float(it[2])
                outf.write('%6d%12.5f%7d     %s\n' % (it[1],it[2],c,it[0]) )
    inf.close()
    outf.close()
    
def formatCont(before,after):
    inf=open(before,'r')
    outf=open(after,'w')
    outf.write(inf.readline()) #fist line
    for line in inf.readlines():
        it=line.split()
        it[0]=int(it[0])  ;  it[1]=int(it[1])
        outf.write( '%5d %5d\n'%(it[0],it[1]) )
    inf.close()
    outf.close()

#gen_dat expects to find the following files: contact.dat, template.dat, seq.dat, seq.fasta
def gen_dat(libdir=yangpdb,blastf='',outd='.',co=0.3):
    outd=addAbsPath(outd)
    #write protein.seq file from seq.fasta
    header,seq=importFastaEntry(open('seq.fasta','r')) #;Bye(header+'\n'+seq)
    writeFasta('protein',seq,filename='protein.seq',n=60)

    #create blast.out
    if blastf: os.system('/bin/cp '+blastf+' blast.out')
    else:
        blastOuts=blastRun('protein.seq') #run psi-blast
        os.system('/bin/cp '+blastOuts['blast']+' blast.out')

    #create protein.frq
    prof={}
    for A in valid_A: prof[A]=0 #initialize generic profile
    log=[]
    for i in range(0,len(seq)): log.append(copy.deepcopy(prof))    
    #go to last round
    fpt=open('blast.out','r') ; line=fpt.readline()
    r=re.compile('\s*Results from round 3|CONVERGED\!')
    while not r.match(line):    line=fpt.readline() #scan until match
    #read all blast alignments
    while(1):
        alg=blastAlignment(fpt,querySeq=seq) #extended alignment to full sequence of "Query"
        if alg.error: break #reading error or just end-of-file reached
        if alg.pd['I']<0.98 and alg.pd['E']<0.001:
            gs=alg.Yalg.project( alg.ix['Sbjct'], alg.ix['Query'] ).gs #gapped sequence of Sbjct onto Query
            for i in range(0,len(gs)):
                if gs[i] in valid_A: log[i][ gs[i] ]+=1
    pout=open('protein.frq','w') ; pout.write(`len(seq)`+'\n')
    for i in range(0,len(seq)):
        pout.write('%3d'%(i+1))
        norm=0 ; prof=log[i]
        for A in valid_A: norm+=prof[A]
        for A in 'CMFILVWYAGTSQNEDHRKP':
            x=1.0*prof[A]
            if norm!=0: x=x/norm
            pout.write('%10.7lf'%x) #precise A order required
        pout.write('\n')
    pout.close()

    #write in.dd
    pout=open('in.dd','w')
    pout.write('protein.seq\nprotein.frq\n')
    pout.write(libdir+'/\n') #notice the trailing '/'
    pout.write('seq.pdb\nseq.dat\n')
    pout.write('%5.3lf'%(co)+'\n') #sequence cutoff
    pout.write('template.dat\ncontact.dat\n')
    pout.write('8   15.5\n') #zscore cutoff for medium and easy classes for PROSPECTOR_3 templates. This two numbers are adjustable when PROSPECTOR_3 is updated
    pout.close()

    #run black box
    os.system('/bin/cp '+zal+' '+dat+' .')
    print 'running zal'  ;  os.system('./zal')
    print 'running dat'  ;  os.system('./dat')
    for file in ['comb.dat','combCA.dat','dist.dat','distL.dat','chain.dat']:
        if not os.path.exists('./'+file): Bye('From mkdat.py: failed to generate '+file)
    if os.getcwd()!=outd: pastry('/bin/mv comb*.dat dist*.dat chain.dat '+outd)
    
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: mkdat.py',
               ' -a _R_header five letter PDB header',
               ' -b _AR_fastaf sequence file in FASTA format',
               ' -c _AR_seqdat secondary structure propensity',
               ' -d _A_templ template.dat-like file',
               ' -e _A_templ2 template file from prospector (in place of -d flag)',
               ' -f _A_pred contact.dat-like file',
               ' -g _A_pred2 predicted contacts from prospector (in place of -d flag)',
               ' -i _A_libdir library of structures (def:/library/yzhang/PDB)',
               ' -k _A_outd output directory (def: current dir)',
               ' -l _A_blastf results from a previous blast run',
               ' -m __co sequence cutoff(def:0.3)'
               )
    ih.parse(locals(),sys.argv)

    if not co:co=0.3
    else: co=float(co)

    print '...running mkdat.py...'
    if not templ and not templ2: ih.abort('ERROR: pass either -d or -e flag')
    if not pred and not pred2  : ih.abort('ERROR: pass either -f or -g flag')
    if not libdir: libdir=yangpdb
    if not outd: outd=os.getcwd()
    for file in (fastaf,seqdat,templ,templ2,pred,pred2,libdir,blastf):
        if file and not os.path.exists(file):ih.abort('ERROR: '+file+' does not exists')

    if not os.path.exists('./seq.fasta'): os.system('/bin/cp '+fastaf+' ./seq.fasta')
    if not os.path.exists('./seq.dat'): os.system('/bin/cp '+seqdat+' ./seq.dat')

    #template file
    if templ2 and not templ:
        templ='./template.dat'  ;  formatTempl(templ2,templ)
    else:
        if not os.path.exists('./template.dat'):
            os.system('/bin/cp '+templ+' ./template.dat 2>/dev/null')

    #contacts file
    if pred2 and not pred:
        pred='./contact.dat'  ;  formatCont(pred2,pred)
    else:
        os.system('/bin/cp '+pred+' ./contact.dat 2>/dev/null')

    gen_dat(libdir=libdir,blastf=blastf,outd=outd,co=co)
    
    #print '...mkdat.py calling mkdat.pl...'
    #pastry(codedir+'/perl/tasser/tasser1.0/mkinput/mkdat.pl -l '+libdir+' -d . -o '+outd)
        
    print '...finished mkdat.py...'
    sys.exit(0)
