#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq,projectsdir
from utilities.small_utilities import chomp
from jobs.job import gen_lf,pastry
from string import upper
inpHand('Usage: number_secondary_elements.py',
        ' -a __header five-letter pdb code',
        ' -b _A_hlist list of headers, mandatory if flag -a not provided',
        ' -c _R_ss letter representing a type of secondary structure (H=helix; B=residue in isolated beta bridge; E=extended beta strand; G=310 helix; I=pi helix; T=hydrogen bonded turn; S=bend)',
        ' -d __minL minimal lenght of one secondary structure element (default H=7, B=1, E=4, G=6, I=1, T=3, S=1)',
        ' -e __maxL maximal length of one secondary structure element(default H=24, )',
        ' -f __minS minimal separation between two consecutive secondary structure elements (def=2)',
        ' -g _AR_outf output file (headers and number of secondary elements of a particular kind',
        ' -h __help this script will read one, or a list of headers, and for each header will calculate the number of secondary elements of one particular kind'
        ).parse(locals(),sys.argv)

#default minimal lengths for particular elements of secondary structure
minlengths={'H':6, 'B':1, 'E':4, 'G':6, 'I':1, 'T':3, 'S':1}
maxlengths={'H':24, 'B':1, 'E':11, 'G':6, 'I':1, 'T':3, 'S':1}
if not header and not hlist:
    sys.stderr.write('ERROR: either -a or -b flag must be passed')
    sys.exit(1)
    
if not minS: minS=2
if not minL:  minL=minlengths[ss] #print hlist,outf
if not maxL:  maxL=maxlengths[ss] #print hlist,outf

if header:
    headers=[header,]
else:
    headers=open(hlist,"r").readlines()
    chomp(headers)
#print headers;sys.exit(1)
#database with secondary structure assignments
ssf=projectsdir+'/db/ss.txt'

ssl=open(ssf,'r')
out=open(outf,'w')

for header in headers:
    print header
    #find corresponding entry in ssf and retrieve a sinle string (seq)
    #containing the secondary structure.
    pattern=upper(header[0:4])+':'+upper(header[-1]); #print pattern
    #retrieve line with the match, and next 49 lines. Typically will contain
    #all secondary structure assignment of the header, plus extra lines
    #correspoding to a different(s) header(s)
    #grep=chomp(os.popen('grep -A 50 "'+pattern+'" '+ssf).readlines())
    grep=os.popen('grep -A 50 "'+pattern+'" '+ssf).readlines()    
    chomp(grep)
    grep.pop(0) #remove first line
    seq='' #single-line string containing secondary assignments
    for line in grep :
        if(re.compile('>').search(line)): break #we reached another header
        seq=seq+line

    #translate the sequence to lengths of different segments, and also
    #translate to '1' if segment of type ss, '0' otherwise
    i=b=c=0 #running sequence indexes
    is_ss=[]
    L=[]
    for item in seq : #go letter by letter in the sequence
        if item==ss :
            b=b+1 #b is length of sec.str.elem of type ss
            if c>0: #we finished going through a sec.str.elem not of type ss
                is_ss.append(0) #signal sec.str.elem 'i' is not of type ss
                L.append(c) #length of the sec.str.elem
                i=i+1 
                c=0
        else :
            c=c+1 #c is length of sec.str.elem of type other than ss
            if b>0:  #we finished going through a sec.str.elem of type ss
              is_ss.append(1) #signal sec.str.elem 'i' is not of type ss
              L.append(b) #length of the sec.str.elem
              i=i+1
              b=0
    #store last segment
    if b>0:
        is_ss.append(1)
        L.append(b)
    else:
        is_ss.append(0)
        L.append(c)
    
    #convert to type ss all sec.str.elem of type other than 'ss' that are
    #smaller than minS and flanked by sec.str.elem of type 'ss'
    i=0 
    end=len(is_ss)-1 #last segment index
    for l in is_ss : #l is either 0 (not of type ss) or 1 (of type ss)
        if l==0 and 0<i<end : #segment not of type ss and not at endpoints
            if L[i]<minS : is_ss[i]=1 #convert to type ss
        i=i+1

    #merge contiguous segments of type ss, if any (due to previous step)
    i=0
    end=len(is_ss)
    while i<end-1:
        if is_ss[i]==1:
            if is_ss[i+1]==1:
                L[i+1]=L[i]+L[i+1]
                L[i]=0     #append info from 'i' to 'i+1', leaving a "blank"
                is_ss[i]=0
        i=i+1

    #number of sec.str.elem of type 'ss' bigger than minL
    nnss=0 #number of elements of secondary structure of type other than ss
    i=0
    end=len(is_ss)
    while i<end:
        if is_ss[i]==1 and L[i]>=minL : nnss=nnss+1
        if L[i]>=minL+maxL : nnss=nnss+1 #extra count for a too long segment
        i=i+1

    out.write("%s %2d\n" % (header,nnss)) #write to output file
    #sys.stdout.write("%s %2d\n" % (header,nnss)) #write to standard output
    
sys.exit(0)

