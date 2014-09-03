#!/usr/bin/python
import os,sys,re
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand


if __name__=='__main__':

    ih=inpHand('Usage: mix_by_column.py [options]',
               ' -a _AR_f1 file 1',
               ' -b _AR_f2 file 2',
               ' -c __c1 column number in file 1 to select (def: first column)',
               ' -d __c2 column number in file 2 to select (def: first column)',
               ' -e __compare compare method: _INBOTH_, _INF1_, _INF2_ (def: _INBOTH_)',
               ' -f __comments symbol denoting beginning of comment line (def: \'#\')',
               ' -g __inc_com include comment lines (y/n) (def:\'n\')',
               ' -i _A_outf output file name (def: standar output)'
               )
    ih.parse(locals(),sys.argv)

    #initiliaze optional input
    if not c1: c1=0  
    else: c1=int(c1)-1 #substract one because of C arrays begin @ zero
    if not c2: c2=0
    else: c2=int(c2)-1
    compare_methods=['_INBOTH_', '_INF1_', '_INF2_']
    if not compare: compare='_INBOTH_'
    if compare not in compare_methods: compare='_INBOTH_'
    if not comments:comments='#'
    inc_coms=['y','Y','yes','n','N','no']
    if not inc_com: inc_com='y'
    if inc_com not in inc_coms: inc_com='y'
    if inc_com=='Y' or inc_com=='yes': inc_com='y'
    if inc_com=='N' or inc_com=='no': inc_com='n'

    #read input
    f1l=chomp(open(f1,'r').readlines())
    f2l=chomp(open(f2,'r').readlines())

    #compare f1l and f2l according to the selected method.
    if not outf: out=sys.stdout
    else:  out=open(outf,'w')
    p=re.compile(comments)

    
    if compare!='_INF2_': #either in both files, or in file 1 but not in file2
        for l1 in f1l: #line by line of file 1
            #print 'l1=',l1
            if p.search(l1):
                continue  #out.write(l1)
            ll1=l1.split() #split its columns
            x1=ll1[c1]    #;print x1 #select column number c1  
            inf2=0        #x1 not in file 2
            for l2 in f2l:
                if p.search(l2):
                    continue   #out.write
                ll2=l2.split()
                x2=ll2[c2]
                if x1==x2:
                    inf2=1 #x1 in file 1
                    if compare=='_INBOTH_': out.write(l1+'  '+l2+'\n')
                    break #exit l2 loop
            if inf2==0 and compare=='_INF1_': out.write(l1+'\n')
    else: #in file 2 but not in file 1
        for l2 in f2l:
            if p.search(l2):
                continue  #out.write(l1)
            ll2=l2.split() #split its columns
            x2=ll2[c2]    #select column number c1
            inf1=0
            for l1 in f1l:
                if p.search(l1):
                    continue   #out.write
                ll1=l1.split()
                x1=ll1[c1]
                if x1==x2:
                    inf1=1
                    break
            if inf1==0 : out.write(l2+'\n')
