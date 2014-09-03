#!/usr/bin/python

"""We do a global sequence alignment between variant and PDB. If
sequence identity of aligned region, after removing singlet states as
spurious identically aligned residues, is bigger than 0.35, then we
extract from this region the identical aligned residues.

                  variant xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                  PDB                 yyyyyyyyyyyyyyyyy

Then we:
A. Extract contacts for CA and SC
B. Filter predicted contacts with the PDB distances between identical aligned residues
C. Substitute predicted secondary structure with that of the identical aligned residues
D. If coverage of identical aligned residues > 0.7, add as first template. If coverage of
   identical aligned residues > coverage of first template, add the identical aligned
   residues as first template
"""

import os,sys
from string import lowercase
from utilities.codedir import scratchdir
from utilities.small_utilities import Bye,junkName,chomp
from seq.align import getYangManageResultsAlignmentFormat
from seq.fastaManager import importFastaEntry
from seq.letters import one2three,three2one
from math import sqrt
from tasser.common import tasserconcuts
from secondary.seqdat2ss import translateSec
from secondary.pdb2ss import outputSec
from secondary.ss2seqdat import genSeqdat

# some global variables
min_templ_cov=0.5 #minimum coverage for the template to be considered global template

#   =================================================

def refresh_seqdat(chain_pdb,id_pairs,ss):
    """Update predicted secondary structure with template info

    chain_pdb: PDB template of template sequence
    id_pairs: correspondence between positions of query and template sequence
              for identically aligned residues
    ss: predicted secondary structure for the query"""
    
    ss0=outputSec(chain_pdb).replace('T','C') #string with secondary structure from the PDB file
    #Bye(`len(ss0)`+' '+`len(chain_pdb)`+' '+ss+'\n'+ss0)
    ssl=list(ss) #cast predicted secondary structure into a list format
    for (i,j) in id_pairs: ssl[i]=ss0[j] #update predicted sec. with that of the PDB file
    ss0=''.join(ssl) #;Bye(ss0)
    return ss0 #cast to a string format and output

#   =================================================

def distance(r1,r2):
    sum=0.0
    for i in range(3):
        d=r1[i]-r2[i]
        sum+=d*d
    return sqrt(sum)

#   =================================================

def refresh_combCA(combCAl,templ_pdb):
    """refresh CA contacts with template-derived info

    combCAl: lines of combCA.dat file, except for first line
    templ_pdb: template in prospector-output format

    We first remove those lines of combCAl containing contacts that are
    incompatible with templ_pdb. Then we retrieve contacts from templ_pdb
    that are not in the filtered combCAl and add them."""
    #print '\n'.join(templ_pdb)
    #read template CA coordinates into a dictionary
    coords={}
    resIDs=[]
    for l in templ_pdb:
        resID=int(l[21:26])
        resIDs.append(resID)
        coords[resID]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]
    #Bye(coords)
    
    #remove contacts involving residues in templ_pdb. They will be replaced later.
    removed=0
    n=0
    remaining=len(combCAl) #;Bye(remaining)
    while n<len(combCAl):
        i,j,x=combCAl[n] ; i=int(i) ; j=int(j) #;print n,remaining
        remaining=remaining-1
        if i in resIDs and j in resIDs:
            del combCAl[n] #;print 'REMARK',i,j,d
            removed+=1
            n-=1
        n+=1
        
    #add contacts from template
    extra=0
    L=len(resIDs)
    for i in range(L-1):
        resi=resIDs[i]   #position along the sequence for residue "i"
        for j in range(i+1,L):
            resj=resIDs[j]  #position along the sequence for residue "i"
            if abs(resj-resi)>2: #not a local-contact
                if distance(coords[resi],coords[resj])<7.0:
                    a=resi ; b=resj
                    if resi>resj: a=resj ; b=resi
                    combCAl.append([`a`,`b`,'1.00000']) #new contact
                    extra+=1
                    #print a,b
    #Bye(`removed`+' '+`extra`)

#   =================================================

def refresh_comb(var_seq,combl,templ_pdb):
    """refresh SC contacts with template-derived info

    combl: lines of comb.dat file, except for first line
    templ_pdb: template in prospector-output format

    We first remove those lines of combl containing contacts that are
    incompatible with templ_pdb. Then we retrieve contacts from templ_pdb
    that are not in the filtered combl and add them."""

    #print '\n'.join(templ_pdb)
    #read  CA coordinates into a dictionary
    coords={}
    resIDs=[]
    for l in templ_pdb:
        resID=int(l[21:26])
        resIDs.append(resID)
        coords[resID]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]

    #create list of residues in three letter code for var_seq
    XXX=[ one2three[X] for X in var_seq ] #;Bye(XXX)
    
    #remove contacts involving residues in templ_pdb. They will be replaced later.
    removed=0
    n=0
    while n<len(combl):
        i,j,x=combl[n] ; i=int(i) ; j=int(j)
        if i in resIDs and j in resIDs:
            del combl[n]
            removed+=1
            n-=1
        n+=1
    
    #add contacts from template
    extra=0
    tc=tasserconcuts()  #;Bye(tc)
    L=len(resIDs)
    for i in range(L-1):
        resi=resIDs[i]
        for j in range(i+1,L):
            resj=resIDs[j]
            if abs(resj-resi)>2: #not a local-contact
                xxxi=XXX[resi-1] ; xxxj=XXX[resj-1] #substract one because of index shift 
                if distance(coords[resi],coords[resj]) < tc[xxxi][xxxj]:
                    a=resi ; b=resj
                    if resi>resj: a=resj ; b=resi
                    combl.append([`a`,`b`,'1.00000']) #new contact from templ_pdb
                    extra+=1
    #Bye(`removed`+' '+`extra`)

#   =================================================

def refresh_dist(distl,templ_pdb):
    """refresh local distance restraints with template-derived info

    distl: lines of dist.dat file, except for first line
    templ_pdb: template in prospector-output format

    We first remove those lines of distl containing restraints that are
    incompatible with templ_pdb. Then we retrieve restraints from templ_pdb
    that are not in the filtered distl and add them."""

    #read  CA coordinates into a dictionary
    coords={}
    resIDs=[]
    for l in templ_pdb:
        resID=int(l[21:26])
        resIDs.append(resID)
        coords[resID]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]

    #remove restraints involving residues in templ_pdb. They will be replaced later.
    removed=0
    n=0
    remaining=len(distl) #;Bye(remaining)
    while n<len(distl):
        i=int(distl[n][0]) ; j=int(distl[n][1])
        remaining=remaining-1
        if i in resIDs and j in resIDs:
            del distl[n]
            n-=1
            removed+=1
        n+=1

    #add local restraints from template
    extra=0
    L=len(resIDs)
    for i in range(L-1):
        resi=resIDs[i]   #position along the sequence for residue "i"
        for j in range(i+1,min(i+6,L)):
            resj=resIDs[j]  #position along the sequence for residue "i"
            if abs(resj-resi)>2 and abs(resj-resi)<6: #local-distance restraint
                a=resi ; b=resj
                extra+=1
                if resi>resj: a=resj ; b=resi
                d=distance(coords[a],coords[b])                
                distl.append([`a`,`b`,`b-a`,'%6.3lf'%d,'0.0']) #new contact
    #Bye(`removed`+' '+`extra`)

#   =================================================

def refresh_distL(distLl,templ_pdb):
    """refresh long-range distance restraints with template-derived info

    distLl: lines of distL.dat file, except for first line
    templ_pdb: template in prospector-output format

    We first remove those lines of distLl containing restraints that are
    incompatible with templ_pdb. Then we retrieve restraints from templ_pdb
    that are not in the filtered distLl and add them."""

    #read  CA coordinates into a dictionary
    coords={}
    resIDs=[]
    for l in templ_pdb:
        resID=int(l[21:26])
        resIDs.append(resID)
        coords[resID]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]

    #remove restraints involving residues in templ_pdb. They will be replaced later.
    removed=0
    n=0
    remaining=len(distLl) #;Bye(remaining)
    while n<len(distLl):
        i=int(distLl[n][0]) ; j=int(distLl[n][1])
        remaining=remaining-1
        if i in resIDs and j in resIDs:
            del distLl[n]
            removed+=1
            n-=1
        n+=1

    #add long-range restraints from template
    extra=0
    L=len(resIDs)
    for i in range(0,L-1,3):
        resi=resIDs[i]   #position along the sequence for residue "i"
        resj_prev=resi+6-2    #6 is the minimum initial sequence jump, 2 is the consecutive jump
        for j in range(i+1,L):
            resj=resIDs[j]  #position along the sequence for residue "i"
            if resj>resj_prev+2:
                resj_prev=resj
                a=resi ; b=resj
                extra+=1
                if resi>resj: a=resj ; b=resi
                d=distance(coords[a],coords[b])                
                distLl.append([`a`,`b`,'%6.3lf'%d]) #new long-range distance

    #Bye(`removed`+' '+`extra`)
#   =================================================

def generate_template(chain_pdb,var_seq,chain_seq,indexPairsIdentical):
    """create a prospector-like template

    chain_pdb: pdb of template with template sequence
    var_seq: query sequence
    chain_seq: template sequence
    indexPairsIdentical: list of pairs (i,j) denoting position of identical
                         residues in the alignment of var_seq and chain_seq

    We return a list of pdb lines that constitue the template of the residues
    of var_seq which were identical to those of chain_seq as given by indexPairsIdentical
    """
    templ=[]  #;Bye(len(chain_pdb))
    for (var_i,chain_i) in indexPairsIdentical:
        var_res=one2three[ var_seq[var_i] ]        #residue identity of variant
        chain_res=one2three[ chain_seq[chain_i] ]  #residue identity of pdb chain
        #print `var_i`+' '+var_res+' '+`chain_i`+' '+chain_res
        l=chain_pdb[chain_i]
        n='%5d'%(1+var_i) #add one because origin of indexes is zero (like C array)
        templ.append('ATOM     1   CA  '+var_res+' '+n+l[26:54]+'%5d'%chain_i+' '+chain_res)
    #Bye('\n'.join(templ))
    return templ

#   =================================================

#Bartosz put an image of the PDB only on cng0002
if os.popen('hostname').readline().strip() != 'cng0002':
    sys.stderr.write('You must be on cng0002 to run this script\n')
    sys.exit(1)

#temporary files to store output of some programs
junkf=junkName()

#list of annotated PDB's for each reference
#pin=open( os.path.join(scratchdir,'spliceVariants/preparing_database/refs300_nonduplicated_withPDB.dat'),'r')
pin=open( os.path.join(scratchdir,'spliceVariants/preparing_database/junk'),'r') #for debugging

#read reference by reference
l0=pin.readline()
while l0 and l0[0:2]=='AC':
    ref=l0[3:11]   #;print l0
    root=ref[0:6]
    #gather all annotated PDB's
    pdbheaders=[]
    l0=pin.readline()
    while l0 and l0[0:2]=='DR':
        pdbheaders.append(l0[10:14].lower())
        l0=pin.readline()  #;print l0
    #gather all variants of the reference
    cmd='for g in `ls -1 '+scratchdir+'/spliceVariants/input/'+root+'*`;do basename $g;done'
    vars=chomp(os.popen(cmd).readlines()) #;print vars
    #work one variant at a time
    for var in vars:
        #import fasta file
        var_fastaf=scratchdir+'/spliceVariants/input/'+var #;Bye(var_fastaf)
        var_header,var_seq=importFastaEntry(open(var_fastaf,'r')) #;Bye(var_seq)
        #create working directory under subdirectory out4 and chdir there
        wd=scratchdir+'/spliceVariants/out4/'+root+'/'+var[-2:]
        os.system('/bin/mkdir -p '+wd)
        #copy and inflate tasser input from out2 in the working directory
        file=os.path.join(scratchdir,'spliceVariants/out2',root,var[-2:],var+'.in.tasser.tar')
        if not os.path.exists(file):
            sys.stderr.write(file+' does not exists\n')
            continue #skip to next variant
        os.system('tar xf '+file+' -C '+wd) #;Bye(file)
        #read comb.dat contacts, which are to be filtered
        combl=[]
        for l in open(wd+'/comb.dat').readlines()[1:]: combl.append(l.split())
        #read combCA.dat contacts, which are to be filtered
        combCAl=[ l.split() for l in open(wd+'/combCA.dat').readlines()[1:] ]
        #read seq.dat file, which is to be filtered
        ss=translateSec(wd+'/seq.dat') #;Bye(ss)
        #read dist.dat distance restraints, which are to be filtered
        distl=[l.split() for l in open(wd+'/dist.dat').readlines()[1:] ]
        #read distL.dat distance restraints, which are to be filtered
        distLl=[l.split() for l in open(wd+'/distL.dat').readlines()[1:] ]
        #variant length
        L=int(open(wd+'/rmsinp').readlines()[1].strip()) #;Bye(`L`+' '+`len(var_seq)`)
        #read coverage of first template in chain.dat
        first_templ_cov=float( os.popen('head -2 '+wd+'/chain.dat|tail -1').readline().strip() )/L
        templ_cov=0.0 #initialize coverage of current template
        best_templ=[] #initialize coverage of best found template
        record_best_templ=False #flag recording of best_templ
        prev_tc=0.0 #coverage of current best_templ
        refresh_flag=False #marks if we do any update in tasser input
        #read each of the annotated PDB. Then update contact information, secondary
        #structure information, and maybe assign as putative first template
        for pdbheader in pdbheaders:
            pdbf='/local/images/pdb-2007062900/pdb/pdb'+pdbheader+'.ent' #;Bye(pdbf)
            if not os.path.exists(pdbf):
                sys.stderr.write(pdbf+' does not exists\n')
                continue #go to next PDB if not in local library
            #extract sequences of all chains from ATOM record including no gaps
            cmd='pdbseq -a 2 -o '+junkf+' '+pdbf+' &>/dev/null'
            os.system(cmd) #;Bye(junkf)
            #work chain by chain
            pin2=open(junkf,'r') ; l2=pin2.readline() #;Bye(l2)
            while l2:
                if l2.find('sequence')>=0:
                    chain_id=l2.strip()[-1] #retrieve which chain in the PDB we're dealing with
                    if chain_id==':': chain_id=' ' #when there's only one chain                    
                    #retrieve CA coordinates of the particular chain (ignore missing CA's)
                    cmd='chnout -c "'+chain_id+'" '+pdbf+' 2>/dev/null|grep " CA "' #;Bye(cmd)
                    chain_pdb=chomp(os.popen(cmd).readlines()) #;Bye('\n'.join(chain_pdb))
                    chain_seq=''.join([three2one[x[17:20]] for x in chain_pdb])  #;Bye(chain_seq)
                    #align the sequence of each chain agains sequence of variant
                    x=getYangManageResultsAlignmentFormat(var_seq,chain_seq)
                    #print x;sys.exit(1)
                    #number of identical res. over the variant length                    
                    id_cov=float( x.identicalL() ) / L #;Bye(id_cov)
                    #number of identical res. over aligned length
                    id_alg=x.seqIdentity() #;Bye(id_alg)
                    if id_alg>0.35: #if clearly evolutionary related
                        refresh_flag=True
                        print var_header,'id_alg=%4.2lf id_cov=%4.2lf'%(id_alg,id_cov),pdbheader,chain_id

                        #generate a prospector-like template from the alignment "x"
                        #and the extracted chain_pdb
                        id_pairs=x.listIndexPairsIdentical() #;print id_pairs;sys.exit(1)
                        templ_pdb=generate_template(chain_pdb,var_seq,chain_seq,id_pairs)
                        #refresh secondary structure
                        ss=refresh_seqdat(chain_pdb,id_pairs,ss)
                        #remove contacts inconsistent with templ_pdb, add contacts from templ_pdb
                        refresh_combCA(combCAl,templ_pdb)
                        refresh_comb(var_seq,combl,templ_pdb)
                        #remove distance restraints inconsistent with templ_pdb, add new ones
                        refresh_dist(distl,templ_pdb)
                        refresh_distL(distLl,templ_pdb)
                        #check if the template is worth to be a global template
                        tc=float( len(templ_pdb) ) / L #;Bye(tc) #coverage of template
                        if tc > prev_tc: #bigger cov. than previous templ
                            best_templ=templ_pdb[:]
                            prev_tc=tc
                        if tc > min_templ_cov or tc > first_templ_cov: record_best_templ=True
                l2=pin2.readline()
        #output changes
        if refresh_flag:
            open(wd+'/combCA.dat','w').write( `len(combCAl)`+'\n'+
                                            '\n'.join( [' '.join(l) for l in combCAl] ) )
            open(wd+'/comb.dat','w').write( `len(combl)`+'\n'+
                                          '\n'.join( [' '.join(l) for l in combl] ) )
            open(wd+'/dist.dat','w').write( `len(distl)`+'\n'+
                                          '\n'.join( [' '.join(l) for l in distl] ) )
            open(wd+'/distL.dat','w').write( `len(distLl)`+'\n'+
                                           '\n'.join( [' '.join(l) for l in distLl] ) )
            open(wd+'/seq.dat','w').write( genSeqdat(var_seq,ss) )
            open(wd+'/best_PDB_templ.pdb','w').write( '\n'.join(best_templ) )
        if record_best_templ:
            chaindatl=open(wd+'/chain.dat','r').readlines() #input chain.dat
            n,diff=chaindatl[0].split() ; n=`1+int(n)`  #increase number of templates by one
            buf ='  '+n+'  '+diff+'\n' #number of templates and difficulty level
            buf+='%5d\n'%len(best_templ)
            buf+='\n'.join( [ l[21:38]+' '+l[38:46]+' '+l[46:54] for l in best_templ ] )+'\n'
            buf+=''.join(chaindatl[1:]) #rest of templates as given in chain.dat
            open(wd+'/chain.dat','w').write(buf)

    
os.system('/bin/rm '+junkf)
sys.exit(0)

