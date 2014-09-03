#!/usr/bin/python
"""

For a pair of variants, let L denote the (long) sequence with the
insertion, and S the (short) sequence without the insertion.

We obtain Top and Best combo models to native, ie, Ltop, Lbest, Stop,
Sbest, and several superpositions with TMalign


"""
import sys,os
from utilities.small_utilities import chomp,junkName
from spicker.spickerYangResultsManager import spickOut,RMSD
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import alignOut
from pdb.create_single_chain import listCAs
from prospector.prospManager import prospOut
from jobs.job import job

myroot='/gpfs1/scratch/jose/spliceVariants/benchmark'
ind=myroot+'/out' #results from spicker are here
tmpldir='/gpfs1/archive/adrian/dat/pdb_may06.apr06.aux'
pdbdir='/gpfs1/archive/adrian/dat/pseudo_splicing/ca' #native pdb files
outd=myroot+'/l1_200_l2_200_next_0_nint_1_lbig_5_seqid_0.6_cov_0.2/superpositions' #output dir
pout=open('compare_variants.dat','w')

algs=myroot+'/preparing_database/pdb.40.300.alns'
algsx=myroot+'/preparing_database/pdb.40.300.alns.idx'

pairs=[]   #will contain the pairs of variant names
lengths=[] #will contain sequence lengths for each variant pair

#typical line of file pdb.40.300.alns_l1_200_l2_200_next_0_nint_1_lbig_5_seqid_0.6_cov_0.2.dat:
#1m8gA 1hybA 167 161 167 161 159 0.964 0.988   0   0   0   1   6
datfile=myroot+'/preparing_database/pdb.40.300.alns_l1_200_l2_200_next_0_nint_1_lbig_5_seqid_0.6_cov_0.2.dat'
#datfile=myroot+'/preparing_database/toy.dat'
for line in os.popen('grep -v -e "#" '+datfile).readlines():
    pairs.append( line.split()[0:2] )
    a,b=line.split()[2:4]
    a=int(a) ; b=int(b) ; lengths.append( (a,b) )
#some explanatory lines
pout.write('# ./compare_variants.py\n')

pout.write('# For a pair of variants, let L denote the (long) sequence with the insertion, and S the (short) sequence without the insertion\n')

pout.write('#We obtain Top and Best combo models to native, ie, Ltop, Lbest, Stop, Sbest, and several superpositions with TMalign. We also output pdb files with the superpositions\n')

pout.write('#1-L 2-TM(Lbest,natL) 3-TM(Ltop,natL) 4-S 5-TM(Sbest,natS) 6-TM(Stop,natS) 7-TM(natL,natS) 8-TM(Lbest,Sbest) 9-TM(Ltop,Stop)  10-RMSD(Sbest,natS) 11-length(natL) 12-RMSD(Lbest,natL) 13-length(insertion) 14-RMSD(Lbest_ins,natL_ins) 15-length(rest=natL-insertion) 16-RMSD(Lbest_rest,natL_rest) 17-length(extended-insertion) 18-RMSD(Lbest_ext,natS_ext) 19-Aligned-in-the-insertion 20-Identical-in-the-insertion 21-#aligned-in-N-terminal 22-#aligned-in-C-terminal\n')

pout.write('#  1     2     3    4      5     6     7    8     9       10  11    12   13   14  15  16     17  18   19  20  21  22\n')

junkN=junkName() ; junkM=junkName() #two temporary files

for i in range(0,len(pairs)):
    if lengths[i][0] > lengths[i][1]:
        L=pairs[i][0] ; S=pairs[i][1] #denotes Long and Short sequence, respectively
    else:
        L=pairs[i][1] ; S=pairs[i][0]
    
    natL=pdbdir+'/'+L+'.pdb' #native file for long sequence
    natS=pdbdir+'/'+S+'.pdb' #native file for short sequence
    
    sL=spickOut(dir=ind+'/'+L[1]+'/'+L,nat=natL) #spickOut object for long sequence
    sS=spickOut(dir=ind+'/'+S[1]+'/'+S,nat=natS)
    
    if sL.readError or sS.readError:
        sys.stderr.write('#ERROR: reading sL or sS\n')
        sys.stderr.write('#      sL.readError='+sL.readError+'\n')
        sys.stderr.write('#      sS.readError='+sS.readError+'\n')
        continue #error reading spicker outputs
    
    rankL=sL.rankIDsByDens()[0:5] #ID's for top five clusters, ranked by density
    rankS=sS.rankIDsByDens()[0:5]

    LtopID=rankL[0] #cluster ID for top cluster
    StopID=rankS[0]

    Ltopcombo=sL.combo[ LtopID ] #filename for top combo
    Stopcombo=sS.combo[ StopID ]

    LbestID=sL.rankIDsByTMtoNat(list=rankL)[0] #ID for cluster with best TM score to native
    SbestID=sS.rankIDsByTMtoNat(list=rankS)[0]

    Lbestcombo=sL.combo[ LbestID ] #filename for best combo
    Sbestcombo=sS.combo[ SbestID ]

    #do several alignments with TMalign
    # TM(Lbest,natL) alignment of best cluster to natice
    prefix=outd+'/'+L+'BEST_'+L+'.TMalign'
    os.system('TMalign '+Lbestcombo+' '+natL+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_Lbest_natL=TMalignOut(prefix+'.dat').tm #TMalignOut object reporting the TM score
    # TM(Ltop,natL)
    prefix=outd+'/'+L+'TOP_'+L+'.TMalign'
    os.system('TMalign '+Ltopcombo+' '+natL+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_Ltop_natL=TMalignOut(prefix+'.dat').tm
    # TM(Sbest,natS)
    prefix=outd+'/'+S+'BEST_'+S+'.TMalign'
    os.system('TMalign '+Sbestcombo+' '+natS+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_Sbest_natS=TMalignOut(prefix+'.dat').tm
    # TM(Stop,natS)
    prefix=outd+'/'+S+'TOP_'+S+'.TMalign'
    os.system('TMalign '+Stopcombo+' '+natS+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_Stop_natS=TMalignOut(prefix+'.dat').tm    
    # TM(natL,natS) note we REQUIRE short sequence the second one for d0 purposes
    prefix=outd+'/'+L+'_'+S+'.TMalign'
    os.system('TMalign '+natL+' '+natS+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_natL_natS=TMalignOut(prefix+'.dat').tm
    # TM(Lbest,Sbest)
    prefix=outd+'/'+L+'BEST_'+S+'BEST.TMalign'
    os.system('TMalign '+Lbestcombo+' '+Sbestcombo+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_Lbest_Sbest=TMalignOut(prefix+'.dat').tm
    # TM(Ltop,Sbest)
    prefix=outd+'/'+L+'TOP_'+S+'TOP.TMalign'
    os.system('TMalign '+Ltopcombo+' '+Stopcombo+' -o '+prefix+'.pdb > '+prefix+'.dat')
    tm_Ltop_Stop=TMalignOut(prefix+'.dat').tm

    #refine best combos,build side-chains, and output them
    Lbest_basename=os.path.basename(Lbestcombo)
    Sbest_basename=os.path.basename(Sbestcombo)
    Lbest_rebuilt=L+'BEST.rebuilt.pdb'
    Sbest_rebuilt=S+'BEST.rebuilt.pdb'
    args='-a "/bin/cp -f '+Lbestcombo+' '+Sbestcombo+' ." -b "dukka.x '+Lbest_basename+' && pulchra best.pdb && /bin/mv best.pdb.rebuilt '+Lbest_rebuilt+' && dukka.x  '+Sbest_basename+' && pulchra best.pdb && /bin/mv best.pdb.rebuilt '+Sbest_rebuilt+'" -c "/bin/mv '+Lbest_rebuilt+' '+Sbest_rebuilt+' '+outd+'" -d '+outd+' -x ref'+L
    #job(exe='generic_job.py', exed='/gpfs1/active/jose/code/python/combo_jobs', args=args).qsub()
    
    #find RMSD values for insertion in the long sequence, as well as
    #rest and global and local-extended (insertion plus five residues
    #on either side) find where the insertion begins and ends in the
    #long sequence. We look at the alignment first pull the alignment
    #among all the alignments in the alignments file    
    index=int(os.popen('grep '+L+' '+algsx+' | grep '+S).readline().split()[0]) #;print index
    palg=open(algs,'r')#open file with all alignments
    palg.seek(index)   #goto alignment between L and S
    alg=alignOut(palg).alg #;print alg;sys.exit(1)#create alignment object
    palg.close()
    #who's the long sequence in the alignment?
    Lid=1 ; Sid=2 #initialization
    if alg.gs[1].sL < alg.gs[2].sL: #we compare lengths of ungapped sequences
        Lid=2 ; Sid=1
    bigGap=alg.extractBiggestInternalGap()    #;print bigGap;sys.exit(1)
    #the big gap must be in the gapped sequence of the short sequence. Find the position
    gsS=alg.gs[Sid].gs #;print gsS #string containing the short sequence including gaps
    #position in the gapped short sequence where gap begins, which is the same as to say position
    #in long gapped sequence where insertion begins
    n=gsS.find(bigGap) #;print n
    #use the correspondence between indexes of the gapped and ungapped sequence
    begin=alg.gs[Lid].gsi2si[n] #position in the long UNGAPPED sequence where insertion begins
    end=begin+len(bigGap)-1     #;print begin,end
    #extended loop includes five residues on each side of the loop
    begin_ext=begin-5
    end_ext=end+5
    if begin_ext<0: begin=1
    if end_ext >= alg.gs[Lid].sL: end_ext=alg.gs[Lid].sL-1 #last index if overflow
    #what are the extension of the N- and C- aligned regions ?
    Naligned=begin
    Caligned=alg.gs[Lid].sL-end-1
    #dump CA's into a list
    casN=listCAs(open(natL,'r')) #;print casN
    casM=listCAs(open(Lbestcombo,'r'))

    #RMSD(Lbest,natL) global rmsd
    open(junkN,'w').write( ''.join(casN) )
    open(junkM,'w').write( ''.join(casM) )
    lall=len(casN)
    rmsd=float( chomp(os.popen('rmsd.x '+junkM+' '+junkN).readline()) )
    
    #RMSD(Lbest_ins,natL_ins) insertion
    open(junkN,'w').write( ''.join(casN[begin+1:end+2]) )
    open(junkM,'w').write( ''.join(casM[begin+1:end+2]) )
    lloop=1+end-begin
    rmsdLoop=float( chomp(os.popen('rmsd.x '+junkM+' '+junkN).readline()) )
    #print rmsdLoop;sys.exit(1)

    #RMSD(Lbest_rest,natL_rest) everything but the insertion
    open(junkN,'w').write( ''.join(casN[0:begin])+''.join(casN[end+1:]) )
    open(junkM,'w').write( ''.join(casM[0:begin])+''.join(casM[end+1:]) )
    lrest=lall-lloop
    rmsdRest=float( chomp(os.popen('rmsd.x '+junkM+' '+junkN).readline()) ) 

    #RMSD(Lbest_loop_extension,natL_loop_extension)
    open(junkN,'w').write( ''.join(casN[begin_ext+1:end_ext+2]) )
    open(junkM,'w').write( ''.join(casM[begin_ext+1:end_ext+2]) )
    lext=1+end_ext-begin_ext
    rmsdExt=float( chomp(os.popen('rmsd.x '+junkM+' '+junkN).readline()) )

    #RMSD(Sbest,natS)
    rmsd_Sbest_Snat=float(RMSD(Sbestcombo,natS))
    #consider now if the first prospector template for L contains the insertion sequence
    dir=myroot+'/out/'+L[1]+'/'+L   #prospector results are under this directory
    p=prospOut(dir=dir,parser='genomesmay06') #load results into object    
    alg=p.tpls[0]['alg'] #alignment between long sequence L and sequence T of first template
    #print alg #; sys.exit(1)
    gsL=alg.gs[1] #gapped sequence for L
    gsT=alg.gs[2] #gapped sequence for T
    #we know where the insertion begins and ends in the ungapped sequence of L, thus retrieve the
    #corresponding indexes for the gapped sequence of L in the prospector alignment. However,
    #it may be that the PROSPECTOR alignment between L and T finish before the insertion begins,
    #or before the insertion ends
    ins_aligned=0 ; ins_identical=0 #initialize
    if gsL.sL>=begin: #the insertion is at least partially included in the prospector alignment
        b=gsL.si2gsi[begin-1] #; print 'b=',b
        if gsL.sL>=end: #the insertion is included in the prospector alignment in its entirety
            e=gsL.si2gsi[end-1]
        else:
            e=gsL.si2gsi[gsL.sL-1] #the insertion is truncated in the PROSPECTOR alignment
        ins=alg.extract(b,e) #; print ins #a portion of the prospector alignment
        ins_aligned=ins.alignedL20() #;print 'ins_aligned=',ins_aligned
        ins_identical=ins.identicalL() #number of aligned residues in the insertion region

    
    pout.write('%s %5.3f %5.3f %s %5.3f %5.3f %5.3f %5.3f %5.3f  %5.2lf %3d %5.2lf %3d %5.2lf %3d %5.2lf %3d %5.2lf %3d %3d %3d %3d\n'%(L,tm_Lbest_natL,tm_Ltop_natL,S,tm_Sbest_natS,tm_Stop_natS,tm_natL_natS,tm_Lbest_Sbest,tm_Ltop_Stop,rmsd_Sbest_Snat,lall,rmsd,lloop,rmsdLoop,lrest,rmsdRest,lext,rmsdExt,ins_aligned,ins_identical,Naligned,Caligned) )

os.system('/bin/rm '+junkN+' '+junkM) #clean-up
pout.close()
sys.exit(0)
