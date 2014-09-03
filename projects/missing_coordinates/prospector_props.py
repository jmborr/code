#!/usr/bin/python
"""
Given:
(1) a query sequence
(2) a set of missing atom chunks (MAC) in the query sequence
(3) a set of templates given by PROSPECTOR

Do:
For every MAC and template, calculate the coverage of the template on
the MAC. Then plot the average coverage for a given template global
Z-score Z* and MAC lenght L*, where average is taken over all
(MAC,template) pairs such that Z<Z* and L<L*
     <cov>= (sum_i \theta(Z_i-Z^*)*\theta(L_i-L^*)*cov_i ) /
                                (sum_i \theta(Z_i-Z^*)*\theta(L_i-L^*))

In addition, we store all intermediate info that shows up when reading
the PROSPECTOR templates in pickle-dump files
"""
import os,sys
from mac import collf,MAC,loadMac,loadMACs,updColletion
from utilities.small_utilities import Bye,junkName,chomp
from prospector.prospManager import templ_file

#global variables
macs=loadMACs(collf=collf) #load MAC's into dictionary
#macs={ '117E2.283_286':loadMac('117E2.283_286',file='/gpfs1/scratch/jose/missing_coordinates/out/117E2/117E2.283_286.mac'), }
currd=os.getcwd() #current directory
wd=os.path.join(currd,junkName()) #working directory
os.system('/bin/mkdir -p '+wd)
maxZ=50 ;maxL=300 ;avcov=[]; cumhis=[]
for Z in range(maxZ):
    avcov.append( [0.0]*maxL )
    cumhis.append( [0]*maxL )

#process each MAC
os.chdir(wd)
remaining=len(macs) ;nerr=0
for mac in macs.values(): #retrieve MAC object
    print remaining, mac.pdbheader
    isProspOut=False
    prospf=''
    prospf_0=os.path.join(mac.dumpdf,mac.pdbheader+'.out.prospector.')
    for suffix in ('tar','tar.bz2','tbz2'):
        prospf=prospf_0+suffix
        if os.path.exists(prospf):
            isProspOut=True
            break
    if not isProspOut:
        sys.stderr.write('ERROR: no prospector output in '+os.path.join(mac.dumpdf)+'\n')
        nerr+=1
        continue
    mac.update('ptor.outf',prospf) #update mac.ptor['outf']
    templfn=mac.pdbheader+'rap3orienrev5s.pdb' #template file base name
    mac.update('ptor.tmpl.fbn',templfn)
    mac.update('ptor.tmpl.list',{}) #initialize list that will contain template-related info
    tmplist=mac.ptor['tmpl']['list'] #handy reference
    ops='xf'
    if '.bz2' in prospf or '.tbz2' in prospf: ops='jxf'
    os.system('tar '+ops+' '+prospf+' -C '+wd+' pdbbpdborienrev && /bin/mv pdbbpdborienrev/* '+wd)
    p=templ_file( templfn ) #template-file object
    #print wd
    p.loadTemplates()       #load templates into memory
    for tmpl in p.tpls.values():
        #print tmpl
        mac.ptor['tmpl']['list'][ tmpl['header'] ]= {'seqid':tmpl['seqID'],
                                                     'Z':tmpl['Z'],
                                                     'nalg':tmpl['nalg']
                                                     }
        mac_tmpl=mac.ptor['tmpl']['list'][ tmpl['header'] ]  #handy reference
        mac_tmpl['maccovn']=0
        mac_tmpl['maccov']=0.0
        if mac.b_e[1] <= tmpl['alg'].gs[1].sL: #alignment may cover the MAC region
            #print mac.b_e[0],mac.b_e[1],tmpl['alg'].gs[1].sL
            chunk=tmpl['alg'].extractII(mac.b_e[0]-1,mac.b_e[1]-1)#extract alignment for MAC region
            mac_tmpl['maccovn']=chunk.nAligned() # number template residues covering the MAC region
            mac_tmpl['maccov']=float(mac_tmpl['maccovn'])/mac.L #coverage
            #update average coverage histogram
            for i in range( min(mac_tmpl['Z'],maxZ) ):
                for j in range( min(mac.L,maxL) ):
                    cumhis[i][j]+=1
                    avcov[i][j]+=mac_tmpl['maccov']
    mac.pickleDump() #update MAC pickle file
    remaining-=1
    
os.chdir(currd)
os.system('/bin/rm -r '+wd)
updColletion() #update the single file containing all pickles
for i in range(maxZ):
    for j in range(maxL):
        if cumhis[i][j]==0: continue
        print '%3d %3d %4.2lf'%(i,j,avcov[i][j]/cumhis[i][j])

if nerr: sys.stderr.write('There were '+nerr+' missing prospector outputs\n')

sys.exit(0)

