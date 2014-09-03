#!/usr/bin/python

#import some modules to ease the pain
import sys,os,re
from inputArgs import inpHand
from small_utilities import chomp,replace_at
from pdbDomains import storeSingleChain,extractSegment,OneLetterSeq
from align_yang_f77 import align

# Function to parse file containing multiple fasta entries. Typical
# entry is of the form:
# >1u6g_B mol:protein length:108 RING-box protein 1
# MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWA
# CNHAFHFHCISR
# unfortunately, '2b9p_2' in 'pdb_seqres.txt' would meand chain id "2"
# in pdb file 2b9p, while in 'seq.pdb' would mean domain number 2 in
# pdb 2b9p. We have to provide for this with variable 'mode'
def get_chains(fastaf,mode,filter=[]):
    chains={} #dictionary to store all pairs (id,seq)
#match '>1cdb_','>1dky_B','> 133l_','> 12e8L2','> 12e8L','> 7taa_1','>1apy.3'
    pat=re.compile('>\s?([\w\.]+)')
    if mode==0: pat2=re.compile('_[A-Z,0-9]$')
    pin=open(fastaf,'r')
    line=pin.readline() #very fist line has to be a fasta header
    while line:
        if line[0]=='>':
            id=pat.match(line).group(1) #retrieve sequence id
            if mode==0:  #transform '1d7p_M' to '1d7pM'
                if filter:
                    if id[0:4] not in filter:
                        line=pin.readline()
                        continue
                if pat2.search(id):
                    id=id.replace('_','')
            if mode==1: #make sure chain ID in upper case
                id.replace('.','_') #correct some stupid dots
                id=replace_at(id,4,id[4].upper())
            seq='' #initialize query sequence
            line=pin.readline()    
            while line and line[0]!='>': #get all the sequence chunks
                seq=seq+chomp(line)
                line=pin.readline()
                chains[id]=seq #add the pair (id,seq) to dictionary chains
    pin.close()
    return chains

def generate_pdb(pdb,chainId,outpdb,segment,endline=''):
    [s,e]=segment.split('-')
    #'s' residue sequence number, 'si' is residue insertion code
    [s,si]=pat2.match(s).groups()
    if not si: si=' '# no insertion code (most prob. case)
    [e,ei]=pat2.match(e).groups()
    if not ei: ei=' '
    #extract segment does that, extract a chain segment from the file
    #and store it in another file. We may run pulchra if the file
    #contains only CA atoms    
    return extractSegment(pdb,chainId,outpdb,int(s),int(e),si,ei,endline)


#command-line arguments
ih=inpHand(
    'Usage: updateMkLib.py ',
    '  -a _A_pdbd pdb directory (def: /mirrors/pdb_mirror/data/structures/all/pdb)',
    '  -b _A_pdbl pdb list (def: /mirrors/pdb_mirror/derived_data/pdb_seqres.txt)',
    '  -c _A_dbd  database of parsed pdb\'s (def: /mirrors/pdb_mirror/jose/PDB)',
    '  -d _A_dbl  list of fasta parsed files (def: dbd/seq.pdb)',
    '  -e _A_homol list of homologs (def: dbd/homologs)',
    '  -f __sic sequence identity cut-off (def: 0.70)',
    '  -g _A_scopf scop file (def:dbd/dir.cla.scop.txt_1.69)',
    '  -i _A_added file added.pdb (def:none)',
    '  -j _A_obsolf file with obsolete entries (def:none)'
    )
ih.parse(locals(),sys.argv)

#get current directory
currd=os.getcwd()

#initialize default variables
if not pdbd : pdbd ='/mirrors/pdb_mirror/data/structures/all/pdb'
if not pdbl : pdbl ='/mirrors/pdb_mirror/derived_data/pdb_seqres.txt'
if not dbd  : dbd  ='/mirrors/pdb_mirror/jose/PDB'
if not dbl  : dbl  =dbd+'/seq.pdb'
if not homol: homol=dbd+'/homologs'
if not scopf: scopf=dbd+'/dir.cla.scop.txt_1.69'
if not sic  : 
    sic  =0.70
else:
    sic=float(sic)

#check that input files and directories exist
for file in [pdbd,pdbl,dbd,dbl,homol,scopf,added,obsolf]:
    if file: #scop is not required, thus may be null string
        if not os.path.exists(file):
            sys.stderr.write('ERROR from updateMkLib.py: '+file+' does not exists!\n')
            sys.exit(1)
        
#retrieve the new pdb id's of 'added'. Typical line in file 'added' is
#-rw-r--r--    1 update   rcsbdata  2817504 Apr  4  2006 pdb1zgu.ent
#if added:
#    addedpdbs=open.(added,'r').readlines() ; chomp(addedpdbs)
#    for i in range(len(addedpdbs)):
#        addedpdbs[i]=addedpdbs[i].split()[-1][3:7] #pdb id
    
#store id's and sequences into two dictionaries
chains_new=get_chains(pdbl,0)
newlist=open('/mirrors/pdb_mirror/jose/PDB/newhomologs','r').readlines()
print len(newlist)
chomp(newlist)
for id in chains_new.keys():
    if id not in newlist:
        del chains_new[id]
newlist=chains_new.keys() ; print len(newlist)
#else: chains_new=get_chains(pdbl,0,filter=addedpdbs)
#chains_parsed=get_chains(dbl,1) #store all seq.pdb

#filter out domain id's of chains_parsed
#print 'library size=',len(chains_parsed.keys())
#for id in chains_parsed.keys():
#    if len(id)>5: del chains_parsed[id] #six-letter id is a domain
#print 'size after pruning domains id\'s=',len(chains_parsed.keys())

#remove from chains_new all in 'list' file (chain_parsed) and all in
#'homologs' file
#print 'size of pdb_seqres.txt=',len(chains_new.keys())
#parsedlist=chains_parsed.keys()
#for id in chains_new.keys():
#    if id in parsedlist: del chains_new[id]
#print 'size after removal of library entries=',len(chains_new.keys())
#hlist=open(homol,'r').readlines() ; chomp(hlist)
#for id in chains_new.keys():
#    if id in hlist: del chains_new[id]
#print 'size after removal of entries in homologs file=',len(chains_new.keys())

#Some new pdbs may be homodimers or similar. Thus we should remove homologous
#within new pdb first
#newlist=chains_new.keys() ; newlist.sort()
#now newlist looks like:
#[...,'101m_','1ckaA','1ckaB','1olm_',...]
#i=0
#while(1):
#    id=newlist[i]
#    seq=chains_new[id]
#    l=len(seq)
#    pdbid=newlist[i][0:4] #pdb id is the first four characters
#    j=i+1    
#    while(newlist[j][0:4]==pdbid): #both chains belong to same pdb file
#        id2=newlist[j]
#        seq2=chains_new[id2]
#        #nid, number of identical residues in the aligned region
#        (nid,nali)=align(seq,seq2)
#        si=float(nid)/l
#        if si>sic: #homologous
            #print id2,'is homologous to',id
#id2 should go to homologous file, irrespective of whether id is
#homologous to any sequence in the parsed library
#            os.system('echo '+id2+' >> '+homol)
#            del chains_new[id2] #remove homologous entry
#if we erase id2 from newlist, then we shift the whole newlist to the
#left from the j position. Thus no need to issue the j=j+1 statement
#            del newlist[j]
#        else:
#            j=j+1
#        if(j==len(newlist)): break #index out of range
#    i=i+1
#    if(i>=len(newlist)-1): break #index out of range
#print 'size after removal of homodimer partners=',len(chains_new.keys())

#if we pass scop file, then for each id in chains_new, pull all id's
#in chains_new and chains_parsed. These are putative homologous, as
#they belong to the same family. This is an extra filter.
#typical scop entry:
#d1s6aa_ 1s6a    A:      a.1.1.1 105306  cl=46456,cf=46457,sf=46458,fa=46459,dm=46460,sp=81667,px=105306
#if not added: # not necessary when proccesing a few pdbs
#    scopfam={}
#    psf=re.compile('sf=(\d+)') #will retrieve the scop family number
#    pin=open(scopf,'r') ; line=pin.readline()
#    while(line): #read scop file line by line
#        if line[0]=='d': #one valid entry
#            id=line[1:5]+line[5].upper() #pdbid plus chain id
#            scopfam[id]=psf.search(line).group(1) #scop family number
#        line=pin.readline()
#    pin.close()
    #now gather all members of the same family in two
    #dictionaries. For a given family number, scopnews[family] will
    #contain a list of all entries in chains_new that belong to the
    #same family. Analogous construction for scopparsed and
    #chains_parsed    
#    scopnews={}
#    scopparsed={}
#    for id in scopfam.keys():
#        family=scopfam[id]
#        if id in parsedlist:
#            if scopparsed.has_key(family): #already some id in same family
#                scopparsed[family].append(id)
#            else: #new family entry
#                scopparsed[family]=[id,]
#        elif id in newlist:
#            if scopnews.has_key(family):
#                scopnews[family].append(id)
#            else:
#                scopnews[family]=[id,]
#    print 'len(scopnews)=',len(scopnews),'len(scopparsed)=',len(scopparsed)
    #For a given family, scopnews[family]+scopparsed[family] contains
    #a list of putative homologous, since all belong to the same
    #family. Let's do sequence alignment to prune this list
    #First prune all the scopnews[family] lists
#    for family in scopnews.keys():
#        putativehom=scopnews[family] #list of putative homologous
#        i=0
#        j=1
#        if len(putativehom)==0: print 'SIZE OF putativehom==0!'
#        while(1):
#            id=putativehom[i]
#            seq=chains_new[id]
#            l=len(seq)
#            while(j<len(putativehom)):
#                id2=putativehom[j]
#                seq2=chains_new[id2]
#                l2=len(seq2)
#                (nid,nali)=align(seq,seq2)
#                if float(nid)/l > sic: #seq homologous to seq2
#                    os.system('echo '+id+' >> '+homol)
#                    del chains_new[id] #remove homologous entry
#                    #remove entry from the family. Since putativehom
#                    #is a reference to scopnews[family], we are
#                    #actually reducing the scopnews[family] list
#                    del putativehom[i]
#                    i=i-1 #needed trick, because we do i=i+1 later
#                    break
#                elif float(nid)/l2 > sic: #seq2 homologous to seq
#                    os.system('echo '+id2+' >> '+homol)
#                    del putativehom[j]
#                else:
#                    j=j+1
#            i=i+1
#            if i>=len(putativehom)-1: break
#    print 'size after filter scop homologs within new entries=',len(chains_new.keys())
#    #now scan chains_new against chains_parsed for homologous
#    for family in scopnews.keys():
#        s_news=scopnews[family]
#        if scopparsed.has_key(family):
#            s_parsed=scopparsed[family]
#            for id in s_news:
#                seq=chains_new[id]
#                l=len(seq)
#                for id2 in s_parsed:
#                    seq2=chains_parsed[id2]
#                    (nid,nali)=align(seq,seq2)
#                    if float(nali)/l>sic:
#                        os.system('echo '+id+' >> '+homol)
#                        del chains_new[id]
#                        break
#    newlist=chains_new.keys() ; newlist.sort()
#    scopfam.clear() ; scopnews.clear() ; scopparsed.clear()  #liberate memory
#    print 'size after filter scop homologs to library=',len(chains_new.keys())
#
##remove from chains_new all above sic sequence identity with
##chains_new. Sequence identity is number of identical residues in the
##alignment divided by sequence length
#i=0
#for id in chains_new.keys():
#    i=i+1
#    seq=chains_new[id] #query sequence
#    linv=1.0/len(seq) #sequence length
#    for id2 in parsedlist: #align query sequence to all in chains_parsed
#        seq2=chains_parsed[id2]
#        l2=len(seq2)
#        #nid: number of identical residues in the alignment
#        #nali: number of aligned residues
#        if l2*linv<sic: #trick to avoid alignment of long query sequences
#            nid=0
#        else:
#            (nid,nali)=align(seq,seq2)
#        si=nid*linv #our definition of sequence identity
#        if si>sic: #new entry is homologous to some template sequence
#            print i,id,'homolog to ',id2,len(seq),len(seq2),nid,si
#            os.system('echo '+id+' >> '+homol) #add entry to 'homologs' file
#            del chains_new[id] #remove the homolog entry
#            break #go to next new entry
#if new entry is not homologous to any template sequence, add it to
#chains_parsed. This way we avoid two new homologous entries to be
#inserted in the library
#    if(si<sic):
#        print i,id,'is not homologous',len(seq)
#        os.system('echo '+id+' >> /mirrors/pdb_mirror/jose/PDB/newhomologs') #delete this later
#        chains_parsed[id]=seq2 #add to library
#        parsedlist.append(id)  #add to the list of library id's


#store domains listed in scop database
scopdom={}
scopseg={} #for each domain, store the chain segments
#matches '134-190', '-134-190', '134A-190' segments
pat=re.compile('\s*-?(\d+[A-Z]?-\d+[A-Z]?)')
pat2=re.compile('([A-Z]):')
pin=open(scopf,'r') ; line=pin.readline()
while(line): #read scop file line by line
    #print line
    if line[0]=='d': #one valid entry
        fields=line.split()
        domainname=fields[0]  #name given by SCOP classification
        pdbID=domainname[1:5]
        chainID=domainname[5].upper()
        domainID=domainname[6]
        #print 'domainname=',domainname
        #print 'pdbID=',pdbID
        #print 'chainID=',chainID
        #print 'domainID=',domainID
        if not chainID=='.': #domain made up of a single chain
            chainname=pdbID+chainID
            id=chainname+domainID #domain name in our clasification
            if chainID=='_': cid=''
            else: cid=chainID+':'
            pattern=cid+'-?(\d+[A-Z]?-\d+[A-Z]?)'
            fragments=re.compile(pattern).findall(fields[2])
            #print fields[2],pattern,fragments
            #d2pcbb_ 2pcb B: => fragments==[], ie, domain span the whole chain
            if fragments:
                scopseg[id]=','.join(fragments)
                if not scopdom.has_key(chainname):
                    scopdom[chainname]=[domainID,]
                else: scopdom[chainname].append(domainID)
        else: #domain made up of two or more chains
              #first find id of each chain
              chainIDs=pat2.findall(fields[2])
              #create a domain for each chain
              for chainID in chainIDs:                      
                  chainname=pdbID+chainID
                  id=chainname+domainID
                  pattern=chainID+':-?(\d+[A-Z]?-\d+[A-Z]?)'
                  fragments=re.compile(pattern).findall(fields[2])
                  #d1lq8.1 1lq8  A:,B: => fragments==[]
                  if fragments:
                      scopseg[id]=','.join(fragments)
                      if not scopdom.has_key(chainname):
                          scopdom[chainname]=[domainID,]
                      else: scopdom[chainname].append(domainID)
    line=pin.readline()
pin.close()

#Save whole chain and domains.If chain not in scop, apply protein domain parser
pat2=re.compile('(\d+)([A-Z]?)') #matches '4526', '4526A'
chains_domains={} #will store pairs (id,seq) for the derived domains

#this just for debugging purposes
#for id in chains_new.keys():
#    if scopdom.has_key(id):
#        print id,scopdom[id]
#        for domName in scopdom[id]:
#            print 'domName=',domName,scopseg[id+domName]

#go through the new entries
i=0
for id in chains_new.keys():
    i=i+1
    print i,id,
    #store chain in pdb file
    pdb=pdbd+'/pdb'+id[0:4]+'.ent'  #whole pdb file
    pdbchain=dbd+'/new/'+id+'.pdb' #; print pdb,pdbchain,
    #extrach the chain from the whole pdb file, and store. We may have to run
    #pulchra if the pdb contains only CA atoms.
    if storeSingleChain(pdb,id[4],pdbchain,'TER\n')==0: #id[4] is chain id
        #we can't import/export the chain for whatever reasons, just
        #don't enter chain in the library        
        del chains_new[id]
    else:
       #because we don't care about missing residues (either at the
       #termini or inside the chain, the sequence of the stored pdb is
       #smaller than the sequence of chains_new[id]. Let's store this
       #reduced sequence    
       chains_domains[id]=OneLetterSeq(pdbchain)
       #If stored in scop database, use the info there to extract the domains
       if scopdom.has_key(id):
           print 'in scop',
           for dom_id in scopdom[id]:
               domName=id+dom_id
               print 'dom_id=',dom_id,       
               tmpfile='./junkdomain12345.pdb'
               pdbdomain=dbd+'/new/'+domName+'.pdb'
               #erase just in case
               os.system('/bin/rm '+pdbdomain+' 2>/dev/null')
               #domain may be made up of more than one chain segment
               segments=scopseg[domName].split(',')
               error=0
               for segment in segments:
                   if generate_pdb(pdb,id[4],tmpfile,segment)==0:
                       error=1
                       break
                   os.system('cat '+tmpfile+' >> '+pdbdomain+' && /bin/rm '+tmpfile+' 2>/dev/null')
               if not error:
                   #armonize residue number among the chunks
                   storeSingleChain(pdbdomain,id[4],tmpfile,'TER\n')
                   os.system('/bin/mv '+tmpfile+' '+pdbdomain)
                   #retrieve domain seq from the stored pdb file of the domain
                   chains_domains[domName]=OneLetterSeq(pdbdomain)
               print '\n'
       #If id not stored in scop database, apply protein domain parser,
       #pdp. Typical output (to standard output) of pdp: 2cvxA
       #76-195/196-218,428-625,686-746/219-427/626-685 the previous
       #protein has three domains in chain A. Second domain made up of
       #three chunks!
       else:
           print 'apply pdp...'
           line=chomp(os.popen('./pdp '+pdb+' '+id[4]).readline())[5:]
           print 'line=',line
           n=1+line.count('/') #number of domains
           if n>1: #multidomain protein
               domains=line.split('/')
               print 'domains=',domains
               for i in range(n): # `i` will be domain id
                   print 'i='+`i`+',',
                   pdbdomain=dbd+'/new/'+id+`i`+'.pdb'
                   #list of chunks for the domains
                   segments=pat.findall(domains[i])
                   print 'segments=',segments,'=>',
                   error=0
                   if len(segments)==1: #domain made up of a single chunk
                       if generate_pdb(pdb,id[4],pdbdomain,segments[0],endline='TER\n')==0: error=1
                   else: #many chunks make up a domain
                       tmpfile='./junkdomain12345.pdb'
                       #erase just in case
                       os.system('/bin/rm '+pdbdomain+' 2>/dev/null')
                       for segment in segments: #add chunk by chunk
                           print segment,
                           if generate_pdb(pdb,id[4],tmpfile,segment)==0:
                               error=1
                               break
                           os.system('cat '+tmpfile+' >> '+pdbdomain+' && /bin/rm '+tmpfile+' 2>/dev/null')
                       if not error:
                           #armonize residue number among the chunks
                           storeSingleChain(pdbdomain,id[4],tmpfile,'TER\n')
                           os.system('/bin/mv '+tmpfile+' '+pdbdomain)
                   #retrieve domain seq from the stored pdb file of
                   #the domain
                   if not error: chains_domains[id+`i`]=OneLetterSeq(pdbdomain)
                   print '\n'
           else: print 'single domain'
#remove the byproducts of 'pdp' executable
os.system('/bin/rm results.pdp details.pdp '+tmpfile+' 2>/dev/null')

#store new sequences in 'new.fasta' file
pout=open('new.fasta','w')
for id in chains_domains.keys():
    pout.write('> '+id+'\n') #write header
    seq=chains_domains[id]       #sequence
    while seq:
        pout.write(seq[0:60]+'\n') #write in chunks of 60 residues
        seq=seq[60:]
    pout.write('\n')#the seq.pdb format requires a new line
pout.close()
#add new sequences to seq.pdb
#os.system('cat new.fasta >> seq.pdb')

#store new id's in 'new.list' file
list=chains_domains.keys()
list.sort() #order the list
open('new.list','w').write('\n'.join(list))
#add new id's in 'list' file
#os.system('cat new.list>>list')

#end of nightmare :)
sys.exit(0);
