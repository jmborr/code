#!/usr/bin/python
'''
We expect a log file name with this structure:
                    root-directory+globbed+specific-dir+/logname
 globbed can be /x/xxxxx or /xxxxx
 
'''
import os,sys,re,copy
from inputArgs.inputArgs import inpHand,deglobb,extractHeader
from utilities.small_utilities import Bye,chomp,junkName
from secondary.seqdat2ss import translateSec
from random import randint

ih=inpHand('Usage: success_fail.py [options]',
           ' -a _A_logn0 globbed log file name (def: ./xxxxx/*log)',
           ' -b _A_outf report file, besides success.list and fail.list (def:root-directory/success_fail.txt))',
           ' -c _A_successf success file name (def: root-directory/success.list)',
           ' -d _A_failf fail file name (def: root-directory/fail.list)',
           ' -e __prefix prefix to append to default report, success, and fail files(def:None)',
           ' -f _A_origListf original list of jobs, to complete fail.list',
           )
ih.parse(locals(),sys.argv)

currd=os.getcwd()
if not logn0: logn0=os.path.join(currd,'xxxxx/*log')
if not prefix: prefix=''


#get the root and especific directory, remember that general log file
#name is  rootd+globbed+specd+/logn
isglobbed=False
rootd=''
specd=''
subdir=''
logn=os.path.basename(logn0) #basename
dirn=os.path.dirname(logn0)  #dirname
pat={'x/xxxxx':re.compile('(\S+)/x/xxxxx'),'xxxxx':re.compile('(\S+)/xxxxx')}
m=logn0.find('/x/xxxxx')
if m>=0:
    isglobbed=True
    rootd=logn0[0:m]#root directory
    subdir=dirn[m:] #globbed directory relative to rootd
    specd=dirn[m+len('/x/xxxxx'):] #everything after globbed and before /logn
else:
    m=logn0.find('/xxxxx')
    if m>=0:
        isglobbed=True
        rootd=logn0[0:m]  #Bye('rootd='+rootd)
        subdir=dirn[m:]
        specd=dirn[m+len('/xxxxx'):] #;Bye('specd='+specd)
    else:
        rootd=os.path.dirname(logn0)

#Obtain log file names    
os.chdir(rootd)
cmd='/usr/bin/find . -name "'+logn+'"|grep "'+specd+'"'#;Bye(cmd)
logfs=chomp(os.popen(cmd).readlines())
if not logfs: Bye('I could not find any log files')
#Bye('\n'.join(logfs))

if not successf: successf=os.path.join(rootd,prefix+'success.list') #success file
successp=open(successf,'w')
if not failf: failf=os.path.join(rootd,prefix+'fail.list')          #fail file
failp=open(failf,'w')
faill=[]
observed=[]
for logf in logfs:    
    if isglobbed: listentry=extractHeader(subdir,os.path.dirname(logf))
    else: listentry=logf
    cmd='tail -1 '+logf+'|cut -d \' \' -f 4' #retrieve the numeric code "X" for the "exit mode = X" line
    exitmode=0
    try:
        exitmode=int( os.popen(cmd).readline() )           
    except ValueError:
        sys.stderr.write(logf+' last line is not "exit mode = X"\n')
        continue
    if exitmode:
        faill.append(logf)
        failp.write(listentry+'\n')  #;print logf
    else:
        successp.write(listentry+'\n')
    observed.append(listentry)

if origListf: #check if there are missing log files
    origList=chomp(open(origListf,'r').readlines())    
    for listentry in origList:
        if listentry not in observed: failp.write(listentry+'\n')
    
failp.close()
successp.close()

buf='' #buffer
#report ten random cases of failed
nfail=len( open(failf,'r').readlines() )
alljobs=len(logfs)
if origListf: alljobs=len(origList)
buf+=`nfail`+' failed jobs out of '+`alljobs`+'\n'
if nfail:
    for i in range(min(nfail,10)):
        n=randint(0,min(nfail,10)-1)
        logf=faill[n]
        logfb=os.path.basename(logf)
        buf+='EXAMPLE %2d ##################  '%(i+1)+logfb+'  ##################\n'
        buf+=''.join(os.popen('cat '+logf).readlines())


#write report
if not outf: open(prefix+'success_fail.txt','w').write(buf)
else:        open(outf,'w').write(buf)
os.chdir(currd)
sys.exit(0)
