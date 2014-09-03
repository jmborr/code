#!/usr/bin/python
import os,sys,re
from utilities.codedir import codedir
from utilities.small_utilities import my_dict_sort,chomp,Bye
from time import sleep,time
from stat import ST_MTIME
from inputArgs.inputArgs import inpHand
from random import shuffle

nnodes=1154 #number of nodes in cluster

def jobs():
    report=''    
    lines=chomp(os.popen('qstat|grep jose').readlines())
    runninglines=os.popen('qstat -n -u jose|grep " R " -A 1').readlines()
    nrl=len(runninglines)
    for i in range(len(lines)):
        items=lines[i].split()
        id=items[0]
        name=items[1]
        runningtime=items[3]
        if runningtime=='0':runningtime='00:00:00'
        status=items[4]
        node=''
        if status=='R':
            for i in range(nrl):
                if id in runninglines[i]: node=runninglines[i+1][3:10]            
        line='%8s %15s %15s %s %8s\n'%(node,id,name,status,runningtime)
        report=report+line
    return report[0:-1] #trick to remove last '\n'


def prunejobIDs(name,state=''):
    pattern=re.compile(name) #;Bye(name)
    if state: state=' '+state+' ' #add blank spaces to each side
    jobIDs=[]
    lines=os.popen('qstat |grep jose|grep -v " C "').readlines()
    for line in lines:
        jobID=line.split()[0]
        if name:
            if pattern.search(line):
                if not state or (state and state in line):
                    jobIDs.append(jobID)
        else:
            jobIDs.append(jobID)
    return jobIDs

    
def py_qdel(name,state=''):

    jobIDs=prunejobIDs(name,state=state)
    if jobIDs:
        nj=20 #chunk of jobID's to be removed in a single qdel call
        #invoque qdel in chunks of nj jobs
        while jobIDs:
            buf=' '.join(jobIDs[0:nj])
            os.system('qdel '+buf)
            sleep(3.0) #rest in between consecutive qdel calls
            jobIDs=jobIDs[nj:]


def py_qalter(name,priority,state='Q'):
    
    if int(priority)<-1024: priority='-1024'
    elif int(priority)>1023: priority='1023'
    
    jobIDs=prunejobIDs(name,state=state)
    if jobIDs:
        nj=20
        while jobIDs:
            buf=' '.join(jobIDs[0:nj])
            cmd='qalter -p '+priority+' '+buf 
            os.system(cmd)
            jobIDs=jobIDs[nj:]


def load_people():
    status_types=('ALL','RUNNING','QUEUED','HELD','COMPLETED','SUSPENDED')
    loads={'RUNNING':{},'ALL':{},'QUEUED':{},'HELD':{},'COMPLETED':{},'SUSPENDED':{}}
    #initialize loads for jose
    for type in status_types: loads[type]['jose']=0
    all=os.popen('qstat|grep -v -e "----------------"|grep -v Name').readlines()
    for line in all:
        pieces=line.split()
        user=pieces[2]
        for type in status_types:
            if user not in loads[type].keys(): loads[type][user]=0

        status=pieces[4]
        if status!='C': loads['ALL'][user]+=1
        for type in status_types[1:]:
            if status==type[0]: loads[type][user]+=1

    #sort by number of running jobs
    items=[(v,k) for k,v in loads['RUNNING'].items()]
    items.sort()
    items.reverse()
    
    for type in status_types: 
        loads[type]['total']=sum(loads[type].values())
        for key in loads[type].keys(): loads[type][key]='%5d'%(loads[type][key])
        
    buf= '           RUNNING    ALL   QUEUED    HELD   COMPLETED  SUSPENDED\n'
    buf+='%10s'%('jose')+'    '+loads['RUNNING']['jose']+' '+loads['ALL']['jose']+'    '+loads['QUEUED']['jose']+'   '+loads['HELD']['jose']+'       '+loads['COMPLETED']['jose']+'      '+loads['SUSPENDED']['jose']+'\n'
    buf+='%10s'%('total')+'    '+loads['RUNNING']['total']+' '+loads['ALL']['total']+'    '+loads['QUEUED']['total']+'   '+loads['HELD']['total']+'       '+loads['COMPLETED']['total']+'      '+loads['SUSPENDED']['total']+'\n'

    for (v,key) in items:
        if key not in ('total','jose'):
            buf+='%10s'%(key)+'    '+loads['RUNNING'][key]+' '+loads['ALL'][key]+'    '+loads['QUEUED'][key]+'   '+loads['HELD'][key]+'       '+loads['COMPLETED'][key]+'      '+loads['SUSPENDED'][key]+'\n'
    return buf
        
def run_command(cmd):
    for i in range(1,nnodes+1):
        node='cng'+'%04d'%(i)
        os.system("/usr/bin/rsh "+node+" '"+cmd+"'")
        sleep(0.5)

def report_load_bnodes():
    bnodes={}
    p=re.compile('load average\:\s(\d\.\d+)')
    nodenumbers=range(1,15)
    shuffle(nodenumbers)
    for i in nodenumbers: bnodes['b'+`i`]=4.00
    for bnode in bnodes.keys():
        cmd='/usr/bin/rsh '+bnode+' "uptime" 2>/dev/null|grep "load average:"' #; print cmd
        line=os.popen(cmd).readline()
        if p.search(line): bnodes[bnode]=float(p.search(line).group(1))
    #returns a list of pairs
    return my_dict_sort(bnodes,'value',mode='ascending')

def rshb():
    os.system('/usr/bin/rsh '+report_load_bnodes()[0][0])

listFunc=(jobs,py_qdel,load_people,run_command,report_load_bnodes,rshb)

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: run_command.py',
               ' -a __fun function to call ("rshb","load_people","report_load_bnodes")',
               ' -b __arglist line containing all arguments, if any, to previous function')
    ih.parse(locals(),sys.argv)

    if fun in ('jobs','rshb','load_people','report_load_bnodes'): print eval(fun+'()')


