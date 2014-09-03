#!/usr/bin/python

import pickle,os,sys,re,copy,pdb
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import inpHand

joink=os.path.join
dirname=os.path.dirname
basename=os.path.basename
isthere=os.path.exists
class Prod:

    def __init__(self,TYP,dir,module_name=None):
        """load Prod directory structure to a dictionary
        All file names with respect to dir
        Prod
        |_dumpf       #file containing all object info
        |_simdir      
        |_TYP
        |_rstf
        |_times
        | |_dt        #elementary time step(picosec)
        | |_ntwx      #number of elementary time steps in between frames
        | |_twx       #time(picosec) in between frames
        |_init
        | |_infile
        | |_topfile
        | |_rstfile
        | |_reffile
        |_prev
          |_crdlist
          |_vellist
          |_rstlist
          |_outlist
          |_first     #first frame, in string representation
          |_last      #last frame,  in string representation
          |_frf       #frame format
       ============Up to here is minimal Prod directory structure
        |_strip       #removed atoms from mask (usually [:WAT,Na+])
          |_maskl
          |_nat       #number of remaning atoms
          |_topfile
          |_reffile
          |_crdfile
          |_velfile
          |_centroid  #centroid by RMS among CA atoms
          |_average   #structure from averaging CA atoms
          """
        self.dumpf='pickle.dump'
        dumpf=joink(dir,self.dumpf)
        if isthere(dumpf): #create object from the dump file
            if module_name:
                self.armonize_module_name(module_name,dumpf)
            obj=pickle.load(open(dumpf,'r'))
            if obj.simdir!=dir: obj.updateSimdir(dir)
            self.__dict__=copy.deepcopy(obj.__dict__)
            self.armonize_prev_outlist()
            self.armonize_rst_outlist()
            self.armonize_strip()
            return None
        self.TYP=TYP
        self.simdir=dir
        self.init={'infile' :'init/'+TYP+'.in',
                   'topfile':'init/'+TYP+'.top',
                   'rstfile':'init/'+TYP+'.rst',
                   'reffile':'init/'+TYP+'.ref',
                   }
        self.prev={'crdlist':[],'vellist':[],'rstlist':[],'outlist':[] }
        #crdlist, vellist, rstlist, outlist
        cmds={'crdlist':'cd '+dir+' && ls -1 prev|grep -P "'+TYP+'.+\.crd.*"',
              'vellist':'cd '+dir+' && ls -1 prev|grep -P "'+TYP+'.+\.vel.*"',
              'outlist':'cd '+dir+' && ls -1 prev|grep -P "'+TYP+'.+\.out.*"',
              'rstlist':'cd '+dir+' && ls -1 prev|grep -P "'+TYP+'.+\.rst.*"',
            }
        for key,cmd in cmds.items():
            #print key,cmd
            self.prev[key]=chomp(os.popen(cmd).readlines())
            junkl=[ 'prev/'+x for x in self.prev[key] ]
            self.prev[key]=junkl
            #print '\t%s %s'%(self.prev[key][0],self.prev[key][-1])

        #rstf,first,last,frf
        rstf='rst.list'
        list=open(joink(dir,rstf),'r').readlines()
        match=re.compile(TYP+'\.(\d+)\.rst').search(list[-1])
        prev=self.prev
        prev['last']=match.group(1)
        prev['frf']='%0'+str(len(prev['last']))+'d'
        prev['first']=prev['frf']%(1)
        self.timeBetweenFrames()
        #import strip, if extant
        strip=joink(dir,'strip')
        if isthere(strip):
            self.strip={}
            cmds={'crdfile':'cd '+dir+' && ls -1 strip/'+TYP+'*.crd*',
                  'velfile':'cd '+dir+' && ls -1 strip/'+TYP+'*.vel*',
                  }
            for key,cmd in cmds.items():
                print cmd
                self.strip[key]=os.popen(cmd).readline().strip()
        self.dump()

    def armonize_strip(self):
        if 'strip' not in self.__dict__.keys():
            sys.stderr.write('ERROR no "strip" attribute.'+\
                             'Do you want to add? (y/n)')
            if sys.stdin.readline()[0].lower()=='y':
                self.strip={}
            else:
                return None
        if 'crdfile' not in self.strip.keys():
            self.strip['crdfile']='strip/'+self.TYP+'.00001_50000.crd.gz'
        if 'velfile' not in self.strip.keys():
            self.strip['velfile']='strip/'+self.TYP+'.00001_50000.vel.gz'
        self.dump()
        return self
        
    def armonize_prev_outlist(self):
        """unzip out files within prev"""
        newlist=[]
        for f in self.prev['outlist']:
            if '.gz' in f:
                newlist.append(f.replace('.gz',''))
                os.system('gunzip '+joink(self.simdir,f))
            else:
                newlist.append(f)
        self.prev['outlist']=newlist
        self.dump()
        
    def armonize_rst_outlist(self):
        """unzip rst files within prev"""
        newlist=[]
        for f in self.prev['rstlist']:
            if '.gz' in f:
                newlist.append(f.replace('.gz',''))
                os.system('gunzip '+joink(self.simdir,f))
            else:
                newlist.append(f)
        self.prev['rstlist']=newlist
        self.dump()
        
    def armonize_module_name(self,module_name,dumpf):
        """replace module name in the pickle dump if necessary"""
        line=open(dumpf,'r').readline().strip()
        match=re.compile('\(c(.*)').match(line)
        if match:
            name=match.group(1)
            if module_name!=name:
                print name,module_name
                buf=open(dumpf,'r').read().replace(name,module_name,1)
                open(dumpf,'w').write(buf)
        else:
            sys.stderr.write('ERROR could not find module name\n')
            sys.exit(1)

    def dump(self):
        """store object in the production directory"""
        dumpf=joink(self.simdir,self.dumpf)
        pickle.dump( self,open(dumpf,'w'), True )

    def timeBetweenFrames(self):
        """print time in between frames, in picoseconds"""
        file=joink(self.simdir,self.init['infile'])
        buf=''.join(open(file,'r').readlines())
        self.times={} ; pt=self.times
        pt['dt']=float(re.compile(' dt=(\d+\.*\d*)').search(buf).group(1))
        pt['ntwx']=int(re.compile('ntwx=(\d+)').search(buf).group(1))
        pt['twx']=pt['dt']*pt['ntwx']
        return {'self':self,'twx':pt['twx']}

    def updateSimdir(self,simdir):
        """update simdir directory"""
        if simdir[-5:]!='/Prod':
            sys.stderr.write('simdir must finish in /Prod\n')
            raise SystemExit,1
        self.simdir=simdir
        self.dump()
        
    def importStripTopCrd(self,topfile,reffile):
        """Import stripped topology and reference file"""
        from amber.amber10 import top
        buf='Could not find %s\n'
        for file in (topfile,reffile):
            if not isthere(file):
                sys.sdterr.write(buf%file)
                sys.exit(1)
        if 'strip' not in self.__dict__.keys():
            self.strip={}
            os.system('/bin/mkdir -p '+joink(self.simdir,'strip'))
        self.strip['topfile']='strip/'+self.TYP+'.top'
        file=joink(self.simdir,self.strip['topfile'])
        os.system('/bin/cp '+topfile+' '+file)
        self.strip['reffile']='strip/'+self.TYP+'.rst'
        file=joink(self.simdir,self.strip['reffile'])
        os.system('/bin/cp '+reffile+' '+file)
        self.strip['nat']=top(topfile).numberAtoms()
        self.dump()

    def count_frames(self,mode='solvated'):
        from amber10 import ptraj
        print 'counting frames...'
        topfile=joink(self.simdir,self.init['topfile'])
        crdlist=self.prev['crdlist']
        nframes=0
        if mode=='unsolvated':
            try:
                topfile=joink(self.simdir,self.strip['topfile'])
                crdlist=[self.strip['crdfile'],]
            except:
                sys.stderr.write('ERROR: could not find topfile/crdlist\n')
                return 0
        if isthere(topfile):
            output=junkName()
            pt=ptraj(topfile,[joink(self.simdir,x) for x in crdlist])
            pt.go(outfile=output)
            if not isthere(output):
                sys.stderr.write('ERROR: no ptraj.output (%s)\n'%(output))
                sys.exit(1)
            pattern=re.compile('Successfully read in (\d+) sets')
            match=pattern.search( open(output).read() )
            if match:
                nframes=int(match.group(1))
                print 'Found %d frames'%(nframes,)
                print 'If number of frames incorrect, it is either\n'+\
                      ' an error in the topology or trajectory files'
            else:
                buf='ERROR counting frames!\nExamine %s.\n'+\
                     'Remove %s\n'%(output,)
                sys.stderr.write(buf)
                sys.exit(1)
            os.system('/bin/rm %s'%(output,))
        else:
            sys.stderr.write(topfile+' not found\n')
            sys.exit(1)
        #sys.stdout.write('Should we abort ?(y/n):')
        #if sys.stdin.readline()[0].lower()=='n': sys.exit(1)
        return nframes
        
    def doStrip(self,maskl=[':WAT',':Na+']):
        """merge trajectories and strip"""
        from amber.amber10 import top,ptraj
        print "merge trajectories and strip of solvent and ions"""
        if 'strip' not in self.__dict__.keys():
            self.strip={}
            os.system('/bin/mkdir -p '+joink(self.simdir,'strip'))
        self.strip['maskl']=maskl
        init=self.init
        prev=self.prev
        dowegzip=False
        if '.gz' not in prev['crdlist'][0]: #trajectories are not zipped!
            #sys.stdout.write('Will gzip the hydrated trajectories?(y/n):')
            #if sys.stdin.readline().strip().lower()[0]=='y': dowegzip=True
            dowegzip=True
        t=joink(self.simdir,init['topfile'])  #topology
        #strip from coordinate file
        crdfile='strip/'+self.TYP+'.'+prev['first']+'_'+prev['last']+'.crd'
        crds=[ joink(self.simdir,crdf) for crdf in prev['crdlist'] ]
        ptr=ptraj(t,crds)
        ptr.strip(maskl)
        ptr.go(outcrd=joink(self.simdir,crdfile))
        os.system('gzip -f '+joink(self.simdir,crdfile))
        self.strip['crdfile']=crdfile+'.gz'
        
        #strip from velocity file        
        velfile='strip/'+self.TYP+'.'+prev['first']+'_'+prev['last']+'.vel'
        vellist=[ joink(self.simdir,crd) for crd in prev['vellist'] ]
        ptr=ptraj(t,vellist)
        ptr.strip(maskl)
        ptr.go(outcrd=joink(self.simdir,velfile))
        os.system('gzip -f '+joink(self.simdir,velfile))
        self.strip['velfile']=velfile+'.gz'
        #tar and compress files        
        if dowegzip:
            for ftype in ('crd','vel','rst','out'):
                tarf='%s.%s.tar'%(self.TYP,ftype)
                cmd='cd %s/prev && '%(self.simdir,tarf)
                cmd+='for f in `ls -1|grep %s`;do tar uf %s $f;done && '%(ftype,tarf)
                cmd+='gzip %s && '%(tarf,)
                cmd+='for f in `ls -1|grep %s|grep -v tar`;do /bin/rm $f;done'%(ftype,tarf)
                print 'DO NOT FORGET TO RUN:\n'+cmd
        self.dump() #update the pickle dump file

    def getAvg(self,condense=True):
        """find averages from *.out files
        condense: concatenate all *.out files into one
        Prod
        |_prev
          |_avg       #averages
            |_EAMBER (non-restraint)
            |_TEMP(K)
            |_EKtot
            |_RESTRAINT
            |_Etot
            |_1-4 NB
            |_DIHED
            |_EHBOND
            |_VDWAALS
            |_EPtot
            |_1-4 EEL
            |_PRESS
            |_ANGLE
            |_EELEC
            |_BOND
        """
        import re
        p=re.compile('\s*=\s*')
        q=re.compile('([\w,\-,\.,\d,\(,\),\s]+)=([\-,\.,\d]+)')
        avg={}
        N=0
        for f in self.prev['outlist']: #cycle over simulations
            if not N%100: print '%s %d'%( f,len(avg.keys()) )
            pf=open(joink(self.simdir,f)) ; l=pf.readline()
            while 'A V E R A G E S   O V E R ' not in l: l=pf.readline()
            l=pf.readline().strip()
            while not l: l=pf.readline().strip() #blank lines
            while '=' in l:
                for key,val in q.findall( p.sub('=',l) ):
                    key=key.strip()
                    if key not in avg.keys():
                        avg[key]=float(val)
                    else:
                        avg[key]+=float(val)
                l=pf.readline().strip()
            N+=1
        N=len(self.prev['outlist'])
        for key in avg.keys(): avg[key]/=N
        self.prev['avg']=avg
        if condense:
            outf=junkName() ; N=0
            for f in self.prev['outlist']:
                os.system('echo %s >> %s'%(self.prev['outlist'],outf))
                ff=joink(self.simdir,f)
                os.system('cat %s >> %s'%(ff,outf))
                N+=1
            f='prev/'+self.TYP+'.out'
            self.prev['outfile']=f
            cmd='/bin/mv '+outf+' '+joink(self.simdir,f)
            os.system(cmd)
        print 'average temperature %f'%(self.prev['avg']['TEMP(K)'])
        self.dump()
        
    def rmsJump(self,outcrd=False):
        """find all-atom rmsd jump between consecutive frames
        _strip
         |_rmsJump
         |_maxrmsJump
         |_fit2prev (optional trajectory)
        """
        from amber.amber10 import ptraj
        print 'rmsJump: all-atom rmsd jump between consecutive frames'
        crdf=self.strip['crdfile']
        rmsdf=crdf.replace('.crd.gz','.rmsJump.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        ptr.rms(args=' previous out '+joink(self.simdir,rmsdf))
        if outcrd:
            print "output fit trajectory"
            crdout=crdf.replace('.crd.gz','.fit2prev.crd')
            ptr.go(outcrd=joink(self.simdir,crdout))
            os.system('gzip -f '+joink(self.simdir,crdout))
            self.strip['fit2prev']=crdout+'.gz'
        else:
            ptr.go()
        badframe=0 ; maxrms=0.0 #find maximum rmsJump
        for line in open(joink(self.simdir,rmsdf),'r').readlines():
            frame,rms=line.split(); rms=float(rms)
            if rms>maxrms:
                maxrms=rms
                badframe=frame
        print badframe,str(maxrms) #report worst case
        self.strip['rmsJump']=rmsdf
        self.strip['maxrmsJump']=maxrms
        self.dump()#update the pickle dump file

    def rmsRef(self):
        """find all-atom rmsd with respect to reference structure
        _strip
         |_rmsRef
         |_maxrmsRef
        """
        print "find all-atom rmsd with respect to reference structure"
        from amber.amber10 import ptraj
        crdf=self.strip['crdfile']
        rmsdf=crdf.replace('.crd.gz','.rmsRef.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        ptr.reference(joink(self.simdir,self.strip['reffile']))
        ptr.rms(args=' reference out '+joink(self.simdir,rmsdf))
        ptr.go()
        badframe=0 ; maxrms=0.0 #find maximum rmsJump
        for line in open(joink(self.simdir,rmsdf),'r').readlines():
            frame,rms=line.split(); rms=float(rms)
            if rms>maxrms:
                maxrms=rms
                badframe=frame
        print badframe,str(maxrms) #report worst case
        self.strip['rmsRef']=rmsdf
        self.strip['maxrmsRef']=maxrms
        self.dump()#update the pickle dump file
        
    def rmsFirst(self):
        """find all-atom rmsd with respect to first structure
        _strip
         |_rmsFirst
         |_maxrmsFirst
        """
        print "find all-atom rmsd with respect to first structure"
        from amber.amber10 import ptraj
        crdf=self.strip['crdfile']
        rmsdf=crdf.replace('.crd.gz','.rmsFirst.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        ptr.reference(joink(self.simdir,self.strip['reffile']))
        ptr.rms(args=' first out '+joink(self.simdir,rmsdf))
        ptr.go()
        badframe=0 ; maxrms=0.0 #find maximum rmsJump
        for line in open(joink(self.simdir,rmsdf),'r').readlines():
            frame,rms=line.split(); rms=float(rms)
            if rms>maxrms:
                maxrms=rms
                badframe=frame
        print badframe,str(maxrms) #report worst case
        self.strip['rmsFirst']=rmsdf
        self.strip['maxrmsFirst']=maxrms
        self.dump()#update the pickle dump file

    def rmsCentr(self):
        """find all-atom rmsd with respect to reference structure
        _strip
         |_rmsCentr
         |_maxrmsCentr
         |_centrNframe #frame number for the centr (counting starts at 1)
        """
        print "find all-atom rmsd with respect to reference structure"
        if 'centroid' not in self.strip.keys(): self.centr_and_avg()
        from amber.amber10 import ptraj
        centrf=junkName() #store temporary restart file for the centroid
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,self.strip['centroid']))
        ptr.trajout(centrf,args=' restart')['self'].go()
        centrf+='.1'
        crdf=self.strip['crdfile']
        rmsdf=crdf.replace('.crd.gz','.rmsCentr.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        ptr.reference(centrf)
        ptr.rms(args=' reference @CA out '+joink(self.simdir,rmsdf))
        ptr.go()
        os.system('/bin/rm '+centrf) #remove the temporary directory
        badframe=0 ; maxrms=0.0 ; minrms=100.0 ; minframe=1 ; N=1
        for line in open(joink(self.simdir,rmsdf),'r').readlines():
            frame,rms=line.split(); rms=float(rms)
            if rms>maxrms:  #find maximum rms to centroid
                maxrms=rms
                badframe=frame
            if rms<minrms:  #find minimum rms to centroid. Should be zero
                minrms=rms
                minframe=N
            N+=1 #next frame
        print badframe,str(maxrms) #report worst case
        self.strip['rmsCentr']=rmsdf
        self.strip['maxrmsCentr']=maxrms
        self.strip['centNframe']=minframe  #;print self.strip['centNframe']
        self.dump()#update the pickle dump file


    def centr_and_avg(self):
        """find centroid and average structures (CA atoms only)
        _strip
         |_centroid
         |_average
        """
        print 'finding centroid and average structures'
        from amber.amber10 import ptraj
        crdf=self.strip['crdfile']
        self.strip['centroid']=crdf.replace('.crd.gz','.centr.crd')
        self.strip['average']=crdf.replace('.crd.gz','.avg.crd')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        prefix=junkName()
        args='out '+prefix+' representative amber '+\
              'average amber means clusters 1 rms sieve 100 @CA'
        ptr.cluster(args=args)
        ptr.go()
        #cleanup
        file=joink(self.simdir,self.strip['centroid'])
        cmd='/bin/mv '+prefix+'.rep.c0 '+file
        os.system(cmd)
        file=joink(self.simdir,self.strip['average'])
        cmd='/bin/mv '+prefix+'.avg.c0 '+file
        os.system(cmd)
        cmd='/bin/rm '+prefix+'*'
        os.system(cmd)
        self.dump() #update the object dump file
        return {'self':self,}

    def fit2centr(self):
        """fit trajectory to centroid
        _strip
         |_fit2centr
        """
        print 'fitting trajectory to centroid'
        if 'centroid' not in self.strip.keys(): self.centr_and_avg()
        from amber.amber10 import ptraj
        centrf=junkName() #store temporary restart file for the centroid
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,self.strip['centroid']))
        ptr.trajout(centrf,args=' restart')['self'].go()
        centrf+='.1'
        crdf=self.strip['crdfile']
        outf=crdf.replace('.crd.gz','.fit2centr.crd')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        ptr.reference(centrf)
        ptr.rms(args='reference @CA')
        ptr.go(outcrd=joink(self.simdir,outf))
        os.system('/bin/rm '+centrf)
        os.system('gzip -f '+joink(self.simdir,outf))  ;  outf+='.gz'
        self.strip['fit2centr']=outf
        self.dump() #update the object dump file
        return {'self':self,}

    def Bfactors(self,args=None):
        """find Bfactors from trajectory fitted to centroid
        _strip
         |_bfact
        """
        print 'calculating Bfactors by atom'
        if 'fit2centr' not in self.strip.keys(): self.fit2centr()
        from amber.amber10 import ptraj
        crdf=self.strip['fit2centr']
        outf=crdf.replace('.crd.gz','.bf.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        if not args: args='out '+joink(self.simdir,outf)+' byatom bfactor'
        ptr.atomicfluct(args=args)
        ptr.go()
        self.strip['bfact']=outf
        self.dump() #update the object dump file
        return {'self':self,}
    
    def BfactorsByRes(self,args=None):
        """find Bfactors from trajectory fitted to centroid
        _strip
         |_bfactres
        """
        print 'calculating Bfactors by residue'
        if 'fit2centr' not in self.strip.keys(): self.fit2centr()
        from amber.amber10 import ptraj
        crdf=self.strip['fit2centr']
        outf=crdf.replace('.crd.gz','.bfres.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        if not args: args='out '+joink(self.simdir,outf)+' byres bfactor'
        ptr.atomicfluct(args=args)
        ptr.go()
        self.strip['bfactres']=outf
        self.dump() #update the object dump file
        return {'self':self,}

    def BfactorsByMask(self,args=None,mask='*'):
        """find Bfactors by mask from trajectory fitted to centroid
        NOTE: this by mask calculates Bfactors for each atom in the mask,
              it is not an average for all atoms in the mask
        _strip
         |_bfactmask
        """
        print 'calculating Bfactors for each atom in mask',mask
        if 'fit2centr' not in self.strip.keys(): self.fit2centr()
        from amber.amber10 import ptraj
        crdf=self.strip['fit2centr']
        outf=crdf.replace('.crd.gz','.bfmask.dat')
        t=joink(self.simdir,self.strip['topfile'])
        ptr=ptraj(t,joink(self.simdir,crdf))
        ptr.strip(entities=['!'+mask,])
        if not args: args='out '+joink(self.simdir,outf)+' byatom bfactor'
        ptr.atomicfluct(args=args)
        ptr.go()
        self.strip['bfactmask']=outf
        self.dump() #update the object dump file
        return {'self':self,}
 
if __name__=='__main__':
    
    ih=inpHand('Usage: Prod2.py',
               ' -a _RA_prodd absolute path to production directory',
               ' -b _R_TYP protein ID',
               )
    ih.parse(locals(),sys.argv)

    pr=Prod(TYP,prodd)
    #pr.count_frames(mode='solvated') #check number of frames
    print 'getAvg' ;  pr.getAvg(condense=True)
    print 'doStrip' ; pr.doStrip()
    #pr.rmsJump()
    #pr.rmsRef()
    #pr.centr_and_avg()
    #pr.rmsCentr()
    #pr.rmsFirst()
    #pr.fit2centr()
    #pr.Bfactors()
    for (key,val) in  pr.strip.items(): print key,val
