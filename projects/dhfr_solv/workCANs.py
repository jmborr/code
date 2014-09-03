#!/usr/bin/python

import argparse
import sys,os
import numpy
from pdb import set_trace as tr
import ContactMapAnalysis.ContactMapAnalysisAPI as CMAPI
import ContactMapAnalysis.ContactMapAnalysis as CMA

def Sassena(*kargs,**kwargs):
    
    jobList = ['generate input files',
               'translate with VMD',
               ]
    
    JOB = kwargs['job']  #what job to do
    
    if JOB == 'generate input files':
        """generate sassena input files for incoherent scattering"""
        template_0="""<root>
        <sample>
          <structure>
            <file>_PDBFILE_</file>
            <format>pdb</format>
          </structure>
          <framesets>
            <frameset>
              <file>_DCDFILE_</file>
              <format>dcd</format>
            </frameset>
          </framesets>
        </sample>
        <stager>
          <target>system</target>
        </stager>
        <scattering>
            <type>self</type>
            <dsp>
              <type>autocorrelate</type>
              <method>fftw</method>
            </dsp>
            <vectors>
                <type>scans</type>
                <scans>
                  <scan>
                    <from>_QINIT_</from>
                    <to>_QEND_</to>
                    <points>_NQ_</points>
                    <base>
                      <x>1</x>
                      <y>0</y>
                      <z>0</z>
                    </base>
                  </scan>
                </scans>
            </vectors>
            <average>
              <orientation>
                <type>vectors</type>
                <vectors>
                  <type>sphere</type>
                  <algorithm>boost_uniform_on_sphere</algorithm>
                  <resolution>_RESOLUTION_</resolution>
                  <seed>5</seed>
                </vectors>
              </orientation>
            </average>
            <signal>
              <file>_SIGNALFILE_</file>
              <fqt>true</fqt>
              <fq0>true</fq0>
              <fq>true</fq>
              <fq2>true</fq2>
            </signal>
        </scattering>
</root>"""
        options={'_PDBFILE_':kwargs['pdbfile'],
                 '_DCDFILE_':kwargs['dcdfile'],
                 '_QINIT_':0.1,
                 '_QEND_':3,
                 '_NQ_':30,
                 '_RESOLUTION_':500,
                 '_SIGNALFILE_':kwargs['signalfile']
                 }
        
        template = template_0
        for key,val in options.items(): template = template.replace(key,str(val))
        open(kwargs['outfile'],'w').write(template)
    
    elif JOB == 'translate with VMD':
        """Sassena can't proccess AMBER dcd files. We use VMD to create new dcd files. 
        These Sassena can handle. This JOB will create the VMD input files"""
        template_0="""mol new dhp.top type parm7 
animate read dcd _DCDHEADER_.dcd beg 0 end -1 waitfor all 
animate write dcd _DCDHEADER__vmd.dcd beg 0 end -1 waitfor all 
quit"""
        template = template_0
        options={'_TOPFILE_':'dhp.top',
                 '_DCDHEADER_':kwargs['dcdheader'],
                 }
        for key,val in options.items(): template = template.replace(key,str(val))
        open(kwargs['outfile'],'w').write(template)
        #run VMD in the command line: vmd -dispdev none  -e kwargs['outfile']
        
    else:
        print 'Job not recognized. Nothing to do'
 
 
def ContactMap(*kargs,**kwargs):
    
    jobList = ['generate residue ContactMaps',
               'number of atomic contacts',
               'occupancy',
               'ISO- ISO+',
               'RMSF- RMSF+',
               'RMSF- RMSF+ 50%',
               'isopropanol residence times',
               'residence times statistics',
               'cluster rows by K-means',
               'extract centroid frames'
               ]
    JOB = kwargs['job']  #what job to do
    def updateValues(values,fileName,nRows,column):
        iRes = 0
        for line in open(fileName).readlines():
            if line[0]=='#': continue
            if iRes >= nRows: break
            values[iRes] += float( line.split()[column-1] )
            iRes += 1
    if JOB == 'generate residue ContactMaps':
        """calculate the isopropanol-protein residue-based contact maps"""
        #variables defining the ContactMapProtocol
        rowSelection = kwargs['rowSelection'] #protein residues
        colSelection = kwargs['colSelection'] #isopropanols
        cutOff=float( kwargs['cutOff'] )    #atomic cutOff
        pargs = {'rowSelection':rowSelection,'colSelection':colSelection,'cutOff':cutOff,'byres':True}
        protocol = CMAPI.ContactMapProtocol(**pargs)
        #generate the ContactMapList with the trajectory and protocol
        PSFile = kwargs['PSFile']      #topology
        trajfile = kwargs['trajfile']  #trajectory
        cml = CMAPI.GenerateContactMapList(PSFile,trajfile,protocol) #contactMapList
        cml.saveToFile(kwargs['outf'],fmt='HDF5') #save to file in HDF5 format 
        
    elif JOB == 'number of atomic contacts':
        """calculate number of isopropanol-protein contacts for each frame"""
        cml = CMA.ContactMapList()
        cml.loadFromFile(kwargs['inFile'])
        outStr = "#frame #contacts\n"
        frame = 1
        for cm in cml:
            outStr += '%5d %3d\n'%(frame,cm.nnz)
            frame += 1
        open(kwargs['outFile'],'w').write(outStr)
        
    elif JOB == 'occupancy':
        """Calculate the average number of isopropanols in contact with each residue"""
        nRes = 161 #number of DHFR residues plus cofactor plus ligand
        occ = numpy.zeros(1+nRes)
        cml = CMA.ContactMapList()
        cml.loadFromFile( kwargs['inFile'],offSet=int(kwargs['offSet']) )
        occ = cml.Occupancy(1+nRes) #residue index begins at 1, not zero
        outStr = '#residue occupancy\n'
        ires = 1
        for x in occ[1:]: #residue index begin at 1, not zero
            outStr += '%3d %5.2f\n'%(ires,x)
            ires += 1
        open(kwargs['outFile'],'w').write(outStr)
    
    elif JOB == 'ISO- ISO+':
        """Calculate the sets iso- and iso+ by:
        1. Rank each residue by its cumulative occupancy
        2. Calculate a cumulative RMSF_M at 0% isopropanol by increasing the number of residues, ranked by 1."""
        rootDir=kwargs['rootDir']
        nRes = 161
        
        #calculate cumulative occupancy
        cumOcc= [0]*nRes
        for isoVol in ('0.05','0.10','0.15','0.20','0.25','0.50'):
            occFile='%s/%s/solv/Prod/isoOccupancy.dat'%(rootDir,isoVol)
            updateValues(cumOcc, occFile, nRes, 2)
        #write cumulative occupancy
        buf = '#cumulative residue occupancy\n#residue CumulativeOccupancy\n'
        for iRes in range( nRes ):
            buf += '%3d %5.2f\n'%(1+iRes, cumOcc[iRes])
        open(rootDir+'/analysis/cumulativeOccupancy.dat','w').write(buf)
        #rank residues by cumulative occupancy
        cumOcc = numpy.array( cumOcc )
        ranking = numpy.argsort(cumOcc)
        #load RMSF at 0% isopropanol
        rms = [0]*nRes
        isoVol = '0.00'
        rmsFile = '%s/%s/unsolv/rmsf.dat'%(rootDir,isoVol)
        updateValues(rms, rmsFile, nRes, 2)
        rms = numpy.array(rms)
        #save the cumulative RMS
        buf='#Cumulative RMS, residues ranked by cumulative occupancy\n#rank residue cumulativeRMS\n'
        rms = rms[ranking] #shuffle according to cumulative occupancy
        cumulativeRms = 0.00
        for iRank in range( nRes ):
            cumulativeRms += rms[iRank]
            buf += '%3d %3d %9.6f\n'%(1+iRank, 1+ranking[iRank], cumulativeRms)
        open(rootDir+'/analysis/cumulativeRms.dat','w').write(buf)
        
    elif JOB == 'RMSF- RMSF+':
        """Calculate the RMSF for the ISO- and ISO+ sets based on JOB=='ISO- ISO+' for each isopropanol concentration"""
        nRes = 161
        #read the file containing the ranking cumulative occupancy and the cumulative RMSF
        rootDir = kwargs['rootDir']
        cumRMSfile = rootDir+'/analysis/cumulativeRms.dat'
        ranking = numpy.zeros(nRes)
        updateValues(ranking, cumRMSfile, nRes, 2)
        ranking = ranking.astype(int) - 1
        cumRMSF = numpy.zeros(nRes)
        updateValues(cumRMSF, cumRMSfile, nRes, 3)
        #find ISO- and ISO+ sets
        nISOplus = int(kwargs['nISOplus']) #number of residues making up the ISO+ set
        RMSFcut = cumRMSF[-1] - cumRMSF[-1-nISOplus] #cumulative RMSF of the ISO+ set
        isoMinus = ranking[ numpy.where( cumRMSF < RMSFcut ) ]  #residue indexes for ISO-
        isoPlus = ranking[ -nISOplus: ] #residue indexes for ISO+
        #go over the different isopropanol concentrations and calculate the cumulative RMSF for ISO- and ISO+
        buf = '#Cumulative RMSF for ISO- and ISO+\n#iso-conc RMSF_ISO- RMSF_ISO+\n'
        for isoVol in ('0.00', '0.05','0.10','0.15','0.20','0.25','0.30','0.50','1.00'):
            RMSF = numpy.zeros(nRes)
            rmsFile = '%s/%s/unsolv/rmsf.dat'%(rootDir,isoVol)
            updateValues(RMSF, rmsFile, nRes, 2)
            rmsfIsoMinus = RMSF[isoMinus].sum()
            rmsfIsoPlus = RMSF[isoPlus].sum()
            buf += '%s%6.2f %6.2f\n'%(isoVol, rmsfIsoMinus, rmsfIsoPlus)
        open(rootDir+'/analysis/cumulativeRMSisoSets.dat','w').write(buf)
    elif JOB == 'isopropanol residence times':
        """calculate for each residue the residence time of each isopropanol that was in contact"""
        cml = CMA.ContactMapList()
        cml.loadFromFile( kwargs['inFile'],offSet=int(kwargs['offSet']) )
        nRes=161 #number of DHFR residues plus cofactor plus ligand
        rtl = cml.ResidenceTimes(1+nRes)    #create a list of residence times for each residue
        rtl.setDt(kwargs['dt']) #set the unit time, in picoseconds
        rtl.saveToFile(kwargs['outFile'],fmt='HDF5') #save the list to file

    elif JOB == 'residence times statistics':
        """calculate for each residue a few statistics of its list of isopropanol residence times"""
        rtl = CMA.loadResidenceTimesListFromFile( kwargs['inFile'] )
        outStr='#residue maximum average stdev\n'
        for ires in range(1,rtl.n):
            m=0 ; a=0; s=0
            if rtl.data[ires]:
                L=numpy.array(rtl.data[ires])  #list of isopropanol residence times to this particular residue
                m=rtl.dt*max(L) ; a=rtl.dt*numpy.average(L) ; s=rtl.dt*numpy.std(L)
            outStr += '%3d %6d %6.1f %6.1f\n'%(ires,m,a,s)
        open(kwargs['outFile'],'w').write(outStr)

    elif JOB == 'accumulate maximum residence time':
      """add the maximum residence time for the given concentrations"""
      template_file='/projects/research/dhfr_solv/out/isopropanol/%s/solv/Prod/mdprod/isoResidenceTimes.dat'
      ccs=kwargs['concentrations'].split('-') # list of concentrations
      nres=161 # 161 residues
      pps=[open(template_file%cc) for cc in ccs] #list of file pointers
      maxtime=[0.]*nres # cumulative maximum residence times
      eof=False # end-of-file
      while not eof:
        for pp in pps:
            line=pp.readline()
            if not line:
              eof=True
              break
            if line[0]!='#':
              maxtime[int(line.split()[0])-1]+=float(line.split()[1])
      buf='#residue maximum\n'
      for ires in range(nres): buf+='%3d %8.1f\n'%(1+ires,maxtime[ires])
      open(kwargs['outfile'],'w').write(buf)

    elif JOB == 'cluster rows by K-means':
        nRes = 161 #number of DHFR residues plus cofactor plus ligand
        cml = CMA.ContactMapList()
        cml.loadFromFile( kwargs['inFile'],offSet=int(kwargs['offSet']) )
        k = int( kwargs['nCentroids'] )
        clusterid,cdata,nearestMember = cml.ClusterRowsByKmeans(k,nrows=nRes,saturateOccupancy=True,jobName=kwargs['jobName'])
        population = numpy.zeros(k)
        for cid in clusterid: population[cid] += 1
        buf = '#centroid population frame\n'
        for i in range(k):
            buf += '%2d %4d %4d\n'%(1+i,population[i],1+nearestMember[i]) #frames begin with 1
        open(kwargs['jobName']+'.txt','w').write(buf)
        print buf
    
    elif JOB == 'extract centroid frames':
        ptrajTemplate_0="trajin dt10ps.dcd _FRAME_ _FRAME_ 1\ntrajout _CENTROID_ restart"
        offSet = int(kwargs['offSet'])
        for line in open(kwargs['directory']+'/clusterRowsByKmeans.txt').readlines():
            if line[0]=='#': continue
            centroid, population, frame = [ int(x) for x in line.split() ]
            frame += offSet
            ptrajTemplate = ptrajTemplate_0
            ptrajTemplate=ptrajTemplate.replace('_FRAME_','%d'%frame)
            centroidName = 'centroidIsobinding.rep.rst'%centroid
            ptrajTemplate=ptrajTemplate.replace('_CENTROID_',centroidName)
            open(kwargs['directory']+'/centroidIsoBinding.in', 'w').write(ptrajTemplate)
            os.system('cd %s; ptraj11 dhp.top < centroidIsoBinding.in'%kwargs['directory'])
            os.system('cd %s; /bin/mv %s.%s %s'%(kwargs['directory'],centroidName,frame,centroidName))
    else:
        print 'Job not recognized. Nothing to do'



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='service provider for dhfr_solv project')
    parser.add_argument('service',help='requested service, the name of a function defined in this module')
    parser.add_argument('--kargs',help='required arguments of service. Ex: "arg1,arg2,arg3"')
    parser.add_argument('--kwargs',help='optional arguments of service. Ex: "arg1=val1,arg2=val2"')
    args = parser.parse_args()
    service=args.service
    reqargs=[]
    if args.kargs: reqargs=args.kargs.split(',')
    optargs={}
    if args.kwargs: optargs=dict( map( lambda x: x.split('='), args.kwargs.split(',') ) )
    exitCode=0
    try:
        locals()[service](*reqargs,**optargs)
    except:
        sys.stderr.write('Error: service '+service+' not found or error in service\n')
        exitCode=1
    sys.exit(exitCode)
