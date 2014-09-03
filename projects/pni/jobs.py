import argparse,sys,os,numpy,tempfile
from pdb import set_trace as trace

import ContactMapAnalysis as CMA
import ContactMapAnalysisAPI as CMAPI

def ContactMap(*kargs,**kwargs):
    
    jobList = ['generate residue ContactMaps',
               'get contact map',
               'number of atomic contacts',
               'contact map movie',
               ]
    JOB = kwargs['job']  #what job to do

    def filter(csr):
        """filter out contacts that are too close along the sequence
        csr: contact map, a csr_matrix object
        """
        cutOff = 0
        i=0
        for irow in range( len( csr.indptr ) - 1 ):
            for icol in csr.indices[ csr.indptr[irow] : csr.indptr[irow+1] ]:
                if abs( irow - icol ) <= cutOff: csr.data[i] = 0
                i += 1
        csr.eliminate_zeros()
        return csr

    if JOB == 'generate residue ContactMaps':
        """calculate the protein-protein residue-based contact maps"""
        #variables defining the ContactMapProtocol
        rowSelection = kwargs['rowSelection'] #protein residues
        colSelection = rowSelection           #self-contact map
        if 'colSelection' in kwargs.keys(): colSelection = kwargs['colSelection']
        cutOff=float( kwargs['cutOff'] )      #atomic cutOff
        pargs = {'rowSelection':rowSelection,'colSelection':colSelection,'cutOff':cutOff,'byres':True,'filter':filter}
        protocol = CMAPI.ContactMapProtocol(**pargs)
        #generate the ContactMapList with the trajectory and protocol
        PSFile = kwargs['PSFile']      #topology
        trajfile = kwargs['trajfile']  #trajectory
        cml = CMAPI.GenerateContactMapList(PSFile,trajfile,protocol) #contactMapList
        cml.saveToFile(kwargs['outf'],fmt='HDF5') #save to file in HDF5 format 

    elif JOB == 'number of atomic contacts':
        """calculate number of protein-protein contacts for each frame"""
        cml = CMA.ContactMapList()
        cml.loadFromFile(kwargs['inFile'])
        outStr = "#frame #contacts\n"
        frame = 1
        for cm in cml:
            outStr += '%5d %3d\n'%(frame,cm.nnz)
            frame += 1
        open(kwargs['outFile'],'w').write(outStr)
    
    elif JOB == 'get contact map':
        """Retrieve a contact map from the HDF5 and save to
           a file suitable to be loaded by CMVIEW"""
        hdfFile = kwargs['hdfFile']
        index = int(kwargs['frameNumber']) - 1
        cm = CMAPI.GetContactMapAsString(hdfFile,fileFormat='HDF5',index=index,stringFormat='CMVIEW')
        open( kwargs['outFileName'], 'w' ).write(cm)
        
    elif JOB=='contact map movie':
        """create animated GIF using the contact maps"""
        mkd = tempfile.mkdtemp(dir='/tmp') #create temporary directory to store the PNG files
        outFile = mkd+'/junk.dat' #store the current contact map in CMVIEW format
        cmList = CMA.ContactMapList()
        cmList.loadFromFile( kwargs['hdfFile'] )  #get the contact maps to a list
        index = 1
        for cm in cmList:
            cmString = cm.asString(stringFormat='CMVIEW')
            open(outFile,'w').write(cmString)
            #generate the PNG with CMVIEW
            os.system('%s -Y -f %s -I %s/cmview%04d.png'%( kwargs['CMVIEW'], outFile,mkd, index))
            index +=1
        #Use 'convert' from ImageMagick to paste all png's into the animated GIF
        loop = 1 #don't do infinite looping of the animation
        if kwargs['loop'] in ('yes','y','Y'): loop = 0
        os.system( 'convert -delay 100 -loop %d %s/cmview*.png %s'%(loop,mkd,kwargs['gifFile']))
        os.system( '/bin/rm -rf %s'%mkd) #clean-up

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
