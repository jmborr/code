#!/usr/bin/python

import argparse
import sys,os
import h5py

from pdb import set_trace as trace

def Sassena(*kargs,**kwargs):
    
    jobList = ['add sassena version',
               'generate mini file',
               ]
    JOB = kwargs['job']  #what job to do
    
    if JOB == 'add sassena version':
        f=h5py.File( kwargs['infile'] )
        f.attrs["sassena_version"] = kwargs['version']
        f.close()
        
    elif JOB == 'generate mini file':
        """generate a small file for testing"""
        f=h5py.File( kwargs['infile'], 'r' )
        g=h5py.File( kwargs['outfile'], 'w' )
        #copy attributes
        for key,val in f.attrs.items(): g.attrs[key]=val
        #copy a reduced version of the datasets
        nkeep = 5
        mkeep = 8
        for key in f.keys():
            if key=='fqt':
                g[key] = f[key][0:nkeep,0:mkeep,:]
            else: 
                g[key] = f[key][0:nkeep]
        f.close()
        g.close()
        
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
