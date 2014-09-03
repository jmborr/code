#!/usr/bin/python

import pdb
import os,sys,re,numpy
from copy import deepcopy
from tempfile import mkstemp,mkdtemp
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
pwrite=sys.stdout.write
joink=os.path.join; jk=os.path.join
isthere=os.path.exists
basename=os.path.basename

my_project='DiffSphere'                      #substitute 'my_project'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')

def Volino(*kargs,**kwargs):
    """Checks Figure 1 of Volino80 paper
    kwargs:
        X[=10.0]
        dx[=0.01]
        task[=
            xnl  #check the xnl by plotting the curves in xmgrace and looking at the roots
            A0l  #check Anl for n=0 and l from 1 to 5
            A1l  #check Anl for n=1 and l from 0 to 2
            lA1l #check (2l+1)*Anl for n=1 and l from 0 to 2 
            jl   #check spherical bessel functions jl for l running from 0 t0 3
               ]
    """
    
    from numpy import array,square
    from scipy.special import sph_jn  #spherical bessel and its derivative up to order n

    xnl={} #some of the coefficients needed 
    xnl[0]={0:0.0,      1:2.081576, 2:3.342094, 3:4.514100, 4:5.646704, 5:6.756456 }
    xnl[1]={0:4.493409, 1:5.940370, 2:7.289932, 3:8.583755, 4:9.840446, 5:11.070207}

    def arrayx(X,dx):
        """generate a 1D grid from dx to X with spacing dx"""
        x0=dx; x=[]
        while x0<X: 
            x.append(x0); x0+=dx
        return array(x)

    def jl(l,x):
        """generate a numpy array of spherical bessel function values given numpy array x"""
        z=[]
        for xi in x: z.append( sph_jn(l,xi)[0][-1] )
        return array(z)

    def diff(l,x):
        "the zeros of the returned curve represents the xnl for a given l and running n"
        if l==0: return jl(1,x)
        return l*jl(l,x)-x*jl(l+1,x)

    def A(n,l,x):
        """Volino's coefficients"""
        z=square(xnl[n][l])
        return 6*z/(z-l*(l+1)) * square( (x*jl(l+1,x) - l*jl(l,x)) / (square(x)-z) )

    #create 1D grid x
    try: X=float(kwargs['X'])
    except: X=10.0 #default value
    try: dx=float(kwargs['dx'])
    except: dx=0.01 #default value
    x=arrayx(X,dx)
    
    try: task=kwargs['task']
    except: sys.stderr.write('ERROR: argument "task" is required\n'); sys.exit(1)
    
    #pdb.set_trace()
    
    def VolinoHelper(start,end,x,y):
        """y is a function of 'l', and returns a numpy array"""
        for l in range(start,end):
            z=y(l)
            for i in range( len(x) ): print x[i],z[i]
            print '&'
       
    if task=='xnl':
        def diffl(l): return diff(l,x)
        VolinoHelper( 0,6,x,diffl )  #check the xnl by plotting the curves in xmgrace and looking at the roots
    elif task=='A0l':
        def A0(l): return A(0,l,x)  #returns a numpy array
        VolinoHelper( 1,6,x,A0 )    #check Anl for n=0 and l from 1 to 5
    elif task=='lA1l':
        def lA1(l):  return (2*l+1)*A(1,l,x)
        VolinoHelper( 0,3,x,lA1 )    #check (2l+1)*Anl for n=1 and l from 0 to 2        
    elif task=='A1l':
        def A1(l):  return A(1,l,x)
        VolinoHelper( 0,3,x,A1 )    #check Anl for n=1 and l from 0 to 2
    elif task=='jl':
        def j(l): return jl(l,x)
        VolinoHelper( 0,4,x,j )     #check spherical bessel functions jl for l running from 0 t0 3
            
    return None

def help(*kargs):

    """list available services"""
    
    services=['datatree','doitall','updateCollection','initCollection','createSimdir',
              'view']
    os.system('/usr/bin/clear')
    if not kargs:
        print 'Services='
        for service in services:
            print '  '+service+' : '+globals()[service].__doc__.split('\n')[0]
        print 'Type "workCANs.py -a help -b service" for more detail on "service"'
    else:
        sv=kargs[0]
        print globals()[sv].__doc__        


if __name__=='__main__':
    inpHand('Usage: workCANs.py [options]',
            ' -a _R_service required service. Type "help" for list of services.',              
            ' -b __servargs service required arguments in python syntax, and separate them with ",". Escape symbols \' and " (def:None)',
            ' -c __servoptargs service optional arguments (Example: "key=val" for one argument and "key1=val1, key2=val2,..." for two or more arguments (def:None)',
            ).parse(locals(),sys.argv)

    #parse required arguments [arg1,arg2,...]
    reqargs=[]
    if servargs:
        reqargs=servargs.split(',')
    #parse optional arguments {key1:val1, key2:val2,...}
    optsargs={}
    if servoptargs:
        for pair in servoptargs.split(','):
            key,val=pair.split('=')
            optsargs[key]=val
    #execute the service
    locals()[service](*reqargs,**optsargs)
