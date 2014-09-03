#!/usr/bin/python

import os,sys,numpy,pdb

chain=['','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
       'P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4',
       '5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
       'P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4',
       '5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
       'P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4',
       '5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
       'P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4',
       '5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
       'P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4',
       '5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
       'P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4',
       '5','6','7','8','9']

templ=((' N  ','ATOM  _N_  N   GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' H1 ','ATOM  _N_  H1  GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' H2 ','ATOM  _N_  H2  GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' H3 ','ATOM  _N_  H3  GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' CA ','ATOM  _N_  CA  GLY _A_ _R_    _XYZ_  1.00  0.00'),
       ('2HA ','ATOM  _N_ 2HA  GLY _A_ _R_    _XYZ_  1.00  0.00'),
       ('3HA ','ATOM  _N_ 3HA  GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' C  ','ATOM  _N_  C   GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' O  ','ATOM  _N_  O   GLY _A_ _R_    _XYZ_  1.00  0.00'),
       (' OXT','ATOM  _N_  OXT GLY _A_ _R_    _XYZ_  1.00  0.00')
       )

#positional parameters as in Table1 of Jonsson72. We rename the atoms
#to PDB naming convention. We need to multiply them by "zoom"
zoom=1E-5
gly={' OXT':numpy.array([ 30494,  9439, 23539]),
     ' O  ':numpy.array([-14722, 14150, 10708]),
     ' N  ':numpy.array([ 30116,  8984,-25904]),
     ' C  ':numpy.array([  7504, 12486,  6619]),
     ' CA ':numpy.array([  6474, 14485,-21308]),
     ' H3 ':numpy.array([ 28972, 10036,-45414]),
     ' H2 ':numpy.array([ 49450, 11929,-13184]),
     ' H1 ':numpy.array([ 29935,   561,-22613]),
     '2HA ':numpy.array([  7688, 23444,-24322]),
     '3HA ':numpy.array([-13322, 11439,-35718])
     }

#lattice parameters. Monoclinic P2_1/n (a.k.a P2_1/c)
a=5.1054; b=11.9688; c=5.4645; beta=3.141592*(111.697/180)
abc=numpy.array([5.1054,11.9688,5.4645])
x=numpy.array([a+c*numpy.cos(beta),b,c*numpy.sin(beta)])
def asym_unit():
    """build the asymmetric unit"""
    buf=''
    N=0
    R=1
    for (name,line) in templ:
        fa,fb,fc=gly[name]*zoom
        xyz=[a*fa+c*fc*numpy.cos(beta),b*fb,c*fc*numpy.sin(beta)]
        XYZ=''
        for r in xyz: XYZ+='%8.3lf'%r
        #pdb.set_trace()
        line=line.replace('_A_',chain[R])
        line=line.replace('_N_','%4d'%N)
        line=line.replace('_R_',str(R))
        buf+=line.replace('_XYZ_',XYZ)+'\n'
        open('mono-GLY.pdb','w').write(buf)
        N+=1
        
def unit(n,m,l,buf,R,N):
    """generate one lattice unit at (n,m,l)"""
    for (name,line) in templ:
        #general coordinate (n+x,m+y,l+z)
        fa,fb,fc=gly[name]*zoom;
        fa=n+fa; fb=m+fb; fc=l+fc
        xyz=[a*fa+c*fc*numpy.cos(beta),b*fb,c*fc*numpy.sin(beta)]
        XYZ=''
        for r in xyz: XYZ+='%8.3lf'%r
        line=line.replace('_N_','%4d'%N)
        line=line.replace('_A_',chain[R])
        line=line.replace('_R_','%4d'%R)
        buf+=line.replace('_XYZ_',XYZ)+'\n'
        N+=1
    buf+='TER\n'
    R+=1
    return (buf,N,R)

def unit2(n,m,l,buf,R,N):
    """generate one lattice unit at (n,m,l)"""
    for (name,line) in templ:
        #general coordinate (n+1/2-x,m+1/2-y,l+1/2-z)
        fa,fb,fc=gly[name]*zoom;
        fa=n-fa; fb=m-fb; fc=l-fc
        xyz=[a*fa+c*fc*numpy.cos(beta),b*fb,c*fc*numpy.sin(beta)]
        XYZ=''
        for r in xyz: XYZ+='%8.3lf'%r
        line=line.replace('_N_','%4d'%N)
        line=line.replace('_A_',chain[R])
        line=line.replace('_R_','%4d'%R)
        buf+=line.replace('_XYZ_',XYZ)+'\n'
        N+=1
    buf+='TER\n'
    R+=1
    return (buf,N,R)

def unit3(n,m,l,buf,R,N):
    """generate one lattice unit at (n,m,l)"""
    for (name,line) in templ:
        #general coordinate (n+1/2-x,m+1/2-y,l+1/2-z)
        fa,fb,fc=gly[name]*zoom;
        fa=n+0.5-fa; fb=m+0.5+fb; fc=l-fc
        xyz=[a*fa+c*fc*numpy.cos(beta),b*fb,c*fc*numpy.sin(beta)]
        XYZ=''
        for r in xyz: XYZ+='%8.3lf'%r
        line=line.replace('_N_','%4d'%N)
        line=line.replace('_A_',chain[R])
        line=line.replace('_R_','%4d'%R)
        buf+=line.replace('_XYZ_',XYZ)+'\n'
        N+=1
    buf+='TER\n'
    R+=1
    return (buf,N,R)

def unit4(n,m,l,buf,R,N):
    """generate one lattice unit at (n,m,l)"""
    for (name,line) in templ:
        #general coordinate (n+1/2-x,m+1/2-y,l+1/2-z)
        fa,fb,fc=gly[name]*zoom;
        fa=n+0.5+fa; fb=m+0.5-fb; fc=l+fc
        xyz=[a*fa+c*fc*numpy.cos(beta),b*fb,c*fc*numpy.sin(beta)]
        XYZ=''
        for r in xyz: XYZ+='%8.3lf'%r
        line=line.replace('_N_','%4d'%N)
        line=line.replace('_A_',chain[R])
        line=line.replace('_R_','%4d'%R)
        buf+=line.replace('_XYZ_',XYZ)+'\n'
        N+=1
    buf+='TER\n'
    R+=1
    return (buf,N,R)


def crystal(L=2):
    buf=''; N=1; R=1
    for n in range(2*L):
        for m in range(L):
            for l in range(2*L):
                buf,N,R=unit(n,m,l,buf,R,1)
                buf,N,R=unit2(n,m,l,buf,R,1)
                buf,N,R=unit3(n,m,l,buf,R,1)
                buf,N,R=unit4(n,m,l,buf,R,1)
    open('poly-GLY.pdb','w').write(buf)

def crystal2(L=2):
    buf=''; N=1; R=1
    for n in range(4):
        for m in range(1):
            for l in range(4):
                buf,N,R=unit(n,m,l,buf,R,1)
                buf,N,R=unit2(n,m,l,buf,R,1)
                buf,N,R=unit3(n,m,l,buf,R,1)
                buf,N,R=unit4(n,m,l,buf,R,1)
    open('poly-GLY.pdb','w').write(buf)
    
#asym_unit()
#crystal2()
#sys.exit(0)

#buf=''; R=1
#buf,N,R=unit(0,0,0,buf,R,1)
#buf,N,R=unit2(0,0,0,buf,R,1)
#buf,N,R=unit3(0,0,0,buf,R,1)
#buf,N,R=unit4(0,0,0,buf,R,1)

#buf,N,R=unit(1,1,1,buf,1,1)
#buf,N,R=unit(1,2,1,buf,1,1)

#buf,N,R=unit(1,0,0,buf,R,1)
#buf,N,R=unit2(1,0,0,buf,R,1)
#buf,N,R=unit(2,0,0,buf,R,1)
#buf,N,R=unit2(2,0,0,buf,R,1)

#buf,N,R=unit(0,0,1,buf,R,1)
#buf,N,R=unit2(0,0,1,buf,R,1)
#buf,N,R=unit(0,0,2,buf,R,1)
#buf,N,R=unit2(0,0,2,buf,R,1)

#buf,N,R=unit(1,0,1,buf,R,1)
#buf,N,R=unit2(1,0,1,buf,R,1)
#buf,N,R=unit(1,0,2,buf,R,1)
#buf,N,R=unit2(1,0,2,buf,R,1)

#buf,N,R=unit(2,0,1,buf,R,1)
#buf,N,R=unit2(2,0,1,buf,R,1)
#buf,N,R=unit(2,0,2,buf,R,1)
#buf,N,R=unit2(2,0,2,buf,R,1)

#open('junk.pdb','w').write(buf)
