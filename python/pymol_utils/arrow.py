from pymol import cmd
from pymol.cgo import *
from math import *

#
# Some functions to allow drawing arrows (vectors) in Pymol
# In need of proper documentation...
#
# Please don't distribute (parts of) this file, without credits
#
# (c)2006 Tsjerk A. Wassenaar, PhD, University of Utrecht
#
# t s j e r k w .at. g m a i l .dot. c o m
# http://nmr.chem.uu.nl/~tsjerk/
#
# EXAMPLE of pml file:
# load pdbfile.pdb, selection-name
# run arrow.py
# cgo_arrow([x1,y1,z1],[x2,y2,z2],name="arrow1",color=[0.5,0,0.5])
# cgo_arrow([x3,y3,z3],[x4,y4,z4],name="arrow2",color=[0,0,1])
# center
# zoom


def t( X ):
   if not X: return X
   Y = []
   for i in range( len( X[0] ) ):
       Y.append( [] )
       for j in X:
           Y[i].append( j[i] )
   return Y

def v_add( a, b ):   return ( a[0]+b[0], a[1]+b[1], a[2]+b[2] )
def v_sub( a, b ):   return ( a[0]-b[0], a[1]-b[1], a[2]-b[2] )
def vecprod( a, b ): return ( a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] )
def inprod( a, b=None ):
   if b: return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
   else: return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]
def svmult( s, a ):  return ( s*a[0], s*a[1], s*a[2] )
def norm( a ):       return svmult( 1/sqrt(inprod( a )), a )

def mvmult( R, x ):
   y = []
   for i in R: y.append( inprod( i, x ) )
   return tuple(y)

def mv_add( X, a ):
   Y = []
   for i in X:
       Y.append( v_add( i, a ) )
   return Y

def mmmult( R, X ):
   Y = []
   for i in X: Y.append( mvmult( R, i ) )
   return Y

def smatrix( v ): return [[ v[0], 0, 0 ], [ 0, v[1], 0 ], [ 0, 0, v[2] ]]
def rmatrix( v ):
   cosx, sinx = cos( v[0] ), sin( v[0] )
   cosy, siny = cos( v[1] ), sin( v[1] )
   cosz, sinz = cos( v[2] ), sin( v[2] )

   return mmmult( mmmult( [[1,0,0],[0,cosx,-sinx],[0,sinx,cosx]],[[cosy,0,-siny],[0,1,0],[siny,0,cosy]] ),[[cosz,-sinz,0],[sinz,cosz,0],[0,0,1]] )

def block( i, dphi ):
   ddphi = 0.25*dphi
   phi0  = i*dphi
   phi1  = phi0+ddphi
   phi2  = phi1+ddphi
   phi3  = phi2+ddphi
   phi4  = phi3+ddphi
   sqrt2 = sqrt(2)
   return [  (-0.5*sqrt2,-0.5*sqrt2*cos(phi2),0.5*sqrt2*sin(phi2)),
             (1,0,0),
             (0,cos(phi0),-sin(phi0)),
             (0,cos(phi1),-sin(phi1)),
             (0,cos(phi2),-sin(phi2)),
             (0,cos(phi3),-sin(phi3)),
             (0,cos(phi4),-sin(phi4))
          ]

def cgo_triangle_fan( X ):
   Y = []
   while ( X ):
       i = X.pop(0)
       Y.extend( [ NORMAL, i[0], i[1], i[2], ] )
       for i in range( 6 ):
           i = X.pop(0)
           Y.extend( [ VERTEX, i[0], i[1], i[2], ] )
   return Y

def cgo_arrow1( S, E, r=0.2, hr=0.4, hl=1.0, color=[1,0,0]):
   P0 = S
   D  = v_sub( E, S )
   DL = inprod( D, D )
   P1 = v_add( S, svmult( (DL-hl)/DL, D ) )
   P2 = E

   # Define a vector orthogonal to P1-P0
   V = v_sub( P1, P0 )
   V = norm( V )
   if V[2] != 0:
       A = ( 1, 1, -(V[0]+V[1])/V[2] )
   elif V[1] != 0:
       A = ( 1, -V[0]/V[1], 0 )
   else:
       A = ( 0, -V[0], 0 )
   A = norm( A )

   B = vecprod( V, A )
   print (inprod(V), inprod(B), inprod(A))
   R = t([ svmult( hl,V ), svmult( hr,A ), svmult( hr,B ) ])

   # Define the transformation matrix (scale and rotation)
   #C = v_sub( P2, P1 )
   #scale  = ( hl, hr, hr )
   #rotate =(0,acos( C[0]/sqrt(C[0]**2+C[2]**2) ), acos(C[0]/sqrt(C[0]**2+C[1]**2)))
   #R = mmmult( smatrix( scale ), rmatrix( rotate ) )
   obj = [CYLINDER, S[0], S[1], S[2], P1[0], P1[1], P1[2], r]+\
         color+color+[COLOR,]+color+[BEGIN, TRIANGLE_FAN ]

   N = 10
   dphi = 2*pi/N
   crds = []
   for i in range(N+1):
       crds.extend( block( i, dphi ) )
   crds = mv_add( mmmult( R, crds ), P1 )
   obj.extend( cgo_triangle_fan( crds ) )
   obj.extend( [ END, ] )

   return obj

def cgo_arrow( S, E, r=0.3, hr=0.5, hl=2.0, color=[1,0,0], name="arrow", state=1):
   obj = cgo_arrow1( S, E, r=r, hr=hr, hl=hl, color=color)
   cmd.load_cgo( obj, name, state )

def cgo_arrows( X, r=0.2, hr=0.4, hl=1.0, name="arrow", state=1 ):
   obj = []
   for i in X:
       obj.extend( cgo_arrow1( (i[0], i[1], i[2]), (i[3], i[4],i[5]), r=r, hr=hr, hl=hl ) )
   cmd.load_cgo( obj, name, state )

def network(selA, inpf, selB=None, col=3, maxR=0.5, minR=0.01,
            hr=0.0, hl=0.0, color=[1,0,0],
            minimum=None, maximum=None, invertSign='no',
            maxCO=None, minCO=None, scaling='linear'):
   """
   Read a file with residue indexes and a value:
     #resi  resj val
       1      3  5.6
       ...
     Will draw arrows (lines for hl=0) between CA atoms with thickness
     scaled between minR and maxR
   """
   #load pairs
   resi=[]; resj=[]; vals=[]
   for l in open(inpf).readlines():
      if l[0]=='#': continue
      items = l.split()
      x = float(items[ col-1])
      if invertSign[0].lower() == 'y': x *= -1
      if minCO and x < minCO: continue #disregard too small values
      if maxCO and x > maxCO: continue #disregard too big values
      resi.append( int(items[0]) )
      resj.append( int(items[1]) )
      vals.append( x )
   #scale values
   if minimum == None: minimum = min(vals)
   if maximum == None: maximum = max(vals)
   f = ( maxR-minR )/( maximum-minimum )
   N = len(vals)
   for i in range(N):
      x = vals[i]
      if x > maximum: x = maximum   #saturate over maximum
      elif x < minimum: x = minimum #saturate below minimum
      if scaling=='linear':
         vals[i] = minR + f * (x - minimum)
   #Retrieve CA atoms from the selections
   mA = cmd.get_model('%s and name CA'%selA)
   mB = mA
   if selB: mB = cmd.get_model('%s and name CA'%selB)
   #store coordinates for the selections
   resA = {} ; resB = {}
   for atom in mA.atom:
      resA[ atom.resi_number ] = atom.coord
   for atom in mB.atom:
      resB[ atom.resi_number ] = atom.coord
   #plot arrows/lines joining the CA
   for i in range(N):
      ri = resi[i]; rj = resj[i]; v = vals[i]
      if ri not in resA.keys() : continue
      if rj not in resB.keys() : continue
      S = resA[ri] #coordinates of CA atom of residue ri
      E = resB[rj] #coordinates of CA atom of residue ri
      cgo_arrow(S, E, r=v, hr=hr, hl=hl, color=color,
                name='%03d-%03d'%(ri,rj) )
