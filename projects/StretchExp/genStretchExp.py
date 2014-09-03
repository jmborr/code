import sys
from math import exp,gamma,sqrt
from pdb import set_trace as trace

buf = '#Stretched Exponential = Height * exp( -(x/Lifetime)^StrechingExponent )\n'

print 'input the Height, Lifetime, and Streching Exponent (separeted by spaces)'
h,t,b = [float(x) for x in sys.stdin.readline().split() ]
buf += '#Height = %f\n#Lifetime = %f\n#Stretching Exponent = %f\n#\n'%(h,t,b)

avg = t*gamma(1./b)/b #mean value
dev = sqrt( t*( t*gamma(2./b)-gamma(1./b) )/b ) #standard deviation
buf += '#average = %f\n#standard deviation=%f\n#\n'%(avg,dev)

extent = 4.0  #data will extend four standard deviations around the average
xmin = max(0,avg-extent*dev)
xmax = avg + extent*dev
xN = 100
dx = (xmax - xmin)/xN

x=xmin
while x < xmax:
    buf += '%f %f\n'%(x, h*exp(-(x/t)**b) )
    x += dx

print 'input the output file name'
filename = sys.stdin.readline().strip()
open(filename,'w').write( buf )

sys.exit(0)
