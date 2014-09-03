import os
from math import log,pow
lines=open('phi_with_mutant_crude.dat','r').readlines()

n=0
x=0.0
x2=0.0
for line in lines:
    list=line.split()
    exper=float(list[6])
    n+=1
    x+=pow((1.2-exper),3)+pow((0.22+exper),3)
    x2+=pow((1.2-exper),5)+pow((0.22+exper),5)
x=x/(1.42*3*n)   
x2=x2/(1.42*5*n) 

y=0
for line in lines:
    theor=float(list[1])
    exper=float(list[6])
    y+=pow((theor-exper),2)
y=(y/n-x)/x2

os.system('./integral.exe '+`y`)

