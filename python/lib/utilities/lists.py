#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

#remove duplicates from a list
def unique(lst):
    d = {}
    for item in lst: d[item] = None
    return d.keys()
    
def multiply(xl,yl):
    v=[]
    for i in range(len(xl)):
        v.append(xl[i]*yl[i])
    return v

def default_comparison(x_1,x_2):
    if x_1 <= x_2: return 1
    return None

#bubble sort according to some two-argument function, which we default to increasing x's
#we can pass an extra list yl, that is changed in parallel to xl
def sort_lists(xl,yl=[],compare_function=default_comparison):
    xl2=xl[:]
    if yl: yl2=yl[:]
    N=len(xl)-2
    i=-1
    while i<N:
        i=i+1  ;  j=i+1
        while j>0:
            if compare_function(xl2[j],xl2[j-1]):
                a=xl2[j]  ;  xl2[j]=xl2[j-1]  ;  xl2[j-1]=a
                if yl:
                    a=yl2[j]  ;  yl2[j]=yl2[j-1]  ;  yl2[j-1]=a
                j=j-1
            else: break
    if yl: return [xl2,yl2]
    return xl2

def read_file_lines(descriptor):
    if isinstance(descriptor,str):
        try: pin=open(descriptor,'r')
        except:
            sys.stderr.write('lists::read_file_lines: Could not open file '+descriptor+'\n')
            return []
    return chomp(pin.readlines())

