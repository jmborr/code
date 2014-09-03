#!/usr/bin/python

import os,sys
from utilities.codedir import codedir


def concut( contf=os.path.join(codedir,'f77/tasser/tasser.1.0/common/contact.comm') ):
    """dictionary of contacts between SC from the contact.comm file"""
    dists={}
    pin=open(contf,'r')
    x=pin.readline().split()
    for i in range(20):
        y=pin.readline().split()
        resy=y[0]
        dists[resy]={}
        for j in range(20):
            dists[resy][x[j]]=float(y[j+1])
    return dists

def devs( contf=os.path.join(codedir,'f77/tasser/tasser.1.0/common/contact.comm') ):
    """dictionary of standard deviations in contacts between SC from the contact.comm file"""
    deviations={}
    pin=open(contf,'r')
    for i in range(21): pin.readline() #skip first 21 lines
    x=pin.readline().split()
    for i in range(20):
        y=pin.readline().split()
        resy=y[0]
        deviations[resy]={}
        for j in range(20):
            deviations[resy][x[j]]=float(y[j+1])
    return deviations

def tasserconcuts( contf=os.path.join(codedir,'f77/tasser/tasser.1.0/common/contact.comm') ):
    """dictionary of contacts between SC from the contact.comm file according to tasser"""
    tc={}
    c=concut(contf=contf)
    d=devs(contf=contf)
    for resi in c.keys():
        tc[resi]={}
        row=c[resi]
        #the contacts in tasser are less restringent than in concut
        for resj in row.keys(): tc[resi][resj]=c[resi][resj]+2.5*d[resi][resj]
    return tc
