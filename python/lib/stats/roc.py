#!/usr/bin/python

def roc(input,ycut):

    #parse input. At end we have a list of [x,y] pairs
    #we can pass the name of a file, a pointer to a file, or a list of [x,y] pairs
    inpt=None
    if isinstance(input,str): inpt=open(input,'r')
    else if isinstance(input,file): inpt=input
    list=[]
    if inpt:
        ll=inpt.readlines()
        for l in ll:
            list.append( l.split() )
    else: list=copy.deepcopy(input) #literal copy

    #sort list by increasing x
    yy={}
    for (x,y) in list: yy[x]=y
    xx=yy.keys()
    xx.sort()
    list=[]
    for x in xx: list.append( [x,yy[x]] )
    yy=[] ; xx=None #liberate memory

    #find the predicted positive rate for every x value
