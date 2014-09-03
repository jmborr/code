#!/usr/bin/python

import sys,os,re,numpy
from utilities.lists import sort_lists
from math import sqrt,exp

def listproduct(xl,yl):
    N=len(xl)
    if N!=len(yl):
        sys.stderr.write('ERROR statistics::listproduct: lists of different length')
        return []
    xyl=[]
    for i in range(N): xyl.append(xl[i]*yl[i])    
    return xyl

def average(z):
    av=0.
    for x in z:av=av+x
    return float(av)/len(z)

def std2(z):
    av=average(z)
    a=0
    for x in z: a=a+x*x
    return float(a)/len(z)-av*av

def standard_dev(z):
    return sqrt(std2(z))

def deviation_of_the_mean(z):
    return sqrt(std2(z)/len(z))

def kurtosis(z):
    av=average(z)
    dev=standard_dev(z)
    kurt=0
    for x in z: kurt+=(x-av)**4
    return kurt/( (len(z)-1)*dev**4 )

def excess_kurtosis(z):
    return kurtosis(z)-3

def skewness(z):
    av=average(z)
    dev=standard_dev(z)
    sk=0
    for x in z: sk+=(x-av)**3
    return sk/( (len(z)-1)*dev**3 )

#we can do a running average in different ways:
#(1)Average all points in bins of width dx. Then we return a list composed of four list:
#   average x and y-coordinates, standard deviation and standard deviation of the mean
#(2)if we treat y as an ordered list (forget about x), then yl[m]+yl[m+1]+..+yl[m+n]
def running_average(xl,yl,dx=0,n=0):
    if dx:
        [xl2,yl2]=sort_lists(xl,yl) #order yl according to increasing x's
        #avyl:average over the bin, stdl:standard deviation, stdml:standard deviation of the mean
        avxl=[]  ;  avyl=[]  ;  stdl=[]  ;  stdml=[]
        xlimit=xl2[0]  ;  N=len(xl2)  ;  i=0
        while i<N:
            x=xl2[i]  ;  y=yl2[i]  ;  a=0  ;  b=0  ;  c=0  ;  n=0  ;  xlimit=xlimit+dx            
            while x<=xlimit:
                a=a+x  ;  b=b+y  ;  c=c+y*y  ;  n=n+1
                i=i+1
                if i==N: break
                x=xl2[i]  ;  y=yl2[i]
            avxl.append(a/n)  ;  avyl.append(b/n)
            a=sqrt(c/n-b*b/(n*n))   ;  b=a/sqrt(n)
            stdl.append(a)  ;  stdml.append(b)
        return [avxl,avyl,stdl,stdml]        


def corrcoeff(xl,yl):
    xyl=listproduct(xl,yl)
    return ( average(xyl)-average(xl)*average(yl) )/( standard_dev(xl)*standard_dev(yl) )


def linreg(xl,yl):

    """Linear regression of y = ax + b"""
    
    if len(xl) != len(yl):  raise ValueError, 'unequal length'
    N = len(xl)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in map(None, xl, yl):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
    meanerror = residual = 0.0
    for x, y in map(None, xl, yl):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2
    RR = 1 - residual/meanerror
    ss = residual / (N-2)
    Var_a, Var_b = ss * N / det, ss * Sxx / det
    return {'slope':a, 'intercept':b, 'corrcoef':RR}

def invLapl(data,trange=[0,10],nt=100):
    """numerically calculate the laplace transform for a data set
    data=[(s_0,L_0),..,(s_{ns-1},L_{ns-1}]
    We fill fit the data to the following equation:
    L(s)=Sum_{it=0}^{it=nt-1} f_it * exp(-s*trange*it/nt)
    returns [(t_0,f_0),..,(t_{nt-1},f_{nt-1}]
    """
    from Scientific.Functions.LeastSquares import leastSquaresFit
    dt=float(trange[-1]-trange[0])/nt; t0=trange[0];
    p=[1.0/nt]*nt #initial guess
    #dynamically parameterize function invLapModel with parameters t0
    #and dt    
    def invLaplModel(param,s):
        sum=0; t=t0
        for it in range(len(param)):
            sum+=param[it]*exp(-s*t); t+=dt
        return sum*dt

    results=leastSquaresFit(invLaplModel,p,data)[0]
    
    invL=[];  t=t0
    for it in range(nt):
        invL.append( (t,results[it]) ); t+=dt
    return invL

def calc_varcov(av,dev,varcov,dim,ndata):
    """calculate the variance-covariance matrix
    av_i=sum(x_i)
    dev_i=sum(x_i*x_i)
    varcov_{ij}=sum(x_i*x_j)"""
    av/=ndata; varcov/=ndata; dev=numpy.sqrt(dev/ndata-av*av)
    varcov-=numpy.dot(av.reshape(dim,1),av.reshape(1,dim))
    return {'av':av,'dev':dev,'varcov':varcov}
    
def calc_cov(av,dev,cov,dim,ndata):
    """calculate covariance matrix"""
    results=calc_varcov(av,dev,cov,dim,ndata)
    dev=results['dev']
    results['cov']=results['varcov']/numpy.dot(dev.reshape(dim,1),dev.reshape(1,dim))
    del results['varcov']
    return results
