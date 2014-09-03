import numpy,pdb

def histogram(a,bins=10,rg=None,normed=False):
    a = numpy.asarray(a).ravel()
    if not numpy.iterable(bins):
        if rg is None:
            rg = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in rg]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = numpy.linspace(mn, mx, bins, endpoint=False)

    # best block size probably depends on processor cache size
    block = 65536
    n = numpy.sort(a[:block]).searchsorted(bins)
    for i in range(block,len(a),block):
        n += numpy.sort(a[i:i+block]).searchsorted(bins)
    n = numpy.concatenate([n, [len(a)]])
    n = n[1:]-n[:-1]
    
    if normed:
        db = bins[1] - bins[0]
        return 1.0/(a.size*db) * n, bins
    else:
        return n, bins

def Integral(xvals,hist):
    """normalize a histogram"""
    I=0; N=len(xvals)
    xprev=xvals[0]; yprev=hist[0]
    for i in range(1,N):
        xnext=xvals[i]; ynext=hist[i]
        I+=(xnext-xprev)*(ynext+yprev)/2.
        xprev=xnext; yprev=ynext
    return I

def chiSquare(xvals,hist,fit='normal',lowcut=None):
    """Calculate Chi-Square fit of a histogram to a fit function
    lowcut: remove from the calculation those bins with population
            lower than lowcut*max(hist)
    """
    N=len(hist)
    if fit=='normal':
        """calculate a normal distribution based on the histogram data"""
        M=hist.max()
        X=[]; Y=[]
        if not lowcut:
            X=xvals; Y=hist
        else:
            M=M*float(lowcut)
            for i in range(N):
                x=xvals[i]; y=hist[i]
                if y>M:
                    X.append(x); Y.append(y)
            X=numpy.array(X); Y=numpy.array(Y)
        I=Integral(X,Y)
        Y=Y/Integral(X,Y) #normalize histogram
        N=len(X);
        av=(X*Y).sum()/Y.sum(); std2=(Y*(X-av)**2).sum()/Y.sum();
        a=1/(2*std2); b=numpy.sqrt(a/numpy.pi)
        chi=0.0
        for i in range(N):
            x=X[i]; y=Y[i]
            theor=b*numpy.exp(-a*x*x)
            chi+=(y-theor)**2/theor
    return chi/Y.sum()
