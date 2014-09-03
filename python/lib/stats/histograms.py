import numpy,pdb,copy

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

def chiSquare(Qs,xvals,hist,fit='normal',lowcut=0.00):
    """Calculate Chi-Square fit of a histogram to a fit function
    xvals: ordinate values for hist
    lowcut: remove from the calculation those bins with population
            lower than lowcut*max(hist)

    THEORY

    We will transform the incoming histogram into a probability
    density with standar deviation equal to one and average of zero.
    
    In biology the most common application for chi-squared is in
    comparing observed counts of particular cases to the expected
    counts. For example, the willow tree (Salix) is dioecious, that
    is, like people (and unlike most plants) a willow tree will have
    just male or female sex organs. One might expect that half of the
    willows are male and half female. If you examine N willow trees
    and count that x1 of them are male and x2 of them are female, you
    will probably not find that exactly x1=0.5*N and x2=0.5*N. Is the
    difference significant enough to rule out the 50/50 hypothesis? We
    could almost calculate the chi-squared, but we don't know the
    standard deviation for each count. Never fear: most counts are
    distributed according to the Poisson distribution, and as such the
    standard deviation equals the square root of the expected
    count. Thus we can calculate X2:

       \sum_i^N  (O_i-E_i)^2/E_i
       O_i: observed probability
       E_i: theoretical probability    
    """
    N=len(hist)
    M=hist.max()*float(lowcut)
    dX=(xvals[-1]-xvals[0])/(N-1)
    if fit=='normal':
        """calculate a normal distribution based on the histogram data"""
        X=copy.copy(xvals); H=copy.copy(hist)
        X=numpy.array(X);   H=numpy.array(H)
        av=(X*H).sum()/H.sum()
        std=numpy.sqrt( (H*(X-av)**2).sum()/H.sum() )
        b=1.0/numpy.sqrt(2*numpy.pi); c=1.0/(H.sum()*(dX/std))
        chi=0.0
        for i in range(N):
            if H[i]<M: continue #too much error in low numbers
            x=(X[i]-av)/std #normalize X-axis
            y=c*H[i]        #histogram to probability density
            theor=b*numpy.exp(-0.5*x*x)
            chi+=(dX/std)*(y-theor)**2/theor
    return chi
