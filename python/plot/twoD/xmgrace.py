import pdb

def binPlot(m,outf,absmin=None,nbin=10):
    """chop a matrix m into ten matrixes ranging from lowest to
    highest value"""

    from utilities.readingWritingFiles import write_from_numpy

    print 'max=',str(m.max())
    minimum=m.min(); delta=(m.max()-minimum)/nbin
    current=minimum+delta

    header='';    buf=['']*nbin
    for ibin in range(nbin):
        lower=minimum+delta*(ibin); upper=lower+delta
        header+='@    s%d legend  "%f %f"\n'%(ibin,lower,upper)
        buf[ibin]+='@target G0.S%d\n'%ibin

    N=len(m); M=len(m[0])
    for i in range(N):
        for j in range(M):
            z=m[i][j];
            if absmin and z<absmin: continue
            index=int( (z-minimum)/delta )
            if index==nbin: index-=1
            buf[index]+='%d %d\n'%(i,j)
    open(outf,'w').write( header+'&\n'.join(buf) )
    
