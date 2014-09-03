import os,numpy,pdb
from tempfile import mkstemp

"""
You should have installed gnuplot

"""


def matrixGreyPlot(m,outf,title=None):
    from utilities.readingWritingFiles import write_from_numpy
    """Plot a 2D color map of the numpy.array 'm' """
    cmd="""set xrange [0:_N_]; set yrange [0:_N_]
set size ratio 1
set pm3d map
set palette gray negative
set pointsize 10
set title "_TITLE_"
set output '_OUTF_'
set terminal postscript enhanced color
splot '_INFILE_' matrix"""
    (a,fn)=mkstemp(); pt=open(fn,'w')
    write_from_numpy(pt,m,ncols=len(m)) #write matrix to temporary file
    if not title: title='2D plot'
    for (key,rpl) in ( (  '_N_',str( len(m) )  ),
                       (  '_MIN_',str( m.min() )  ),
                       (  '_MAX_',str( m.max() )  ),
                       ('_TITLE_',title),
                       ('_OUTF_',outf), ('_INFILE_',fn) ):
        cmd=cmd.replace(key,rpl)
    #pdb.set_trace()
    (a,gn)=mkstemp(); pt=open(gn,'w').write(cmd) #commands to temp file
    os.system('gnuplot '+gn) #make plot
    os.system('/bin/rm %s %s'%(fn,gn)) #remove temp files

def matrixPlot(m,outf,title=None):
    from utilities.readingWritingFiles import write_from_numpy
    """Plot a 2D color map of the numpy.array 'm' """
    cmd="""set xrange [ 0 : _N_ ]
set yrange [ 0 : _N_ ]
set cbrange [ _MIN_ : _MAX_ ]
set size ratio 1
set pm3d map
set palette rgbformulae 22,13,-31
set terminal postscript enhanced color
set title "_TITLE_"
set output '_OUTF_'
splot '_INFILE_' matrix"""

#set size ratio 1
    (a,fn)=mkstemp(); pt=open(fn,'w')
    write_from_numpy(pt,m,ncols=len(m)) #write matrix to temporary file
    if not title: title='grey map'
    for (key,rpl) in ( (  '_N_',str( len(m) )  ),
                       (  '_MIN_',str( m.min() )  ),
                       (  '_MAX_',str( m.max() )  ),
                       ('_TITLE_',title),
                       ('_OUTF_',outf), ('_INFILE_',fn) ):
        cmd=cmd.replace(key,rpl)
    #pdb.set_trace()
    (a,gn)=mkstemp(); pt=open(gn,'w').write(cmd) #commands to temp file
    #pdb.set_trace()
    os.system('gnuplot '+gn) #make plot
    os.system('/bin/rm %s %s'%(fn,gn)) #remove temp files
