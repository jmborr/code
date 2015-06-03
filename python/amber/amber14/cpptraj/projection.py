import os
import numpy
import scipy

class projection(object):
    '''Manage output from projection command
For mass-weighted covariant matrix, the transformation between cartesian (X)
and mode coordinates (Q) is:
                    Q = E^t * M^{1/2} * X
For mode 'i', projection Q_i = e_i^t * M^{1/2} * X, where e_i is
the ith quasi-harmonic eigenmode stored in the evecs file, here
represented as the transpose of a vector of length 3*nparticles.
'''
    
    def __init__(self):
        self.handle = None  # handle to projection file
        self.nmode = 0      # number of modes
        self.nframe = 0    # number of frames
        self.frameidx = [] # file indexes for each frame

    def setIndex(self, filename):
        '''Prepare certain attributes for fast reading of file
Argurments:
 filename: (str) self-evident

Raise:
 IOError: file does not exists
'''
        if not os.path.exists(filename):
            raise IOError('File {0} does not exist'.format(filename))
        self.handle = open(filename,'r')
        line = self.handle.readline() # header
        self.nmode = len(header.split())
        while line:
            self.frameidx.append(self.handle.tell())
            line = self.handle.readline()
        self.nframe = len(self.frameidx)
        
    def is_frameidx_loaded(self):
        if not self.frameidx:
            raise Exception('Filename not loaded')

    def getFrame(self, framenumber):
        '''Find and return selected frame
Input:
 framenumber: (int) Notice the first frame must be index "1"

Returns:
 numpy array containing mode projections for the selected frame
'''
        self.is_frameidx_loaded() # check if filename loaded
        index = framenumber - 1
        if index < 0:
            raise Exception('Framenumber must be bigger than zero')
        if framnumber > self.nframe:
            raise Exception('Frame number bigger than number of frames ({0})'.format(self.nframe))
        self.handle.seek(index)
        frame = [float(x) for x in self.handle.readline().split()]
        return numpy.array(frame)

    def getProjection(self, modenumber):
        '''Return the projection of the trajectory along a particular mode
Input:
 imode: (int) Notice the mode indexes begin with one, not zero

Returns:
 numpy array containing the projections for selected mode along the trajectory
'''
        self.is_frameidx_loaded() # check if filename loaded
        index = modenumber - 1
        if index < 0:
            raise Exception('Framenumber must be bigger than zero')
        if modenumber > self.nmode:
            raise Exception('modenumber bigger than number of modes ({0})'.format(self.nmode))
        self.handle.seek(0)
        projections = []
        for iframe in range(self.nframe):
            line = self.handle.readline()
            projections.append(float( line.split()[index] ))
        return numpy.array(projections)

    def loadToMemory(self, firstframe=1, lastframe=self.nframe, firstmode=1, lastmode=self.nmode):
        '''Load projection file to memory
Input:
 [firstframe]: (int) (def:0) First frame to read. Index begin with one, not zero
 [lastframe]: (int) (def: last frame in trajectory) lastframe will also be read
 [firsmode]: (int) (def: first mode). Index begin with one, not zero
 [lastmode]: (int) (def: last mode) lastmode will also be read

Returns:
 bidimensional numpy.array, first index in frame index, second index in mode index
 array shape = (lastframe-firstframe+1, lastmode-firstmode+1)
'''
        self.is_frameidx_loaded() # check if filename loaded
        if firstframe < 1 or lastframe > self.nframe:
            raise Exception('requested frames out of bounds (1 to {0})'.format(self.nframe))
        if firstmode < 1 or lastmode > self.nmode:
            raise Exception('requested modes out of bounds (1 to {0})'.format(self.nmode))
        first_mode_index = firstmode - 1
        last_mode_index = lastmode - 1
        n_mode_read = lastmode - firstmode + 1
        frame_index_start = firstframe - 1
        self.handle.seek(frame_index_start)
        n_frame_read = lastframe - firstframe + 1
        projections = numpy.zeros(n_frame_read, n_mode_read)
        for iframe in range(n_frame_read):
            values = [ float(x) for x in self.handle.read().split() ]
            projections[iframe] = numpy.array(values[first_mode_index, last_mode_index]])
        return projections

    def fit2Normal(self):
        '''Test the values of each model projection to a normal distribution
Returns:
 fits: (list) a list with as many entries as modes. Each entry is a
       dictionary with the following keys:
       'avg': average of the values
       'std': standard deviation
       'chi': Chi-square of fit to a normal distribution
               mean "avg" and standard deviation "std"
'''
        from scipy.stats.mstats import normaltest
        if scipy.__version__ < '0.14':
            raise Exception('fit2Normal requires scipy version 0.14 and above')
        self.is_frameidx_loaded() # check if filename loaded
        for imode in range(1, 1+self.nmode):
            proj = self.getProjection(imode)
            normaltest(proj)
            

