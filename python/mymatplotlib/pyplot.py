import matplotlib.pyplot as plt
import numpy

from pdb import set_trace as tr

def imshowfromfile(file, nx, ny, usecols=None, xlabel='X-axis', ylabel='Y-axis', zlabel='Z-axis',
                    cmaptype='gist_rainbow', vmax=None):
    '''
    '''
    data = numpy.loadtxt(file, usecols=usecols)
    data = data.transpose()
    #tr()
    x = data[0]
    y = data[1]
    z = data[2]
    z = z.reshape(nx,ny)
    z = z.transpose()
    cmap =  plt.get_cmap(cmaptype)
    im = plt.imshow(z[::-1,:],  vmax=vmax, extent=(min(x), max(x), min(y), max(y)), cmap=cmap)
    plt.xlabel(xlabel, fontsize='xx-large')
    plt.xticks(numpy.arange(min(x),max(x)+1,20), fontsize='large')
    plt.ylabel(ylabel, fontsize='xx-large')
    plt.yticks(numpy.arange(min(y),max(y)+1,20), fontsize='large')
    plt.grid()
    cl=plt.colorbar(shrink=0.9)
    cl.set_label(zlabel, fontsize='xx-large')
    cl.ax.tick_params(labelsize=13)


    
    
