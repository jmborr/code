import pdb
from copy import deepcopy
from Scientific.IO.NetCDF import NetCDFFile

def cloneNetCDF(outfile,templatefile):
    "return NetCDF object based on a templatefile"""
    pin=NetCDFFile(templatefile, 'r')
    pout=NetCDFFile(outfile, 'w')
    #clone dimensions
    for name,size in pin.dimensions.items(): pout.createDimension(name,size)
    #clone variables
    for name,var in pin.variables.items():
        datatype=var.typecode()
        dimensions=var.dimensions
        pout.createVariable(name, datatype, dimensions)
        pout.variables[name][:]=deepcopy(pin.variables[name][:])
    #clone global attributes
    for attributeName in dir(pin):
        if attributeName not in dir(pout):
            setattr(pout,attributeName,getattr(pin,attributeName) )
    return pout
