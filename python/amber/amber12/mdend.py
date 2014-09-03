import re
import os
from pdb import set_trace as tr

class MdendException(Exception):
    """Base class for Mdend exceptions"""
    def __init__(self,value):
        self.value=value
    def __str__(self):
        return repr(self.value)

class FormatException(MdendException):
    """File format error exception"""
    def __str__(self):
        return '%s does not conform to expected file format'%repr(self.value)

    
class Mdend(object):

    def __init__(self):
        pass

    def load_file(self, filename):
        """Load file to object"""
        try:
            self.validate_mdend_file(filename)
            self.load_file_header(filename)
            self.load_file_body(filename)
            self.filename=filename
        except Exception:
            self.filename=None
    
    def validate_mdend_file(self, filename):
        """check format conforms to expected format"""
        pattern=re.compile('^L\d+\s')
        if not os.path.exists(filename):
            raise IOError(filename)
        for line in open(filename).readlines():
            if not pattern.match(line):
                raise FormatException(filename)

    def load_file_header(self, filename):
        """load field names"""
        pfile=open(filename)
        prev_line_index=-1
        self.fields=[]
        while True:
            line=pfile.readline()
            curr_line_index = int( line.split()[0][1:] )
            if curr_line_index < prev_line_index:
                break # finished reading the header
            prev_line_index = curr_line_index
            self.fields+=line.split()[1:]
        self.lines_per_record = 1+prev_line_index
        self.contents={}
        for field in self.fields: self.contents[field]=[]
        return self.fields

    def skip_header(self,pfile):
        for iline in range(self.lines_per_record):
            pfile.readline()

    def load_next_record(self, pfile):
        values=[]
        for iline in range(self.lines_per_record):
            line=pfile.readline()
            if not line:
                return []
            values += [ float(value) for value in line.split()[1:] ]
        for i in range(len(values)):
            field=self.fields[i]
            value=values[i]
            self.contents[field].append(value)
        return values

    def load_file_body(self,filename):
        pfile=open(filename)
        self.skip_header(pfile)
        while self.load_next_record(pfile):
            pass
        return pfile

    def plot(self, fields, x='Nsteps'):
        """Make a simple Figure of one or more fields versus some common axis

        Creates a figure containing as many plots as fields in the fields argument

        Args:
          fields [str or list]: one field or a list of fields, the Y-axes
          x [str]: X-axis

        Returns:
          matplotlib.figure object containing the plots
        """
        from matplotlib.pyplot import figure
        if isinstance(fields,str):
            fields=[fields,]
        nplots=len(fields)
        fig=figure( figsize=(7*nplots,4) )
        fig.subplots_adjust(wspace=0.5)
        for ifield in range(nplots):
            field=fields[ifield]
            ifig=fig.add_subplot(1, nplots, ifield, xlabel=x, ylabel=field)
            ifig.plot(self.contents[x], self.contents[field] )
        return fig

