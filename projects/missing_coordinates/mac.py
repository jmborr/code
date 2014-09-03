#!/usr/bin/python
"""Program that:
(1) library for class MAC
(2) provides several services related to MAC objects when executed as
    standalone
(3) initialize global variables pertaining to the missing_coordinates
    project"""
import pickle,os,sys
from utilities.codedir import scratchdir
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,Bye

#global variables
scroot=os.path.join(scratchdir,'missing_coordinates')
scout=scroot+'/out'
filteredlist=scroot+'/preparing/filterout500.list'
maclist=scroot+'/preparing/mac500.list'
collf=scout+'/collection.mac' #single file containing all mac objects.

class MAC:
    """missing atoms chunk (MAC) object implementation
    a MAC is a container to hold ever increasing info

    Recommended attributes:
    'pdbheader': string
    'b_e': list [b,e] with begin and end positions. b==1 for first residue
    'L', L=e-b+1 length of the MAC region
    'macheader': string
    'dumpf': picke dump file
    'dumpbf': basename of pickle dump file
    'dumpdf': dirname of pickle dump file
    'pdbp': dictionary for pdb parent properties
    'chunk': dictionary for MAC properties
    'ptor': dictionary for PROSPECTOR properties
      |_'outf': file with PROSPECTOR output ( (un)compressed tar file )
      |_'tmpl': dictionary for templates properties
          |_:'fbn',  file base name
          |_:'list', list of templates
              |_:header,  template header
                  |_:'seqid',   sequence id. between whole query and template
                  |_:'nalg', number of residues aligned by template on query sequence
                  |_:'maccovn', number of residues of the MAC covered by the template
                  |_:'maccov', float(maccovn)/self.L
                  |_:'Z', Z-score
    """
    def __init__(self,**kwargs):
        self.__dict__=kwargs.copy() #shallow copy

    def genMacHeader(self):
        self.macheader=self.pdbheader+'.'+'%03d_%03d'%(self.b_e[0],self.b_e[1])
        
    def genDumpf(self,dir):
        """generate name of picke dump file, as well as dirname and basename"""
        if 'macheader' not in self.__dict__.keys(): self.genMacHeader()
        self.dumpdf=dir
        self.dumpbf=self.macheader+'.mac'
        self.dumpf=os.path.join(dir,self.dumpbf)

    def pickleDump(self,dir=dir):
        """store object in file"""
        if dir and 'dumpf' not in self.__dict__.keys(): self.genDumpf(dir)
        pickle.dump( self,open(self.dumpf,'w'), True )  

    def for_write_report(self,attribute):
        buf=' '+`attribute`
        return buf

    def update(self,attribute,value):
        """initialize or update the attribute with the passed value"""
        #create attributes. Go into dictionary levels if neccessary
        levels=attribute.split('.')
        ref=self.__dict__
        for level in levels[0:-1]:
            if level not in ref.keys(): ref[level]={}
            ref=ref[level]
        #assign value
        ref[levels[-1]]=value

def loadMac(macheader,collf=collf,file=''):
    """loadMac(macheader,collf=collf,file='')
    load a MAC object either from the collection-file or from a file
    containing the object"""
    if file:
        return pickle.load( open(file,'r') )
    #find index in indexfile
    index=0
    for line in open(collf+'.idx','r').readlines():
        if macheader in line:
           index=int( line.split()[1] )
           break
    #fetch the MAC from collf
    pin=open(collf,'r')
    pin.seek(index)
    return pickle.load(pin)

def loadMACs(root=scout,list=maclist,collf=''):
    """loadMACs(root=scout,list=maclist,collf='')    
    handy function to memory-load all MAC objects given a file-list of MAC headers"""
    macs={}
    #read from collection file if passed
    if collf:
        collp=open(collf,'r')
        while True:
            try:
                mac=pickle.load(collp)
                macs[mac.macheader]=mac
            except: break
        return macs
    #read from individual files
    headers=list
    if isinstance(list,str):
        headers=chomp(open(list,'r').readlines())
    for header in headers:
        macfile=os.path.join(scout,header.split('.')[0],header+'.mac')
        if os.path.exists(macfile):
            macs[header]=pickle.load( open(macfile,'r') )
        else:
            sys.stderr('ERROR in mac.loadMACs: file "'+macfile+'" does not exists\n')
    Bye(len(macs))
    return macs

def updColletion(collf=collf,list=maclist):
    """updColletion(coll=collf,list=maclist)    
    update the single file containing all pickles"""
    buf='' #store indexes
    macs=loadMACs(list=list) #load all MACs into memory
    pout=open(collf,'w')
    for macheader in macs.keys():
        buf+=macheader+' '+`pout.tell()`+'\n'
        pickle.dump(macs[macheader],pout,True)
    open(collf+'.idx','w').write(buf) #write index file
    
def write_report(attrl,outf='',collf=collf):
    """write_report(attrl,outf)
    For each MAC in coll, write to file outf the list of attributes of
    attrl on a single line. By default, every line begins with
    macheader attribute, irrespective of whether we list it in
    attrl"""
    #write first line naming the attributes
    buf='# 1-macheader'
    attr_i=2
    for attr in attrl:
        buf+=' '+`attr_i`+'-'+attr
        attr_i+=1
    buf+='\n'
    #write info for each MAC
    macs=loadMACs(collf=collf)
    for mac in macs.values():
        for attr in attrl:
            buf+=mac.for_write_report(attr)
        buf+='\n'

            
def help(*kargs):
    """help(sv)
    list of services"""
    services=['updColletion','write_report']
    if not kargs:
        print 'services=',services
        print 'Type "mac.py -a help -b sv=service" for more detail on "service"'
    else:
        sv=kargs[0]
        print globals()[sv].__doc__        


if __name__=='__main__':
    inpHand('Usage: mac.py [options]',
            ' -a _R_service required service. Type "help" for list of services.',              
            ' -b __servargs service arguments in python syntax. Escape symbols \' and " (def:None)',
            ).parse(locals(),sys.argv)

    if not servargs:
        locals()[service]()
    else:
        locals()[service](servargs)
