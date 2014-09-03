#!/usr/bin/python
"""Program that:
(1) library for class CAN, a container of whatever list of data capsules we may want to put in
(2) provides several services related to CAN objects when executed as
    standalone
"""
import pickle,os,sys
from utilities.codedir import scratchdir
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,Bye
from jobs.job import pastry

#global variables
joink=os.path.join #a handy reference to this method

def addtotree(parent,dict,treestore):
    """iteratively add nested dictionaries to gtk.TreeStore object"""
    for key in dict.keys():
        val=dict[key]
        thetype=type(val).__name__
        if thetype=='dict':
            child=treestore.append(parent,[key])
            addtotree(child,val,treestore) #iterate
        elif thetype in ('bool','int','float','str','list','tuple',
                         'NoneType','ndarray'):
            if thetype in ('float'):
                child=treestore.append(parent,[repr(key)+' : '+repr(str(val))])
            elif thetype in ('list'): #print first and last elements only
                if val:
                    reduced=repr([repr(val[0]),'..',repr(val[-1])])
                else:
                    reduced=repr(val)
                chile=treestore.append(parent,[repr(key)+' : '+reduced])
            else:
                child=treestore.append(parent,[repr(key)+' : '+repr(val)])
        else:
            child=treestore.append(parent,[key])
            addtotree(child,val.__dict__,treestore) #iterate

def delete_event(widget, event, data=None):
    import gtk
    """exit the window viewer"""
    gtk.main_quit()
    return False
        

class CAN:
    """missing atoms chunk (CAN) object implementation
    a CAN is a container to hold ever increasing info

    Required attributes:
    'id': string identifying the CAN object
    'dumpf': picke dump file
    'dumpbf': basename of pickle dump file
    'dumpdf': dirname of pickle dump file
    Recommended attributes:
    'capsule': dictionary holding all info
      |_'yadayada': one property with literal name 'yadayada'
          |_'otheryadayada'
      |_header: one property with unknown literal name
    """
    def __init__(self,**kwargs):
        self.__dict__=kwargs.copy() #shallow copy


    def genDumpf(self,dir):

        """generate name of picke dump file, as well as dirname and
        basename"""

        if 'id' not in self.__dict__.keys(): self.genId()
        self.dumpdf=dir
        self.dumpbf=self.id+'.can'
        self.dumpf=joink(dir,self.dumpbf)

    def pickleDump(self,dumpfile=''):

        """store object in file"""

    def pickleDump(self,dumpfile=''):

        """store object in file"""

        if not dumpfile:
            if 'dumpf' not in self.__dict__.keys(): self.genDumpf(dir)
            dumpfile=self.dumpf;
            try:
                pickle.dump( self,open(dumpfile,'w'), True )
            except:
                os.system('/bin/cp %s.bak %s'%(dumpfile,dumpfile))
                Bye('ERROR: could not pickle dump '+self.id+'\n')
            try:
                junkobj=pickle.load( open(dumpfile) )
            except:
                os.system('/bin/cp %s.bak %s'%(dumpfile,dumpfile))
                Bye('ERROR: failed load %s after pickle dump. '+\
                    'Reverting to old can file\n'%(self.id,))
            pastry('/bin/cp '+dumpfile+' '+dumpfile+'.bak') #backup


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


    def tview(self):
        """representation of object output to terminal"""
        from pprint import PrettyPrinter
        pp=PrettyPrinter(indent=2)
        pp.pprint(self.__dict__)
        
    def view(self):

        """primitive viewer of object properties"""
        
        import pygtk
        pygtk.require('2.0')
        import gtk

        #Store object info in a gtk.TreeStore object
        treestore = gtk.TreeStore(str)
        treeview = gtk.TreeView(treestore) # create the TreeView
        tvcolumn = gtk.TreeViewColumn('properties') #create the TreeViewColumn to display the data
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText() #create a CellRendererText to render the data
        tvcolumn.pack_start(cell, True) # add the cell to the tvcolumn and allow it to expand
        tvcolumn.add_attribute(cell, 'text', 0) #set the cell "text" attribute to column 0
        treeview.set_search_column(0) #make it searchable
        tvcolumn.set_sort_column_id(0) #Allow sorting on the column
        addtotree(None,self.__dict__,treestore)
        #Insert info into a gtk.Window object
        window =gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_title(self.id)
        window.set_default_size(900,700)
        window.set_position(gtk.WIN_POS_CENTER)
        window.connect('delete_event',delete_event)
        vbox =gtk.VBox(False, 8)
	window.add(vbox)
        scwin=gtk.ScrolledWindow()
        scwin.add(treeview)
        vbox.pack_start(scwin)
        window.show_all()
        gtk.main()
        
        
class CANl:
    
    def __init__(self,collf,canlistf,repod,other={}):
        """__init__(self,collf,canlist,repod)
        collf: collection file generated by pickle module containing all CAN objects.
        canlist: list of CAN's. For example, list of query sequences
        repod: CAN objects are also stored inside a directory, each in a dump
               file. This allows one to update one CAN object instead of whole
               collf file"""
        self.collf=collf
        self.canlistf=canlistf
        self.canlist=chomp(open(canlistf,'r').readlines())
        self.repod=repod
        if other: #additional attributes
            for key,value in other.items(): self.__dict__[key]=value


    def loadCAN(self,id,collf='',file=''):

        """loadCAN(id,collf='',file='') load a CAN object either from
        the dump file (default behavior), the collection-file, or from
        a specified file containing the object"""
        
        if file:
            try:
                object=pickle.load( open(file,'r') )
                msg='ERROR (loadCAN) while loading '+\
                     file+'\n      resorting to backup file...\n'
                sys.stderr.write(msg)
            except:
                backupfile=file+'.bak' #resort to backup file
                if not os.path.exists(backupfile):
                    msg='ERROR (loadCAN) backup does not exists!\n'
                    sys.stderr.write(msg)
                    return None
                os.system('/bin/cp '+file+'.bak '+file) #restore backup
                object= pickle.load( open(file,'r') )
            return object
        if collf:
            index=0 #find index in indexfile
            for line in open(self.collf+'.idx','r').readlines():
                if id in line:
                    index=int( line.split()[1] )
                    break
            pin=open(collf,'r') #fetch the CAN from collf
            pin.seek(index)
            return pickle.load(pin)
        file=joink(self.repod,id+'.can') #;print file
        if not os.path.exists(file):
            sys.stderr.write('ERROR file '+file+' does not exists\n')
            return None
        try:
            #print 'loading '+file
            return pickle.load( open(file,'r') ) #dump file
        except:
            msg='ERROR while loading '+os.path.basename(file)+\
                 '\nloading the backup file\n'
            sys.stderr.write(msg)
            try:
                file+='.bak'
                return pickle.load( open(file,'r') )
            except:
                msg='could not load backup file '+\
                     os.path.basename(file)+' either!'
                Bye(msg)


    def loadCANs(self,repod='',list=[],collf=''):

        """loadCANs(repod='',list=[],collf='') handy function to
        memory-load all CAN objects given a list of CAN
        headers. Default behaviour is to read all dump files"""

        if not repod: repod=self.repod
        if not list: list=self.canlist
        cans={}
        #read from collection file if passed
        if collf:
            collp=open(collf,'r')
            while True:
                try:
                    can=pickle.load(collp)
                    cans[can.id]=can
                except: break
            return cans
        #read from individual files
        headers=list
        if isinstance(list,str):
            headers=chomp(open(list,'r').readlines())
        for header in headers:
            canfile=joink(repod,header+'.can')
            if os.path.exists(canfile):
                cans[header]=pickle.load( open(canfile,'r') )
            else:
                sys.stderr('ERROR in CANl.loadCANs: file "'+canfile+'" does not exists\n')
                Bye(len(cans))
        return cans


    def updColletion(self,collf='',list=[]):

        """updColletion(coll='',list='')    
        update the single file containing all pickles"""

        if not collf: collf=self.collf
        if not list: list=self.canlist
        buf='' #store indexes
        cans=self.loadCANs(list=list) #load all CANs into memory
        pout=open(collf,'w')
        for id in cans.keys():
            buf+=id+' '+`pout.tell()`+'\n'
            pickle.dump(cans[id],pout,True)
        open(collf+'.idx','w').write(buf) #write index file

    
    def write_report(self,attrl,outf='',collf=''):

        """write_report(attrl,outf='',collf='')
        For each CAN in coll, write to file outf the list of attributes of
        attrl on a single line. By default, every line begins with
        id attribute, irrespective of whether we list it in
        attrl"""

        if not collf: collf=self.collf
        #write first line naming the attributes
        buf='# 1-id'
        attr_i=2
        for attr in attrl:
            buf+=' '+`attr_i`+'-'+attr
            attr_i+=1
        buf+='\n'
        #write info for each CAN
        cans=self.loadCANs(collf)
        for can in cans.values():
            for attr in attrl:
                buf+=can.for_write_report(attr)
            buf+='\n'


    def iterate(self,list='',collf=''):

        """iterator over CAN list

        iterate(self,list='',collf='') iterates over all dump files or
        over a subset of them, given by list. If we want to iterate
        over the collection file instead of the dump files, then we
        pass collf """

        if not list:list=self.canlist
        for header in self.canlist:
            yield self.loadCAN(header,collf=collf)


    def cycleOverList(self,methodname,*kargs,**kwargs):
        """cycleOverList(methodname,*kwargs,**kwargs)
        
        Execute the CAN method with name 'methodname' and arguments kwargs, kwargs
        iteratively on each varset. Parse some arguments withing kwargs that are not part
        of the arguments that should be passed to methodname
        By default we iterate over all ID's. We can restrict which ID's with optional parameters:
        id: pass only one ID
        idlist: pass a list of ID's
        idlistf: pass file name that is a list of ID's
        The optional parameter we pass them inside **kwargs"""
        from time import sleep
        silent=False
        if 'silent' in kwargs.keys():
            if str(kwargs['silent'])=='True': silent=True
            del kwargs['silent']
        pause=0.0
        if 'pause' in kwargs.keys():
            pause=float(kwargs['pause'])
            del kwargs['pause']
        idlist=self.canlist
        if 'id' in kwargs.keys():
            idlist=(kwargs['id'],)
            del kwargs['id']
        elif 'idlist' in kwargs.keys():
            idlist=kwargs['idlist'].split() #split by blank space
            del kwargs['idlist']
        elif 'idlistf' in kwargs.keys():
            idlist=chomp(open(kwargs['idlistf'],'r').readlines())
            del kwargs['idlistf']
        remaining=len(idlist)
        retVals={}
        for id in idlist: #cycle through every CAN object
            if not silent: print remaining,id
            c=self.loadCAN(id) #load from corresponding dump file
            retVal=getattr(c,methodname)(*kargs,**kwargs)
            if type(retVal)==type({}):
                dump=retVal['dump']
                del retVal['dump']
                retVals[c.id]=retVal
            else:
                dump=retVal
            if dump==True:
                #print 'update dump file'
                c.pickleDump()
            remaining-=1 ; sleep(pause)
        return retVals


    def cycleOverList2(self,**kwargs):
        """
        kwargs={'local':{ 'kargs'=[], 'kwargs'={'silent'=True,} },
                'method'={ 'name':'foo', 'kargs'=[], 'kwargs'={} }
                }
        """
        silent=kwargs['local']['kwargs']['silent'] #;Bye(silent)
        methodname=kwargs['method']['name'] #;Bye(methodname)
        try:    largs=kwargs['method']['kargs']
        except: largs=[]
        try:    lwargs=kwargs['method']['kwargs']
        except: lwargs={}
        remaining=len(self.canlist)
        for header in self.canlist: #cycle through every CAN object
            if not silent: print remaining,header
            c=self.loadCAN(header) #load from corresponding dump file
            x=getattr(c,methodname)(*largs,**lwargs)
            if x!=False: c.pickleDump() #update dump file only if returned anything but False
            remaining-=1


    def view(self):

        """primitive viewer of object properties"""
        
        import pygtk
        pygtk.require('2.0')
        import gtk

        #Store object info in a gtk.TreeStore object
        treestore = gtk.TreeStore(str)
        treeview = gtk.TreeView(treestore) # create the TreeView
        tvcolumn = gtk.TreeViewColumn('properties') #create the TreeViewColumn to display the data
        treeview.append_column(tvcolumn)
        cell = gtk.CellRendererText() #create a CellRendererText to render the data
        tvcolumn.pack_start(cell, True) # add the cell to the tvcolumn and allow it to expand
        tvcolumn.add_attribute(cell, 'text', 0) #set the cell "text" attribute to column 0
        treeview.set_search_column(0) #make it searchable
        tvcolumn.set_sort_column_id(0) #Allow sorting on the column
        for c in self.iterate(): addtotree(None,self.__dict__,treestore)
        #Insert info into a gtk.Window object
        window =gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_title(self.id)
        window.set_default_size(900,700)
        window.set_position(gtk.WIN_POS_CENTER)
        window.connect('delete_event',delete_event)
        vbox =gtk.VBox(False, 8)
	window.add(vbox)
        scwin=gtk.ScrolledWindow()
        scwin.add(treeview)
        vbox.pack_start(scwin)
        window.show_all()
        gtk.main()
        
        
