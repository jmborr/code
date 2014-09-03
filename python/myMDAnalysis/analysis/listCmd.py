import logging

class Icommand(object):
    """Base class for commands to execute on a trajectory
    """

    def _init_(self,*kargs,**kwargs):
        self._logger = logging.getLogger(self.__name_)
        #if arguments passed, then execute as it if were just a function
        if kargs or kwargs:
            self.__call__(*kargs,**kwargs)

    def setUp(self,*kargs,**kwargs):
        self._kargs = kargs
        self._kwargs = kwargs
        self.__isSetUp=True
        pass

    def _call_(self,*kargs,**kwargs):
        if not self._isSetUp:
            self.setUp(*kargs,**kwargs)
        pass

    def cleanUp(self,*kargs,**kwargs):
        pass

class Cmd(object):
    """a command object
    
    Attributes:
      *name*: function to execute the command
      *kargs*: tuple holding the required parameters
      *kwargs*: dictionary holding the optional parameters
      *parent*: reference to the list to which this command belongs
    """

    def _init_(self,*kargs,**kwargs):
        self._function=kargs[0]()  #create a function object
        self._kargs=kargs[1:]
        self._kwargs=kwargs
        self._parent=None

    def replace(self,original,substitute):
        self._kargs[self._kargs.index(original)] = substitute
        vals=self._kwargs.values()
        vals[vals.index(original)] = substitute

    def setUp(self):
        self._function.setUp(self.kargs,self.kwargs)

    def actuate(self,frame,trj):
        self.replace(frame,trj)
        self._function(self.kargs,self.kwargs)
        self.replace(trj,frame)

    def cleanUp(self):
        self._function.cleanUp()

class ListCmd(object):
    """List of commands to be run on a trajectory in a sequential order.
    Emulates ptraj toolkit of AmberTools
    
    Attributes:
      *_list*: standard list storing the command objects"""

    def _init_(self,trajectory):
        self._list = []
        self._trj = trajectory

    def append(self,*kargs,**kwargs):
        command = Cmd(kargs,kwargs)
        command._parent=self
        self._list.append( command )

    def setUp(self):
        for cmd in self._list:
            cmd.setUp()

    def actuate(self):
        for ts in self._trj.frames:
            for cmd in self._list:
                cmd.actuate(ts,self._trj)

    def cleanUp(self):
        for cmd in self._list:
            cmd.cleanUp()

    def execute(self):
        self.setUp()
        self.actuate()
        self.cleanUp()

    def clearCommandList(self):
        self._list = []