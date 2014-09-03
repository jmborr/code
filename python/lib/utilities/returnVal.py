#!/usr/bin/python

import sys

class returnVal:

    def __init__(self,st=0,msg='',errMsg='',val=None):
        self.st=st         #error status
        self.msg=msg       #message
        self.errMsg=errMsg #error message
        self.val=val       #returned object

    def clearErr(self): #set error status to OK, clear the error message
        self.st=0
        self.errMsg=''
        
    def hasErr(self,): #return 1 if error status is set
        if self.st!=0: return 1
        return None

    def setErr(self,status=1,msg=''): #set the error status and include an error message
        self.st=status
        self.errMsg=msg

    def abort(self):
        sys.stdout.write(self.msg+'\n')
        sys.stderr.write(self.errMsg+'\n')
        sys.exit(self.st)
    
