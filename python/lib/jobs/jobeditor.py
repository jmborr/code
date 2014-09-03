from Tkinter import *
import tkMessageBox
import string,os,pickle,sys,re
from string import join
from utilities.small_utilities import chomp,todayDate
from utilities.codedir import projectsdir,scratchdir
from random import randint

class jobeditor:
    def __init__(self,name,jobclass,storagedir):
        self.storagedir=storagedir  #stash away some references
        self.jobclass=jobclass
        self.row=0
        self.current=None
        self.root=root=Tk()  #create window and size it
        root.minsize(1000,700)

        root.rowconfigure(0,weight=1)    #define how columns and rows scale
        root.columnconfigure(0,weight=0) #when the window is resized
        root.columnconfigure(1,weight=1)
        
        self.frame=Frame(root)
        self.frame.grid(columnspan=2,sticky=E+W+N+S)
        self.scrollbar = Scrollbar(self.frame, orient=VERTICAL)
        self.listbox=Listbox(self.frame, yscrollcommand=self.scrollbar.set,
                             selectmode=SINGLE) #main listbox
        self.listbox.grid(row=self.row,column=0,sticky=E+W+N+S)
        self.scrollbar.config(command=self.listbox.yview)
        self.scrollbar.grid(row=self.row,column=1,sticky=E)
        self.scrollbar.pack(side=RIGHT, fill=Y)
        self.listbox.pack(side=LEFT, fill=BOTH, expand=1)
        #action self.select upon clicking with the mouse on the listbox
        self.listbox.bind('<ButtonRelease-1>',self.select)
        self.row=self.row+1

        for field in jobclass.fields: #creat attributes
            setattr(self,field,self.add_variable(root,field))

        #buttons in a separate window
        commW=Toplevel()

        self.title=Label(commW,text='ACTIONS', font=("Helvetica",12),height=1)
        self.title.grid(row=0,columnspan=4,sticky=E+W)

        self.add_button(commW,1,0,'Search',self.doSearch)
        self.add_button(commW,1,1,'Del entry',self.delentry)
        self.add_button(commW,1,2,'Reload All',self.load_data)
        self.add_button(commW,1,3,'Del File',self.delfile)

        self.add_button(commW,2,0,'args help',self.args_help)
        self.add_button(commW,2,1,'change',self.change)
        self.add_button(commW,2,2,'rewrite',self.rewrite)
        self.add_button(commW,2,3,'save as new',self.newsave)

        self.add_button(commW,3,0,'qsub',self.qsub)
        self.add_button(commW,3,1,'bsub',self.bsub)
        self.add_button(commW,3,2,'sub' ,self.sub )
        self.add_button(commW,3,3,'rsub',self.rsub)

        self.add_button(commW,4,0,'quit',self.quit)
                
        self.load_data()
        
    def add_variable(self,root,varname):
        Label(root,text=varname,font=("Times", 16, "bold italic") ).grid(row=self.row,column=0,sticky=E)
        value=Label(root,text='',font=("Helvetica", 16),
                    background='gray90',relief=SUNKEN,anchor=W,
                    justify=LEFT,wraplength=800)
        value.grid(row=self.row,column=1,sticky=E+W)
        self.row=self.row+1
        return value

    def add_button(self,root,row,column,text,command):
        button=Button(root,text=text,command=command)
        button.grid(row=row,column=column,sticky=E+W,padx=5,pady=5)

    def load_data(self):
        self.listbox.delete(0,END) #delete is a method of a Listbox object
        #will store all instantianted jobclass objects stored in storagedir
        self.items=[]
        #for filename in os.listdir(self.storagedir):
        for filename in chomp(os.popen('ls -1 -t '+self.storagedir).readlines()):
            #item is an object of jobclass
            item=pickle.load(open(os.path.join(self.storagedir,filename)))
            self.items.append(item)
            #`item` is the string representation of object item. This is
            #what we get to see
            self.listbox.insert(END,`item`)#insert is a Listbox method
            #select_set is a Listbox method, highlight the first item in
            #self.listbox
            self.listbox.select_set(0) 
            self.select(None)

    def refresh(self):
        self.listbox.delete(0,END) #delete is a method of a Listbox object
        for item in self.items:
            self.listbox.insert('end',`item`)
            self.listbox.select_set(0)
            self.select(None)

    def select(self,event):
        #selection is an object of an unknown class to us
        selection=self.listbox.curselection()
        #int(selection[0]) is the index of the corresponding entry in
        #self.listbox, thus self.items[..] returns the appropriate
        #jobclass jobject
        self.selection=self.items[int(selection[0])]
        #go through all the names of attributes
        self.refresh_selection()

    def refresh_selection(self):
        for field in self.jobclass.fields:
            #reference to the Label object wich is also referenced by
            #the attribute of self whose name is stored in field
            label=getattr(self,field)
            #reference to the string object which is also referenced by
            #the attribute of self.selection whose name is stored in field
            try:
                labelstr=getattr(self.selection,field)
            except:
                continue
            #windows operating system does not handle correctly '\r'
            labelstr=string.replace(labelstr,'\r','')
            #a method of Label object to refresh its 'text' attribute
            label.config(text=labelstr)

    def delfile(self):
        if tkMessageBox.askyesno("Print", "Delete File?"):
            os.remove(self.selection.file)
            self.delentry()
    
    def qsub(self):
        rn='a'+`randint(0,randint(0,9999999))` #7 digits (log name 8 characters long)
        self.qop={'queue':'batch','tmpd':'/tmp/jose','mem(mb)':'300',
                  'jobname':rn,'outdir':projectsdir+'/tmp','logname':rn+'.log',
                  'outputmove':'','extraflags':''}                
        newW=self.modifyDict('qsub options','modify')
        self.root.wait_window(newW)
        if tkMessageBox.askyesno("Print", "submit to queueing system?"):
            #for key in self.qop.keys(): print 'qop['+key+']='+self.qop[key]
            qs=self.selection.qsub
            qop=self.qop
            qs(queue=qop['queue'],tmpd=qop['tmpd'],mem_limit=qop['mem(mb)'],jobname=qop['jobname'],
               outdir=qop['outdir'],logname=qop['logname'],outputmove=qop['outputmove'],
               extraflags=qop['extraflags'])

    def bsub(self):        
        if tkMessageBox.askyesno("Print", "submit to queueing system?"):
            self.selection.bsub()

    def sub(self):
        if tkMessageBox.askyesno("Print", "run job locally?"):
            jobname='gui'+chomp(os.popen('date +%H%M%S').readline())
            outdir=scratchdir+'/qsub/'+todayDate()
            self.selection.sub(jobname,outdir,libRev='',redirect=0)

    def rsub(self):
        if tkMessageBox.askyesno("Print", "submit directly to some node?"):
            self.selection.rsub()

    def modifyDict(self,title,buttontext,width=600,height=200):
        newW=Toplevel()
        self.focus=newW
        row=0
        self.entries={}
        newW.minsize(width,height)
        newW.rowconfigure(0,weight=1)    #define how columns and rows scale
        newW.columnconfigure(0,weight=0) #when the window is resized
        newW.columnconfigure(1,weight=1)
        Label(newW,text=title,font=("Helvetica",16)).grid(columnspan=2) 
        row=row+1
        for key in self.qop.keys():
            self.entries[key]=self.add_entry(newW,key,row,self.qop[key])
            row=row+1
        button=Button(newW,text=buttontext,command=self.doModifyDict)
        button.grid(row=row,column=0,sticky=E+W,padx=5,pady=5)
        button2=Button(newW,text='cancel',command=self.cancel)
        button2.grid(row=row,column=1,sticky=W,padx=5,pady=5)        
        return newW
    
    def newEntryW(self,title,command,buttontext,mode,width=300,height=200):
        newW=Toplevel()
        self.focus=newW
        row=0
        self.entries={}
        newW.minsize(width,height)
        newW.rowconfigure(0,weight=1)    #define how columns and rows scale
        newW.columnconfigure(0,weight=0) #when the window is resized
        newW.columnconfigure(1,weight=1)
        Label(newW,text=title,font=("Helvetica",16)).grid(columnspan=2) 
        row=row+1
        for field in self.jobclass.fields: #create attributes
            if mode==1:   #do a search
                self.entries[field]=self.add_entry(newW,field,row,'')
            elif mode==2: #do changes
                label=getattr(self.selection,field)
                self.entries[field]=self.add_entry(newW,field,row,label)
            row=row+1
        button=Button(newW,text=buttontext,command=command)
        button.grid(row=row,column=0,sticky=E+W,padx=5,pady=5)
        button2=Button(newW,text='cancel',command=self.cancel)
        button2.grid(row=row,column=1,sticky=W,padx=5,pady=5)

        
        row=row+1
        return newW
    
    def doSearch(self):
        newW=self.newEntryW('search by REGEXP',self.sendSearch,'submit',1)

    def change(self):
        newW=self.newEntryW('change entry',self.dochange,
                            'do changes',2,width=1265)
    
    def add_entry(self,newW,varname,row,text):
        Label(newW,text=varname).grid(row=row,column=0,sticky=E)
        value=Entry(newW,background='gray90',font=("Helvetica",16),relief=SUNKEN)
        value.grid(row=row,column=1,sticky=E+W)
        value.insert(INSERT,text)
        return value

    def dochange(self):
        for field in self.jobclass.fields:
            setattr(self.selection,field,self.entries[field].get())
        self.refresh_selection()

    def doModifyDict(self):
        self.qop={}
        for key in self.entries.keys():
            self.qop[key]=self.entries[key].get()
        self.focus.destroy()
        
    def sendSearch(self):
        for varname in self.jobclass.fields:
            value=self.entries[varname].get()
            p=re.compile(value)
            x=[]
            if value:
                l=len(self.items)
                for i in range(l): #search item by item
                    field=getattr(self.items[i],varname)
                    if not p.search(field):
                        x.append(i)
            for i in x: #remove those with no matches
                del self.items[i]
        self.refresh()

    def delentry(self):
        selection=self.listbox.curselection()
        index=int(selection[0])
        del self.items[index]
        self.listbox.delete(ACTIVE)
        self.listbox.select_set(index)

    def rewrite(self):
        pickle.dump(self.selection,open(self.selection.file,'w'))
    
    def newsave(self):
        self.selection.store(self.storagedir)

    def args_help(self):
        newW=Toplevel()
        newW.minsize(300,200)
        newW.rowconfigure(0)
        newW.columnconfigure(0) #when the window is resized
        Label(newW,text=join(self.selection.help(),''),justify='left',
              font=("Helvetica", 16)).grid(columnspan=1)

    def quit(self):
        self.root.destroy()

    def cancel(self):
        self.focus.destroy()

