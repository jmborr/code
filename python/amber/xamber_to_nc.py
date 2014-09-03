#!/usr/bin/python

import os
from Tkinter import *
from tkFont import Font
import FileDialog
from amber.amber2MMTK import getBox,fixPDB
from amber.amber_to_nc import usage,convertAmberToMMTK
from nMoldyn.gui import *

isthere=os.path.exists

class xamber_to_nc(Frame):

    def __init__(self,master):

        Frame.__init__(self,master)
        self.pack()
        self.master=master
        self.Default()
        self.CreateMenu()
        self.createChoiceWindow()

    def Default(self):
        self.pdb=None
        self.fixedpdb=None
        self.natoms=None
        self.prmcrd=None
        self.prmvel=None
        self.box=[0.0,0.0,0.0]
        self.dtime=None

    def CreateMenu(self):
        menu_bar=Menu(self)
        self.createFileMenu(menu_bar)
        self.createHelpMenu(menu_bar)
        self.master.config(menu=menu_bar)

    def createFileMenu(self,menu_bar):
        menu=Menu(menu_bar,tearoff=0)
        menu_bar.add_cascade(label="File", menu=menu)
        menu.add_command(label='Import PDB',
                         command=self.import_pdb)
        menu.add_command(label='Import prmcrd',
                         command=self.import_prmcrd)
        menu.add_command(label='Import prmvel',
                         command=self.import_prmvel)
        menu.add_separator()
        menu.add_command(label='Quit',command=self.quit)
        self.file_menu=menu

    def createHelpMenu(self,menu_bar):
        menu=Menu(menu_bar,tearoff=0)
        menu_bar.add_cascade(label="Help", menu=menu)
        menu.add_command(label='Contents',
                         command=self.HelpContents)
        menu.add_command(label='About Amber to MMTK',
                         command=self.aboutInfo)
        self.help_menu=menu

    def ConvertWindow(self):
        """open a new frame"""
        
    def aboutInfo(self):
        aboutmesg="""xamber_to_nc v1.0 written by:
Jose Borreguero.
Special thanks to:
Konrad Hinsen,
Paolo Calligari,
Krzysztof Murzyn
"""        
        about=Toplevel(self)
        about.title('About xamber_to_nc')
        f0=Frame(about,bd=2,relief='groove')
        f0.pack(side=TOP,padx=3,pady=3,fill=BOTH)
        Label(f0,text=aboutmesg).grid(column=0,row=0)
        f2 = Frame(about,bd=2,relief='groove')
        f2.pack(side=TOP,fill=X,padx=3,pady=3)
        Button1=Button(f2,text='Close',command=about.destroy,underline=0)
        Button1.pack(padx=1,pady=1,side=RIGHT)
        about.resizable(width=NO,height=NO)
        about.initial_focus = about
        about.grab_set()
        about.initial_focus.focus_set()
        about.wait_window(about)

    def HelpContents(self):
        contents=Toplevel(self)
        contents.title('Contents')
        f0=Frame(contents,bd=2,relief='groove')
        f0.pack(side=TOP,padx=3,pady=3,fill=BOTH)
        Label(f0,text=usage).grid(column=0,row=0)
        f2 = Frame(contents,bd=2,relief='groove')
        f2.pack(side=TOP,fill=X,padx=3,pady=3)
        Button1=Button(f2,text='Close',
                       command=contents.destroy,underline=0)
        Button1.pack(padx=1,pady=1,side=RIGHT)
        contents.resizable(width=NO,height=NO)
        contents.initial_focus = contents
        contents.grab_set()
        contents.initial_focus.focus_set()
        contents.wait_window(contents)
        
    def createChoiceWindow(self):
        frame=Frame(self, bd=2)
        frame.pack(side=TOP, fill=Y)
        f11=Frame(frame,bd=2,relief='groove')
        f11.pack(side=TOP, padx=7, pady=7)
        
        f13=Frame(f11)
        f13.pack(side=TOP, padx=7, pady=7)

        l0 = Label(f13,text='Required:',anchor=W)
        l0.grid(column=0,row=0,sticky='news',pady=7)
        self.myFont = Font(font=l0["font"]).copy()
                
        l1 = Label(f13,text='prmpdb',anchor=W)
        l1.grid(column=0,row=1,sticky='news',pady=7)
        self.p41=FilenameEntry(f13,type_file='PDB',browse_pattern='*.pdb')
        self.p41.grid(column=1,row=1,pady=7)
        
        l2 = Label(f13,text='prmcrd',anchor=W)
        l2.grid(column=0,row=2,sticky='news',pady=7)
        self.p42=FilenameEntry(f13,type_file='Amber',
                               browse_pattern='*.crd')
        self.p42.grid(column=1,row=2,pady=7)
        
        l3 = Label(f13,text='nc_file',anchor=W)
        l3.grid(column=0,row=3,sticky='news',pady=7)
        self.p43=FilenameEntry(f13,type_file='MMTK',
                               browse_pattern='*.nc')
        self.p43.grid(column=1,row=3,pady=7)

        l4 = Label(f13,text='Options:',anchor=W)
        l4.grid(column=0,row=4,sticky='news',pady=7)
        
        l5 = Label(f13,text='prmvel',anchor=W)
        l5.grid(column=0,row=5,sticky='news',pady=7)
        self.e44=FilenameEntry(f13,type_file='Amber',browse_pattern='*.vel')
        self.e44.grid(column=1,row=5,pady=7)

        l6 = Label(f13,text='Box',anchor=W)
        l6.grid(column=0,row=6,sticky='news',pady=7)

        self.e41 = FloatEntry(f13,'', self.box[0])
        self.e41.grid(column=1,row=6,pady=7)
        self.e41.bind('<Return>',lambda event,f=frame,s=self: s._updBox(f))
        self.t41 = Label(f13, text='(x, Ang)',anchor=W)
        self.t41.grid(column=2, row=6, pady=7)

        self.e42 = FloatEntry(f13,'', self.box[1])
        self.e42.grid(column=1,row=7,pady=7)
        self.e42.bind('<Return>',lambda event,f=frame,s=self: s._updBox(f))
        self.t42 = Label(f13, text='(y, Ang)',anchor=W)
        self.t42.grid(column=2, row=7, pady=7)

        self.e43 = FloatEntry(f13,'', self.box[2])
        self.e43.grid(column=1,row=8,pady=7)
        self.e43.bind('<Return>',lambda event,f=frame,s=self: s._updBox(f))
        self.t43 = Label(f13, text='(z, Ang)',anchor=W)
        self.t43.grid(column=2, row=8, pady=7)

        l9 = Label(f13,text='dtime:',anchor=W)
        l9.grid(column=0,row=9,sticky='news',pady=7)
        self.e44 = FloatEntry(f13,'', self.dtime)
        self.e44.grid(column=1,row=9,pady=7)
        self.t44 = Label(f13, text='(ps)',anchor=W)
        self.t44.grid(column=2, row=9, pady=7)

                
    def import_pdb(self):
        fd = FileDialog.LoadFileDialog(self)
        pdb=fd.go(key='pdb')
        if not isthere(pdb): return
        self.pdb=pdb
        #pop prompting there will be a try to fix the PDB
        self.fixedPDB,self.fixes,self.natoms=fixPDB(pdb)
        print 'Created'+self.fixedPDB
        
    def import_prmcrd(self):
        fd = FileDialog.LoadFileDialog(self)
        prmcrd=fd.go(key='prmcrd')
        if not isthere(prmcrd): return
        self.prmcrd=prmcrd
        
    def import_prmvel(self):
        fd = FileDialog.LoadFileDialog(self)
        prmvel=fd.go(key='prmvel')
        if not isthere(prmvel): return
        self.prmvel=prmvel



root=Tk()
root.title('Amber to MMTK')
x=xamber_to_nc(root)
x.mainloop()
