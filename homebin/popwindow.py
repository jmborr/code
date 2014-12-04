#!/usr/bin/python

import sys
import os
import argparse
from Tkinter import *
from pdb import set_trace as tr

class App:

    def __init__(self, master, message):

        frame = Frame(master)
        frame.pack()

        self.button = Button(frame, height=16, width=16, text="QUIT",
                             fg="red", bg="yellow", command=frame.quit)
        self.button.pack(side=LEFT)
        
        self.hi_there = Button(frame, text=message,command=self.say_hi)
        self.hi_there.pack(side=LEFT)
        
    def say_hi(self):
        print 'Don\'t tickle me!'

parser = argparse.ArgumentParser(description='Pops a window with a message')
parser.add_argument('sleep',help='argument to Unix sleep command, for a delayed message')
parser.add_argument('message',help='Quote message if it contains more than one word')

args=parser.parse_args()
if args.sleep: os.system('sleep %s'%args.sleep)
root = Tk()
root.geometry("%dx%d%+d%+d" % (600, 400, 0, 0))
app = App(root, args.message)
root.mainloop()

