#!/usr/bin/python2
from jobs.jobeditor import jobeditor
from jobs.job import job
from utilities.codedir import codedir,libseq
from Tkinter import mainloop
db=codedir+'/python/db/job'
jobeditor("Jobs",job,db)
mainloop()
