def upl(fname): 
	f=open(fname,'r')
	i=1
	upl=f.readline()
	upls=f.readlines()
	for upl in upls:		
		cns=upl.split()
		print 'cmd.distance('+'\'restr'+str(i),'\',\'resi '+cns[0]+' and name ca\',\'resi '+cns[1]+' and name ca\''+')'
#		cmd.distance('\'restr'+str(i),'\',\'resi '+cns[0]+' and name ca\',\'resi '+cns[1]+' and name ca\'')
		cmd.distance('restr'+str(i),'resi '+cns[0]+' and name ca','resi '+cns[1]+' and name ca')
		i=i+1
	f.close()
	cmd.hide('labels')
	cmd.set ('dash_gap', 0.05)
	#cmd.do ("orient")
	#cmd.set ('movie_delay', 1500)
 
