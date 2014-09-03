#!/usr/bin/python

import os,sys
import MySQLdb


#mimic database layout
tables={'author':{'first':'','last':'','username':''},
        'presentation':{'date':'','title':'','filename':''},
#        'dbdir':{'rootd':''}
        }

#connect to the database
try:
    conn=MySQLdb.connect(user='jose',passwd='password',host='localhost',db='research_presentations')
except MySQLdb.Error, e:
     print "Error %d: %s" % (e.args[0], e.args[1])
     sys.exit (1)

cursor = conn.cursor ()

#insert new entry
for tablename in tables.keys():
    table=tables[tablename]
    cols='INSERT INTO '+tablename+' ('
    vals='VALUES ('
    for column in table.keys():
        cols+=column+','
        #print 'insert '+column+':'
        sys.stdout.write( 'insert '+column+':')
        #sys.stdout.flush()
        table[column]=sys.stdin.readline().strip()
        vals+='\''+table[column]+'\','
    print cols[:-1]+') '+vals[:-1]+');'
    cursor.execute(cols[:-1]+') '+vals[:-1]+');')

cursor.close ()
conn.close ()                     
