'''
Created on Dec 20, 2012

@author: jmborr
'''

import xml.etree.ElementTree as xml

root=xml.Element('root')

# Append a child element to root
child=xml.Element('child')
root.append(child)

# Set an attribute
child.attrib['name']='Charlie'

# Write to file
file=open('/tmp/junk.xml', 'w')
xml.ElementTree(root).write(file) 
file.close()