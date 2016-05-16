#!/usr/bin/python
import hashlib
import os,sys
 
def md5sum(filename, blocksize=65536):
    """Calculate md5 hash sum of a file provided """
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()
def save_sum(filename,md_sum):
    """Save hash sum into file with appropriate filename"""
    md_fname = filename+'.md5'
    f = open(md_fname,'w')
    f.write(md_sum)
    f.close()
 
if __name__ == '__main__':
 
    if len(sys.argv)<2 or not os.path.isfile(sys.argv[1]):
        print "Create md5 has for a data file to be used in testing"
        print "See http://www.mantidproject.org/Data_Files_in_Mantid"
        print "Usage: hash_file.py file_name"
        exit(1)
 
    filename = sys.argv[1]
 
    path,fname = os.path.split(filename)
    hash_sum = md5sum(filename)
    print "MD SUM FOR FILE: {0} is {1}".format(fname,hash_sum)
 
   # save hash sum in file with original file name and extension  .md5
    save_sum(os.path.join(path,fname),hash_sum)
 
   # Rename hashed file into hash sum name. 
    hash_file = os.path.join(path,hash_sum)
    if os.path.isfile(hash_file):
        print "file: {0} already exist".format(hash_sum)
    else:
        os.rename(filename,hash_file)
