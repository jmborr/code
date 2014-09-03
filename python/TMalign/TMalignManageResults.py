#!/usr/bin/python

import sys,os,re,math
from utilities.small_utilities import chomp,Bye
from seq.alignYangManageResults import alignment,gappedSeq

rs={}
rs['Chain1']=re.compile('Chain 1:(\S+)\s')
rs['Chain2']=re.compile('Chain 2:(\S+)\s')
rs['size1']=re.compile('Chain 1.*Size=\s*(\d+)')
rs['size2']=re.compile('Chain 2.*Size=\s*(\d+)')
rs['Aligned length']=re.compile('Aligned length=\s*(\d+)')
rs['tm']=re.compile('TM-score=(\d\.\d+)')
rs['seq1']=re.compile('Angstrom\)\s*\n([^\n]+)\n')
rs['identities']=re.compile('Angstrom\)\s*\n[^\n]+\n([^\n]+)')
rs['seq2']=re.compile('Angstrom\)\s*\n[^\n]+\n[^\n]+\n([^\n]+)\n')
rs['gaps']=re.compile('(-+)')
rs['init']=re.compile('[A-Z]')
rs['rmsd']=re.compile('RMSD=\s*(\d+\.\d+)')

def TMalign(pdb1,pdb2):
	return float(os.popen('TMalign '+pdb1+' '+pdb2+'|grep "TM-score="|cut -d\'=\' -f 4|cut -d\',\' -f 1').readline())
	
def get_d0(N):
	return 1.24*(N-15)**(1.0/3.0)-1.8
	
def reverse(s): 
	"""Return the sequence string in reverse order.""" 
	letters = list(s) 
	letters.reverse() 
	return ''.join(letters)

def loadAlg(inpf):
	alg=[]
	if isinstance(inpf,file):	
		#search for the first line of the entry
		while(True):
			line=inpf.readline()
			if not line: return [] #reached end of file
			if line.find('Chain 1:')>=0:
				alg.append(line) ; break
	        #keep on reading until we find other entry or reach end of file
		while(True):
			line=inpf.readline()
			if not line: break #end-of-file reached
			if line.find('Angstrom')>=0:
				alg.append(line)
				break
			alg.append(line)
		line=inpf.readline() ; alg.append(line) #first gapped sequence
		line=inpf.readline() ; alg.append(line) #"identities" lines
		line=inpf.readline() ; alg.append(line) #second gapped sequence
	return alg

class TMalignOut:
	
	def __init__(self,inpf):		
		'''inpf can be the following:
		(1) handle to input file containing results of some Yang\'s TMalign program
		(2) string containing name of previous file
		'''
		self.readError=''
		self.inpf=None
		self.all=None
		self.alg=None
		
		#type checking
		if isinstance(inpf,str): #name of the alignment file
			self.inpf=inpf
			lines=open(inpf,'r').readlines()
			self.all=''.join(lines) #join all lines of the input file
		elif isinstance(inpf,file): #file handle to the alignment file
			lines=loadAlg(inpf)
			if lines==[]:
				self.readError='no alignment found'
				return
			self.all=''.join(lines)

		self.Chain={}
		self.size={}
		self.seq={}
		self.ngaps={}
		
		self.Chain[1]=rs['Chain1'].search(self.all).group(1)
		self.Chain[2]=rs['Chain2'].search(self.all).group(1)
		self.size[1]=int(rs['size1'].search(self.all).group(1))
		self.size[2]=int(rs['size2'].search(self.all).group(1))
		self.AlignedLength=int(rs['Aligned length'].search(self.all).group(1))
		self.tm=float(rs['tm'].search(self.all).group(1))
		self.seq[1]=rs['seq1'].search(self.all).group(1)
		self.identityLine=rs['identities'].search(self.all).group(1)
		self.nAligned=len( re.findall(':',self.identityLine) ) #number of aligned residues
		self.seq[2]=rs['seq2'].search(self.all).group(1)
		self.rmsd=float(rs['rmsd'].search(self.all).group(1))
		self.isAligned={1:[],2:[]} #tell if residue is aligned
		self.alg=alignment(gappedSeq(self.seq[1]),gappedSeq(self.seq[2]))
		#print self.Chain[1],self.Chain[2],self.size[1],self.size[2],self.tm
		#print self.seq[1]
		#print self.seq[2]

	def init_ngaps(self):
		for i in [1,2]:
			insertions=rs['gaps'].findall(self.seq[i])
			#print insertions
			begin=0 ; end=len(insertions)
			#avoid beginning gaps
			if self.seq[i][0]=='-':	begin=1
			#avoid final gaps
			if self.seq[i][len(self.seq[i])-1]=='-': end=-1
			n=0
			for insertion in insertions[begin:end]: n+=len(insertion)
			#print 'n=',n
			self.ngaps[i]=n

	def gapDensity(self):
		if not self.ngaps: self.init_ngaps()
		#find beginning of the alignment
		init=rs['init'].search(self.seq[1]).start()
		n=rs['init'].search(self.seq[2]).start()
		if n>init: init=n
		#print 'init=',init		
		#find ending of the alignment
		seq=reverse(self.seq[1])
		end=len(seq)-rs['init'].search(seq).start()-1
		seq=reverse(self.seq[2])
		n=len(seq)-rs['init'].search(seq).start()-1
		if n<end: end=n
		#print 'end=',end
		#find gap density
		alignL=1+end-init #length of the alignment
		numgaps=self.ngaps[1]+self.ngaps[2] #number of gaps within the alignment
		#print alignL,numgaps
		return (1.*numgaps)/alignL

	def initIsAlignedArray(self):
		for i in range( len(self.identityLine) ):
			isId=0
			if self.identityLine[i]==':': isId=1 #corresponding residues are aligned
			for chainID in (1,2):
				if self.seq[chainID][i]!='-': #position corresponds to a residue
					self.isAligned[chainID].append(isId)
				
