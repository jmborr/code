import os,numpy,pickle,pdb
from tempfile import mkstemp

class mopin:
    """obtain Z-matrix from cartesian coordinates"""

    def __init__(self,**kwargs):
        """self.pdb: template pdb
        self.cn: connectivity map
        """
        self.__dict__=kwargs

    def pickleDump(self,pout):
        """save object with picke"""
        pickle.dump(self,pout,True)

    def pickleLoad(self,pin):
        """load object with pickle"""
        tmp=pickle.load(pin)
        self.__dict__.update(tmp.__dict__)

    def loadPDB(self,pdbf):
        self.pdb=open(pdbf).readlines()

    def topmap(self):
        """find the connectivity map self.cn[i] outputs indexes of the
        three atoms to which atom 'i' is connected through a bond, a
        bend angle, and a dihedral angle, respectively. Remember atom
        indexing here begins at ZERO.
        dihedrals are not defined for the first three atoms, thus a '-1'
        is placed for the index of the non-existent atom."""
        tmpf=mkstemp()[-1]; open(tmpf,'w').write(''.join(self.pdb))
        pt=os.popen('babel -ipdb '+tmpf+' -omopin')
        cn=[]; iat=0
        for l in pt.readlines()[3:]:
            l=l.strip()
            cn+=[ int(l[48:53])-1, int(l[53:57])-1, int(l[57:])-1 ]
            iat+=1
        self.cn=numpy.array(cn,dtype=int).reshape(iat,3)
        os.system('/bin/rm -f '+tmpf)
        return self.cn
    
    def zm(self,xyz):
        iat=0; tmpf=mkstemp()[-1]; tmpp=open(tmpf,'w')
        z=numpy.zeros(len(xyz)*3).reshape(len(xyz)*3)
        #create temporary pdb file with coordinates of xyz
        for l in self.pdb:
            if l[0:5]=='ATOM ':
                tmpp.write(l[0:30]+'%8.3f%8.3f%8.3f'%xyz[iat]+'\n')
                iat+=1
            else:
                pbdpt.write(l)
        tmpp.close()
        #find Z-Matrix with babel
        for l in os.popen('babel -ipdb '+tmpf+' -omopin').readlines()[3:]:
            items=l.split();
            d=float(items[1]); bend=float(items[3]); dihed=float(items[5])
            z[iat]=[d,bend,dihed]
        #cleaning of the temporary file
        os.system('/bin/rm '+tmpf) 
        return z #return the z-matrix

    def zmTraj(self,iterator,outp):
        """calculate Z-Matrix for a trajectory of coordinates
        iterator: return numpy.array(nat,3) coordinates"""
        xyz=iterator(); nat=len(xyz); iframe=0; buf=''; maxSize=50E6 #50MB
        tmpf=mkstemp()[-1];
        while xyz.any():
            tmpp=open(tmpf,'w')
            iat=0
            #create temporary pdb file with coordinates of xyz
            for l in self.pdb:
                if l[0:5]=='ATOM ':
                    a,b,c=xyz[iat]
                    tmpp.write(l[0:30]+'%8.3f%8.3f%8.3f'%(a,b,c)+'\n')
                    iat+=1
                else:
                    tmpp.write(l)

            tmpp.close()
            #find Z-Matrix with babel
            pt=os.popen('babel -ipdb '+tmpf+' -omopin')
            for l in pt.readlines()[3:]:
                its=l.split();
                buf+=' %s %s %s'%(its[1],its[3],its[5])
            buf+='\n'
            #check if we need to flush the buffer
            if len(buf)>maxSize:
                outp.write(buf); buf=''
            xyz=iterator(); iframe+=1
            print 'iframe=',iframe
                
        #flushing remaining of buf
        outp.write(buf); 
        #cleaning of the temporary file
        os.system('/bin/rm '+tmpf) 

    def dihedrals(self,iterator,outp):
        """calculate dihedrals for a trajectory of coordinates
        iterator: return numpy.array(nat,3) coordinates"""
        #dihedrals are not defined for first three atoms
        xyz=iterator(); nat=len(xyz); iframe=0; factor=180./numpy.pi
        buf=''; maxSize=50E6 #50MB
        hists=numpy.zeros(nat*360).reshape(nat,360)
        while xyz.any():
            dihedrals=numpy.zeros(nat)
            for i in range(3,nat):
                one=i;two,three,four=self.cn[i]
                v1=xyz[two  ]-xyz[one]
                v2=xyz[three]-xyz[two]
                v3=xyz[four ]-xyz[three]
                v12=numpy.cross(v1,v2); v12n=(v12*v12).sum()
                v23=numpy.cross(v2,v3); v23n=(v23*v23).sum()
                cosPhi=numpy.inner(v12,v23)/numpy.sqrt(v12n*v23n)
                phi=int(factor*numpy.arccos(cosPhi))
                hists[i][phi]+=1
                dihedrals[i]=phi
            buf+=' '.join([str(x) for x in dihedrals])+'\n'
            #check if we need to flush the buffer
            if len(buf)>maxSize:
                outp.write(buf); buf=''
            xyz=iterator(); iframe+=1
            if not iframe%500: print 'iframe=',iframe
        #flushing remaining of buf
        outp.write(buf); 
        #cleaning of the temporary file
        return hists

    def torsions(self,iterator,outp):
        """calculate torsion angles for a trajectory of coordinates
        iterator: return numpy.array(nat,3) coordinates"""
        #dihedrals are not defined for first three atoms
        xyz=iterator(); nat=len(xyz); iframe=0; factor=180./numpy.pi
        buf=''; maxSize=50E6 #50MB
        hists=numpy.zeros(nat*360).reshape(nat,360)
        pdb.set_trace()
        while xyz.any():
            torsions=numpy.zeros(nat)
            for i in range(3,nat):
                one=i;two,three,four=self.cn[i]
                v1=xyz[two  ]-xyz[one]
                v2=xyz[three]-xyz[two]; v2n=(v2*v2).sum()
                v3=xyz[four ]-xyz[three]
                #calculate x-y-z system of reference
                z=-v2/numpy.sqrt(v2n)
                v23=numpy.cross(v2,v3); v23n=(v23*v23).sum()
                y=-v23/numpy.sqrt(v23n)
                x=numpy.cross(y,z)
                #find torsion angle
                w=numpy.dot(v1,z)*z-v1; wn=(w*w).sum()
                cosPhi=numpy.dot(w,x)/numpy.sqrt(wn)
                phi=int(factor*numpy.arccos(cosPhi))
                if numpy.dot(w,y)<0: phi-=360
                hists[i][phi]+=1
                torsions[i]=phi
            buf+=' '.join([str(x) for x in torsions])+'\n'
            #check if we need to flush the buffer
            if len(buf)>maxSize:
                outp.write(buf); buf=''
            xyz=iterator(); iframe+=1
            print 'iframe=',iframe
            if not iframe%1000: print 'iframe=',iframe
            #only for debugging
            if iframe==5000:
                print 'only for debugging'; break
            #end of only for debugging
        #flushing remaining of buf
        outp.write(buf); 
        #cleaning of the temporary file
        return hists




            
