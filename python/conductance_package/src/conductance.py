#!/usr/bin/python
"""
Type:
python conductance.py -a prmtop -b mdcrd -c mdvel -d pdbfile -e indexes.dat -f outDirectory
"""
import os,sys,numpy,pdb

isthere=os.path.exists
joink=os.path.join

def Ktraj(topC,velf,pdbf,outd):
    """create trajectory of kinetic energy"""
    from amber10 import crd
    from smallThings import insertBfact
    from readingWritingFiles import write_from_numpy

    nat=topC.ITITL['NATOM']
    Kcrd=joink(outd,'K.crd');
    if isthere(Kcrd): os.system('/bin/rm '+Kcrd)
    Kpt=open(Kcrd,'a')
    crdobj=crd(velf,nat)
    Kav,Kdev=crdobj.kinetic_energy(topC.AMASS,outptf=Kpt)
    os.system('gzip -f '+Kcrd); Kcrd+='.gz'  #gzip trajectory
    print '#save average Kinetic energy as .pdb'
    buf=insertBfact(pdbf,Kav)
    Kpdb=joink(outd,'Kav.pdb')
    open(Kpdb,'w').write(buf)
    print '#save average Kinetic energy as .dat'
    buf='\n'.join([str(x) for x in Kav])
    Kdat=joink(outd,'Kav.dat')
    open(joink(outd,Kdat),'w').write(buf)
    return Kcrd,Kpdb,Kdat

def gnmodes(topC,crdf,outd):
    """calculate quasiharmonic modes in the isotropic gaussian
    network approximation"""
    from gn3D import gn3D
    from amber10 import crd

    nat=topC.ITITL['NATOM']
    gn=gn3D(nat)
    gnfile=joink(outd,'gn.dump')
    if isthere(gnfile):
        gn.pickleLoad(open(gnfile))
    else:
        crdobj=crd(crdf,nat)#coordinate trajectory object
        gn.loadM(topC.AMASS) #load masses
        gn.findVarCovar(crdobj.next)#mass weighted varcov matr
        gn.sd=outd
        gn.pickleDump(open(gnfile,'w'))
    return gn

def find_neighborlist(pdbf):
    """find neighbor atom contact list"""
    from smallThings import neighborList
    neig,dist,xyz=neighborList(open(pdbf),co=6.0)
    return {'all':{'co':6.0,'map':neig,'distances':dist}}

def get_VdW(pdbf):
    """retrieve VdW radii"""
    radii={'H':1.20, 'C':1.70, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.80,
          'S':1.80, 'Cl':1.75, 'Cu':1.4, 'Na':2.27}
    VdW=[]
    for l in open(pdbf).readlines():
        if l[0:5]=='ATOM ':
            name=l[12:16]
            if    'Na' in name: VdW.append( radii['Na'] )
            elif  'Cl' in name: VdW.append( radii['Cl'] )
            else: VdW.append( radii[ l[13] ] )
    return numpy.array(VdW)

def findResIndexes(pdbf,indexf):
    """for each residue, find the list of atom indexes

    NOTE: first atom index is 0, NOT 1

    NOTE: we create and extra residue with the cold-spot atoms, thus
    removing these atoms from whatever residues in the PDB belong to."""

    nres=0; iline=0
    #search for cold-atoms
    catalytic=[]
    for line in open(indexf).readlines():
        if line.strip()=='1': catalytic.append(True)
        else: catalytic.append(False)
        
    ptin=open(pdbf); l=ptin.readline()
    resindexes=[]; currires=None
    group=[]; catgroup=[]
    iat=0 #atom index
    while l:
        if l[0:5]=='ATOM ':
            ires=l[22:26];      #residue number
            resname=l[17:20]    #residue name
            if not currires: currires=ires #initialize residue index
            if ires!=currires:  #new residue
                nres+=1
                resindexes.append(group)
                currires=ires
                group=[]        #empty bucket
            if catalytic[iat]: catgroup.append(iat)
            else: group.append(iat)
            iat+=1
        l=ptin.readline()
    resindexes.append(group); nres+=1 #enter last group
    resindexes.append(catgroup); nres+=1 #catalytic atoms as last group
    return {'indexes':resindexes,'nres':nres}
    
def meanFieldFlux(topf,crdf,velf,pdbf,indexf,outd):
    """calculate mean flux according to
            J_{i,j}=-w_{i,j}*k_B*(T_i-T_j)
            w_{i,j}^2=k_B*(T_i+T_j)/(2*m_i*m_j*<r_i*r_j>)
            r_i: vector displacement of atom 'i'
    """
    from readingWritingFiles import write_from_numpy
    from smallThings import insertBfact
    from math import log10
    from amber10 import top
    import fluxGraph

    topC=top(topf)        #instantiate topology object
    nat=topC.ITITL['NATOM']
    
    print '#load av. kinetic energies, transform to Temperatures'
    Kcrd,Kpdb,Kdat=Ktraj(topC,velf,pdbf,outd)
    Kdat=joink(outd,'Kav.dat')
    pt=open(Kdat)
    Ts=numpy.array([float(x) for x in pt.readlines()])


    print '#load variance covariance matrix, find displacements'
    gn=gnmodes(topC,crdf,outd)
    vc=gn.varcov
    msd=numpy.sqrt(vc.diagonal()/gn.M) #remember is mass weighted
    std=vc.ravel().compress(vc.ravel()>0.0).std() #std of posit cov
    v_cutoff=std*0.1

    print '#load native contact map'
    contmap=find_neighborlist(pdbf)
    cmap=contmap['all']['map']
    ds=contmap['all']['distances']
    
    print '#filter map by cutoff: VdW_i+msd[i]+VdW_j+msd[j]'
    VdW=get_VdW(pdbf)
    cnts=[]
    for i in range(nat):
        ds_i=ds[i]; map_i=cmap[i]; cnt=[]; a=VdW[i]+0.5*msd[i]
        for j in range(len(ds_i)):
            if ds_i[j]<a+VdW[j]+0.5*msd[j]:
                cnt.append( map_i[j] )
        cnts.append(cnt)

    print '#load atom symbols'
    atomSymbol=[]
    for x in open(pdbf).readlines():
        if x[0:5]=='ATOM ': atomSymbol.append( x[13] )
    
    print '#calculate interatomic flux'
    f=numpy.zeros(nat*nat).reshape(nat,nat)
    for iat in range(nat):
        cnt_i=cnts[iat]; T_i=Ts[iat]
        #Below line when SHAKE performed on Hydrogen Bonds
        #if atomSymbol[iat]=='H': T_i*=1.409#SHAKE :( Debugging
        for j in range(len(cnt_i)):
            jat=cnt_i[j]; v=vc[iat][jat]
            if v<=0: continue #negative force constant
            T_j=Ts[jat]
            #Below line when SHAKE performed on Hydrogen Bonds
            #if atomSymbol[jat]=='H': T_j*=1.409#Debugging
            avT=(T_i+T_j)/2.0; w=numpy.sqrt(avT/v)
            #f[i][j]: flow arriving to 'i' from 'j'
            f[iat][jat]=w*(T_j-T_i); f[jat][iat]=-f[iat][jat]
    outf=joink(outd,'meanFieldFlux.dat'); pt=open(outf,'w')
    comment='#Mean Field Flux'
    write_from_numpy(pt,f,format=' %9.5f',ncols=8,comment=comment)
    pt.close()
    
    print '#calculate interatomic absorvances'
    avs=f.sum(axis=1)
    outf=joink(outd,'meanFieldAbsov.dat'); pt=open(outf,'w')
    comment='#Mean Field Absorvances'
    write_from_numpy(pt,avs,format=' %f',ncols=1,comment=comment)
    pt.close()
    outf=joink(outd,'meanFieldAbsov.pdb'); pt=open(outf,'w')
    buf=insertBfact(pdbf,(-1)*avs)
    open(outf,'w').write(buf)
    print '#calculate interatomic Conductance'
    trf=numpy.abs(f).sum(axis=1)
    outf=joink(outd,'meanFieldCond.dat'); pt=open(outf,'w')
    comment='#Mean Field Cond'
    write_from_numpy(pt,trf,format=' %f',ncols=1,comment=comment)
    pt.close()
    outf=joink(outd,'meanFieldCond.pdb'); pt=open(outf,'w')
    buf=insertBfact(pdbf,avs)
    open(outf,'w').write(buf)
    
    print '#calculate interresidue flux'
    resinfo=findResIndexes(pdbf,indexf)
    nres=resinfo['nres']
    rf=numpy.zeros(nres*nres).reshape(nres,nres) #flux between residues
    for ires in range(nres-1):
        idx=resinfo['indexes'][ires] #list of atoms in ires
        for jres in range(ires+1,nres):
            jdx=resinfo['indexes'][jres] #list of atoms in jres
            a=0; #add atomic fluxes
            for iat in idx:
                for jat in jdx: a+=f[iat][jat]
            #rf[ires][jres] is the flow arriving to ires from jres
            rf[ires][jres]=a; rf[jres][ires]=-a
    outf=joink(outd,'meanFieldFluxRes.dat'); pt=open(outf,'w')
    comment='#Mean Field Flux'
    write_from_numpy(pt,rf,format=' %9.5f',ncols=8,comment=comment)
    pt.close()
    outf=joink(outd,'meanFieldFluxRes.pml')
    
    fluxGraph.fluxGraph(rf,pdbf,indexf,outf)
    print 'copy arrow.py from src/ to the output directory, then run "pymol meanFieldFluxRes.pml'

    print '#calculate interresidue absorvances'
    ravs=rf.sum(axis=1) #net flux through a residue
    outf=joink(outd,'meanFieldAbsovRes.dat'); pt=open(outf,'w')
    comment='#Mean Field Absorvances'
    write_from_numpy(pt,ravs,format=' %f',ncols=1,comment=comment)
    pt.close()
    avs=numpy.zeros(nat)
    onezero=numpy.zeros(nat)
    for ires in range(nres):
        for iat in resinfo['indexes'][ires]:
            avs[iat]=log10(abs(ravs[ires])) #note the LOG!!!
    outf=joink(outd,'meanFieldAbsovRes.pdb'); pt=open(outf,'w')
    buf=insertBfact(pdbf,(-1)*avs,)
    open(outf,'w').write(buf)

    print '#calculate interresidue Cond (rf:interresidue flux)'
    rtrf=numpy.abs(rf).sum(axis=1)
    outf=joink(outd,'meanFieldCondRes.dat');
    pt=open(outf,'w')
    comment='#Mean Field Cond'
    write_from_numpy(pt,rtrf,format=' %f',ncols=1,comment=comment)
    pt.close()
    avs=numpy.zeros(nat)
    for ires in range(nres):
        for iat in resinfo['indexes'][ires]:
            avs[iat]=log10(rtrf[ires])
    outf=joink(outd,'meanFieldCondRes.pdb');
    pt=open(outf,'w')
    buf=insertBfact(pdbf,avs,)
    open(outf,'w').write(buf)            
    return False #True

if __name__=='__main__':
    from inputArgs import inpHand    
    ih=inpHand('Usage: conductance.py [Required args] [Optional args]',
               ' -a _RA_topf topology file without water and ions',
               ' -b _RA_crdf coordinates file',
               ' -c _RA_velf velocities file',
               ' -d _RA_pdbf unsolvated pdb file',
               ' -e _RA_indexf file with one-number code to find cold-spot atoms',
               ' -f _A_outd directory to output results (default: current dir)',
               )
    ih.parse(locals(),sys.argv)

    #Resolve optional arguments
    if not outd: outd=os.getcwd()

    #Check existence of files and directories
    for f in (topf,crdf,velf,pdbf,indexf):
        assert isthere(f), 'file '+f+' not found!'
    if not isthere(outd): os.system('/bin/mkdir -p '+outd)

    meanFieldFlux(topf,crdf,velf,pdbf,indexf,outd)
    
    print 'Are there any junk files left over?'
    sys.exit(0)
