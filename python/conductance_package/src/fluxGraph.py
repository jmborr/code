#!/usr/bin/python

import os,sys,numpy,pdb

def rainbow(val):
    if val<0 or val>1:
        sys.stderr.write('ERROR: val should be in the [0,1] range')
        return None
    
    interval=int(val*6)
    residue=int(255*(val*6.0-float(interval)))

    if   interval==0: return [255        ,residue    ,0          ]
    elif interval==1: return [255-residue,255        ,0          ]
    elif interval==2: return [0          ,255        ,residue    ]
    elif interval==3: return [0          ,255-residue,255        ]
    elif interval==4: return [residue    ,0          ,255        ]
    elif interval==5: return [255        ,0          ,255-residue]
    elif interval==6: return [255        ,0          ,0          ]

def red2blue(val):
    return rainbow(val*0.66666666666)

def fluxGraph(rf,pdbf,indexf,outf,co=0.00,coM=1E+06):
    """pdb containing one line per residue, with a representative atom"""
    valid_AAA=('ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLN', 'GLU',
               'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
               'PRO', 'SER', 'THR', 'TRP', 'TYR', 'UNK', 'VAL')
    nres=len(rf) #number of residues
    #search for cold-atoms
    ix=0
    catalytic=[]
    for line in open(indexf).readlines():
        if line.strip()=='1': catalytic.append(True)
        else: catalytic.append(False)
        ix+=1
    
    #group lines of the PDB in terms of residues.
    ix=0; currires=None; name=None
    groups=[]; currgroup=[]; catgroup=[]; resnames=[]; iress=[]; mres=0
    ptin=open(pdbf); l=ptin.readline(); natom=0
    while l:
        if l[0:5]=='ATOM ':
            if catalytic[ix]:
                catgroup.append(l)
            else:    
                ires=l[22:26];      #residue number
                resname=l[17:20]    #residue name
                if not currires:
                    currires=ires #very first residue
                    name=resname
                if ires!=currires:  #new residue
                    groups.append(currgroup)
                    currgroup=[]
                    resnames.append(name)
                    name=resname
                    iress.append( int(currires) )
                    currires=ires
                    mres+=1
                else: currgroup.append(l)
            ix+=1
        if mres==nres-1: break #we finished
        l=ptin.readline()
    groups.append(currgroup)
    resnames.append(name)
    iress.append( int(currires) )
    groups.append(catgroup);
    resnames.append('CLD')
    iress.append( int(iress[-1]+1) )
    #gather representative coordinates for each group
    rxyz=[]
    for i in range(nres):
        isValid=False
        group=groups[i]; name=resnames[i]
        avxyz=numpy.zeros(3)
        for l in group:
            atomname=l[12:16]
            xyz=[float(x) for x in [l[30:38],l[38:46],l[46:54]] ]
            avxyz+=numpy.array(xyz)
            if name in valid_AAA and atomname==' CA ': #standard residue
                isValid=True
                rxyz.append( xyz )
                break
        if not isValid: rxyz.append( avxyz/len(group) ) #geometric center
    #create pymol pml file
    maxval=numpy.max(rf)
    print 'max_co=',str(maxval)
    buf='load '+pdbf+'\nrun arrow.py\n'
    buf+='hide lines\n'
    buf+='show cartoon\n'
    buf+='select CA, name CA\n'
    buf+='show spheres, CA\n'
    buf+='color white, CA\n'
    buf+='set sphere_scale=0.2\n'
    #pdb.set_trace()
    for i in range(nres):
        ri=rxyz[i]
        xi='[%8.3f,%8.3f,%8.3f]'%(ri[0],ri[1],ri[2])
        for j in range(nres):
            #if j==nres-1 and i==5: pdb.set_trace()
            #eliminate flux between residues further away than 7.5A
            dd = 0.0
            for k in range(3):
                dd += (rxyz[i][k]-rxyz[j][k])**2
            #rf[i][j] is flux arriving to i from j
            if rf[i][j] > co and rf[i][j] < coM and abs(i-j)>2\
                   and dd < 7.5**2:
                rj=rxyz[j]
                #retain only flux going towards the cold atoms
                p = 0.; a = 0.; b = 0.
                for k in range(3):
                    rc = rxyz[-1] #coordinates of cold atoms
                    p += (ri[k]-rj[k])*(rc[k]-rj[k])
                    a += (ri[k]-rj[k])*(ri[k]-rj[k])
                    b += (rc[k]-rj[k])*(rc[k]-rj[k])
                p = p / numpy.sqrt(a*b)
                if p < 0.5: continue #flux not towards the cold atoms
                xj='[%8.3f,%8.3f,%8.3f]'%(rj[0],rj[1],rj[2])
                name=resnames[i]+`iress[i]`+'-'+resnames[j]+`iress[j]`
                color=numpy.array(red2blue( (rf[i][j]-co)/(maxval-co) ))
                color=(1.0*color)/255
                color='[%4.2f,%4.2f,%4.2f]'%(color[0],color[1],color[2])
                buf+='cgo_arrow(%s,%s,r=0.15,hr=0.5,hl=2.0,name="%s",color=%s)\n'%(xj,xi,name,color)
    buf+='center\nzoom\n'
    open(outf,'w').write(buf)
      
if __name__=='__main__':
    from inputArgs import inpHand    
    ih=inpHand('Usage: fluxGraph.py [Required args] [Optional args]',
               ' -a _RA_topf topology file with not water and ions',
               ' -b _RA_fluxf meanFieldFluxRes.dat file',
               ' -c _RA_pdbf pdb file',
               ' -d _RA_indexf file with one-number code to find cold-spot atoms',
               ' -e _A_outf output file (default: meanFieldFluxRes.pml)',
               ' -f __co cut-off flux (default=0.00)',
               ' -g __coM maximum cut-off flux (default=1E6)'
               )
    ih.parse(locals(),sys.argv)
    
    #Resolve optional arguments
    if not outf: outf='meanFieldFluxRes.pml'
    if co: co=float(co)
    else:  co=0.00
    if coM: coM=float(coM)
    else:  coM=1E+06

    from amber10 import top
    topC=top(topf)        #instantiate topology object
    nres=topC.ITITL['NRES']+1 #extra 1 for cold atoms
    
    #load meanFieldFlux.dat as matrix
    from readingWritingFiles import read_to_numpy
    #pdb.set_trace()
    rf=read_to_numpy( open(fluxf),nres*nres,shape=(nres,nres))

    fluxGraph(rf,pdbf,indexf,outf,co=co,coM=coM)
