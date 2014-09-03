import pymol.cmd as pcmd

def bondthem(sel,co=4.0):
    """
    connect atoms in the selection closer than cut-off co with a bond
    """
    co2 = co*co
    m = cmd.get_model("%s"%sel)
    nat = len( m.atom )
    for i in range(nat-1):
        iat = m.atom[i]
        ist = '/%s//%s/%s/%s'%(sel,iat.chain,iat.resi,iat.name)
        for j in range( 1+i, nat ):
            jat = m.atom[j]
            d2 = 0.0
            for k in range(3):
                d = iat.coord[k] - jat.coord[k]
                d2 += d*d
            if d2 < co2:
                jst = '/%s//%s/%s/%s'%(sel,jat.chain,jat.resi,jat.name)
                pcmd.bond( ist, jst )

cmd.extend("bondthem",bondthem)

def bondCAatoms(file,cols=(1,2)):
    """place a bond between CA atoms given by cols containing the
    corresponding residue numbers
    """
    icol = cols[0] - 1
    jcol = cols[1] - 1
    pin=open(file)
    for l in pin.readlines():
        print l.strip()
        x = l.split()
        i = int( x[icol] )
        j = int( x[jcol] )
        iat = "resid %d and name CA"%i
        jat = "resid %d and name CA"%j
        pcmd.bond(iat,jat)
        pcmd.show( "lines", "%s or %s"%(iat,jat) )

cmd.extend("bondCAatoms",bondCAatoms)
