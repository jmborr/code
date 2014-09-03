#!/usr/bin/python

import os,sys

def red2blue(val):
    return rainbow(val*0.66666666666)

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
    
if __name__=='__main__':
    for i in range(0,21):
        val=float(i)/20
        print 'color',red2blue(val)
