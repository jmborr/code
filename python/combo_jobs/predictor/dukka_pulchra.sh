#!/bin/bash
/bin/mv $2/$1 . #transfer the combo file to current directory
dukka.x $1
pulchra best.pdb
/bin/mv best.pdb $2/${1}.rebuilt #transfer the rebuilt combo file
exit 0
