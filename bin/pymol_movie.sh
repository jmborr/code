#!/bin/bash
#In pymol, type:
# mset 1 x64
# util.mroll(1,64,1)
# mpng mov
#this will output 64 frames mov0000.png to mov064.png
#Now in a terminal:
wd=$1 #working directory
cd $wd
for f in `ls -1 mov*.png`;do echo $f && convert $f mov${f:3:4}.gif;done >/dev/null
gifsicle --delay=20 --colors 256 --loop mov????.gif > anim.gif
/bin/rm -f mov????.png mov????.gif
echo "produced anim.gif"