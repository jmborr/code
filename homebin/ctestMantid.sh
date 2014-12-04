#!/bin/bash

cd /home/jbq/code/mantidproject/mantid/Code/debug

make -j8 AllTests

echo -e "\n\nFinished building Mantid. Testing now. Logfile is /tmp/ctestMantid.log\n"
ctest -j8 &> '/tmp/ctestMantid.log'

popwindow.py "finished testing of Mantid.\nLog file /tmp/ctestMantid.log "
