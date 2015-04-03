import os
from tempfile import mkstemp

def exec_cpptraj(topfile, script, delscriptfile=True):
    handle, scriptfile = mkstemp(prefix='junk', dir='/tmp')
    open(scriptfile,'w').write(script)
    os.system('module load amber/amber14; cpptraj -p {0} -i {1}'.format(topfile, scriptfile))
    if delscriptfile:
        os.system('/bin/rm {0}'.format(scriptfile))
    return scriptfile
