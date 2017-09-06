import os
from tempfile import mkstemp
from pdb import set_trace as tr

def exec_cpptraj(topfile, script, delscriptfile=True):
    handle, scriptfile = mkstemp(prefix='junk', dir='/tmp')
    handle, stdoutputfile = mkstemp(prefix='junk', dir='/tmp')
    open(scriptfile,'w').write(script)
    if topfile:
        os.system('cpptraj -p {0} -i {1} > {2}'.format(topfile, scriptfile,stdoutputfile))
    else:
        os.system('cpptraj -i {0} > {2}'.format(scriptfile,stdoutputfile)) 
    if delscriptfile:
        os.system('/bin/rm {0}'.format(scriptfile))
    stdoutput = open(stdoutputfile).readlines()
    os.system('/bin/rm {0}'.format(stdoutputfile))
    return scriptfile, stdoutput
