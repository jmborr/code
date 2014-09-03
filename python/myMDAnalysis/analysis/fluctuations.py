import numpy

from copy import deepcopy
import MDAnalysis

from myMDAnalysis.analysis.listCmd import Icommand

class avg_conf(Icommand):
    """Average conformation

    :Arguments:
        *traj*
            trajectory, :class:`MDAnalysis.Universe` object. If it contains
            a trajectory then it is used.
        *average*
            average conformation, :class:`MDAnalysis.Universe` object. Stores
            the output average conformation coordinates
    """
    #override Icommand.setUp
    def setUp(self,*kargs,**kwargs):
        self._logger.info("Calculating average conformation")
        self._counter = 0

    #override Icommand.__call__
    def __call__(self,traj,average):
        if len(traj.frames) > 1:
            self.setup()
            for ts in traj.frames:
                self.__call__(ts,average)
            self.cleanup()
        else:
            average.coordinates += ts.coordinates
            self._counter += 1

    #override Icommand.cleanUp
    def cleanUp(self,average):
        average.coordinates /= self._counter
        self._logger.info("Finished calculating average over "+`self._counter`+" conformations")
        
def avg_conf(traj):
    """Average conformation

    :Arguments:
      *traj*
         trajectory, :class:`MDAnalysis.Universe` object. If it contains
         a trajectory then it is used.

    :Returns:
      *avg_conf*
         average conformation, :class:`MDAnalysis.Universe` object. A copy of
         argument *traj* with atom coordinates being the average coordinates
    """
    average = deepcopy( traj )
    avg_coord = None

    logger.info("Calculating average conformation")

    for k,ts in enumerate( traj.trajectory ):
        if not avg_coord: avg_coord = ts.coordinates().copy()
        avg_coord += ts.coordinates()
    average.atoms.coordinates = avg_coord / k

    return average