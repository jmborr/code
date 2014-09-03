import numpy

import MDAnalysis
import MDAnalysis.core.qcprot as qcp
from MDAnalysis import SelectionError
from MDAnalysis.core.log import ProgressMeter

from MDAnalysis.analysis.align import _process_selection
from myMDAnalysis.analysis.listCmd import Icommand

class rms_fit_trj(Icommand):
    """RMS-fit trajectory to a reference structure using a selection.

    :Arguments:
      *traj*
         trajectory, :class:`MDAnalysis.Universe` object
      *reference*
         reference coordinates; :class:`MDAnalysis.Universe` object
         (uses the current time step of the object)
      *select*
         1. any valid selection string for
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.selectAtoms` that produces identical
            selections in *mobile* and *reference*; or
         2. a dictionary ``{'mobile':sel1, 'reference':sel2}`` (the
            :func:`fasta2select` function returns such a
            dictionary based on a ClustalW_ or STAMP_ sequence alignment); or
         3. a tuple ``(sel1, sel2)``

         When using 2. or 3. with *sel1* and *sel2* then these selections can also each be
         a list of selection strings (to generate a AtomGroup with defined atom order as
         described under :ref:`ordered-selections-label`).
      *filename*
         file name for the RMS-fitted trajectory or pdb; Not created if not passed
      *rmsdfile*
         file name for writing the RMSD timeseries [``None``]
      *mass_weighted*
         do a mass-weighted RMSD fit
      *tol_mass*
         Reject match if the atomic masses for matched atoms differ by more than
         *tol_mass* [0.1]

    Both reference and trajectory must be :class:`MDAnalysis.Universe`
    instances. If they contain a trajectory then it is used. The
    output file format is the same as the input *traj*.
    """

    #override Icommand.setUp
    def setUp(self,traj, reference, select='all', filename=None, rmsdfile=None,
                 mass_weighted=False, tol_mass=0.1):
        if filename:
            self._writer=traj.trajectory.Writer(filename, remarks='RMS fitted trajectory to reference')
            self._select = _process_selection(select)

        self._ref_atoms = reference.selectAtoms(*select['reference'])
        self._traj_atoms = traj.selectAtoms(*select['mobile'])
        self._natoms = self._traj_atoms.numberOfAtoms()

        if len(self._ref_atoms) != len(self._traj_atoms):
            raise SelectionError("Reference and trajectory atom selections do not contain "+
                                 "the same number of atoms: N_ref=%d, N_traj=%d" % \
                                 (len(self._ref_atoms), len(self._traj_atoms)))

        self._logger.info("RMS-fitting on %d atoms." % len(self._ref_atoms))

        mass_mismatches = (numpy.absolute(self._ref_atoms.masses() - self._traj_atoms.masses()) > tol_mass)
        if numpy.any(mass_mismatches):
            # diagnostic output:
            self._logger.error("Atoms: reference | trajectory")
            for ar,at in zip(self._ref_atoms,self._traj_atoms):
                if ar.name != at.name:
                    self._logger.error("%4s %3d %3s %3s %6.3f  |  %4s %3d %3s %3s %6.3f" % (ar.segid, ar.resid, ar.resname, ar.name, ar.mass,at.segid, at.resid, at.resname, at.name, at.mass,))
            errmsg = "Inconsistent selections, masses differ by more than %f; mis-matching atoms are shown above." % tol_mass
            self._logger.error(errmsg)
            raise SelectionError(errmsg)
        del mass_mismatches

        if mass_weighted:
            # if performing a mass-weighted alignment/rmsd calculation
            self._weight = self._ref_atoms.masses()/self._ref_atoms.masses().mean()
        else:
            self._weight = None

        # reference centre of mass system
        # (compatibility with pre 1.0 numpy: explicitly cast coords to float32)
        self._ref_com = self._ref_atoms.centerOfMass().astype(numpy.float32)
        self._ref_coordinates = self._ref_atoms.coordinates() - self._ref_com

        # allocate the array for selection atom coords
        self._traj_coordinates = self._traj_atoms.coordinates().copy()

        # RMSD timeseries
        self._nframes = len(traj.trajectory)
        self._rmsd = numpy.zeros((self._nframes,))

        # R: rotation matrix that aligns r-r_com, x~-x~com
        #    (x~: selected coordinates, x: all coordinates)
        # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
        self._rot = numpy.zeros(9,dtype=numpy.float64)      # allocate space for calculation
        self._R = numpy.matrix(self._rot.reshape(3,3))

        self._percentage = ProgressMeter(self._nframes, interval=10,format="Fitted frame %(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")

        self._iframe = 0



    #override Icommand.__call__
    def __call__(self,traj, reference, select='all', filename=None, rmsdfile=None,
                 mass_weighted=False, tol_mass=0.1):
        if not self.__isSetUp:
            self.setUp(traj, reference, select=select, filename=filename, rmsdfile=rmsdfile,
                 mass_weighted=mass_weighted, tol_mass=tol_mass)
        if len(traj.frames > 1):
            for ts in traj.frames:
                self.__call__(ts,reference, select='all', filename=None, rmsdfile=None,
                              mass_weighted=False, tol_mass=0.1)
            self.cleanUp()
        else:
            # shift coordinates for rotation fitting
            # selection is updated with the time frame
            x_com = self._traj_atoms.centerOfMass().astype(numpy.float32)
            self._traj_coordinates[:] = self._traj_atoms.coordinates() - x_com

            # Need to transpose coordinates such that the coordinate array is
            # 3xN instead of Nx3. Also qcp requires that the dtype be float64
            # (I think we swapped the position of ref and traj in CalcRMSDRotationalMatrix
            # so that R acts **to the left** and can be broadcasted; we're saving
            # one transpose. [orbeckst])
            k = self._iframe
            self._rmsd[k] = qcp.CalcRMSDRotationalMatrix(self._ref_coordinates.T.astype(numpy.float64),
                                                         self._traj_coordinates.T.astype(numpy.float64),
                                                         self._natoms, self._rot, self._weight)
            self._R[:,:] = self._rot.reshape(3,3)

            # Transform each atom in the trajectory (use inplace ops to avoid copying arrays)
            # (Marginally (~3%) faster than "ts._pos[:] = (ts._pos - x_com) * R + ref_com".)
            ts._pos   -= x_com
            ts._pos[:] = ts._pos * self._R # R acts to the left & is broadcasted N times.
            ts._pos   += self._ref_com

            if filename: 
                self._writer.write(self._traj.atoms) # write whole input trajectory system
                self._percentage.echo(ts.frame)
                self._iframe += 1



    def cleanUp(self):
        self._logger.info("Wrote %d RMS-fitted coordinate frames to file %r",1+self._iframe, self._filename)
        if not self._rmsdfile is None:
            numpy.savetxt(self._rmsdfile,self._rmsd)
            self._logger.info("Wrote RMSD timeseries  to file %r", self._rmsdfile)
