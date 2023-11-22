"""
SASAAnalysis --- :mod:`mdakit_sasa.analysis.SASAAnalysis`
===========================================================

This module contains the :class:`SASAAnalysis` class, this class follow
the standarised :class:`AnalysisBase` that will be familiar for MDAnalisys users.

This module is the entry point for this MDA Kit.

"""
from typing import Union, TYPE_CHECKING

from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
import freesasa
import os
import logging

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe, AtomGroup

logger = logging.getLogger(__name__)

class SASAAnalysis(AnalysisBase):
    """SASAAnalysis class.

    This class is used to compute the solvant accessible area of a trajectory.
.

    Parameters
    ----------
    universe_or_atomgroup: :class:`~MDAnalysis.core.universe.Universe` or :class:`~MDAnalysis.core.groups.AtomGroup`
        Universe or group of atoms to apply this analysis to.
        If a trajectory is associated with the atoms,
        then the computation iterates over the trajectory.
    select: str
        Selection string for atoms to extract from the input Universe or
        AtomGroup

    Attributes
    ----------
    universe: :class:`~MDAnalysis.core.universe.Universe`
        The universe to which this analysis is applied
    atomgroup: :class:`~MDAnalysis.core.groups.AtomGroup`
        The atoms to which this analysis is applied
    results: :class:`~MDAnalysis.analysis.base.Results`
        results of calculation are stored here, after calling
        :meth:`SASAAnalysis.run`
    start: Optional[int]
        The first frame of the trajectory used to compute the analysis
    stop: Optional[int]
        The frame to stop at for the analysis
    step: Optional[int]
        Number of frames to skip between each analyzed frame
    n_frames: int
        Number of frames analysed in the trajectory
    times: numpy.ndarray
        array of Timestep times. Only exists after calling
        :meth:`SASAAnalysis.run`
    frames: numpy.ndarray
        array of Timestep frame indices. Only exists after calling
        :meth:`SASAAnalysis.run`
    """

    def __init__(
        self,
        universe_or_atomgroup: Union["Universe", "AtomGroup"],
        select: str = "all",
        **kwargs
    ):
        super().__init__(universe_or_atomgroup.trajectory, **kwargs)
        self.universe = universe_or_atomgroup.universe
        self.atomgroup = universe_or_atomgroup.select_atoms(select)

    def _prepare(self):
        self.results.total_area = np.zeros(
            self.n_frames,
            dtype=float,
        )
        self.results.residue_area = np.zeros(
            (self.n_frames, len(self.universe.residues.resids)),
            dtype=float,
        )


    def _single_frame(self):
        """Calculate data from a single frame of trajectory"""
        
        structure = freesasa.Structure()  
        # FreeSasa structure accepts PDBS if not available requires to reconstruct the structure using `addAtom`
        for a in self.atomgroup:
            x,y,z = a.position
            structure.addAtom(a.type.rjust(2), a.resname, a.resnum.item(), a.segid, x, y, z)
        
        # Define 1 cpu for windows avoid freesasa code to calculate it.
        parametes =  freesasa.Parameters()
        if self._is_windows():
            parametes.setNThreads(1)

        result = freesasa.calc(structure, parametes)

        residue_areas = [result.residueAreas()[s][r] for s in sorted(list(result.residueAreas().keys())) for r in sorted(list(result.residueAreas()[s].keys()))]
        
        self.results.total_area[self._frame_index] = result.totalArea()
        
        # Defend agains residue counts mismatch
        if  len(self.universe.residues.resids)!= len(residue_areas):
            logger.error(f'Residude count do not match the expectation, residue SASA not in results { len(self.universe.residues.resids)} != {len(residue_areas)}')
        else:
            self.results.residue_area[self._frame_index] = [r.total for r in residue_areas]
        
    def _conclude(self):
        self.results.mean_total_area= self.results.total_area.mean()

    def _is_windows(self):
        return os.name == 'nt'