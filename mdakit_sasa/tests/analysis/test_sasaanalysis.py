import pytest
import numpy as np

from numpy.testing import assert_allclose

from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
from mdakit_sasa.tests.utils import make_Universe
from MDAnalysis.core.topologyattrs import Atomnames, Resnames, Resids, Resnums, Segids, Atomtypes

class TestSASAAnalysis:

    # fixtures are helpful functions that set up a test
    # See more at https://docs.pytest.org/en/stable/how-to/fixtures.html
    @pytest.fixture
    def universe(self):
        u = make_Universe(
            n_frames=3,
        )
        
        for ts in u.trajectory:
            ts.positions[:ts.frame] *= -1
            
        u.add_TopologyAttr(Atomnames(["H"] * len(u.atoms)))
        u.add_TopologyAttr(Resnames(["GLY"] * len(u.residues)))
        u.add_TopologyAttr(Resids(list(range(0, len(u.residues)))))
        u.add_TopologyAttr(Resnums(list(range(0, len(u.residues)))))
        u.add_TopologyAttr(Segids(["A"] * len(u.segments)))
        u.add_TopologyAttr(Atomtypes(["O"]* len(u.atoms)))
        return u

    @pytest.fixture
    def analysis(self, universe):
        return SASAAnalysis(universe)

    @pytest.mark.parametrize(
        "select, n_atoms",  # argument names
        [  # argument values in a tuple, in order
            ("all", 125),
            ("index 0:9", 10),
            ("segindex 3:4", 50),
        ]
    )
    def test_atom_selection(self, universe, select, n_atoms):
        # `universe` here is the fixture defined above
        analysis = SASAAnalysis(
            universe, select=select)
        assert analysis.atomgroup.n_atoms == n_atoms

    def test_total_sasa_calculation(self, analysis):
        analysis.run(stop=3)
        assert analysis.n_frames == 3
        
    def test_total_sasa_calculation_results(self, analysis):
        analysis.run(stop=3)
        assert analysis.n_frames == 3
        assert analysis.results['total_area'].dtype ==  np.dtype('float64')
        assert np.all(analysis.results['total_area'] >= 0)

    def test_residue_sasa_calculation_results(self, analysis):
        analysis.run(stop=3)
        assert analysis.n_frames == 3
        assert analysis.results['residue_area'].dtype ==  np.dtype('float64')
        assert np.all(analysis.results['residue_area'] >= 0)
        assert analysis.results['residue_area'].shape == (3,25)