import pytest
import numpy as np

from mdakit_sasa.analysis.parallel_sasa import ParallelSASAAnalysis
from mdakit_sasa.tests.utils import make_Universe
from MDAnalysis.core.topologyattrs import Atomnames, Resnames, Resids, Resnums, Segids, Atomtypes


class TestParallelSASAAnalysis:

    @pytest.fixture
    def universe(self):
        u = make_Universe(n_frames=4)
        u.add_TopologyAttr(Atomnames(["H"] * len(u.atoms)))
        u.add_TopologyAttr(Resnames(["GLY"] * len(u.residues)))
        u.add_TopologyAttr(Resids(list(range(0, len(u.residues)))))
        u.add_TopologyAttr(Resnums(list(range(0, len(u.residues)))))
        u.add_TopologyAttr(Segids(["A"] * len(u.segments)))
        u.add_TopologyAttr(Atomtypes(["O"] * len(u.atoms)))
        return u

    def test_import(self):
        """ParallelSASAAnalysis can be imported from the analysis package."""
        from mdakit_sasa.analysis import ParallelSASAAnalysis as PSA
        assert PSA is not None

    def test_atom_selection(self, universe):
        """Atom selection is applied correctly."""
        analysis = ParallelSASAAnalysis(universe, select="index 0:9")
        assert analysis.atomgroup.n_atoms == 10

    def test_run_returns_self(self, universe):
        """run() returns self for chaining."""
        analysis = ParallelSASAAnalysis(universe, n_jobs=1)
        result = analysis.run()
        assert result is analysis

    def test_total_area_shape(self, universe):
        """results.total_area has one entry per frame."""
        analysis = ParallelSASAAnalysis(universe, n_jobs=1).run()
        assert analysis.results.total_area.shape == (len(universe.trajectory),)

    def test_residue_area_shape(self, universe):
        """results.residue_area has shape (n_frames, n_residues)."""
        analysis = ParallelSASAAnalysis(universe, n_jobs=1).run()
        n_frames = len(universe.trajectory)
        n_residues = len(universe.residues)
        assert analysis.results.residue_area.shape == (n_frames, n_residues)

    def test_mean_total_area(self, universe):
        """results.mean_total_area equals the mean of results.total_area."""
        analysis = ParallelSASAAnalysis(universe, n_jobs=1).run()
        assert np.isclose(
            analysis.results.mean_total_area,
            analysis.results.total_area.mean()
        )

    def test_frame_range(self, universe):
        """start/stop arguments limit the frames processed."""
        analysis = ParallelSASAAnalysis(universe, n_jobs=1).run(start=0, stop=2)
        assert analysis.results.total_area.shape == (2,)

    def test_n_jobs_parameter_accepted(self, universe):
        """n_jobs parameter is accepted and stored correctly without running."""
        # We only verify parameter storage here. Running calcStructuresParallel
        # inside a pytest process that has MDAnalysis OpenMP extensions loaded
        # causes a segfault on some platforms. Parallel execution is verified
        # separately in test_compat.py.
        analysis = ParallelSASAAnalysis(universe, n_jobs=4)
        assert analysis.n_jobs == 4

    def test_fallback_flag(self, universe):
        """_has_native_parallel reflects whether calcStructuresParallel is available."""
        import freesasa
        analysis = ParallelSASAAnalysis(universe, n_jobs=2)
        expected = hasattr(freesasa, "calcStructuresParallel")
        assert analysis._has_native_parallel == expected

    def test_with_atomgroup(self, universe):
        """Accepts an AtomGroup as input."""
        analysis = ParallelSASAAnalysis(universe.atoms, n_jobs=1).run()
        assert analysis.results.total_area.shape == (len(universe.trajectory),)
