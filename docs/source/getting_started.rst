Getting Started
===============

This software is available and tested on Linux/Mac/Windows

Installation
------------

This sofware can be installed from source using the documentation under `github`_,  alterantive using a python package manager.


Using `Pip`_ ::

    pip install mdakit-sasa

Using `Conda`_ ::
    
    TBC    



SASA analisys
---------------------

If you already familiar with MDAnaysis the first steps are simple::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC

    u = mda.Universe(PSF, DCD)

We can now use the universe and the trajectory to compute SASA::

    from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis

    analysis = SASAAnalysis(u)
    analysis.run()

    print(analysis.results.total_area)

Where the total area of the selection is provided for each individual frame::

    [10418.26874946 10254.94074755 10285.61490128 10190.08823062
    10129.67104807 10161.28481992 10193.53742475 10080.24520417

    [...]
    11417.53736883 11525.57789122]



.. _github: https://github.com/pegerto/mdakit_sasa/blob/main/README.md
.. _pip: https://pypi.org/project/pip/
.. _conda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html


