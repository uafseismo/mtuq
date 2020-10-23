Library reference
=================

.. warning::

   MTUQ is still in under development.  Some documentation may be incomplete or describe planned rather than current functionality.


Seismic data structures
-----------------------

============================================================================================================  ============================================================================================================
`mtuq.Dataset <generated/mtuq.Dataset.html>`_                                                                  Container for observed or synthetic seismic data
`mtuq.GreenTensor <generated/mtuq.GreensTensor.html>`_                                                         Green's tensor object
`mtuq.GreenTensorList <generated/mtuq.GreensTensorList.html>`_                                                 Container for Green's tensor objects
============================================================================================================  ============================================================================================================


Data I/O
--------

============================================================================================================  ============================================================================================================
`mtuq.download_greens_tensors <generated/mtuq.download_greens_tensors.html>`_                                  Downloads Green's tensors from syngine
`mtuq.open_db <generated/mtuq.open_db.html>`_                                                                  Opens Green's tensor database from disk
`mtuq.read <generated/mtuq.io.readers.SAC.read.html>`_                                                         Reads seismic data from disk
============================================================================================================  ============================================================================================================



Data processing and inversion
-----------------------------

============================================================================================================  ============================================================================================================
`mtuq.ProcessData <generated/mtuq.ProcessData.html>`_                                                          Data processing function constructor
`mtuq.Misfit <generated/mtuq.Misfit.html>`_                                                                    Misfit function constructor
`mtuq.grid_search <generated/mtuq.grid_search.grid_search.html>`_                                              Evaluates misfit over grids
`mtuq.MTUQDataArray <generated/mtuq.grid_search.MTUQDataArray.html>`_                                          Data structure for storing misfit values on regularly-spaced grids
`mtuq.MTUQDataFrame <generated/mtuq.grid_search.MTUQDataFrame.html>`_                                          Data structure for storing misfit values on irregularly-spaced grids
============================================================================================================  ============================================================================================================


Station and event structures
----------------------------

============================================================================================================  ============================================================================================================
`mtuq.Station <generated/mtuq.Station.html>`_                                                                  Station metdata object
`mtuq.Origin <generated/mtuq.Origin.html>`_                                                                    Origin metdata object
`mtuq.Wavelet <generated/mtuq.Wavelet.html>`_                                                                  Source wavelet or source-time function object
`mtuq.Force <generated/mtuq.Force.html>`_                                                                      Force source
`mtuq.MomentTensor <generated/mtuq.MomentTensor.html>`_                                                        Moment tensor source
============================================================================================================  ============================================================================================================


Data visualization
------------------

============================================================================================================  ============================================================================================================
`mtuq.graphics.plot_data_synthetics <generated/mtuq.graphics.plot_data_synthetics.html>`_                      Compares data and synthetics
`mtuq.graphics.plot_beachball <generated/mtuq.graphics.plot_beachball.html>`_                                  Plots focal mechansim (requires GMT)
============================================================================================================  ============================================================================================================


Uncertainty visualization
-------------------------

============================================================================================================  ============================================================================================================
`mtuq.graphics.plot_likelihood_lune <generated/mtuq.graphics.plot_likelihood_lune.html>`_                      Plots maximum likelihoods on lune (requires GMT)
`mtuq.graphics.plot_marginal_lune <generated/mtuq.graphics.plot_marginal_lune.html>`_                          Plots marginal likelihoods on lune (requires GMT)
`mtuq.graphics.plot_misfit_lune <generated/mtuq.graphics.plot_misfit_lune.html>`_                              Plots misfit on lune (requires GMT)
`mtuq.graphics.plot_likelihood_vw <generated/mtuq.graphics.plot_likelihood_vw.html>`_                          Plots maximum likelihoods on `v-w` rectangle
`mtuq.graphics.plot_marginal_vw <generated/mtuq.graphics.plot_marginal_vw.html>`_                              Plots marginal likelihoods on `v-w` rectangle
`mtuq.graphics.plot_misfit_vw <generated/mtuq.graphics.plot_misfit_vw.html>`_                                  Plots misfit on `v-w` rectangle
`mtuq.graphics.plot_likelihood_force <generated/mtuq.graphics.plot_likelihood_force.html>`_                    Plots force likelihoods on unit sphere (requires GMT)
`mtuq.graphics.plot_marginal_force <generated/mtuq.graphics.plot_marginal_force.html>`_                        Plots force marginals on unit sphere (requires GMT)
`mtuq.graphics.plot_misfit_force <generated/mtuq.graphics.plot_misfit_force.html>`_                            Plots force misfit on unit sphere (requires GMT)
============================================================================================================  ============================================================================================================



Moment tensor and force grids
-----------------------------

============================================================================================================  ============================================================================================================
`mtuq.grid.FullMomentTensorGridRandom <generated/mtuq.grid.FullMomentTensorGridRandom.html>`_                  Grid with randomly-drawn full moment tensors
`mtuq.grid.FullMomentTensorGridSemiregular <generated/mtuq.grid.FullMomentTensorGridSemiregular.html>`_        Grid with semiregularly-spaced full moment tensors
`mtuq.grid.DoubleCoupleGridRandom <generated/mtuq.grid.DoubleCoupleGridRandom.html>`_                          Grid with randomly-drawn double couples
`mtuq.grid.DoubleCoupleGridRegular <generated/mtuq.grid.DoubleCoupleGridRegular.html>`_                        Grid with regularly-spaced double couples
`mtuq.grid.ForceGridRandom <generated/mtuq.grid.ForceGridRandom.html>`_                                        Grid with randomly-drawn forces
`mtuq.grid.ForceGridRegular <generated/mtuq.grid.ForceGridRegular.html>`_                                      Grid with regularly-spaced forces
============================================================================================================  ============================================================================================================


Grid base classes
-----------------
.. note::
   These base classes can be used for specialized uses cases not covered by any of the above grids

============================================================================================================  ============================================================================================================
`mtuq.grid.Grid <generated/mtuq.grid.Grid.html>`_                                                              Grid defined by values on axes
`mtuq.grid.UnstructuredGrid <generated/mtuq.grid.UnstructuredGrid.html>`_                                      Grid defined by lists of coordinates points
============================================================================================================  ============================================================================================================



Database and web service clients
--------------------------------

.. note::
    These classes provide low-level database and web service access.  Normally it is easier to use high-level functions ``mtuq.download_greens_tensors`` or ``mtuq.open_db`` instead.

============================================================================================================  ============================================================================================================ 
`mtuq.io.clients.AxiSEM_NetCDF.Client <generated/mtuq.io.clients.AxiSEM_NetCDF.Client.html>`_                  AxiSEM NetCDF database client (based on instaseis)
`mtuq.io.clients.FK_SAC.Client <generated/mtuq.io.clients.FK_SAC.Client.html>`_                                FK database client
`mtuq.io.clients.syngine.Client <generated/mtuq.io.clients.syngine.Client.html>`_                              Syngine web service client
============================================================================================================  ============================================================================================================ 



Low-level misfit modules
------------------------

.. note::
    Misfit evaluation is the most computationally intensive task performed by this package.  These classes provide low-level access to non-optimized and optimized versions of the misfit function  (the former exists to check correctness of the latter).  Normally users can avoid this complexity by simply calling ``mtuq.misfit.Misfit``.

============================================================================================================  ============================================================================================================
`mtuq.misfit.level0.misfit <generated/mtuq.misfit.level0.html>`_                                               Easy-to-read pure Python implementation
`mtuq.misfit.level1.misfit <generated/mtuq.misfit.level1.html>`_                                               Fast pure Python implementation
`mtuq.misfit.level2.misfit <generated/mtuq.misfit.level2.html>`_                                               Fast Python/C implementation
============================================================================================================  ============================================================================================================

