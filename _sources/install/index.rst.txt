Installation
============

If you're familiar with Python virtual environments, feel free to install MTUQ as follows within an environment of your choice::

   git clone https://github.com/uafgeotools/mtuq.git
   cd mtuq
   pip install -e .


Otherwise, follow these instructions for installation under `conda` (preferred) or `virtualenv`:

.. toctree::

   env_conda
   env_virtualenv



Graphics dependencies
---------------------

.. warning::

   To avoid graphics problems, we recommend installing GMT (>v6.0.0). 

MTUQ uses `GMT <https://www.soest.hawaii.edu/gmt>`_ graphics utilities if present on the system path.  Otherwise, MTUQ falls back to ObsPy, which does not always attempt to plot full moment tensors correctly (see this `issue <https://github.com/obspy/obspy/issues/2388>`_).  For installing GMT under `conda` or `virtualenv`, see above links.

