Installation under virtualenv
=============================

Install `virtualenv`, if not already present, following these `instructions <https://virtualenv.pypa.io/en/latest/installation.html>`_.

Create a virtual environment:

.. code::

   virtualenv --python=python3 mtuq_env


Download MTUQ source code:

.. code::

   git clone https://github.com/uafgeotools/mtuq.git mtuq_src


Install MTUQ:

.. code::

   source mtuq_env/bin/activate
   cd mtuq_src
   pip install numpy
   pip install obspy instaseis
   pip install -e .


If `instaseis` installation fails, try:

.. code::

   CC=gfortran; pip install instaseis 



Unpack seismic waveforms used by examples:

.. code::

    bash ./data/examples/unpack.bash
    bash ./data/tests/unpack.bash


Finally, install GMT following these `instructions <https://github.com/GenericMappingTools/gmt/blob/master/INSTALL.md>`_.  

(For an automated GMT installation, consider `installing MTUQ under conda <https://uafgeotools.github.io/mtuq/install/env_conda.html>`_ instead.)

