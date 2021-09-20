Installation
============

TIES MD
-----------

``TIES_MD`` can be installed with ``Conda`` on ``Linux-64`` and ``Linux-ppc64le`` machines. The steps to do this are as follows,
starting with an install of ``Miniconda``::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh -b -p $prefix
    export PATH="$prefix/bin:$PATH"


Ensure the ``Miniconda`` used matches the platform on which you are running, for example use ``Miniconda3-latest-Linux-ppc64le.sh``
for ``Linux-ppc64le`` machines, and that ``$prefix`` is set to some directory with read write permissions. We provide ``Conda`` installs for
``Python 3.7`` on ``Linux-64`` and ``Python 3.7`` on ``Linux-ppc64le``. If required you can change your python version in ``Conda``
using this command::

    conda install python=3.7.3

Finally add ``conda-forge`` to your available conda channels and install ``TIES MD``::

    conda config --add channels conda-forge
    conda install -c adw62 ties_md

These install instructions are all that is needed to run ``TIES MD`` with ``NAMD``. If you do not need ``OpenMM`` there is no need to read further.

TIES OpenMM
-----------


.. note::
    Crashes have been observed when using ``TIES_MD`` with ``OpenMM`` ``7.5`` and ``7.6``. For now we recommend using
    ``OpenMM`` ``7.4`` which can be installed as conda install -c omnia openmm==7.4.2

To use the OpenMM protocol in ``TIES_MD`` a custom version of ``OpenMMTools`` is needed to perform the alchemical transformations
of the system and allow for thermodynamic integration calculations. In order to install the custom version of ``OpenMMTools`` run::

    mkdir openmmtools_install
    cd openmmtools_install
    git clone -b adw62-PowerPC https://github.com/adw62/openmmtools.git
    pip install ./openmmtools --use-feature=in-tree-build

The install of ``OpenMM`` which was installed with ``TIES_MD`` can be verified by running::

    python -m openmm.testInstallation

for older ``OpenMM`` versions ``< 7.6`` this command was::

    python -m simtk.testInstallation

In some instances the wrong version of ``OpenMM`` and ``CUDAtoolkit`` could be installed this can be corrected by running::

    conda install -c conda-forge openmm cudatoolkit=10.0

where ``10.0`` should be replaced with the particular ``CUDA`` version you want to target. See the
``OpenMM`` `user guild <http://docs.openmm.org/latest/userguide/index.html>`_ for more details on this.

The install of ``TIES MD`` can be tested by downloading and running (:ref:`Tutorial`) any of the examples
provided `here <https://github.com/adw62/TIES_MD/tree/master/TIES_MD/examples>`_. These examples can be download by running::

    git clone https://github.com/adw62/TIES_MD.git


