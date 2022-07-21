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
``Python 3.7/3.8/3.9`` on ``Linux-64`` and ``Python 3.7`` on ``Linux-ppc64le``. If required you can change your python version in ``Conda``
using this command::

    conda install python=3.7.3

Finally add ``conda-forge`` to your available conda channels and install ``TIES MD``::

    conda config --add channels conda-forge
    conda install -c ucl-ccs ties_md

The install of ``OpenMM`` which was installed with ``TIES_MD`` can be verified by running::

    python -m openmm.testInstallation

for older ``OpenMM`` versions ``< 7.6`` this command was::

    python -m simtk.testInstallation

In some instances the wrong version of ``OpenMM`` and ``CUDAtoolkit`` could be installed. If this has happend the above
test of ``OpenMM`` will produce an output which looks like::

    OpenMM Version: 7.7
    Git Revision: 130124a3f9277b054ec40927360a6ad20c8f5fa6

    There are 4 Platforms available:

    1 Reference - Successfully computed forces
    2 CPU - Successfully computed forces
    3 CUDA - Error computing forces with CUDA platform
    4 OpenCL - Successfully computed forces

    CUDA platform error: Error loading CUDA module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)

    Median difference in forces between platforms:

    Reference vs. CPU: 6.30571e-06
    Reference vs. OpenCL: 6.76359e-06
    CPU vs. OpenCL: 8.05194e-07

    All differences are within tolerance.

the critical information here is ``3 CUDA - Error computing forces with CUDA platform`` and
``CUDA platform error: Error loading CUDA module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`` this can be corrected by
changing the install version of ``CUDAtoolkit`` like so::

    conda install -c conda-forge openmm cudatoolkit=10.0

where ``10.0`` should be replaced with the particular ``CUDA`` version you want to target. One can determine an
appropriate value for this by running ``nvidia-smi`` in a terminal which yields::

    +-----------------------------------------------------------------------------+
    | NVIDIA-SMI 460.80       Driver Version: 460.80       CUDA Version: 11.2     |
    |-------------------------------+----------------------+----------------------+
    | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
    |                               |                      |               MIG M. |
    |===============================+======================+======================|
    |   0  Quadro M1000M       Off  | 00000000:01:00.0 Off |                  N/A |
    | N/A   50C    P5    N/A /  N/A |    435MiB /  2002MiB |      4%      Default |
    |                               |                      |                  N/A |
    +-------------------------------+----------------------+----------------------+

The top right value here ``11.2`` can be used as the version of ``CUDA`` you wish to target. If ``nvidia-smi`` does not
return the above output your GPU and drivers my not be configured correctly.

The install of ``TIES MD`` can be tested by downloading and running (:ref:`Tutorial`) any of the examples
provided `here <https://github.com/UCL-CCS/TIES_MD/tree/main/TIES_MD/examples>`_. These examples can be download by running::

    git clone https://github.com/UCL-CCS/TIES_MD.git

TIES OpenMM linux-ppc64le
--------------------------

To use the OpenMM protocol in ``TIES_MD`` on ``linux-ppc64le`` a custom version of ``OpenMMTools`` is needed to perform the alchemical transformations
of the system and allow for thermodynamic integration calculations. In order to install the custom version of ``OpenMMTools`` run::

    mkdir openmmtools_install
    cd openmmtools_install
    git clone -b adw62-PowerPC https://github.com/adw62/openmmtools.git
    pip install ./openmmtools --use-feature=in-tree-build




