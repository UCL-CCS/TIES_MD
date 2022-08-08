Parallelization
================

Alchemical free energy calculations can be parallelized over numerous domains. Some domains of parallelization can be used in
any kind of molecular dynamics simulation, such as the spatial domain were a simulation box is decomposed into smaller cells
all run in parallel. These domains are, in general, more difficult to achieve parallelization across than the ones we discuss here which
focus on alchemical calculations. The two domains we focus on here are repeat/ensemble simulations and alchemical windows.
Ensemble simulation are critical to control the aleatoric error inherent in chaotic molecular dynamics simulation. Each simulation
in an ensemble has no communication with the other simulations and so this is an embarrassingly parallel problem, or a problem for which
parallelization is easy to implement. Likewise there is no communication between individual alchemical windows of the simulation
and so parallelizing these windows is also easy. The remainder of this page will explore how to achieve this parallelization
using ``OpenMM`` and ``NAMD`` with ``TIES``.

TIES-OpenMM
-----------

For `reference <https://github.com/UCL-CCS/TIES_MD/tree/master/TIES_MD/examples/ethane/zero_sum/leg1>`_ we will consider
running an example system from our ``TIES MD`` ``Github`` page. This example can be run without parallelization using this line::

    ties_md --exp_name=sys_solv

This would use 1 available GPU to execute all 8 alchemical windows and the 3 repeat specified in the config file ``TIES.cfg``
If we wanted to parallelize 3 repeats over 3 GPUs on one node we would run::

    ties_md --exp_name=sys_solv --devices=0,1,2

Each ``CUDA`` device will then run 8 windows of the 1 replica. Equally ths could be spit into to separate runs of ``TIES MD``
masked to only see one device::

    ties_md --exp_name=sys_solv --devices=0 --rep_id=0&
    ties_md --exp_name=sys_solv --devices=1 --rep_id=1&
    ties_md --exp_name=sys_solv --devices=2 --rep_id=2&

To run in this configuration the options ``total_reps=3`` and ``split_run=1`` are set in TIES.cfg to tell ``TIES MD`` that
there are a total of 3 replicas being run and that each execution of ``TIES MD`` should run only one. ``--rep_id``
determines which replica each instance will run. ``--rep_id`` only needs to be set when using ``split_run=1``.

If we need further parallelization over alchemical windows we can use the command line option ``--windows_mask``
this option takes a ``Python`` range (start inclusive and end exclusive) of the windows which that instance of
``TIES MD`` should run.::

    ties_md --exp_name=sys_solv --windows_mask=0,1 --devices=0&
    ties_md --exp_name=sys_solv --windows_mask=1,2 --devices=1&
    ties_md --exp_name=sys_solv --windows_mask=2,3 --devices=2&
    ties_md --exp_name=sys_solv --windows_mask=3,4 --devices=3&
    ties_md --exp_name=sys_solv --windows_mask=4,5 --devices=4&
    ties_md --exp_name=sys_solv --windows_mask=5,6 --devices=5&
    ties_md --exp_name=sys_solv --windows_mask=6,7 --devices=6&
    ties_md --exp_name=sys_solv --windows_mask=7,8 --devices=7&

Now sing the configuration options ``total_reps=3`` and ``split_run=0`` the above runs 3 replica of each alchemical
window on a different GPU.

For maximum parallelism we combine parallelizing over replicas and alchemical windows. For clarity we now consider the
same example as above but now with 6 alchemical windows, 2 replica simulations and one simulation per GPU, so in
TIES.cfg ``global_lambdas=0.0, 0.1, 0.4, 0.6, 0.9, 1.0``, ``total_reps=2`` and ``split_run=1``. To scale over multiple node
we could use the resource allocator of the HPC for example `jsrun <https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=SSWRJV_10.1.0/jsm/jsrun.html>`_
on `Summit <https://www.olcf.ornl.gov/summit/>`_. would allow us to run with 2 replicas of 6 windows as follows::

    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=0,1 --rep_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=1,2 --rep_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=2,3 --rep_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=3,4 --rep_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=4,5 --rep_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=5,6 --rep_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=0,1 --rep_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=1,2 --rep_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=2,3 --rep_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=3,4 --rep_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=4,5 --rep_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=5,6 --rep_id=1&

Note here we do not set ``--devices`` as the masking of GPUs is handled by the resource allocator, this is not the general case.
If a resource allocator is not available an alternative method to run multiple simulations across nodes is to use a message passing interface
(``MPI``). The use of ``MPI`` can vary from system to system and there is no universal solution to running across many node
for all HPC systems, however we provide an example (:ref:`NAMD 3`) which would work with
`ThetaGPU <https://www.alcf.anl.gov/support-center/theta/theta-thetagpu-overview>`_.

TIES-NAMD
---------

The parallelization of TIES in ``NAMD2`` follows the same ideas as ``OpenMM`` above. We want to run independent simulations
for all alchemical window and replica simulations. If in TIES.cfg ``split_run=0`` the submission script that
``TIES_MD`` writes will use the ``NAMD`` option ``+replicas X`` this makes each ``NAMD`` run ``X`` replicas and the
run lines in sub.sh will look something like::

    for stage in {0..3}; do
    for lambda in 0.00 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0; do
            srun -N $nodes_per_namd -n $cpus_per_namd namd2 +replicas 5 --tclmain run$stage-replicas.conf $lambda&
            sleep 1
    done
    wait
    done

Alternatively if ``split_run=1`` the run lines will look like::

    for stage in {0..3}; do
    for lambda in 0.00 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0; do
    for i in {{0..4}}; do
        srun -N $nodes_per_namd -n $cpus_per_namd namd2 --tclmain run$stage.conf $lambda $i &
        sleep 1
    done
    done
    wait
    done

Notice now the additional loop over ``$i``. So these run line are creating 65 different instances of ``NAMD`` each
running 1 replica and one alchemical window. Anecdotally using the ``+replicas`` results in less crashes and
we have tested up to ``+replicas 135`` on `ARCHER 2 <https://www.archer2.ac.uk/>`_ with no crashes. In the two above
examples the parallelism over alchemical windows is achieved in the loop over lambda.

Using ``NAMD3`` parallelization can be achieved like so (:ref:`NAMD 3`). ``NAMD`` in general has extensive options to provision
hardware and achieve parallelism, what have outlined here is not exhaustive and we would suggest consulting
the `documentation <https://www.ks.uiuc.edu/Research/namd/2.14/ug/>`_ for more a more comprehensive information.

