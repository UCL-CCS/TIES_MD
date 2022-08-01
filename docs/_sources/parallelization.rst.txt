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

This would use 1 available GPU to execute all 6 alchemical windows and the 1 repeat specified in the config file ``TIES.cfg``
If we wanted to parallelize 2 repeats over 2 GPUs on one node we would change the option repeats in ``TIES.cfg`` to equal 2
and run::

    ties_md --exp_name=sys_solv --devices=0,1

The ``CUDA`` device 0 will then run 6 windows of the first repeat and CUDA device 1 will run 6 windows of the second repeat.
Equally ths could be spit into to runs of ``TIES MD`` masked to only see one device::

    ties_md --exp_name=sys_solv --devices=0 --node_id=0
    ties_md --exp_name=sys_solv --devices=1 --node_id=1

To run in this configuration the options ``total_reps=2`` and ``reps_per_exec=1`` are set in TIES.cfg to tell ``TIES MD`` that
there are a total of 2 replicas being run and that each execution of ``TIES MD`` should run only one. Also note we have set
``--node_id`` to some different values for otherwise identical run lines and this ensures these parallel runs write output
to unique locations. ``--node_id`` only needs to be set when identical replicas of a simulation are run in separate executions
of ``TIES MD``.

If we need further parallelization over alchemical windows we can use the command line option ``--windows_mask``
this option takes a ``Python`` range (start inclusive and end exclusive) of the windows which that instance of
``TIES MD`` should run. So returning to the original reference example with one repeat we would run::

    ties_md --exp_name=sys_solv --windows_mask=0,1 --devices=0&
    ties_md --exp_name=sys_solv --windows_mask=1,2 --devices=1&
    ties_md --exp_name=sys_solv --windows_mask=2,3 --devices=2&
    ties_md --exp_name=sys_solv --windows_mask=3,4 --devices=3&
    ties_md --exp_name=sys_solv --windows_mask=4,5 --devices=4&
    ties_md --exp_name=sys_solv --windows_mask=5,6 --devices=5&

These commands run submitted to a node with 6 GPUS would run one window on each GPU. To scale over multiple node
we could use the resource allocator of the HPC for example `jsrun <https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=SSWRJV_10.1.0/jsm/jsrun.html>`_
on Summit ``ORNL`` would allow us to run with 2 replicas of 6 windows as follows::

    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=0,1 --node_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=1,2 --node_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=2,3 --node_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=3,4 --node_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=4,5 --node_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=5,6 --node_id=0&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=0,1 --node_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=1,2 --node_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=2,3 --node_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=3,4 --node_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=4,5 --node_id=1&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=5,6 --node_id=1&

Note here we do not set ``--devices`` as the masking of GPUs is handled by the resource allocator. If a resource allocator
is not available an alternative method to run multiple simulations across nodes is to use a message passing interface
(``MPI``). The use of ``MPI`` can vary from system to system and there is no universal solution to running across many node
for all HPC systems, however we provide an example here for reference which would work with
`ThetaGPU <https://www.alcf.anl.gov/support-center/theta/theta-thetagpu-overview>`_::

    mpirun -host $node1 -np 1 ties_md --devices=0 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=0,1 --node_id=0&
    mpirun -host $node1 -np 1 ties_md --devices=1 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=1,2 --node_id=0&
    mpirun -host $node1 -np 1 ties_md --devices=2 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=2,3 --node_id=0&
    mpirun -host $node1 -np 1 ties_md --devices=3 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=3,4 --node_id=0&
    mpirun -host $node1 -np 1 ties_md --devices=4 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=4,5 --node_id=0&
    mpirun -host $node1 -np 1 ties_md --devices=5 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=5,6 --node_id=0&
    mpirun -host $node2 -np 1 ties_md --devices=0 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=0,1 --node_id=1&
    mpirun -host $node2 -np 1 ties_md --devices=1 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=1,2 --node_id=1&
    mpirun -host $node2 -np 1 ties_md --devices=2 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=2,3 --node_id=1&
    mpirun -host $node2 -np 1 ties_md --devices=3 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=3,4 --node_id=1&
    mpirun -host $node2 -np 1 ties_md --devices=4 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=4,5 --node_id=1&
    mpirun -host $node2 -np 1 ties_md --devices=5 --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv' --windows_mask=5,6 --node_id=1&


TIES-NAMD
---------

Under construction