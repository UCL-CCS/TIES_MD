HPC Submission scripts
======================

Here we provide some example submission scripts for various HPC systems.

NAMD
----

Here is an example of a submission script for a large system (≈100k atoms) running on
`SuperMUC-NG <https://doku.lrz.de/display/PUBLIC/SuperMUC-NG>`_::

    #!/bin/bash
    #SBATCH --job-name=LIGPAIR
    #SBATCH -o ./%x.%j.out
    #SBATCH -e ./%x.%j.err
    #SBATCH -D ./
    #SBATCH --nodes=130
    #SBATCH --tasks-per-node=48
    #SBATCH --no-requeue
    #SBATCH --export=NONE
    #SBATCH --get-user-env
    #SBATCH --account=XXX
    #SBATCH --partition=general
    #SBATCH --time=10:00:00

    module load slurm_setup
    module load namd/2.14-gcc8-impi

    nodes_per_namd=10
    cpus_per_namd=480

    echo $nodes_per_namd
    echo $cpus_per_namd

    #change this line to point to your project
    ties_dir=/hppfs/work/pn98ve/di67rov/test_TIES/study/prot/ties-l2-l1/com

    for stage in {0..3}; do
    for lambda in 0.00 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0;
    do
            cd $ties_dir/replica-confs
            srun -N $nodes_per_namd -n $cpus_per_namd namd2 +replicas 5 --tclmain sim$stage-replicas.conf $lambda&
            sleep 1
    done
    wait
    done

The first 20 lines of this script could be adapted for a smaller system (≈10k atoms) as follows::

    #!/bin/bash
    #SBATCH --job-name=LIGPAIR
    #SBATCH -o ./%x.%j.out
    #SBATCH -e ./%x.%j.err
    #SBATCH -D ./
    #SBATCH --nodes=13
    #SBATCH --tasks-per-node=45
    #SBATCH --no-requeue
    #SBATCH --export=NONE
    #SBATCH --get-user-env
    #SBATCH --account=XXX
    #SBATCH --partition=micro
    #SBATCH --time=10:00:00

    module load slurm_setup
    module load namd/2.14-gcc8-impi

    #--nodes and nodes_per_namd can be scaled up for large simulations
    nodes_per_namd=1
    cpus_per_namd=45


OpenMM
------

Here we provide an example of ``TIES MD`` running with ``OpenMM`` on `Summit <https://www.olcf.ornl.gov/summit/>`_::

    #!/bin/bash
    #BSUB -P XXX
    #BSUB -W 20
    #BSUB -nnodes 1
    #BSUB -alloc_flags "gpudefault smt1"
    #BSUB -J test
    #BSUB -o otest.%J
    #BSUB -e etest.%J
    cd $LS_SUBCWD
    export PATH="/gpfs/alpine/scratch/adw62/chm155/TIES_test/miniconda/bin:$PATH"
    export ties_dir="/gpfs/alpine/scratch/adw62/chm155/TIES_test/TIES_MD/TIES_MD/examples/ethane/zero_sum/leg1"
    module load cuda/10.1.168
    date
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv'  --windows_mask=0,1 --node_id="0" > $ties_dir/0.out&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv'  --windows_mask=1,2 --node_id="0" > $ties_dir/1.out&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv'  --windows_mask=2,3 --node_id="0" > $ties_dir/2.out&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv'  --windows_mask=3,4 --node_id="0" > $ties_dir/3.out&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv'  --windows_mask=4,5 --node_id="0" > $ties_dir/4.out&
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 ties_md --config_file=$ties_dir/TIES.cfg --exp_name='sys_solv'  --windows_mask=5,6 --node_id="0" > $ties_dir/5.out&
    wait

NAMD 3
------

Here we provide an example of ``TIES MD`` running with ``NAMD3`` on `ThetaGPU <https://www.alcf.anl.gov/alcf-resources/theta>`_::









