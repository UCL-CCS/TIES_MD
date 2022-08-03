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

    #!/bin/bash
    #COBALT -A XXX
    #COBALT -t 100
    #COBALT -n 2
    #COBALT -q full-node
    export mpirun="/lus/theta-fs0/software/thetagpu/openmpi-4.0.5/bin/mpirun"
    export namd3="/lus/theta-fs0/projects/CompBioAffin/awade/NAMD3/NAMD_3.0alpha9_Linux-x86_64-multicore-CUDA/namd3"
    node1=$(sed "1q;d" $COBALT_NODEFILE)
    node2=$(sed "2q;d" $COBALT_NODEFILE)

    cd /lus/theta-fs0/projects/CompBioAffin/awade/many_reps/mcl1/l18-l39/com/replica-confs
    for stage in {0..3}; do
      $mpirun -host $node1 --cpu-set 0 --bind-to core -np 1 $namd3 +devices 0 --tclmain sim$stage.conf 0.00 0&
      $mpirun -host $node1 --cpu-set 1 --bind-to core -np 1 $namd3 +devices 1 --tclmain sim$stage.conf 0.05 0&
      $mpirun -host $node1 --cpu-set 2 --bind-to core -np 1 $namd3 +devices 2 --tclmain sim$stage.conf 0.10 0&
      $mpirun -host $node1 --cpu-set 3 --bind-to core -np 1 $namd3 +devices 3 --tclmain sim$stage.conf 0.20 0&
      $mpirun -host $node1 --cpu-set 4 --bind-to core -np 1 $namd3 +devices 4 --tclmain sim$stage.conf 0.30 0&
      $mpirun -host $node1 --cpu-set 5 --bind-to core -np 1 $namd3 +devices 5 --tclmain sim$stage.conf 0.40 0&
      $mpirun -host $node1 --cpu-set 6 --bind-to core -np 1 $namd3 +devices 6 --tclmain sim$stage.conf 0.50 0&
      $mpirun -host $node1 --cpu-set 7 --bind-to core -np 1 $namd3 +devices 7 --tclmain sim$stage.conf 0.60 0&
      $mpirun -host $node2 --cpu-set 0 --bind-to core -np 1 $namd3 +devices 0 --tclmain sim$stage.conf 0.70 0&
      $mpirun -host $node2 --cpu-set 1 --bind-to core -np 1 $namd3 +devices 1 --tclmain sim$stage.conf 0.80 0&
      $mpirun -host $node2 --cpu-set 2 --bind-to core -np 1 $namd3 +devices 2 --tclmain sim$stage.conf 0.90 0&
      $mpirun -host $node2 --cpu-set 3 --bind-to core -np 1 $namd3 +devices 3 --tclmain sim$stage.conf 0.95 0&
      $mpirun -host $node2 --cpu-set 4 --bind-to core -np 1 $namd3 +devices 4 --tclmain sim$stage.conf 1.00 0&
    wait
    done

This script is running 13 alchemical windows using only 1 replica simulation in each window. Additionally 3 GPUs are idle
on node2. For real world application this script needs to be scaled up. Currently ``TIES MD`` will not attempt to build
``NAMD3`` HPC scripts automatically. For creating general scripts a ``Python`` script can be very helpful the following
script would allow us to scale up on ThetaGPU::

    import os

    if __name__ == "__main__":

        ###OPTIONS###

        #account name
        acc_name = 'XXX'
        #how many nodes do we want
        nodes = 9
        #what thermodynamic leg to run (these may have different wall times)
        leg = 'com'
        #Where is the namd3 binary
        namd3_exe = '/lus/theta-fs0/projects/CompBioAffin/awade/NAMD3/NAMD_3.0alpha9_Linux-x86_64-multicore-CUDA/namd3'

        #############

        cwd = os.getcwd()
        #give com and lig simulations differnt wall times if needed
        if leg == 'com':
            wall_time = 100
        else:
            wall_time = 60
        with open(os.path.join(cwd, 'thetagpu_{}.sub'.format(leg)), 'w') as f:

            #Writing a header
            f.write('#!/bin/bash\n')
            f.write('#COBALT -A {}\n'.format(acc_name))
            f.write('#COBALT -t {}\n'.format(wall_time))
            f.write('#COBALT -n {}\n'.format(nodes))
            f.write('#COBALT -q full-node\n')

            #exporting mpirun and namd3 install locations
            f.write('export mpirun=\"/lus/theta-fs0/software/thetagpu/openmpi-4.0.5/bin/mpirun\"\n')
            f.write('export namd3=\"/lus/theta-fs0/projects/CompBioAffin/awade/NAMD3/NAMD_3.0alpha9_Linux-x86_64-multicore-CUDA/namd3\"\n')

            #writing line to read node file
            for node in range(nodes):
                f.write('node{0}=$(sed \"{1}q;d\" $COBALT_NODEFILE)\n'.format(node+1, node+1))

            #move to ties directory
            f.write('cd {}\n'.format(os.path.join(cwd, 'replica-confs')))

            #iterate over minimization, NVT eq, NPT eq and production
            for stage in ['sim0', 'sim1', 'sim2', 'sim3']:
                count = 0
                node = 1
                #iterate over alchemical windows
                for lam in [0.00, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00]:
                    #iterate over replica siulations
                    for rep in [0, 1, 2, 3, 4]:
                        #write the run line
                        f.write('$mpirun -host $node{} --cpu-set {} --bind-to core -np 1 $namd3 +devices {} --tclmain {} {} {}&\n'.format(node, count%8, count%8, '{}.conf'.format(stage), lam, rep))
                        # count the number of gpus move to next node when gpus all filled
                        count += 1
                        if count%8 == 0:
                            node += 1
                #make sure we wait between simulation stages for all sims to finish
                f.write('wait\n')











