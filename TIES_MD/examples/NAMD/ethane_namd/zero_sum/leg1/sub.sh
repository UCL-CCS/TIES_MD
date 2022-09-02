#!/bin/bash
#Example script for ARCHER2 NAMD2
#SBATCH --job-name=LIGPAIR
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --account=XXX
#SBATCH --partition=standard
#SBATCH --qos=standard

module load namd/2.14-nosmp

#--nodes and nodes_per_namd can be scaled up for large simulations
nodes_per_namd=1
cpus_per_namd=128

#change this line to point to your project
ties_dir=/home/a/post_qa/TIES_MD/TIES_MD/examples/ethane_namd/zero_sum/leg1
cd $ties_dir/replica-confs

# Looping over minimization, equilibration and production stages
for stage in {0..3}; do
for lambda in 0.00 1.00; do
 for i in {0..1}; do
        /home/a/NAMD/NAMD_3.0alpha13_Linux-x86_64-multicore-CUDA/namd3 --tclmain run$stage.conf $lambda $i &
        sleep 1
done
done
wait
done
