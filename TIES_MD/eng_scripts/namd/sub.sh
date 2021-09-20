#!/bin/bash
#Example script fot ARCHER2
#SBATCH --job-name=LIGPAIR
#SBATCH --nodes={nodes}
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --account=XXX
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env
module load namd/2.14-nosmp-gcc10

#--nodes and nodes_per_namd can be scaled up for large simulations
nodes_per_namd=1
cpus_per_namd=128

#change this line to point to your project
ties_dir=/work/PROJECT/STUDY/SYSTEM/LIGAND/THERMODYNAMIC_LEG/

for stage in {{0..2}}; do
win_id=0
for lambda in {lambs};
do
        cd $ties_dir/replica-confs
        srun -N $nodes_per_namd -n $cpus_per_namd --distribution=block:block --hint=nomultithread namd2 +replicas {reps} --tclmain eq$stage-replicas.conf $lambda $win_id&
        (( win_id++ ))
        sleep 1
done
wait
done

for stage in {{1..1}}; do
win_id=0
for lambda in {lambs};
do
        cd $ties_dir/replica-confs
        srun -N $nodes_per_namd -n $cpus_per_namd --distribution=block:block --hint=nomultithread namd2 +replicas {reps} --tclmain sim$stage-replicas.conf $lambda $win_id&
        (( win_id++ ))
        sleep 1
done
wait
done