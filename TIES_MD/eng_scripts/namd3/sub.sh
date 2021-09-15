#!/bin/bash
#COBALT -A CompBioAffin
#COBALT -t 720
#COBALT -n 9
#COBALT -q full-node
export mpirun="/lus/theta-fs0/software/thetagpu/openmpi-4.0.5/bin/mpirun"
namd3=/path/to/namd3

node1=$(sed "1q;d" $COBALT_NODEFILE)
node2=$(sed "2q;d" $COBALT_NODEFILE)
node3=$(sed "3q;d" $COBALT_NODEFILE)
node4=$(sed "4q;d" $COBALT_NODEFILE)
node5=$(sed "5q;d" $COBALT_NODEFILE)
node6=$(sed "6q;d" $COBALT_NODEFILE)
node7=$(sed "7q;d" $COBALT_NODEFILE)
node8=$(sed "8q;d" $COBALT_NODEFILE)
node9=$(sed "9q;d" $COBALT_NODEFILE)

#change this line to point to your project
ties_dir=/work/PROJECT/STUDY/SYSTEM/LIGAND/THERMODYNAMIC_LEG/

for stage in {{0..2}}; do
win_id=0
for lambda in {lambs};
for rep in {reps};
do
do
        cd $ties_dir/replica-confs
        $mpirun -host $node1 -np 1 $namd3 +devices {} $rep --tclmain eq$stage-replicas.conf $lambda $win_id&
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