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


#TO BE ADDED