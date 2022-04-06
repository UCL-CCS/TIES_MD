#!/bin/bash
{header}

#--nodes and nodes_per_namd can be scaled up for large simulations
nodes_per_namd={nodes_per_namd}
cpus_per_namd= 128 \* $nodes_per_namd #128 cpus per node on ARCHER2

#change this line to point to your project
ties_dir={root}
cd $ties_dir/replica-confs

for stage in {{0..2}}; do
win_id=0
for lambda in {lambs}; do
 for i in {{1..{reps}}}; do
        {run_line}
        (( win_id++ ))
        sleep 1
done
done
wait
done

for stage in {{1..1}}; do
win_id=0
for lambda in {lambs}; do
for i in {{1..{reps}}}; do
        {run_line}
        (( win_id++ ))
        sleep 1
done
done
wait
done