#!/bin/bash
{header}

#change this line to point to your project
ties_dir={root}
cd $ties_dir/replica-confs

# these loops are for equilibration
for stage in {{0..2}}; do
for lambda in {lambs};
do
        {run_line0}
        sleep 1
done
wait
done

#this is for production
for stage in {{1..1}}; do
for lambda in {lambs};
do
        {run_line1}
        sleep 1
done
wait
done