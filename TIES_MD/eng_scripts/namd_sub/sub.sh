#!/bin/bash
{header}

#change this line to point to your project
ties_dir={root}
cd $ties_dir/replica-confs

# Looping over minimization, equilibration and production stages
for stage in {{0..3}}; do
for lambda in {lambs};
do
        {run_line}
        sleep 1
done
wait
done