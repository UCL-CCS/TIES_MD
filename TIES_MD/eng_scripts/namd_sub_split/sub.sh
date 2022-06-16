#!/bin/bash
{header}

#change this line to point to your project
ties_dir={root}
cd $ties_dir/replica-confs

for stage in {{0..2}}; do
for lambda in {lambs}; do
 for i in {{1..{reps}}}; do
        {run_line}
        sleep 1
done
done
wait
done

for stage in {{1..1}}; do
for lambda in {lambs}; do
for i in {{1..{reps}}}; do
        {run_line}
        sleep 1
done
done
wait
done