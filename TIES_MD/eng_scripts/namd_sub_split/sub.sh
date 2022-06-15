#!/bin/bash
{header}

#change this line to point to your project
ties_dir={root}
cd $ties_dir/replica-confs

for stage in {{0..2}}; do
win_id=0
for lambda in {lambs}; do
 for i in {{1..{reps}}}; do
        {run_line}
        sleep 1
done
(( win_id++ ))
done
wait
done

for stage in {{1..1}}; do
win_id=0
for lambda in {lambs}; do
for i in {{1..{reps}}}; do
        {run_line}
        sleep 1
done
(( win_id++ ))
done
wait
done