#!/bin/bash
{header}

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