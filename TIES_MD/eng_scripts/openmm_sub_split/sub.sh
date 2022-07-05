#!/bin/bash
{header}

export PATH="{py_bin}:$PATH"
export ties_dir="{root}"
cd $ties_dir

for lambda in {lambs}; do
  for i in {{0..{reps}}}; do
        {run_line}
    done
    done
wait
