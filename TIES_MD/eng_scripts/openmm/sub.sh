#!/bin/bash


export PATH="/ccs/home/adw62/OpenMM/miniconda3/bin:$PATH"
export ties_dir="{root}"
cd $ties_dir

module load cuda/10.1.168
date

for lambda in {lambs}; do
  for i in {{1..{reps}}}; do
        {run_line}
    done
    done
wait
