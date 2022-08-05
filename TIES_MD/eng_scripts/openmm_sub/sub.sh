#!/bin/bash
{header}

export ties_dir="{root}"
cd $ties_dir

for lambda in {lambs}; do
    {run_line}
    done
wait
