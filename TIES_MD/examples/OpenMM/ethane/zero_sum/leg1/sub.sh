#!/bin/bash
#Example script for Summit OpenMM
#BSUB -P XXX
#BSUB -W 240
#BSUB -nnodes 4
#BSUB -alloc_flags "gpudefault smt1"
#BSUB -J LIGPAIR
#BSUB -o oLIGPAIR.%J
#BSUB -e eLIGPAIR.%J

export ties_dir="/home/a/post_qa/TIES_MD/TIES_MD/examples/ethane/zero_sum/leg1"
cd $ties_dir

for lambda in 0 1 2 3 4 5 6 7; do
    jsrun --smpiargs="off" -n 1 -a 1 -c 1 -g 1 -b packed:1 TIES_MD --config_file=$ties_dir/TIES.cfg --exp_name=complex --windows_mask=$lambda,$(expr $lambda + 1) --devices=0 > $ties_dir/$lambda_$i.out&
    done
wait
