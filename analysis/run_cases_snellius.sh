#!/bin/bash

INITIAL_CASE=1
FINAL_CASE=27
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    sbatch run_case_snellius.sh $i
done