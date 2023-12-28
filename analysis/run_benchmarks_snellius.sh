#!/bin/bash

INITIAL_CASE=1
FINAL_CASE=8
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $(( 2 * $i))"
    sbatch run_benchmark_snellius_$(( 2 * $i)).sh
done