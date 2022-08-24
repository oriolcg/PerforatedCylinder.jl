#!/bin/sh

for i in {1..10}
do
    export INITIAL_CASE=$1 
    export LOCAL_CASE=$i
    sbatch run_CNN_NS_distributed.sh
done
