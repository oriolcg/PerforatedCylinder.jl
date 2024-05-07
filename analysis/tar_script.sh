#!/bin/sh
#
#
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=thin
#SBATCH --time=1-00:00:00
#SBATCH -n 588
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh

directory=/gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/VTKs
currdir=$PWD
cd $directory


# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' not found."
    exit 1
fi

# Create an empty array to hold the directory names
folders=()

# Loop through directories in the specified directory
for folder in "$directory"/serial*/; do
    if [ -d "$folder" ]; then
        folders+=("${folder%/}")  # Add the directory name to the array
    fi
done

echo "First 16 directories starting with 'serial' in $directory:"
count=0
for ((i = 0; i < ${#folders[@]} && count < 588; i++)); do
    echo "${folders[i]}"
    cmd="srun -N1 -n1 -c1 --exact tar cvzf ${folders[i]}.tgz ${folders[i]} &"
    echo $cmd
    eval $cmd
    ((count++))
done
wait


# declare -a dirs
# path=/gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/VTKs
# list=($(ls "${path}"))
# echo $list
# # short_list=${list:0:15}
# # echo ${short_list}
# currdir=$PWD
# cd $path
# INITIAL_CASE=1
# FINAL_CASE=16
# for i in $(seq $INITIAL_CASE $FINAL_CASE)
# do
#     echo $i
#     # echo "  \n"
#     # echo "$list[${i}]"
#     # echo $d
#     # cmd="srun -N1 -n1 -c1 --exact tar cvzf ${d}.tgz ${d} &"
#     # echo $cmd
#     # eval $cmd
# done
# wait
cd $currdir
