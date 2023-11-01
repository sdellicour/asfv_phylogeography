#!/bin/bash
#SBATCH --job-name=ASFV_1_wo_4_tips
#SBATCH --time=120:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10240
#SBATCH --qos=gpu_prio
 
module load releases/2019b
module load beagle-lib/3.1.2-gcccuda-2019b

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

java -jar beast_dev_1105_290920.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite ASFV_1_wo_4_tips.xml
