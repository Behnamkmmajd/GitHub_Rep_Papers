#!/bin/bash
#SBATCH --job-name=AE2DGUV
#SBATCH --account=share-nv-ikp
#SBATCH --partition=GPUQ
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint="gpu80g"
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=Log_%j.log
#SBATCH --exclude=idun-09-[01-11,19-20]

set -euo pipefail

module purge
module load Python/3.10.4-GCCcore-11.3.0
source ~/venvs/bubbleAE_gpu/bin/activate

cd /cluster/work/behnamk/non_linear_data_driven/2D_test

srun --ntasks=1 --gpus-per-task=1 python -u AE_2D_edited_multiple_channel.py
