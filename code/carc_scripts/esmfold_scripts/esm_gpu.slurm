#!/bin/bash
#SBATCH --job-name=esmfold_test    # create a short name for your job
#SBATCH --partition=singleGPU
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem 22G
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --time=00:30:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/esmfold_gpu_benchmarking.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=ALL


module purge
module load miniconda3/latest
module load gcc/8.3.0-wbma
module load cuda/11.6.2-jeoh
source activate esmfold
 
python esmfold_inference.py -i data/test/esmfold/benchmarking_seqs.fasta -o data/esmfold_outs/benchmarking/ --max-tokens-per-batch 0
