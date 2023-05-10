#!/bin/bash
#SBATCH --job-name=esm_test    # create a short name for your job
#SBATCH --partition=debug
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem 22G
#SBATCH --time=00:05:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/test_cpu.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=ALL


module purge
module load miniconda3/latest
module load gcc/8.3.0-wbma
module load cuda/11.6.2-jeoh
source activate esmfold

python esmfold_inference.py -i data/test/esmfold/sample.fasta -o data/esmfold_outs/ --cpu-only
