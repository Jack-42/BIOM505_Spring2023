#!/bin/bash
#SBATCH --job-name=esmfold_muts    # create a short name for your job
#SBATCH --partition=singleGPU
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem 22G
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --time=46:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/esmfold_mutant_batch.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=ALL


module purge
module load miniconda3/latest
module load gcc/8.3.0-wbma
module load cuda/11.6.2-jeoh
source activate esmfold

INP_DIR="data/esm_fastas/mutants"
OUT_DIR="data/esmfold_outs/mutants"
for fasta_file in $INP_DIR/*; do
	base_fname=${fasta_file##*/}
	base_fname=${base_fname%.*}
	out=$OUT_DIR/${base_fname}/
	python esmfold_inference.py -i $fasta_file -o $out --max-tokens-per-batch 0
done 
