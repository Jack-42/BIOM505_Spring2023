#!/bin/bash
#SBATCH --job-name=ddgun_profs    # create a short name for your job
#SBATCH --partition=general
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem 32G
#SBATCH --time=04:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/mk_profs_mut.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=END


module purge
module load miniconda3/latest
source activate acdcnn_fix

# Iterate over fasta files in INP_DIR and output ddgun profiles for each
INP_DIR="data/hh_suite/outs/mutants"
OUT_DIR="data/hh_suite/profiles/mutants"
for psi_file in $INP_DIR/*.psi; do
	if [ -f "$psi_file" ]; then
		echo $psi_file
		base_fname=${psi_file##*/}
        	base_fname=${base_fname%.*}
		base_fname=${base_fname/.msa/} # remove ".msa" extension
        	out_prof_file="${base_fname}.prof"
        	out_prof_path=$OUT_DIR/${out_prof_file}
		echo $out_prof_path
        	ddgun mkprof $psi_file > $out_prof_path
	fi
done	
