#!/bin/bash
#SBATCH --job-name=hhblits    # create a short name for your job
#SBATCH --partition=bigmem-1TB
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem 256G
#SBATCH --time=48:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output="outs/%j.out"
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=END

# HH-Blits is expensive to run, so will create profiles for files in only this range from INP_DIR
FILE_I_LB=4000
FILE_I_UB=4200

module purge
module load miniconda3/latest
source activate hh_suite


# Create CSV file and write header
RUNTIMES_FILE="runtimes/${FILE_I_LB}_${FILE_I_UB}_mutant_runtimes.csv"
echo "file,runtime" > $RUNTIMES_FILE

# Iterate over fasta files in INP_DIR and output profiles for each
DB_PATH="data/uniclust30_2018/uniclust30_2018_08/uniclust30_2018_08"
INP_DIR="data/hh_suite/inps/hhsuite_fastas/mutants"
OUT_DIR="data/hh_suite/outs/mutants"
for fasta_file in $(ls $INP_DIR | sed -n "${FILE_I_LB},${FILE_I_UB}p"); do
	start=`date +%s`
        inp_path=$INP_DIR/${fasta_file}
	base_fname=${fasta_file##*/}
        base_fname=${base_fname%.*}
	out_hhr_file="${base_fname}.hhr"
	out_opsi_file="${base_fname}.msa.psi"
	out_hhr_path=$OUT_DIR/${out_hhr_file}
        out_opsi_path=$OUT_DIR/${out_opsi_file}
	hhblits -d $DB_PATH -i $inp_path -o $out_hhr_path -opsi $out_opsi_path -n 2 -cpu $SLURM_CPUS_PER_TASK
	end=`date +%s`
	runtime=$((end-start))
	
	# Write file name and runtime to console and CSV file
    	echo "$base_fname,$runtime" >> $RUNTIMES_FILE
done	
