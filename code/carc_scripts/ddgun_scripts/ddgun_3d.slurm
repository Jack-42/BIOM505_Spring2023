#!/bin/bash
#SBATCH --job-name=ddgun # create a short name for your job
#SBATCH --partition=bigmem-1TB
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem 256G
#SBATCH --time=48:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/ddgun3d_alphafold.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=ALL


module purge
module load miniconda3/latest
source activate ddgun

DB_PATH="data/uniclust30_2018/uniclust30_2018_08/uniclust30_2018_08"
DDGUN_SCRIPT="ddgun/ddgun_3d.py"
WILDTYPE_DIR="data/experimental_structures/wildtypes_A_chain"
MUTFILE_DIR="data/ddgun/mut_files"
OUT_DIR="data/ddgun/out_files/experimental_structures"
for pdb in $WILDTYPE_DIR/*.pdb
do
    base=$(basename $pdb .pdb)
    mutfile=$MUTFILE_DIR/$base.muts
    out=$OUT_DIR/$base.out
    python $DDGUN_SCRIPT $pdb 'A' $mutfile -o $out -d $DB_PATH
done
