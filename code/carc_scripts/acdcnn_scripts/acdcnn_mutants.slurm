#!/bin/bash
#SBATCH --job-name=acdcnn  # create a short name for your job
#SBATCH --partition=general
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem 32G
#SBATCH --time=12:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/mutant_esm_acdcnn.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=ALL


module purge
module load miniconda3/latest
source activate acdcnn_fix

WILDTYPE_PROF_DIR="data/hh_suite/profiles/wildtypes"
MUTANT_PROF_DIR="data/hh_suite/profiles/mutants"
WILDTYPE_PDB_DIR="data/esmfold_outs/wildtypes_all"
MUTANT_PDB_DIR="data/esmfold_outs/mutants_all"
DDG_CSV_PATH="data/small_combined_cleaned.csv"
OUT_CSV_PATH="data/acdcnn/esmfold/mutants/ddg_predictions.csv"
python run_acdcnn_mutants.py -wild_prof_dir $WILDTYPE_PROF_DIR -mut_prof_dir $MUTANT_PROF_DIR -wild_pdb_dir $WILDTYPE_PDB_DIR -mut_pdb_dir $MUTANT_PDB_DIR -ddg_csv $DDG_CSV_PATH -o $OUT_CSV_PATH
