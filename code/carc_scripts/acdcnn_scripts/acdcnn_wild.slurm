#!/bin/bash
#SBATCH --job-name=acdcnn_wild  # create a short name for your job
#SBATCH --partition=general
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem 32G
#SBATCH --time=06:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=outs/wild__alpha_acdcnn.out
#SBATCH --mail-user="your_email@unm.edu"
#SBATCH --mail-type=ALL


module purge
module load miniconda3/latest
source activate acdcnn_fix

PROF_DIR="data/hh_suite/profiles/wildtypes"
WILDTYPE_DIR="data/alphafold_outs/wildtypes"
DDG_CSV_PATH="data/small_combined_cleaned.csv"
OUT_CSV_PATH="data/acdcnn/alphafold/wildtypes/ddg_predictions.csv"
python run_acdcnn.py -prof_dir $PROF_DIR -pdb_dir $WILDTYPE_DIR -ddg_csv $DDG_CSV_PATH -o $OUT_CSV_PATH
