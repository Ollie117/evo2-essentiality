#!/bin/csh
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --error=evo-error.%j.log
#SBATCH --output=evo-run.%j.log
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --partition=DL

python evo_scoring_clean.py genes_evo_strategy_peturbed.csv essentiality_scores_test.csv cuda:0
