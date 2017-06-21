#!/bin/bash
#
#SBATCH --array=0-1
#SBATCH --job-name=1525f39f9f618
#SBATCH --output=slurm_%a.out
/Library/Frameworks/R.framework/Resources/bin/Rscript --vanilla slurm_run.R
