#!/usr/bin/bash
#$ -N sim_nuc_cyt
#$ -t 1-8
#$ -q long
#$ -cwd
#$ -l virtual_free=16G,h_rt=06:00:00
#$ -pe smp 2
#$ -V



# Parameter file and working directory for the analysis
par_file=$(find . -name '*.par' | grep -v 'whole_cell' | sed -n "${SGE_TASK_ID}"p)

par="${par_file##*/}"
work_dir="${par_file%/*}"

# Simulate RNAseq data
cd "${work_dir}"
flux --threads 2 -l -s -p "$par" &> flux.log
