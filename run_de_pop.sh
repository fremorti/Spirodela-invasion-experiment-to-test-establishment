#!/bin/sh

#PBS -N invasions_diff_eq_pop
#PBS -l walltime=23:30:00
#PBS -l nodes=1:ppn=4
#PBS -e error_pop.file
#PBS -o output_pop.file

cd $PBS_O_WORKDIR

module load R-bundle-CRAN/2024.06-foss-2023b
chmod 770 hpc_invasion_de.R
Rscript --vanilla hpc_invasion_de_pop.R