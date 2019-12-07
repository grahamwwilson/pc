#!/bin/sh

sbatch job_launch.sh testdata
sbatch job_launch.sh testmc
sbatch job_launch.sh testmc_R1

exit
