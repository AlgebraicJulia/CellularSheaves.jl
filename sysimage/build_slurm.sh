#!/bin/bash
# Build the dependency image on a CPU node of HiPerGator, inside the session
# sandbox:
#   sbatch -A fairbanksj -p hpg-default -c 8 --mem=32G -t 2:00:00 -o sysimage-%j.out sysimage/build_slurm.sh
# The MPI module provides the system libmpi baked into the image; the CUDA
# toolkit is picked up from the cuda module at run time.
in-sandbox bash -lc "source /etc/profile.d/modules.sh && module load gcc/14.2.0 && module load openmpi/5.0.7 && \
    module load julia/1.12.6 && julia --startup-file=no sysimage/build.jl"
