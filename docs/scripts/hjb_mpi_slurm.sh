#!/bin/bash
# Multi-node run of docs/scripts/hjb_mpi_benchmarks.jl on HiPerGator, one MPI
# rank per Slurm task, with the cluster's OpenMPI 5 (UCX over InfiniBand) and
# Slurm's PMIx, every rank inside the session sandbox (in-sandbox). Submit as
#
#   sbatch -A fairbanksj -p hpg-default -N <nodes> --ntasks-per-node=<ranks> -c <threads> \
#          --mem-per-cpu=4G -t 4:00:00 docs/scripts/hjb_mpi_slurm.sh
#
# with HJB_SIZES etc. exported (see hjb_mpi_benchmarks.jl). One-time setup of
# the Julia environment $HJB_MPI_ENV (default $HOME/bench-sysmpi), inside a job
# with the modules below loaded:
#   julia --project=$HJB_MPI_ENV -e 'using Pkg; Pkg.develop(path="<this repo>"); Pkg.add(["MPI", "MPIPreferences", "Statistics"])'
#   julia --project=$HJB_MPI_ENV -e 'using MPIPreferences; MPIPreferences.use_system_binary(;
#       library_names = ["/apps/mpi/gcc/14.2.0/openmpi/5.0.7_el97/lib/libmpi.so"], mpiexec = "srun")'
# Without the sandbox, drop `in-sandbox`.
env=${HJB_MPI_ENV:-$HOME/bench-sysmpi}
# Submitted with sbatch --export=..., the job also gets the submitting login
# environment (--get-user-env), including XDG_RUNTIME_DIR of the submit host.
# Apptainer keeps the state of in-sandbox's per-node MPI container instance
# there; the directory doesn't exist on the compute node, so the instance fails
# and every rank is killed (SIGTERM) within seconds. Drop it, and lmod's shell
# hooks, before starting the tasks.
unset BASH_ENV XDG_RUNTIME_DIR
unset -f module ml which 2>/dev/null || true
# in-sandbox passes only Slurm/PMI/MPI variables into the container: carry the
# benchmark settings explicitly.
settings=$(env | grep -E '^HJB_' | sed 's/^/export /; s/=/="/; s/$/";/' | tr '\n' ' ')
srun --mpi=pmix in-sandbox bash -lc "source /etc/profile.d/modules.sh && module load gcc/14.2.0 && \
    module load openmpi/5.0.7 && module load julia/1.12.6 && \
    export UCX_WARN_UNUSED_ENV_VARS=n JULIA_THREAD_SLEEP_THRESHOLD=infinite && $settings \
    julia -t \${SLURM_CPUS_PER_TASK:-1} --project=$env --startup-file=no docs/scripts/hjb_mpi_benchmarks.jl"
