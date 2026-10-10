# Dependency system image

Loading CellularSheaves means loading roughly 250 package images. From `/blue` that takes about 40 s, and several minutes when the filesystem is busy. Every MPI rank and every CI step pays it again. A Julia system image built from the dependencies loads in seconds.

The image holds every dependency of CellularSheaves, plus MPI.jl, CUDA.jl and the test dependencies. It does **not** hold the packages under development:

- CellularSheaves itself;
- CliqueTrees and Mumblebee, together with anything that depends on them.

Those load on top of the image as ordinary packages, so you can edit them without rebuilding.

## Build

On a CPU node; it takes 20–40 minutes and a few GB of RAM:

```bash
sbatch -A fairbanksj -p hpg-default -c 8 --mem=32G -t 2:00:00 -o sysimage-%j.out sysimage/build_slurm.sh
```

This creates `$HOME/sysimage/`, or `$CS_SYSIMAGE_DIR` if set. It contains:

- `env/`: the environment the image is built from and run with. It develops this checkout, adds MPI, CUDA and the test dependencies, and holds the MPI and CUDA preferences in `LocalPreferences.toml`.
- `builder/`: PackageCompiler, kept outside the image.
- `cellularsheaves-deps-<julia>-<manifest hash>.so`: the image, plus `current.so`, a link to the latest one.

The build bakes in two compile-time choices, and other environments must not change them:

- **MPI:** the system OpenMPI 5.0.7. Set `CS_MPI_LIBRARY=jll` to keep MPI.jl's bundled MPICH instead.
- **CUDA:** the cluster toolkit 12.9, with CUDA loaded from the `cuda/12.9.1` module at run time. Set `CS_CUDA_VERSION=artifact` to let CUDA.jl download its own runtime instead.

## Run

```bash
sysimage/julia.sh script.jl                     # the image, its environment, its CPU target
sysimage/julia.sh --project=other script.jl     # another environment on the same image
```

Inside the sandbox, MPI ranks are launched as

```bash
srun --mpi=pmix in-sandbox bash -lc "module load gcc/14.2.0 openmpi/5.0.7 julia/1.12.6 && sysimage/julia.sh run.jl"
```

and GPU jobs add `module load cuda/12.9.1`.

## Editing CliqueTrees or Mumblebee

Develop the local checkout into the image's environment:

```bash
CS_DEVELOP=/path/to/CliqueTrees.jl julia sysimage/build.jl
```

This normally does not change the image: CliqueTrees isn't in it, so the build only re-resolves and confirms the image is up to date. Edits to the checkout take effect the next time it loads.

**Rebuild** when a change needs a *different version* of a package that's in the image (the manifest hash changes, and so does the image name). A package compiled into the image can't be replaced in a running session.

## CI

Build once per manifest, then run the tests on the image:

```bash
sysimage/julia.sh test/runtests.jl
```

The image name is keyed by the manifest hash, so a CI cache can reuse it until the dependencies change.
