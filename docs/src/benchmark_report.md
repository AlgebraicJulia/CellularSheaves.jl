# Benchmark report

- **Profile:** `small`
- **Expected shards:** `assembly-large, assembly-small, extension-large, extension-small, solver-large, solver-small`
- **Completed shards:** `assembly-large, assembly-small, extension-large, extension-small, solver-large, solver-small`
- **Missing shards:** `none`

| Benchmark ID | Shard | Median time | Memory | Allocs | Runner | Host |
| --- | --- | ---: | ---: | ---: | --- | --- |
| `coboundary_map/cycle/n100` | `assembly-small` | 120.155 us | 83.89 KiB | 1146 | `slurm` | `c0705a-s11.ufhpc` |
| `coboundary_map/cycle/n20` | `assembly-small` | 29.004 us | 14.44 KiB | 257 | `slurm` | `c0705a-s11.ufhpc` |
| `coboundary_map/cycle/n500` | `assembly-large` | 558.243 us | 367.43 KiB | 5564 | `slurm` | `c0705a-s11.ufhpc` |
| `coboundary_map/path/n100` | `assembly-small` | 115.887 us | 83.41 KiB | 1135 | `slurm` | `c0705a-s11.ufhpc` |
| `coboundary_map/path/n20` | `assembly-small` | 26.279 us | 13.72 KiB | 245 | `slurm` | `c0705a-s11.ufhpc` |
| `coboundary_map/path/n500` | `assembly-large` | 557.642 us | 366.88 KiB | 5553 | `slurm` | `c0705a-s11.ufhpc` |
| `harmonic_extension/cycle/n100` | `extension-small` | 349.323 us | 431.70 KiB | 1861 | `slurm` | `c0705a-s11.ufhpc` |
| `harmonic_extension/cycle/n20` | `extension-small` | 90.629 us | 161.41 KiB | 556 | `slurm` | `c0705a-s11.ufhpc` |
| `harmonic_extension/cycle/n500` | `extension-large` | 1.673 ms | 1.78 MiB | 8339 | `slurm` | `c0705a-s11.ufhpc` |
| `harmonic_extension/path/n100` | `extension-small` | 352.529 us | 430.48 KiB | 1849 | `slurm` | `c0705a-s11.ufhpc` |
| `harmonic_extension/path/n20` | `extension-small` | 92.894 us | 160.31 KiB | 543 | `slurm` | `c0705a-s11.ufhpc` |
| `harmonic_extension/path/n500` | `extension-large` | 1.673 ms | 1.78 MiB | 8327 | `slurm` | `c0705a-s11.ufhpc` |
| `laplacian/matrix_direct/cycle/n100` | `assembly-small` | 126.898 us | 177.20 KiB | 1348 | `slurm` | `c0705a-s11.ufhpc` |
| `laplacian/matrix_direct/cycle/n20` | `assembly-small` | 32.280 us | 45.68 KiB | 301 | `slurm` | `c0705a-s11.ufhpc` |
| `laplacian/matrix_direct/cycle/n500` | `assembly-large` | 585.525 us | 921.38 KiB | 6569 | `slurm` | `c0705a-s11.ufhpc` |
| `laplacian/matrix_direct/path/n100` | `assembly-small` | 116.318 us | 176.22 KiB | 1335 | `slurm` | `c0705a-s11.ufhpc` |
| `laplacian/matrix_direct/path/n20` | `assembly-small` | 31.759 us | 44.79 KiB | 287 | `slurm` | `c0705a-s11.ufhpc` |
| `laplacian/matrix_direct/path/n500` | `assembly-large` | 590.233 us | 920.40 KiB | 6556 | `slurm` | `c0705a-s11.ufhpc` |
| `nearest_global_section/ldl/cycle/n100` | `solver-small` | 12.925 ms | 3.49 MiB | 1404 | `slurm` | `c0705a-s11.ufhpc` |
| `nearest_global_section/ldl/cycle/n20` | `solver-small` | 441.315 us | 163.09 KiB | 355 | `slurm` | `c0705a-s11.ufhpc` |
| `nearest_global_section/ldl/cycle/n500` | `solver-large` | 377.764 ms | 84.51 MiB | 6633 | `slurm` | `c0705a-s11.ufhpc` |
| `nearest_global_section/ldl/path/n100` | `solver-small` | 13.557 ms | 3.42 MiB | 1393 | `slurm` | `c0705a-s11.ufhpc` |
| `nearest_global_section/ldl/path/n20` | `solver-small` | 409.095 us | 148.84 KiB | 343 | `slurm` | `c0705a-s11.ufhpc` |
| `nearest_global_section/ldl/path/n500` | `solver-large` | 408.121 ms | 84.18 MiB | 6622 | `slurm` | `c0705a-s11.ufhpc` |
