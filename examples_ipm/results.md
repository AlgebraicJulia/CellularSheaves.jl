# IPM Benchmark Results

Generated: 2026-07-02

**Configuration:** Runs = 5, Warmup = 2, Tolerance = 1.0e-8

---

## Spline Benchmarks

### Nonnegative Spline LP

#### Mosek
```
Shape          raug splines degree   |V|   |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
nonneg        1e+05      50      6    50    49       3.8      12.5      16.7   3.32x   4.41x
nonneg        1e+05     100      6   100    99       7.4      24.0      27.0   3.24x   3.64x
monotone      1e+05      50      6   100    99       5.2      10.8      13.7   2.09x   2.65x
monotone      1e+05     100      6   200   199      10.3      19.8      23.6   1.93x   2.30x
convex        1e+05      50      6   100    99       6.6       9.7      12.0   1.47x   1.82x
convex        1e+05     100      6   200   199      12.6      19.2      23.3   1.52x   1.85x
```

#### Clarabel
| Shape    | raug  | splines | degree | \|V\| | \|E\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|----------|-------|---------|--------|-------|-------|---------|----------|------------|-------|-------|
| nonneg   | 1e+05 | 50      | 6      | 50    | 49    | 3.9     | 8.2      | 7.8        | 2.11x | 2.00x |
| nonneg   | 1e+05 | 100     | 6      | 100   | 99    | 7.7     | 16.5     | 14.4       | 2.16x | 1.88x |
| monotone | 1e+05 | 50      | 6      | 100   | 99    | 5.2     | 3.6      | 11.6       | 0.70x | 2.22x |
| monotone | 1e+05 | 100     | 6      | 200   | 199   | 10.4    | 7.2      | 22.4       | 0.69x | 2.16x |
| convex   | 1e+05 | 50      | 6      | 100   | 99    | 6.7     | 3.5      | 12.2       | 0.53x | 1.84x |
| convex   | 1e+05 | 100     | 6      | 200   | 199   | 13.0    | 7.3      | 24.5       | 0.56x | 1.89x |

---

### Nonnegative Spline Exact (SDP)

#### Mosek
```
Shape          raug splines degree   |V|   |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
nonneg        1e+05      50      6   150    99      24.6      32.7      39.5   1.33x   1.61x
nonneg        1e+05     100      6   300   199      58.9     102.8      69.6   1.75x   1.18x
monotone      1e+05      50      6   150    99      25.8     128.1      34.0   4.97x   1.32x
monotone      1e+05     100      6   300   199      49.8     735.1     208.0  14.76x   4.18x
convex        1e+05      50      6   150    99      19.4      33.5      25.6   1.72x   1.32x
convex        1e+05     100      6   300   199      46.0     118.2      52.0   2.57x   1.13x
```

#### Clarabel
| Shape    | raug  | splines | degree | \|V\| | \|E\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|----------|-------|---------|--------|-------|-------|---------|----------|------------|-------|-------|
| nonneg   | 1e+05 | 50      | 6      | 150   | 99    | 24.6    | 20.6     | 24.8       | 0.84x | 1.01x |
| nonneg   | 1e+05 | 100     | 6      | 300   | 199   | 60.6    | 43.3     | 44.7       | 0.71x | 0.74x |
| monotone | 1e+05 | 50      | 6      | 150   | 99    | 26.0    | 25.4     | 21.3       | 0.98x | 0.82x |
| monotone | 1e+05 | 100     | 6      | 300   | 199   | 55.4    | 54.5     | 42.5       | 0.98x | 0.77x |
| convex   | 1e+05 | 50      | 6      | 150   | 99    | 18.1    | 14.7     | 20.0       | 0.82x | 1.11x |
| convex   | 1e+05 | 100     | 6      | 300   | 199   | 42.3    | 34.2     | 38.3       | 0.81x | 0.90x |

---

### T-Junction Spline (Adaptive LP)

#### Mosek
```
Config         raug   DOF   |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
3x3 center    1e+05   300    12    81       7.1      20.2      30.8   2.85x   4.33x
4x4 diag      1e+05   550    22   171      23.9      49.2      54.1   2.06x   2.26x
4x4 block     1e+05   700    28   225      41.4      71.3      68.7   1.72x   1.66x
6x6 block     1e+05  1200    48   369      67.8     119.3     129.3   1.76x   1.91x
8x8 4x4blk    1e+05  2800   112   945     276.6     260.4     377.6   0.94x   1.37x
```

#### Clarabel
| Config     | raug  | DOF  | \|V\| | H1  | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|------------|-------|------|-------|-----|---------|----------|------------|-------|-------|
| 3x3 center | 1e+05 | 300  | 12    | 81  | 7.4     | 20.0     | 21.1       | 2.69x | 2.83x |
| 4x4 diag   | 1e+05 | 550  | 22    | 171 | 24.6    | 41.1     | 43.8       | 1.67x | 1.78x |
| 4x4 block  | 1e+05 | 700  | 28    | 225 | 43.7    | 64.1     | 62.4       | 1.47x | 1.43x |
| 6x6 block  | 1e+05 | 1200 | 48    | 369 | 70.7    | 123.6    | 119.6      | 1.75x | 1.69x |
| 8x8 4x4blk | 1e+05 | 2800 | 112   | 945 | 325.6   | 493.2    | 415.5      | 1.51x | 1.28x |

---

### T-Junction SDP

#### Mosek
```
Config           raug   DOF   |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
NN 2x2+1        1e+06   854    35    36      26.0      40.0      49.0   1.54x   1.88x
NN 3x3 ctr      1e+01  1464    60    81      68.2      67.8      76.6   0.99x   1.12x
NN 4x4 diag     1e+06  2684   110   171      91.3     119.8     133.5   1.31x   1.46x
CVX 2x2+1       1e+05  2716    35    36     210.8     159.8     215.1   0.76x   1.02x
CVX 3x3 ctr     1e+03  4656    60    81     416.3     226.5     296.3   0.54x   0.71x
CVX 4x4 diag    1e+04  8536   110   171     850.1     470.2     469.2   0.55x   0.55x
```

#### Clarabel
| Config      | raug  | DOF  | \|V\| | H1  | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|-------------|-------|------|-------|-----|---------|----------|------------|-------|-------|
| NN 2x2+1    | 1e+06 | 854  | 35    | 36  | 25.8    | 30.8     | 37.3       | 1.19x | 1.44x |
| NN 3x3 ctr  | 1e+01 | 1464 | 60    | 81  | 69.4    | 57.1     | 57.2       | 0.82x | 0.82x |
| NN 4x4 diag | 1e+06 | 2684 | 110   | 171 | 116.1   | 105.5    | 105.6      | 0.91x | 0.91x |
| CVX 2x2+1   | 1e+05 | 2716 | 35    | 36  | 287.0   | 628.0    | 453.2      | 2.19x | 1.58x |
| CVX 3x3 ctr | 1e+03 | 4656 | 60    | 81  | 577.5   | 1063.1   | 770.5      | 1.84x | 1.33x |
| CVX 4x4 diag| 1e+04 | 8536 | 110   | 171 | 925.8   | 1577.4   | 1661.8     | 1.70x | 1.80x |

---

### MLE Spline (Exp-Cone Leaves)

#### Mosek
```
Mode          raug    M   n      N leaves    |V|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
density      1e+05    8   4    200    200    208      48.5       8.2       9.3   0.17x   0.19x
density      1e+05   16   4   1000   1000   1016     627.6      34.1      38.8   0.05x   0.06x
density      1e+05   16   4   1000     58     74      12.8       7.5       5.9   0.59x   0.46x
intensity    1e+03    8   4    200    201    209      22.9       8.1       9.3   0.36x   0.41x
intensity    1e+04   16   4   1000    949    965     125.8      37.8      40.9   0.30x   0.33x
intensity    1e+01   16   4   1000     58     74       7.5       5.5       5.6   0.74x   0.75x
```

#### Clarabel
| Mode      | raug  | M  | n | N    | leaves | \|V\| | IPM(ms) | Clarabel | speedup |
|-----------|-------|----|---|------|--------|-------|---------|----------|---------|
| density   | 1e+05 | 8  | 4 | 200  | 200    | 208   | 48.2    | 7.0      | 0.14x   |
| density   | 1e+05 | 16 | 4 | 1000 | 1000   | 1016  | 617.6   | 26.1     | 0.04x   |
| density   | 1e+05 | 16 | 4 | 1000 | 58     | 74    | 12.5    | 3.1      | 0.25x   |
| intensity | 1e+03 | 8  | 4 | 200  | 201    | 209   | 22.5    | 7.3      | 0.33x   |
| intensity | 1e+04 | 16 | 4 | 1000 | 949    | 965   | 124.4   | 24.1     | 0.19x   |
| intensity | 1e+01 | 16 | 4 | 1000 | 58     | 74    | 7.3     | 2.9      | 0.40x   |

---

### Tensor Spline (2D Grid, Lossy Variants)

#### Mosek
```
Mode       Shape         raug   grid   |V|    H1     n   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
regress    nonneg       1e+05   4x4      16    81     4       9.0      21.5      35.9   2.40x   3.99x
regress    nonneg       1e+05   8x8      64   441     4      63.5     105.8     176.2   1.67x   2.78x
regress    nonneg       1e+05  16x16    256  2025     4     502.5     791.7    1497.9   1.58x   2.98x
regress    monotone_x   1e+05   8x8     128   441     4      91.0     107.3     195.5   1.18x   2.15x
regress    nonneg       1e+05   8x8      62   378     4      77.7      98.3     189.0   1.26x   2.43x
intensity  nonneg       1e+01   4x4     805    81     4     323.0      63.3     117.3   0.20x   0.36x
```

#### Clarabel
| Mode      | Shape      | raug  | grid  | \|V\| | H1   | n | IPM(ms) | Clarabel | speedup |
|-----------|------------|-------|-------|-------|------|---|---------|----------|---------|
| regress   | nonneg     | 1e+05 | 4x4   | 16    | 81   | 4 | 8.9     | 25.9     | 2.91x   |
| regress   | nonneg     | 1e+05 | 8x8   | 64    | 441  | 4 | 63.7    | 151.7    | 2.38x   |
| regress   | nonneg     | 1e+05 | 16x16 | 256   | 2025 | 4 | 474.8   | 1418.3   | 2.99x   |
| regress   | monotone_x | 1e+05 | 8x8   | 128   | 441  | 4 | 78.2    | 89.0     | 1.14x   |
| regress   | nonneg     | 1e+05 | 8x8   | 62    | 378  | 4 | 75.0    | 138.3    | 1.85x   |
| intensity | nonneg     | 1e+01 | 4x4   | 805   | 81   | 4 | 333.1   | 81.7     | 0.25x   |

---

## Obstacle Benchmarks

### Obstacle/American Option

#### Mosek
```
Config               raug   DOF  |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
perpetual 4-patch   1e-02   413    4   399      18.0       8.6      10.8   0.48x   0.60x
perpetual 801       1e-02   813    4   799      64.2      48.0      20.8   0.75x   0.32x
perpetual 8-patch   1e-02   829    8   799      36.0      50.6      20.7   1.40x   0.57x
perpetual 1601      1e-02  1657    8  1599     198.9      33.3      42.4   0.17x   0.21x
game δ=0.02         1e-02  2439   12   799     533.4      61.1      21.0   0.11x   0.04x

Time Stepping (200 QPs):
American put chain: IPM 7.14s   Mosek 3.24s   speedup 0.45x
```

#### Clarabel
| Config            | raug  | DOF  | \|V\| | H1  | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|-------------------|-------|------|-------|-----|---------|----------|------------|-------|-------|
| perpetual 4-patch | 1e-02 | 413  | 4     | 399 | 17.4    | 4.8      | 5.6        | 0.27x | 0.32x |
| perpetual 801     | 1e-02 | 813  | 4     | 799 | 55.4    | 9.9      | 10.4       | 0.18x | 0.19x |
| perpetual 8-patch | 1e-02 | 829  | 8     | 799 | 37.8    | 10.2     | 10.4       | 0.27x | 0.28x |
| perpetual 1601    | 1e-02 | 1657 | 8     | 1599| 214.1   | 21.2     | 19.3       | 0.10x | 0.09x |
| game delta=0.02   | 1e-02 | 2439 | 12    | 799 | 515.8   | 10.0     | 10.7       | 0.02x | 0.02x |

Time Stepping (200 obstacle QPs, same Q):
| Benchmark         | IPM     | Clarabel | speedup |
|-------------------|---------|----------|---------|
| American put chain| 7.62s   | 3.08s    | 0.40x   |

---

## Map Benchmarks

### SparseMAP LP

#### Mosek (rho=1.0)
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 50×50          1e+06    2500    4900     120.5     682.9     694.3   5.67x   5.76x
cycle               1e+06    4000    4000     114.5     527.0     570.2   4.60x   4.98x
complete            1e+06      95    4465      75.6     539.6     518.6   7.14x   6.86x
ladder              1e+06    2400    3598      85.7     416.5     440.8   4.86x   5.14x
triangular 40×40    1e+06    1600    4641     109.5     638.7     606.6   5.84x   5.54x
```

#### Clarabel (rho=1.0)
| Graph            | raug  | \|V\| | \|E\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|------------------|-------|-------|-------|---------|----------|------------|-------|-------|
| grid 50x50       | 1e+06 | 2500  | 4900  | 112.0   | 378.9    | 248.2      | 3.38x | 2.22x |
| cycle            | 1e+06 | 4000  | 4000  | 100.6   | 329.0    | 208.3      | 3.27x | 2.07x |
| complete         | 1e+06 | 95    | 4465  | 73.2    | 285.9    | 188.4      | 3.91x | 2.58x |
| ladder           | 1e+06 | 2400  | 3598  | 80.1    | 285.3    | 170.1      | 3.56x | 2.12x |
| triangular 40x40 | 1e+06 | 1600  | 4641  | 101.5   | 343.3    | 217.5      | 3.38x | 2.14x |

---

### Dense Stochastic Channels

#### Mosek (rho=1.0)
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 20×20          1e+06     400     760      88.8     244.5     267.7   2.75x   3.01x
cycle               1e+06     400     400     148.8     186.3     216.5   1.25x   1.46x
complete            1e+06      40     780      18.7      82.2      85.0   4.39x   4.54x
ladder              1e+06     400     598     156.7     212.1     226.0   1.35x   1.44x
triangular 15×15    1e+06     225     630      36.0     150.0     163.7   4.16x   4.54x
```

#### Clarabel (rho=1.0)
| Graph            | raug  | \|V\| | \|E\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM  | D/IPM |
|------------------|-------|-------|-------|---------|----------|------------|--------|-------|
| grid 20x20       | 1e+06 | 400   | 760   | 76.8    | 294.1    | 180.1      | 3.83x  | 2.35x |
| cycle            | 1e+06 | 400   | 400   | 151.8   | 140.9    | 59.4       | 0.93x  | 0.39x |
| complete         | 1e+06 | 40    | 780   | 18.3    | 211.8    | 115.7      | 11.55x | 6.31x |
| ladder           | 1e+06 | 400   | 598   | 159.3   | 206.6    | 75.0       | 1.30x  | 0.47x |
| triangular 15x15 | 1e+06 | 225   | 630   | 35.3    | 205.8    | 124.3      | 5.83x  | 3.52x |

---

### Moment SDP

#### Mosek (rho=1.0)
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 20×20          1e+06     400     760     115.4     533.5     153.1   4.62x   1.33x
cycle               1e+06     800     800      61.6     370.5     289.0   6.01x   4.69x
complete            1e+06      12      66      13.2     327.8      21.9  24.79x   1.66x
ladder              1e+06     500     748      55.2     620.0     303.9  11.24x   5.51x
```

#### Clarabel (rho=1.0)
| Graph      | raug  | \|V\| | \|E\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|------------|-------|-------|-------|---------|----------|------------|-------|-------|
| grid 20x20 | 1e+06 | 400   | 760   | 90.6    | 86.6     | 73.6       | 0.96x | 0.81x |
| cycle      | 1e+06 | 800   | 800   | 55.2    | 50.9     | 46.6       | 0.92x | 0.84x |
| complete   | 1e+06 | 12    | 66    | 10.9    | 15.2     | 13.6       | 1.40x | 1.25x |
| ladder     | 1e+06 | 500   | 748   | 66.0    | 57.1     | 53.4       | 0.87x | 0.81x |

---

### Unital CPTP Channels

#### Mosek (rho=1.0)
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 20×20          1e+06     400     760      59.6     323.0     200.0   5.42x   3.36x
cycle               1e+06     400     400      42.1     340.7     140.1   8.09x   3.33x
complete            1e+06      18     153       3.0      84.0      25.7  27.99x   8.57x
ladder              1e+06     400     598      44.3     227.8     165.2   5.14x   3.73x
```

#### Clarabel (rho=1.0)
| Graph      | raug  | \|V\| | \|E\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|------------|-------|-------|-------|---------|----------|------------|-------|-------|
| grid 20x20 | 1e+06 | 400   | 760   | 65.5    | 73.7     | 71.2       | 1.13x | 1.09x |
| cycle      | 1e+06 | 400   | 400   | 40.7    | 39.3     | 39.4       | 0.97x | 0.97x |
| complete   | 1e+06 | 18    | 153   | 3.0     | 7.4      | 9.0        | 2.49x | 3.04x |
| ladder     | 1e+06 | 400   | 598   | 43.9    | 47.2     | 44.3       | 1.07x | 1.01x |

---

### Wasserstein Gaussian

#### Mosek
```
Config           raug   DOF   |V|   IPM(ms)   Mosek     Mosek-D    P/IPM   D/IPM
transport d=3     1e+01   21     1       0.8       0.9       2.9   1.07x   3.63x
Higham C4 d=5     1e+01  220     4      15.4      25.0      65.1   1.63x   4.24x
Higham C4 d=8     1e+01  544     4      71.5     107.6     239.3   1.50x   3.35x
Higham grid 6x6   1e+01  696    92      36.4      49.3     100.3   1.36x   2.75x
```

#### Clarabel
| Config          | raug  | DOF | \|V\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|-----------------|-------|-----|-------|---------|----------|------------|-------|-------|
| transport d=3   | 1e+01 | 21  | 1     | 0.8     | 0.5      | 2.5        | 0.68x | 3.39x |
| Higham C4 d=5   | 1e+01 | 220 | 4     | 15.7    | 12.6     | 26.5       | 0.81x | 1.69x |
| Higham C4 d=8   | 1e+01 | 544 | 4     | 70.0    | 91.0     | 297.8      | 1.30x | 4.26x |
| Higham grid 6x6 | 1e+01 | 696 | 92    | 36.2    | 17.2     | 25.3       | 0.48x | 0.70x |

---

### Wasserstein Graph

#### Mosek
```
Config                 raug   DOF   |V|    H1   IPM(ms)     Mosek   Mosek-D   P/IPM   D/IPM
tri LP                1e+05    12     3     3       0.1       1.8       2.8  32.99x  51.70x
4-cyc LP              1e+04   900     4     4      19.6       9.9      12.4   0.50x   0.63x
smooth bary           1e+05  1220     4     2      54.4      29.2      18.0   0.54x   0.33x
quadreg 25 (lifted)   1e+03   625   625     1      37.9       9.4       8.8   0.25x   0.23x
laplace 15 (lifted)   1e+02  1215  1215     1      14.2      17.9      21.3   1.26x   1.50x
laplace 20 (lifted)   1e+02  2300  2300     1      38.7      46.9      50.9   1.21x   1.31x
```

#### Clarabel
| Config              | raug  | DOF  | \|V\| | H1 | IPM(ms) | Clarabel | Clarabel-D | P/IPM  | D/IPM  |
|---------------------|-------|------|-------|----|---------|----------|------------|--------|--------|
| tri LP              | 1e+05 | 12   | 3     | 3  | 0.1     | 0.3      | 1.6        | 5.04x  | 31.28x |
| 4-cyc LP            | 1e+04 | 900  | 4     | 4  | 21.7    | 10.0     | 9.6        | 0.46x  | 0.44x  |
| smooth bary         | 1e+05 | 1220 | 4     | 2  | 52.8    | 14.6     | 13.0       | 0.28x  | 0.25x  |
| quadreg 25 (lifted) | 1e+03 | 625  | 625   | 1  | 36.2    | 7.8      | 7.1        | 0.21x  | 0.20x  |
| laplace 15 (lifted) | 1e+02 | 1215 | 1215  | 1  | 14.2    | 13.5     | 14.9       | 0.95x  | 1.04x  |
| laplace 20 (lifted) | 1e+02 | 2300 | 2300  | 1  | 34.9    | 40.4     | 39.2       | 1.16x  | 1.12x  |

---

### Quantum Marginal (SDP)

#### Mosek
```
Config           raug   DOF   |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
N6  ℓ2 ring     1e+03    60     6     1       0.6       3.4       4.2   5.41x   6.75x
N6  ℓ3 ring     1e+06   216     6     1       2.5      16.8      16.0   6.63x   6.32x
N10 ℓ2 ring     1e+01   100    10     1       1.7       4.8       5.8   2.83x   3.43x
N10 ℓ3 ring     1e+06   360    10     1       4.2      23.4      17.0   5.58x   4.06x
N10 ℓ4 ring     1e+03  1360    10     1      63.1      75.0      61.3   1.19x   0.97x
N12 ℓ3 ring     1e+05   432    12     1       5.3      27.1      36.5   5.13x   6.92x
N12 ℓ4 ring     1e+04  1632    12     1      77.6      92.0      76.2   1.19x   0.98x
```

#### Clarabel
| Config       | raug  | DOF  | \|V\| | H1 | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|--------------|-------|------|-------|----|---------|----------|------------|-------|-------|
| N6 l2 ring   | 1e+03 | 60   | 6     | 1  | 0.6     | 0.6      | 2.5        | 1.04x | 4.01x |
| N6 l3 ring   | 1e+06 | 216  | 6     | 1  | 2.5     | 3.5      | 6.7        | 1.38x | 2.64x |
| N10 l2 ring  | 1e+01 | 100  | 10    | 1  | 1.6     | 1.1      | 3.3        | 0.68x | 2.09x |
| N10 l3 ring  | 1e+06 | 360  | 10    | 1  | 4.2     | 6.0      | 11.8       | 1.45x | 2.83x |
| N10 l4 ring  | 1e+03 | 1360 | 10    | 1  | 56.2    | 183.7    | 173.4      | 3.27x | 3.08x |
| N12 l3 ring  | 1e+05 | 432  | 12    | 1  | 5.2     | 7.3      | 13.8       | 1.40x | 2.65x |
| N12 l4 ring  | 1e+04 | 1632 | 12    | 1  | 83.9    | 228.7    | 179.8      | 2.73x | 2.14x |

---

## NUM Benchmarks

### Network Utility Maximization

#### Mosek
```
Config               raug   DOF  |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
kelly L=3 tput      1e+03     7    7     0       0.1         —         —       —       —
kelly L=3 prop      1e+01    19   11     0       0.6       1.0       1.3   1.73x   2.40x
kelly L=3 delay     1e+02    23   15     0       0.2         —         —       —       —
kelly L=3 maxmin    1e+03    12   12     0       0.1         —         —       —       —
kelly L=10 prop     1e+01    54   32     0       2.8       1.0       1.6   0.37x   0.58x
kelly L=20 prop     1e+01   104   62     0       7.7       1.4       2.0   0.18x   0.25x
mesh 20×10 prop     1e+02    90   50     0       5.8       1.4       1.8   0.24x   0.32x
mesh 50×25 prop     1e+02   225  125     0       7.7       1.9       3.7   0.24x   0.48x
mesh 100×50 prop    1e+01   450  250     0      22.7       3.7       7.6   0.16x   0.34x
```

#### Clarabel
| Config          | raug  | DOF | \|V\| | H1 | IPM(ms) | Clarabel | speedup |
|-----------------|-------|-----|-------|----|---------|----------|---------|
| kelly L=3 tput  | 1e+03 | 7   | 7     | 0  | 0.1     | -        | -       |
| kelly L=3 prop  | 1e+01 | 19  | 11    | 0  | 0.6     | 0.3      | 0.52x   |
| kelly L=3 delay | 1e+02 | 23  | 15    | 0  | 0.2     | -        | -       |
| kelly L=3 maxmin| 1e+03 | 12  | 12    | 0  | 0.1     | -        | -       |
| kelly L=10 prop | 1e+01 | 54  | 32    | 0  | 2.8     | 0.5      | 0.18x   |
| kelly L=20 prop | 1e+01 | 104 | 62    | 0  | 7.7     | 0.8      | 0.11x   |
| mesh 20x10 prop | 1e+02 | 90  | 50    | 0  | 6.1     | 0.8      | 0.13x   |
| mesh 50x25 prop | 1e+02 | 225 | 125   | 0  | 7.7     | 1.9      | 0.25x   |
| mesh 100x50 prop| 1e+01 | 450 | 250   | 0  | 23.1    | 3.9      | 0.17x   |

---

## OPF Benchmarks

### Basic OPF (AC-OPF with SOCP + Cycle Cells)

#### Mosek
```
Config           raug   DOF   |V|   IPM(ms)   Mosek     Mosek-D    P/IPM   D/IPM
tri linear        1e+05   27    17       0.3       0.7       0.8   2.09x   2.62x
tri quad          1e+05   27    17       0.3       0.6       1.2   1.76x   3.41x
C4 linear         1e+05   34    21       0.4       0.8       1.1   1.85x   2.52x
C4 quad           1e+05   34    21       0.7       0.7       1.1   1.00x   1.51x
radial linear     1e+05   23    16       0.2       0.5       1.0   2.18x   4.24x
radial quad       1e+05   23    16       0.4       0.6       1.0   1.57x   2.46x
```

#### Clarabel
| Config        | raug  | DOF | \|V\| | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|---------------|-------|-----|-------|---------|----------|------------|-------|-------|
| tri linear    | 1e+05 | 27  | 17    | 0.3     | 0.5      | 2.2        | 1.72x | 8.25x |
| tri quad      | 1e+05 | 27  | 17    | 0.3     | 0.4      | 2.1        | 1.21x | 6.55x |
| C4 linear     | 1e+05 | 34  | 21    | 0.4     | 0.4      | 2.1        | 1.09x | 5.30x |
| C4 quad       | 1e+05 | 34  | 21    | 0.6     | 0.5      | 2.1        | 0.71x | 3.26x |
| radial linear | 1e+05 | 23  | 16    | 0.2     | 0.4      | 1.8        | 1.54x | 7.95x |
| radial quad   | 1e+05 | 23  | 16    | 0.4     | 0.4      | 1.9        | 0.96x | 4.99x |

---

### Extended OPF Benchmark

#### Mosek
```
Config                raug   DOF   IPM(ms)  Mosek     Mosek-D    P/IPM   D/IPM
tri linear            1e+05   27     0.29     0.63     0.81    2.19x    2.84x
tri quad              1e+05   27     0.37     0.52     0.85    1.40x    2.31x
C4 quad               1e+05   34     0.71     0.57     0.89    0.81x    1.25x
tri +cycle            1e+05   48     5.57     1.51     1.87    0.27x    0.34x
C4 +cycle             1e+05   70    11.29     5.00     5.83    0.44x    0.52x
bowtie 2-tri          1e+05   80     5.09     7.00     3.37    1.38x    0.66x
bowtie full           1e+05   74     3.07     6.78     5.56    2.21x    1.81x
SE 3x3 r1             1e+03   57     1.28     2.15     3.52    1.68x    2.75x
SE 5x5 r1             1e+03  185     3.25     4.87     8.93    1.50x    2.75x
SE 7x7 r1             1e+03  385     6.98    11.53    21.23    1.65x    3.04x
SE triangle r2        1e+02   36     2.11     4.83     6.33    2.28x    3.00x
SE c4 r2              1e+02   56     8.09    19.65    23.60    2.43x    2.92x
```

#### Clarabel

**Rung 1 (SOCP):**
| Config     | raug  | DOF | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|------------|-------|-----|---------|----------|------------|-------|-------|
| tri linear | 1e+05 | 27  | 0.28    | 0.46     | 2.11       | 1.67x | 7.60x |
| tri quad   | 1e+05 | 27  | 0.33    | 0.36     | 1.95       | 1.11x | 5.92x |
| C4 quad    | 1e+05 | 34  | 0.63    | 0.44     | 2.04       | 0.70x | 3.22x |

**Rung 2 (SOCP + SDP):**
| Config    | raug  | DOF | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|-----------|-------|-----|---------|----------|------------|-------|-------|
| tri +cycle| 1e+05 | 48  | 5.42    | 0.92     | 1.37       | 0.17x | 0.25x |
| C4 +cycle | 1e+05 | 70  | 11.06   | 1.40     | 2.00       | 0.13x | 0.18x |

**Bowtie (Chordal):**
| Config       | raug  | DOF | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|--------------|-------|-----|---------|----------|------------|-------|-------|
| bowtie 2-tri | 1e+05 | 80  | 9.98    | 1.71     | 2.16       | 0.17x | 0.22x |
| bowtie full  | 1e+05 | 74  | 3.27    | 1.66     | 2.26       | 0.51x | 0.69x |

**SE Rung 1 (SOCP + Dense Q):**
| Config    | raug  | DOF | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|-----------|-------|-----|---------|----------|------------|-------|-------|
| SE 3x3 r1 | 1e+03 | 57  | 1.22    | 0.47     | 2.25       | 0.38x | 1.85x |
| SE 5x5 r1 | 1e+03 | 185 | 3.13    | 1.33     | 3.81       | 0.42x | 1.22x |
| SE 7x7 r1 | 1e+03 | 385 | 6.66    | 2.80     | 6.23       | 0.42x | 0.94x |

**SE Rung 2 (SOCP + SDP + Dense Q):**
| Config        | raug  | DOF | IPM(ms) | Clarabel | Clarabel-D | P/IPM | D/IPM |
|---------------|-------|-----|---------|----------|------------|-------|-------|
| SE triangle r2| 1e+02 | 36  | 2.05    | 0.82     | 1.10       | 0.40x | 0.53x |
| SE c4 r2      | 1e+02 | 56  | 12.69   | 1.43     | 1.64       | 0.11x | 0.13x |

---

## Legend

- **raug**: Augmentation parameter
- **DOF**: Degrees of freedom
- **|V|**: Number of vertices in the sheaf graph
- **|E|**: Number of edges in the sheaf graph
- **H1**: First cohomology dimension
- **IPM(ms)**: Sheaf IPM solve time in milliseconds
- **Mosek / Clarabel**: Primal solver time in milliseconds
- **Mosek-D / Clarabel-D**: Dualized solver time in milliseconds
- **P/IPM**: Primal solver time / IPM time (speedup ratio)
- **D/IPM**: Dualized solver time / IPM time (speedup ratio)

Values > 1.0x indicate Sheaf IPM is faster than the reference solver.

---

## Summary

### IPM wins vs Mosek (P/IPM > 1)

| Benchmark | Best Speedup | Notes |
|-----------|--------------|-------|
| map_lp SparseMAP | 5-7x | Grid, complete, triangular graphs |
| map_sdp complete CPTP | 28x | Small complete graph |
| map_sdp Moment SDP | 6-25x | Complete graph best |
| quantum_marginal | 2-7x | Smaller locality best |
| nonnegative_spline_exact monotone | 15x | 100 splines |
| tjunction_spline | 2-4x | Smaller configs |
| tensor_spline regress | 2-4x | All grid sizes |
| opf_sheaf SE | 1.5-3x | State estimation problems |

### IPM wins vs Clarabel (P/IPM > 1)

| Benchmark | Best Speedup | Notes |
|-----------|--------------|-------|
| map_lp SparseMAP | 3-4x | All graph types |
| map_lp Dense complete | 11x | Complete graph |
| tjunction_spline | 1.5-2.7x | All configs |
| tjunction_sdp CVX | 1.7-2.2x | Convex constraints |
| tensor_spline regress | 2-3x | All grid sizes |
| quantum_marginal ℓ4 | 2.7-3.3x | High locality |

### Reference solvers win

| Benchmark | Speedup | Notes |
|-----------|---------|-------|
| obstacle_option | Clarabel 5-50x faster | Sparse Q exploitation |
| num_sheaf large mesh | Both 3-10x faster | Pure hyperedge coupling |
| mle_spline density | Both 4-25x faster | Many exp-cone leaves |
| lifted diagonal Q | Clarabel 2-5x faster | No dense Q advantage |
