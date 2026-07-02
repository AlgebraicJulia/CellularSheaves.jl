# IPM Benchmark Results

Generated: 2026-07-02

## wasserstein_graph.jl

### Mosek
```
Config                 raug   DOF   |V|    H1   IPM(ms)     Mosek   Mosek-D   P/IPM   D/IPM
tri LP                1e+05    12     3     3       0.1       1.8       2.8  32.99x  51.70x
4-cyc LP              1e+04   900     4     4      19.6       9.9      12.4   0.50x   0.63x
smooth bary           1e+05  1220     4     2      54.4      29.2      18.0   0.54x   0.33x
quadreg 25 (lifted)   1e+03   625   625     1      37.9       9.4       8.8   0.25x   0.23x
laplace 15 (lifted)   1e+02  1215  1215     1      14.2      17.9      21.3   1.26x   1.50x
laplace 20 (lifted)   1e+02  2300  2300     1      38.7      46.9      50.9   1.21x   1.31x
```

### OSQP
```
Config                 raug   DOF   |V|    H1   IPM(ms)      OSQP    OSQP-D   P/IPM   D/IPM
tri LP                1e+05    12     3     3       0.1       2.0       2.2  38.76x  41.69x
4-cyc LP              1e+04   900     4     4      24.0       8.6      10.3   0.36x   0.43x
smooth bary           1e+05  1220     4     2      51.4      11.1      15.4   0.22x   0.30x
quadreg 25 (lifted)   1e+03   625   625     1      39.9       6.5       6.2   0.16x   0.16x
laplace 15 (lifted)   1e+02  1215  1215     1      14.3       9.4      12.6   0.66x   0.88x
laplace 20 (lifted)   1e+02  2300  2300     1      36.1      25.0      36.1   0.69x   1.00x
```

---

## map_lp.jl

### Mosek

**SparseMAP (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 50×50          1e+06    2500    4900     120.5     682.9     694.3   5.67x   5.76x
cycle               1e+06    4000    4000     114.5     527.0     570.2   4.60x   4.98x
complete            1e+06      95    4465      75.6     539.6     518.6   7.14x   6.86x
ladder              1e+06    2400    3598      85.7     416.5     440.8   4.86x   5.14x
triangular 40×40    1e+06    1600    4641     109.5     638.7     606.6   5.84x   5.54x
```

**Dense Stochastic Channels (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 20×20          1e+06     400     760      88.8     244.5     267.7   2.75x   3.01x
cycle               1e+06     400     400     148.8     186.3     216.5   1.25x   1.46x
complete            1e+06      40     780      18.7      82.2      85.0   4.39x   4.54x
ladder              1e+06     400     598     156.7     212.1     226.0   1.35x   1.44x
triangular 15×15    1e+06     225     630      36.0     150.0     163.7   4.16x   4.54x
```

### OSQP

**SparseMAP (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) OSQP      OSQP-D      P/IPM   D/IPM
grid 50×50          1e+06    2500    4900     105.2     609.1     471.5   5.79x   4.48x
cycle               1e+06    4000    4000      97.9     340.1     318.5   3.47x   3.25x
complete            1e+06      95    4465      74.6     306.5     413.8   4.11x   5.55x
ladder              1e+06    2400    3598      93.0     274.9     297.1   2.96x   3.20x
triangular 40×40    1e+06    1600    4641      97.9     370.1     336.9   3.78x   3.44x
```

**Dense Stochastic Channels (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) OSQP      OSQP-D      P/IPM   D/IPM
grid 20×20          1e+06     400     760      78.1     179.5      85.3   2.30x   1.09x
cycle               1e+06     400     400     147.7     121.6      56.8   0.82x   0.38x
complete            1e+06      40     780      18.0      74.3      55.6   4.12x   3.08x
ladder              1e+06     400     598     155.4     124.1      80.2   0.80x   0.52x
triangular 15×15    1e+06     225     630      36.2     115.4      59.9   3.19x   1.66x
```

---

## map_sdp.jl

### Mosek

**Moment SDP (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 20×20          1e+06     400     760     115.4     533.5     153.1   4.62x   1.33x
cycle               1e+06     800     800      61.6     370.5     289.0   6.01x   4.69x
complete            1e+06      12      66      13.2     327.8      21.9  24.79x   1.66x
ladder              1e+06     500     748      55.2     620.0     303.9  11.24x   5.51x
```

**Unital CPTP Channels (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
grid 20×20          1e+06     400     760      59.6     323.0     200.0   5.42x   3.36x
cycle               1e+06     400     400      42.1     340.7     140.1   8.09x   3.33x
complete            1e+06      18     153       3.0      84.0      25.7  27.99x   8.57x
ladder              1e+06     400     598      44.3     227.8     165.2   5.14x   3.73x
```

### Clarabel

**Moment SDP (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) Clarabel  Clarabel-D   P/IPM   D/IPM
grid 20×20          1e+06     400     760     107.9     164.9     168.9   1.53x   1.57x
cycle               1e+06     800     800      56.3     106.2     120.0   1.89x   2.13x
complete            1e+06      12      66      13.5      24.3      17.1   1.80x   1.27x
ladder              1e+06     500     748      66.3     129.6     135.2   1.96x   2.04x
```

**Unital CPTP Channels (rho=1.0):**
```
Graph                raug     |V|     |E|   IPM(ms) Clarabel  Clarabel-D   P/IPM   D/IPM
grid 20×20          1e+06     400     760      60.1     199.8     191.6   3.33x   3.19x
cycle               1e+06     400     400      43.8     115.6     110.9   2.64x   2.53x
complete            1e+06      18     153       3.4      28.8      29.9   8.38x   8.68x
ladder              1e+06     400     598      45.4     147.4     149.3   3.25x   3.29x
```

---

## quantum_marginal.jl

### Mosek
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

### Clarabel
```
Config           raug   DOF   |V|    H1   IPM(ms) Clarabel  Clarabel-D   P/IPM   D/IPM
N6  ℓ2 ring     1e+03    60     6     1       0.6       2.7       3.3   4.27x   5.28x
N6  ℓ3 ring     1e+06   216     6     1       2.5       6.9       7.9   2.76x   3.19x
N10 ℓ2 ring     1e+01   100    10     1       1.5       3.5       4.5   2.31x   3.00x
N10 ℓ3 ring     1e+06   360    10     1       4.4       9.7      13.6   2.22x   3.11x
N10 ℓ4 ring     1e+03  1360    10     1      69.2     190.8     176.5   2.76x   2.55x
N12 ℓ3 ring     1e+05   432    12     1       5.3      11.3      16.1   2.14x   3.06x
N12 ℓ4 ring     1e+04  1632    12     1      80.3     233.2     182.0   2.91x   2.27x
```

---

## wasserstein_gauss.jl

### Mosek
```
Config           raug   DOF   |V|   IPM(ms)   Mosek     Mosek-D    P/IPM   D/IPM
transport d=3     1e+01   21     1       0.8       0.9       2.9   1.07x   3.63x
Higham C4 d=5     1e+01  220     4      15.4      25.0      65.1   1.63x   4.24x
Higham C4 d=8     1e+01  544     4      71.5     107.6     239.3   1.50x   3.35x
Higham grid 6x6   1e+01  696    92      36.4      49.3     100.3   1.36x   2.75x
```

### Clarabel
```
Config           raug   DOF   |V|   IPM(ms)   Clarabel  Clarabel-D  P/IPM   D/IPM
transport d=3     1e+01   21     1       0.8       0.6       2.6   0.76x   3.50x
Higham C4 d=5     1e+01  220     4      13.8      13.3      25.1   0.96x   1.82x
Higham C4 d=8     1e+01  544     4      56.3      88.1     323.5   1.57x   5.75x
Higham grid 6x6   1e+01  696    92      33.9      17.1      17.2   0.51x   0.51x
```

---

## num_sheaf.jl

### Mosek
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

### Clarabel
```
Config               raug   DOF  |V|    H1   IPM(ms) Clarabel   speedup
kelly L=3 tput      1e+03     7    7     0       0.1         —       —
kelly L=3 prop      1e+01    19   11     0       0.6       2.0    3.43x
kelly L=3 delay     1e+02    23   15     0       0.3         —       —
kelly L=3 maxmin    1e+03    12   12     0       0.1         —       —
kelly L=10 prop     1e+01    54   32     0       2.9       2.2    0.77x
kelly L=20 prop     1e+01   104   62     0       7.9       2.6    0.32x
mesh 20×10 prop     1e+02    90   50     0       5.9       2.4    0.41x
mesh 50×25 prop     1e+02   225  125     0       7.7       3.8    0.49x
mesh 100×50 prop    1e+01   450  250     0      23.0       5.7    0.25x
```

---

## obstacle_option.jl

### Mosek
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

### OSQP
```
Config               raug   DOF  |V|    H1   IPM(ms) OSQP      OSQP-D      P/IPM   D/IPM
perpetual 4-patch   1e-02   413    4   399      17.5       5.2       5.2   0.30x   0.30x
perpetual 801       1e-02   813    4   799      53.3       8.8      10.9   0.16x   0.20x
perpetual 8-patch   1e-02   829    8   799      49.9       7.7       8.9   0.15x   0.18x
perpetual 1601      1e-02  1657    8  1599     206.0      16.7      18.6   0.08x   0.09x
game δ=0.02         1e-02  2439   12   799     524.2       8.9      11.0   0.02x   0.02x

Time Stepping (200 QPs):
American put chain: IPM 8.06s   OSQP 1.83s   speedup 0.23x
```

---

## nonnegative_spline_lp.jl

### Mosek
```
Shape          raug splines degree   |V|   |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
nonneg        1e+05      50      6    50    49       3.8      12.5      16.7   3.32x   4.41x
nonneg        1e+05     100      6   100    99       7.4      24.0      27.0   3.24x   3.64x
monotone      1e+05      50      6   100    99       5.2      10.8      13.7   2.09x   2.65x
monotone      1e+05     100      6   200   199      10.3      19.8      23.6   1.93x   2.30x
convex        1e+05      50      6   100    99       6.6       9.7      12.0   1.47x   1.82x
convex        1e+05     100      6   200   199      12.6      19.2      23.3   1.52x   1.85x
```

### OSQP
```
Shape          raug splines degree   |V|   |E|   IPM(ms) OSQP      OSQP-D      P/IPM   D/IPM
nonneg        1e+05      50      6    50    49       3.8       6.9       8.2   1.84x   2.17x
nonneg        1e+05     100      6   100    99       7.4      13.5      17.0   1.81x   2.29x
monotone      1e+05      50      6   100    99       5.1       4.2       8.7   0.82x   1.70x
monotone      1e+05     100      6   200   199       9.9       8.9      16.5   0.90x   1.66x
convex        1e+05      50      6   100    99       6.4       4.6       8.0   0.71x   1.24x
convex        1e+05     100      6   200   199      12.7      14.1      15.8   1.11x   1.24x
```

---

## nonnegative_spline_exact.jl

### Mosek
```
Shape          raug splines degree   |V|   |E|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
nonneg        1e+05      50      6   150    99      24.6      32.7      39.5   1.33x   1.61x
nonneg        1e+05     100      6   300   199      58.9     102.8      69.6   1.75x   1.18x
monotone      1e+05      50      6   150    99      25.8     128.1      34.0   4.97x   1.32x
monotone      1e+05     100      6   300   199      49.8     735.1     208.0  14.76x   4.18x
convex        1e+05      50      6   150    99      19.4      33.5      25.6   1.72x   1.32x
convex        1e+05     100      6   300   199      46.0     118.2      52.0   2.57x   1.13x
```

### Clarabel
```
Shape          raug splines degree   |V|   |E|   IPM(ms) Clarabel  Clarabel-D   P/IPM   D/IPM
nonneg        1e+05      50      6   150    99      23.9      25.8      27.4   1.08x   1.14x
nonneg        1e+05     100      6   300   199      59.6      54.5      50.9   0.91x   0.85x
monotone      1e+05      50      6   150    99      29.2      30.8      25.9   1.05x   0.88x
monotone      1e+05     100      6   300   199      55.0      61.4      48.3   1.12x   0.88x
convex        1e+05      50      6   150    99      19.6      19.4      22.5   0.99x   1.15x
convex        1e+05     100      6   300   199      44.7      42.3      46.9   0.95x   1.05x
```

---

## tjunction_spline.jl

### Mosek
```
Config         raug   DOF   |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
3x3 center    1e+05   300    12    81       7.1      20.2      30.8   2.85x   4.33x
4x4 diag      1e+05   550    22   171      23.9      49.2      54.1   2.06x   2.26x
4x4 block     1e+05   700    28   225      41.4      71.3      68.7   1.72x   1.66x
6x6 block     1e+05  1200    48   369      67.8     119.3     129.3   1.76x   1.91x
8x8 4x4blk    1e+05  2800   112   945     276.6     260.4     377.6   0.94x   1.37x
```

### OSQP
```
Config         raug   DOF   |V|    H1   IPM(ms) OSQP      OSQP-D      P/IPM   D/IPM
3x3 center    1e+05   300    12    81       7.2      10.8      30.8   1.51x   4.30x
4x4 diag      1e+05   550    22   171      25.1      19.4      46.1   0.77x   1.84x
4x4 block     1e+05   700    28   225      43.3      26.2      58.6   0.61x   1.35x
6x6 block     1e+05  1200    48   369      66.4      53.3      79.1   0.80x   1.19x
8x8 4x4blk    1e+05  2800   112   945     301.1     221.9     486.7   0.74x   1.62x
```

---

## tjunction_sdp.jl

### Mosek
```
Config           raug   DOF   |V|    H1   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
NN 2x2+1        1e+06   854    35    36      26.0      40.0      49.0   1.54x   1.88x
NN 3x3 ctr      1e+01  1464    60    81      68.2      67.8      76.6   0.99x   1.12x
NN 4x4 diag     1e+06  2684   110   171      91.3     119.8     133.5   1.31x   1.46x
CVX 2x2+1       1e+05  2716    35    36     210.8     159.8     215.1   0.76x   1.02x
CVX 3x3 ctr     1e+03  4656    60    81     416.3     226.5     296.3   0.54x   0.71x
CVX 4x4 diag    1e+04  8536   110   171     850.1     470.2     469.2   0.55x   0.55x
```

### Clarabel
```
Config           raug   DOF   |V|    H1   IPM(ms) Clarabel  Clarabel-D   P/IPM   D/IPM
NN 2x2+1        1e+06   854    35    36      26.6      36.5      39.0   1.37x   1.47x
NN 3x3 ctr      1e+01  1464    60    81      66.2      67.5      61.1   1.02x   0.92x
NN 4x4 diag     1e+06  2684   110   171      89.2     114.7     116.2   1.29x   1.30x
CVX 2x2+1       1e+05  2716    35    36     202.4     626.6     452.4   3.10x   2.24x
CVX 3x3 ctr     1e+03  4656    60    81     401.7    1053.4     766.2   2.62x   1.91x
CVX 4x4 diag    1e+04  8536   110   171     816.4    1610.2    1720.0   1.97x   2.11x
```

---

## mle_spline.jl

### Mosek
```
Mode          raug    M   n      N leaves    |V|   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
density      1e+05    8   4    200    200    208      48.5       8.2       9.3   0.17x   0.19x
density      1e+05   16   4   1000   1000   1016     627.6      34.1      38.8   0.05x   0.06x
density      1e+05   16   4   1000     58     74      12.8       7.5       5.9   0.59x   0.46x
intensity    1e+03    8   4    200    201    209      22.9       8.1       9.3   0.36x   0.41x
intensity    1e+04   16   4   1000    949    965     125.8      37.8      40.9   0.30x   0.33x
intensity    1e+01   16   4   1000     58     74       7.5       5.5       5.6   0.74x   0.75x
```

### Clarabel
```
Mode          raug    M   n      N leaves    |V|   IPM(ms) Clarabel   speedup
density      1e+05    8   4    200    200    208      48.2      10.6    0.22x
density      1e+05   16   4   1000   1000   1016     618.9      44.3    0.07x
density      1e+05   16   4   1000     58     74      12.6       6.0    0.48x
intensity    1e+03    8   4    200    201    209      22.6      11.1    0.49x
intensity    1e+04   16   4   1000    949    965     123.9      38.7    0.31x
intensity    1e+01   16   4   1000     58     74       7.3       6.2    0.85x
```

---

## tensor_spline.jl

### Mosek
```
Mode       Shape         raug   grid   |V|    H1     n   IPM(ms) Mosek     Mosek-D     P/IPM   D/IPM
regress    nonneg       1e+05   4x4      16    81     4       9.0      21.5      35.9   2.40x   3.99x
regress    nonneg       1e+05   8x8      64   441     4      63.5     105.8     176.2   1.67x   2.78x
regress    nonneg       1e+05  16x16    256  2025     4     502.5     791.7    1497.9   1.58x   2.98x
regress    monotone_x   1e+05   8x8     128   441     4      91.0     107.3     195.5   1.18x   2.15x
regress    nonneg       1e+05   8x8      62   378     4      77.7      98.3     189.0   1.26x   2.43x
intensity  nonneg       1e+01   4x4     805    81     4     323.0      63.3     117.3   0.20x   0.36x
```

### Clarabel
```
Mode       Shape         raug   grid   |V|    H1     n   IPM(ms) Clarabel   speedup
regress    nonneg       1e+05   4x4      16    81     4       8.9      29.7    3.36x
regress    nonneg       1e+05   8x8      64   441     4      65.2     171.1    2.62x
regress    nonneg       1e+05  16x16    256  2025     4     483.2    1888.4    3.91x
regress    monotone_x   1e+05   8x8     128   441     4      75.9     111.3    1.47x
regress    nonneg       1e+05   8x8      62   378     4      75.7     152.6    2.01x
intensity  nonneg       1e+01   4x4     805    81     4     330.7      88.0    0.27x
```

---

## opf_sheaf.jl

### Mosek

**Basic OPF:**
```
Config           raug   DOF   |V|   IPM(ms)   Mosek     Mosek-D    P/IPM   D/IPM
tri linear        1e+05   27    17       0.3       0.7       0.8   2.09x   2.62x
tri quad          1e+05   27    17       0.3       0.6       1.2   1.76x   3.41x
C4 linear         1e+05   34    21       0.4       0.8       1.1   1.85x   2.52x
C4 quad           1e+05   34    21       0.7       0.7       1.1   1.00x   1.51x
radial linear     1e+05   23    16       0.2       0.5       1.0   2.18x   4.24x
radial quad       1e+05   23    16       0.4       0.6       1.0   1.57x   2.46x
```

**Extended OPF:**
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

### Clarabel

**Basic OPF:**
```
Config           raug   DOF   |V|   IPM(ms)   Clarabel  Clarabel-D  P/IPM   D/IPM
tri linear        1e+05   27    17       0.3       0.6       2.5   1.72x   7.66x
tri quad          1e+05   27    17       0.3       0.5       2.1   1.38x   6.27x
C4 linear         1e+05   34    21       0.4       0.5       2.1   1.37x   5.35x
C4 quad           1e+05   34    21       0.7       0.6       2.2   0.84x   3.36x
radial linear     1e+05   23    16       0.2       0.4       2.0   1.86x   8.85x
radial quad       1e+05   23    16       0.4       0.4       2.0   1.05x   5.41x
```

**Extended OPF:**
```
Config                raug   DOF   IPM(ms)  Clarabel  Clarabel-D  P/IPM   D/IPM
tri linear            1e+05   27     0.28     0.49     2.18    1.76x    7.85x
tri quad              1e+05   27     0.34     0.45     2.21    1.34x    6.57x
C4 quad               1e+05   34     0.74     0.59     2.49    0.79x    3.35x
tri +cycle            1e+05   48     5.47     1.07     1.46    0.20x    0.27x
C4 +cycle             1e+05   70    10.49     1.59     2.34    0.15x    0.22x
bowtie 2-tri          1e+05   80     5.05     1.77     2.23    0.35x    0.44x
bowtie full           1e+05   74     5.28     1.94     2.33    0.37x    0.44x
SE 3x3 r1             1e+03   57     1.28     0.55     2.66    0.43x    2.07x
SE 5x5 r1             1e+03  185     3.30     1.46     4.14    0.44x    1.25x
SE 7x7 r1             1e+03  385     6.85     3.14     6.86    0.46x    1.00x
SE triangle r2        1e+02   36     2.09     0.89     1.05    0.42x    0.50x
SE c4 r2              1e+02   56     8.28     1.38     1.59    0.17x    0.19x
```

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

### Reference solvers win

| Benchmark | Speedup | Notes |
|-----------|---------|-------|
| obstacle_option | OSQP 10-50x faster | Sparse Q exploitation |
| num_sheaf large mesh | Mosek 3-6x faster | Pure hyperedge coupling |
| mle_spline density | Mosek 5-20x faster | Many exp-cone leaves |
| lifted diagonal Q | OSQP 2-6x faster | No dense Q advantage |
