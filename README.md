# CS422-HW04
Topology Optimizations


Ensure you have gurobipy package via `python -m pip install gurobipy`
and a gurobi license ig

```
    case                                  ILP lam* (multi)  per-pair LP  direct bound  lam* (simple)    time
    Restricted license - for non-production use only - expires 2027-11-29
    uniform                                       0.700000     0.700000      0.700000       0.700000    0.8s
    topology rows: [[0, 1, 1, 0, 1, 1, 0, 0], [0, 0, 1, 0, 1, 0, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [1, 0, 0, 0, 0, 1, 1, 1], [1, 0, 1, 1, 0, 0, 0, 1], [0, 1, 0, 1, 1, 0, 1, 0], [1, 1, 1, 1, 0, 0, 0, 0], [1, 1, 0, 1, 0, 1, 0, 0]]
    hop profile (#d1,#d2,#d3) per node: {(4, 3, 0)}
    circulant C8(1,2,3,4) on uniform T: 0.625
    4 x cyclic permutation                        1.000000     1.000000      1.000000       0.611111    0.0s
    two groups of 4                               0.875000     0.875000      0.900000       0.812500    0.1s
    integer (4 random derangements)               1.000000     1.000000      1.000000       0.857143    0.0s
    4/3 x (3 random derangements)                 0.857143     0.857143      0.921429       0.778846    0.7s
    Sinkhorn(rand^2), seed 0                      0.822614     0.822614      0.842522       0.800190    0.5s
    Sinkhorn(rand^6), seed 0                      0.871063     0.871063      0.908779       0.766250    0.2s
```