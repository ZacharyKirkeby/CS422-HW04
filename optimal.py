"""Part 2: best topology for a given traffic matrix T Gurobi."""
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def best_topology(T, d, multigraph=True, time_limit=600, verbose=False):
    """
    Maximize concurrent flow for a fixed traffic matrix T.

    T[i,j]   : demand from node i to node j
    x[i,j]   : number of unit-capacity links i -> j
    f[s,i,j] : flow originating at source s on edge i -> j
    lam      : common scaling factor applied to every demand

    multigraph=False limits each directed node pair to one link.
    Returns (lam, X, status, mip_gap).
    """

    T = np.asarray(T, dtype=float)

    if T.ndim != 2 or T.shape[0] != T.shape[1]:
        raise ValueError("T must be square")

    n = T.shape[0]
    V = range(n)

    if n < 2 or d < 1:
        raise ValueError("Invalid network size or degree")
    if not multigraph and d > n - 1:
        raise ValueError("Degree impossible without parallel links")
    if not np.isfinite(T).all() or np.any(T < 0):
        raise ValueError("T must contain finite, nonnegative demands")
    if not np.allclose(np.diag(T), 0):
        raise ValueError("Self demand must be zero")

    row_load = T.sum(axis=1)
    col_load = T.sum(axis=0)

    if np.any(row_load > d + 1e-9) or np.any(col_load > d + 1e-9):
        raise ValueError("T violates hose-model bounds")
    if T.sum() == 0:
        raise ValueError("Zero demand has unbounded concurrent flow")

    # Every source and destination has capacity d.
    # This bounds lambda even when T is not fully saturated.
    lam_ub = d / max(row_load.max(), col_load.max())

    m = gp.Model("best_topology")
    m.Params.OutputFlag = int(verbose)
    m.Params.TimeLimit = time_limit
    m.Params.MIPGap = 0.0

    lam = m.addVar(lb=0, ub=lam_ub, name="lam")
    x = m.addVars(V, V, vtype=GRB.INTEGER, lb=0, ub=(d if multigraph else 1), name="x")
    f = m.addVars(V, V, V, lb=0, name="f")

    for i in V:
        # C1: there cannot be any self loops
        m.addConstr(x[i, i] == 0)

        # C2: Exactly d outgoing and d incoming unit links.
        m.addConstr(gp.quicksum(x[i, j] for j in V) == d)
        m.addConstr(gp.quicksum(x[j, i] for j in V) == d)

    # C3: Aggregate traffic cannot exceed installed link capacity
    # Parallel links are represented by x[i,j] > 1
    for i in V:
        for j in V:
            m.addConstr(gp.quicksum(f[s, i, j] for s in V) <= x[i, j])

    for s in V:
        # C4: Flow returning to its source is unnecessary.
        # Cycles can be removed without reducing delivered flow.
        for i in V:
            m.addConstr(f[s, i, s] == 0)

        # C5: Non-source nodes conserve flow or consume demand.
        # Source conservation follows by summing these equations
        # and using C4, so it need not be added separately.
        for v in V:
            if v == s:
                continue
            inflow = gp.quicksum(f[s, u, v] for u in V)
            outflow = gp.quicksum(f[s, v, w] for w in V)
            m.addConstr(inflow - outflow == lam * T[s, v])

    m.setObjective(lam, GRB.MAXIMIZE)
    m.optimize()

    if m.SolCount == 0:
        raise RuntimeError(f"No feasible solution; status={m.Status}")

    X = np.array([[int(round(x[i, j].X)) for j in V]for i in V])

    return lam.X, X, m.Status, m.MIPGap