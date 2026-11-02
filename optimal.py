"""Part 2: best topology for a given traffic matrix T Gurobi."""
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def best_topology(T, d, multigraph=True, time_limit=600, verbose=False):
    """
    x[i,j]   : number of unit links i -> j (topology, integer)
    f[s,i,j] : traffic that starts at s and uses link i -> j (routing)
    lam      : factor by which every demand is scaled (maximized)
    multigraph=False allows at most one link per ordered pair
    Returns (lam, X, status, mip_gap)
    C1 through C5 corrosponding to the doc
    """

    T = np.asarray(T, dtype=float)
    n = T.shape[0]
    V = range(n)
    m = gp.Model("best_topology")
    m.Params.OutputFlag = int(verbose)
    m.Params.TimeLimit = time_limit
    m.Params.MIPGap = 0.0
    lam = m.addVar(lb=0.0, name="lam")
    x = m.addVars(V, V, vtype=GRB.INTEGER, lb=0,
                  ub=(d if multigraph else 1), name="x")
    f = m.addVars(V, V, V, lb=0.0, name="f")
    for i in V:
        m.addConstr(x[i, i] == 0)                           # (C1)
        m.addConstr(gp.quicksum(x[i, j] for j in V) == d)   # (C2) out
        m.addConstr(gp.quicksum(x[j, i] for j in V) == d)   # (C2) in
    for i in V:
        for j in V:
            m.addConstr(gp.quicksum(f[s, i, j] for s in V)
                        <= x[i, j])                         # (C3)
    for s in V:
        for i in V:
            m.addConstr(f[s, i, s] == 0)                    # (C4)
    for s in V:
        for v in V:
            if v != s:                                      # (C5)
                inflow = gp.quicksum(f[s, u, v] for u in V)
                outflow = gp.quicksum(f[s, v, w] for w in V)
                m.addConstr(inflow - outflow == lam * T[s, v])
    m.setObjective(lam, GRB.MAXIMIZE)
    m.optimize()
    X = np.array([[int(round(x[i, j].X)) for j in V] for i in V])
    return lam.X, X, m.Status, m.MIPGap