import numpy as np, time, gurobipy as gp
from gurobipy import GRB
from optimal import best_topology
n, d = 8, 4
rng = np.random.default_rng(0)

def pair_lp(T, X):
    """Max concurrent flow on the topology X using one commodity per (s,t) pair"""
    V = range(n); E = [(i, j) for i in V for j in V if X[i][j] > 0]
    pairs = [(s, t) for s in V for t in V if s != t and T[s][t] > 1e-12]
    m = gp.Model(); m.Params.OutputFlag = 0
    lam = m.addVar(lb=0); f = {(p, e): m.addVar(lb=0) for p in pairs for e in E}
    for e in E:
        m.addConstr(gp.quicksum(f[p, e] for p in pairs) <= X[e[0]][e[1]])
    for (s, t) in pairs:
        for v in V:
            out = gp.quicksum(f[(s, t), e] for e in E if e[0] == v)
            inn = gp.quicksum(f[(s, t), e] for e in E if e[1] == v)
            m.addConstr(out - inn == (lam * T[s][t] if v == s else (-lam * T[s][t] if v == t else 0)))
    m.setObjective(lam, GRB.MAXIMIZE); m.optimize(); return lam.X

def direct_link_bound(T):
    """B(T) = max_{x,lam} lam  s.t.  lam*sum(T) + sum_ij (lam*T_ij - x_ij)^+ <= n*d  (valid upper bound)."""
    V = range(n); m = gp.Model(); m.Params.OutputFlag = 0; m.Params.MIPGap = 0
    lam = m.addVar(lb=0); x = m.addVars(V, V, vtype=GRB.INTEGER, lb=0, ub=d); y = m.addVars(V, V, lb=0)
    for i in V:
        m.addConstr(x[i, i] == 0)
        m.addConstr(gp.quicksum(x[i, j] for j in V) == d); m.addConstr(gp.quicksum(x[j, i] for j in V) == d)
        for j in V: m.addConstr(y[i, j] >= lam * float(T[i][j]) - x[i, j])
    m.addConstr(lam * float(T.sum()) + gp.quicksum(y[i, j] for i in V for j in V) <= n * d)
    m.setObjective(lam, GRB.MAXIMIZE); m.optimize(); return lam.X

def hop_profile(X):
    """BFS distance histogram from each node"""
    out = []
    for s in range(n):
        dist = {s: 0}; q = [s]
        for u in q:
            for v in range(n):
                if X[u][v] > 0 and v not in dist: dist[v] = dist[u] + 1; q.append(v)
        out.append(tuple(sum(1 for v in dist if dist[v] == k) for k in (1, 2, 3)))
    return out

def derangement():
    while True:
        p = rng.permutation(n)
        if all(p[i] != i for i in range(n)): return np.eye(n)[p]

def sinkhorn(A):
    A = A.copy(); np.fill_diagonal(A, 0)
    for _ in range(3000): A /= A.sum(1, keepdims=True); A /= A.sum(0, keepdims=True)
    return 4 * A


# The magic
cases = {}
cases["uniform"] = 4 / 7 * (np.ones((n, n)) - np.eye(n))
P = np.zeros((n, n))

for i in range(n): P[i, (i + 1) % n] = 4
cases["4 x cyclic permutation"] = P

C = np.zeros((n, n))

for g in (range(0, 4), range(4, 8)):
    for i in g:
        for j in g:
            if i != j: C[i, j] = 4 / 3

cases["two groups of 4"] = C
cases["integer (4 random derangements)"] = sum(derangement() for _ in range(4))
cases["4/3 x (3 random derangements)"] = 4 / 3 * sum(derangement() for _ in range(3))
cases["Sinkhorn(rand^2), seed 0"] = sinkhorn(rng.random((n, n)) ** 2)
cases["Sinkhorn(rand^6), seed 0"] = sinkhorn(rng.random((n, n)) ** 6)

print(f"{'case':36s}{'ILP lam* (multi)':>18s}{'per-pair LP':>13s}{'direct bound':>14s}{'lam* (simple)':>15s}{'time':>8s}")
for name, T in cases.items():
    assert np.allclose(T.sum(1), 4) and np.allclose(T.sum(0), 4)
    t = time.time(); lam, X, st, gap = best_topology(T, d, multigraph=True); dt = time.time() - t
    assert st == GRB.OPTIMAL
    lam_s, _, st_s, _ = best_topology(T, d, multigraph=False)

    print(f"{name:36s}{lam:18.6f}{pair_lp(T, X):13.6f}{direct_link_bound(T):14.6f}{lam_s:15.6f}{dt:7.1f}s")
    
    if name == "uniform":
        print("   topology rows:", X.tolist()); print("   hop profile (#d1,#d2,#d3) per node:", set(hop_profile(X)))
        Xc = np.zeros((n, n), int)
        for i in range(n):
            for k in (1, 2, 3, 4): Xc[i, (i + k) % n] = 1
        print("   circulant C8(1,2,3,4) on uniform T:", round(pair_lp(T, Xc), 6))