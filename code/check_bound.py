# Checks, on random connected non-bipartite regular graphs, the three steps of the proof:
# (1) Var S = (n/2)(1 + mean edge correlation); (2) diag((I+P)^{-1}) >= (2n-3)/(2(n-2)), strict unless K_n;
# (3) Var S <= n(n-2)/(2n-3).
import numpy as np, networkx as nx
from antivoter_var import var_S
rng = np.random.default_rng(20261107)
worst = -1e9; cnt = 0
for n in range(5, 15):
    for d in range(3, n):
        if (n * d) % 2: continue
        for trial in range(6):
            G = nx.random_regular_graph(d, n, seed=int(rng.integers(1 << 30)))
            if not nx.is_connected(G) or nx.is_bipartite(G): continue
            A = nx.to_numpy_array(G); P = A / d
            Minv = np.linalg.inv(np.eye(n) + P)
            star = (2 * n - 3) / (2 * (n - 2))
            v = var_S(G); Kn = n * (n - 2) / (2 * n - 3)
            assert Minv.diagonal().min() >= star - 1e-12
            if d < n - 1: assert Minv.diagonal().min() > star + 1e-12
            assert v <= Kn + 1e-9, (n, d, v, Kn)
            worst = max(worst, v - Kn); cnt += 1
print("graphs", cnt, "max Var - K_n value", worst)
