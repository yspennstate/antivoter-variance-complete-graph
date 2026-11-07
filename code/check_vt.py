# Corollary 1.3 check: on vertex-transitive graphs Var S = n^2 / (2 sum_i 1/(1+lambda_i)).
import numpy as np, networkx as nx
from antivoter_var import var_S
gs = {"Petersen": nx.petersen_graph(), "K_{2,2,2}": nx.complete_multipartite_graph(2, 2, 2),
      "prism C3xK2": nx.circular_ladder_graph(3), "C5xK2": nx.circular_ladder_graph(5),
      "C3xC3": nx.cartesian_product(nx.cycle_graph(3), nx.cycle_graph(3)),
      "circulant(11;1,3)": nx.circulant_graph(11, [1, 3]), "circulant(13;1,5)": nx.circulant_graph(13, [1, 5]),
      "Moebius-Kantor": nx.LCF_graph(16, [5, -5], 8), "C7": nx.cycle_graph(7)}
worst = 0
for name, G in gs.items():
    d = G.degree(next(iter(G.nodes()))); n = G.number_of_nodes()
    if nx.is_bipartite(G): print(name, "bipartite, skipped"); continue
    lam = np.linalg.eigvalsh(nx.to_numpy_array(G) / d)
    f = n * n / (2 * np.sum(1 / (1 + lam))); v = var_S(G)
    worst = max(worst, abs(f - v))
    print(f"{name:20s} n={n:2d} d={d} VarS={v:.10f} formula={f:.10f} Kn={n*(n-2)/(2*n-3):.4f}")
print("max |formula - exact|", worst)
