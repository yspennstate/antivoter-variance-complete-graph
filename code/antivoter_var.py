# Exact stationary variance of S = sum of antivoter spins on a connected regular graph.
# Pair correlations c(v,w) = E[eta_v eta_w], v != w, solve
#   2 c(v,w) + (1/d) sum_{u~v} c(u,w) + (1/d) sum_{u~w} c(v,u) = 0,  c(u,u) = 1.
import itertools, sys
import numpy as np
import networkx as nx
from networkx.generators.atlas import graph_atlas_g

def var_S(G):
    n = G.number_of_nodes(); nodes = list(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}
    d = G.degree(nodes[0])
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}
    A = np.zeros((len(pairs), len(pairs))); rhs = np.zeros(len(pairs))
    nb = [[idx[u] for u in G.neighbors(v)] for v in nodes]
    def add(k, a, b, coef):
        if a == b: rhs[k] -= coef
        else: A[k, pid[(min(a, b), max(a, b))]] += coef
    for k, (a, b) in enumerate(pairs):
        A[k, k] += 2.0
        for u in nb[a]: add(k, u, b, 1.0 / d)
        for u in nb[b]: add(k, a, u, 1.0 / d)
    c = np.linalg.solve(A, rhs)
    return n + 2 * c.sum()

def admissible(G):
    return nx.is_connected(G) and not nx.is_bipartite(G) and not (max(dict(G.degree()).values()) == 2)

if __name__ == "__main__":
    best = {}
    for G in graph_atlas_g()[1:]:
        n = G.number_of_nodes()
        if n < 3: continue
        degs = set(dict(G.degree()).values())
        if len(degs) != 1 or not admissible(G): continue
        v = var_S(G)
        best.setdefault(n, []).append((v, G.degree(0), nx.is_isomorphic(G, nx.complete_graph(n))))
    for n in sorted(best):
        L = sorted(best[n], reverse=True)
        print(n, "graphs", len(L), "top:", [(round(v, 6), d, K) for v, d, K in L[:4]], "K_n formula", round(n*(n-2)/(2*n-3), 6))
