# Exact (rational arithmetic) versions of antivoter_var.py and check_bound.py.
# For each graph: solve the pair-correlation equations over the rationals, and compute the
# diagonal of (I+P)^{-1} over the rationals. Checks Var S <= n(n-2)/(2n-3) with equality only for K_n,
# and diag((I+P)^{-1}) >= (2n-3)/(2(n-2)) with equality only for K_n.
from fractions import Fraction as Fr
import numpy as np
import networkx as nx
from networkx.generators.atlas import graph_atlas_g


def solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


def var_S(G):
    nodes = list(G.nodes()); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    d = G.degree(nodes[0])
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}
    A = [[Fr(0)] * len(pairs) for _ in pairs]; rhs = [Fr(0)] * len(pairs)
    nb = [[idx[u] for u in G.neighbors(v)] for v in nodes]
    w = Fr(1, d)
    for k, (a, b) in enumerate(pairs):
        A[k][k] += 2
        for u in nb[a]:
            if u == b: rhs[k] -= w
            else: A[k][pid[(min(u, b), max(u, b))]] += w
        for u in nb[b]:
            if u == a: rhs[k] -= w
            else: A[k][pid[(min(a, u), max(a, u))]] += w
    return n + 2 * sum(solve(A, rhs))


def min_diag(G):
    nodes = list(G.nodes()); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    d = G.degree(nodes[0])
    Q = [[Fr(0)] * n for _ in range(n)]
    for v in nodes:
        Q[idx[v]][idx[v]] += 1
        for u in G.neighbors(v): Q[idx[v]][idx[u]] += Fr(1, d)
    return min(solve(Q, [Fr(int(i == j)) for j in range(n)])[i] for i in range(n))


def check(G):
    n = G.number_of_nodes(); d = G.degree(next(iter(G.nodes())))
    v = var_S(G); Kn = Fr(n * (n - 2), 2 * n - 3); star = Fr(2 * n - 3, 2 * (n - 2))
    m = min_diag(G)
    complete = d == n - 1
    assert (v == Kn) if complete else (v < Kn), (n, d, v)
    assert (m == star) if complete else (m > star), (n, d, m)
    return v


if __name__ == "__main__":
    for G in graph_atlas_g()[1:]:
        n = G.number_of_nodes()
        degs = set(dict(G.degree()).values())
        if n < 3 or len(degs) != 1 or not nx.is_connected(G) or nx.is_bipartite(G) or degs == {2}: continue
        print("atlas n", n, "d", degs.pop(), "Var S", check(G))
    rng = np.random.default_rng(20261107); cnt = 0; gap = None
    for n in range(5, 15):
        for d in range(3, n):
            if (n * d) % 2: continue
            for trial in range(6):
                G = nx.random_regular_graph(d, n, seed=int(rng.integers(1 << 30)))
                if not nx.is_connected(G) or nx.is_bipartite(G): continue
                v = check(G); cnt += 1
                if d < n - 1:
                    g = Fr(n * (n - 2), 2 * n - 3) - v
                    gap = g if gap is None else min(gap, g)
    print("random graphs", cnt, "all checks exact; least gap to the K_n value among non-complete graphs", gap, float(gap))
