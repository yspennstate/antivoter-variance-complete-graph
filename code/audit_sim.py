"""Second-route check of P0253: estimate Var S by simulating the antivoter dynamics directly,
without the pair-correlation equations, and compare with the closed forms in the paper."""
import numpy as np
import networkx as nx

rng = np.random.default_rng(20261108)


def sim_var(G, chains=4000, burn=400, steps=4000):
    A = [list(G.neighbors(v)) for v in G.nodes()]
    n = len(A)
    deg = len(A[0])
    nb = np.array(A)                      # regular graph: n x d
    x = rng.choice([-1, 1], size=(chains, n))
    acc = []
    idx = np.arange(chains)
    for t in range(burn * n + steps * n):
        v = rng.integers(0, n, size=chains)
        w = nb[v, rng.integers(0, deg, size=chains)]
        x[idx, v] = -x[idx, w]
        if t >= burn * n and (t % n == 0):
            acc.append(x.sum(1).astype(float))
    s = np.array(acc)                     # samples x chains
    per = (s ** 2).mean(0)                # E S = 0 by symmetry
    return per.mean(), per.std() / np.sqrt(chains)


def vt_formula(G):
    P = nx.to_numpy_array(G) / G.degree(0)
    lam = np.linalg.eigvalsh(P)
    return G.number_of_nodes() ** 2 / (2 * np.sum(1 / (1 + lam)))


cases = [("K5", nx.complete_graph(5)), ("K7", nx.complete_graph(7)),
         ("Petersen", nx.petersen_graph()), ("C7", nx.cycle_graph(7)),
         ("C9(1,3)", nx.circulant_graph(9, [1, 3]))]
for name, G in cases:
    n = G.number_of_nodes()
    est, se = sim_var(G)
    print("%-9s n=%2d  sim %.4f +- %.4f   VT formula %.4f   K_n bound %.4f"
          % (name, n, est, se, vt_formula(G), n * (n - 2) / (2 * n - 3)), flush=True)
