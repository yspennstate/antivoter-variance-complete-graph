# The complete graph maximizes the variance of the antivoter model

Yitzchak Shmalo

Paper: [paper/main.pdf](paper/main.pdf) (source [paper/main.tex](paper/main.tex)).

In the antivoter model on a finite graph each vertex, at rate one, adopts the opinion opposite to that of a uniformly
chosen neighbour. Let $S$ be the sum of the $\pm1$ opinions in stationarity. Aldous and Fill (Reversible Markov Chains
and Random Walks on Graphs, Section 14.4) asked whether, among regular graphs on $n$ vertices, $\mathrm{Var}\,S$ is
largest for the complete graph, where it equals $n(n-2)/(2n-3)$. The paper proves that it is: for every connected
non-bipartite regular graph on $n\ge3$ vertices and every stationary distribution, $\mathrm{Var}\,S\le n(n-2)/(2n-3)$,
with equality only for $K_n$. For vertex-transitive graphs it gives
$\mathrm{Var}\,S=n^2/\bigl(2\sum_i(1+\lambda_i)^{-1}\bigr)$, where $\lambda_i$ are the eigenvalues of the
simple random walk.

## Code

Python 3 with numpy and networkx. Each script solves the linear equations for the pair correlations
$c(v,w)=\mathbb E[\eta_v\eta_w]$ exactly.

- `antivoter_var.py`: $\mathrm{Var}\,S$ on every connected non-bipartite regular graph with at most seven vertices
  (cycles excluded), from the networkx graph atlas. Output `antivoter_var.out` (Table 1 of the paper).
- `check_bound.py`: 297 random connected non-bipartite regular graphs with $5\le n\le14$; checks the diagonal bound on
  $(I+P)^{-1}$ and the variance bound. Output `check_bound.out`.
- `check_vt.py`: the vertex-transitive formula against the exact solution on eight graphs. Output `check_vt.out`.

## Funding

The research presented in this paper was supported by the European Research Council (ERC) under the European
Union's Horizon Europe research and innovation programme (grant agreement No. 101041711), by the Simons Foundation
as part of the Collaboration on the Mathematical and Scientific Foundations of Deep Learning, by Heights Labs, by
Convex Nexus Capital, by the Israel Science Foundation (grant number 2258/19), and by the Israel Science Foundation
(ISF Grant 4101/25). Computational resources were provided in part by Google Cloud through its education credits
program (credit 530550103).
