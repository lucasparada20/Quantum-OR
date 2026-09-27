# Depth-First Search for the Register Design Problem

This program implements the depth-first search (DFS). It assigns positions to atoms so that their pairwise interactions approximate the coefficients of a QUBO matrix $Q$.

## Problem notation

The notation is as follows:

- $A=(a_1,\ldots,a_n)$: sorted atom indices
- $\mathcal{G}$: finite grid of candidate positions
- $Q$: target QUBO matrix
- $C_6$: interaction constant
- $R_{\mathrm{reg}}$: register radius
- $d_{ij}$: minimum separation between atoms $i$ and $j$
- $S$: set of already placed atoms
- $p_i$: assigned position of atom $i$
- $C_k$: feasible candidate positions for atom $a_k$
- `bestError`: error of the best complete placement found
- `bestPlacement`: best complete placement found

The continuous register is discretized into $\mathcal{G}$. A grid point
$c\in\mathcal{G}$ belongs to $C_k$ when

```text
||c||₂ <= R_reg
||c - p_j||₂ >= d_(a_k,j)    for every j in S.
```

The C++ example uses one common value, `minimumSeparation_`, for every
$d_{ij}$.

## Objective function

DFS calculates the incremental error associated with placing atom $a_k$
at $c\in C_k$ as

$$
\Delta(a_k,c)=\sum_{j\in S}
\left(\frac{C_6}{\lVert c-p_j\rVert_2^6}-Q_{a_kj}\right)^2.
$$

The error of the resulting partial placement is

```text
newError = error + Delta(a_k, c).
```

Thus, $\Delta(a_k,c)$ is the additional error introduced by the new
placement, while `newError` is the accumulated error of the partial placement. Once all $n$ atoms are placed, `error` equals the total objective
value

$$
\sum_{(i,j)\in E}
\left(\frac{C_6}{\lVert p_i-p_j\rVert_2^6}-Q_{ij}\right)^2.
$$

## DFS heuristic

The search begins with

```text
DFS(1, empty set, 0).
```

For the current atom $a_k$, DFS performs the following operations:

1. Stop if the search limit has been reached.
2. Prune the branch if `error >= bestError`.
3. If $k>n$, save `error` and the current positions as `bestError` and
   `bestPlacement`.
4. Construct $C_k\subseteq\mathcal{G}$ using the register and separation
   constraints.
5. Compute $\Delta(a_k,c)$ for each $c\in C_k$.
6. Sort $C_k$ by non-decreasing $\Delta(a_k,c)$.
7. For each candidate, set `newError = error + Delta(a_k, c)` and do a recursive call when
   `newError < bestError`.
8. Remove $a_k$ from $S$ when backtracking.

Because every term in the objective is squared, each
$\Delta(a_k,c)\geq 0$. The accumulated error cannot decrease as atoms are
added, so pruning a partial placement when `newError >= bestError` is valid.

To construct the atom order, define
$P=\{(i,j):i,j\in N,\ i<j\}$ and sort its pairs by non-increasing
$\lvert Q_{ij}\rvert$. Starting with an empty $A$, traverse the sorted
pairs. For every $(i,j)\in P$, append $i$ if it is not already in $A$,
then append $j$ if it is not already in $A$. Finally, append any remaining
atoms. This places atoms belonging to the strongest target interactions earlier
in the DFS search.

## Search limits

The search stops when either `searchLimit_` DFS nodes have been visited or the
60-second time limit has been reached. When a limit is detected at the start of
a DFS call, the program prints the visited-node count and elapsed time.

## Building and running

Compile with a C++11-compatible compiler:

```bash
g++ -O3 -Wall -Wextra dfs_qubo.cpp -o dfs_qubo
```

Call the executable:

```bash
./dfs_qubo
```

The example $Q$, $C_6$, $R_{\mathrm{reg}}$, grid spacing, common
$d_{ij}$, and search limit are defined in `main()`
