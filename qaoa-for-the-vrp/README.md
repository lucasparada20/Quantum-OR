# A-Quantum-Approximate-Optimization-Algorithm-for-the-Vehicle-Routing-Problem

## Model Formulation

We model the Vehicle Routing Problem (VRP) as a Quadratic Unconstrained Binary Optimization (QUBO) following the approach of these [authors](https://ieeexplore.ieee.org/document/9774961). Let graph $G=(V,A)$, with $V = N \cup \{0\}$ be the set of nodes including the depot $\{0\}$, and $A(0)$ be the set of arcs that connect the depot. We are given a fleet of $k$ vehicles. The model follows:

**(1) Objective (total travel cost)**
$\min \ \sum_{(i,j)\in A} c_{ij}x_{ij}$

Where:
- $x_{ij}\in $\{0,1\}$ is 1 if arc $\(i,j\)$ is used.
- $c_{ij}$ is the travel cost on arc $\(i,j\)$.

**(2) Outgoing arc constraint**  
$\sum_{j \in V} x_{ij} = 1 \quad \forall i \in N$

**(3) Incoming arc constraint**  
$\sum_{j \in V} x_{ji} = 1 \quad \forall i \in N$

**(4) Number of vehicles leaving the depot**  
$\sum_{j \in A(0)} x_{0j} = k$

**(5) Number of vehicles returning to the depot**  
$\sum_{j \in A(0)} x_{j0} = k$

In practice, this is not a VRP as the subtour elimination constraints are lacking. The problem in (1)-(5) resembles an assignment problem, but it is a simple enough problem to embed in a QAOA for operations research practitioners transitioning into quantum information technologies.

## QUBO Formulation for the VRP

```math
\begin{aligned}
\min_{x \in \{0,1\}} \quad 
& \underbrace{\sum_{i \in V} \sum_{\substack{j \in V \\ j \neq i}} c_{ij}\, x_{ij}}_{\text{(1) travel cost}} \\[1em]
&+ A \underbrace{\sum_{i \in N} \left( \sum_{\substack{j \in V \\ j \neq i}} x_{ij} - 1 \right)^{2}}_{\text{(2) outgoing arc constraints}} \\[1em]
&+ A \underbrace{\sum_{i \in N} \left( \sum_{\substack{j \in V \\ j \neq i}} x_{ji} - 1 \right)^{2}}_{\text{(3) incoming arc constraints}} \\[1em]
&+ A \underbrace{\left( \sum_{j \in V \setminus \{0\}} x_{0j} - k \right)^{2}}_{\text{(4) depot out-degree}} \\[1em]
&+ A \underbrace{\left( \sum_{j \in V \setminus \{0\}} x_{j0} - k \right)^{2}}_{\text{(5) depot in-degree}}
\end{aligned}
```

There is a C++ `main.cpp` file that will build and solve the QUBO formulation using the Cplex library. The code requires a hardcoded distance matrix and the parameter $k$. In the current implementation, there are two distance matrices hardcoded, so you can uncomment one to test it. To compile and run the code, you first need to give the path to your Cplex library in the `Makefile` here:

```bash
CPLEX_DIR = /some/path/to/Cplex
```
Next, compile, link, and run:

```bash
make
./vrp_qubo
```

## The Python Source Code

A modular implementation of QAOA is provided in the `src` directory. To call the code, use the following command inside the `src` directory:

```bash
./main
```

The modules are as follows.

* `Qubo.py`: Encodes the problem in (1)-(5) as a QUBO model and provides utilities to map it into an Ising Hamiltonian for use in quantum circuits.
  * **`build_vrp_qubo(D, k, A, print_matrices=False)`** – Constructs the QUBO formulation of the distance matrix $D$, number of vehicles $k$, and penalty weight $A$. Returns linear and quadratic QUBO coefficients, the constant term, and a variable-to-qubit index map.  
  * **`qubo_to_hamiltonian(qubo_lin, qubo_quad, const, var_index_map)`** – Converts the QUBO coefficients into a PennyLane Hamiltonian, mapping binary variables to Pauli-Z operators in the Ising form.  
  * **`compute_qubo_value(sample, Q, h, constant)`** – Evaluates the QUBO objective value for a given binary solution vector using $x^T Q x + h^T x + \text{constant}$.
   
* `Qaoa.py`: Uses PennyLane's `default.qubit` simulator backend and a SciPy optimizer to minimize the QUBO Hamiltonian. It also performs multiple optimization restarts to help avoid poor local minima.
  * **`qaoa_layer(gamma, beta, H)`** – Applies one QAOA layer by evolving under the problem Hamiltonian for time $\gamma$, followed by RX rotations with angle $2*\beta$ on all qubits.  
  * **`qaoa_circuit(params, H, n_wires)`** – Builds the complete QAOA circuit with $p$ alternating layers of cost and mixing unitaries, starting from the uniform superposition state $\|+\rangle^{\otimes n}\$.  
  * **`run(qubo_lin, qubo_quad, const_term, var_map, p=1, n_shots=10000)`** – Converts a QUBO into a Hamiltonian, executes the QAOA circuit on PennyLane’s `default.qubit` simulator, and uses classical optimization with multiple restarts to minimize the expected cost, returning samples, the best energy, and the optimal parameters.

* `GraphOutils.py`: Provides utility functions for decoding, evaluating, and visualizing solutions produced by QAOA on the problem.
  * **`decode_solution(bitstring, var_map, D=None)`** – Decodes a sampled bitstring into selected edges $(i,j)$, and optionally computes the total route cost if a distance matrix is provided.  
  * **`compute_arc_cost(sample, edge_map, D)`** – Computes the travel cost of a QAOA solution by summing the costs of all arcs selected in the bitstring.  
  * **`plot_solution_graph(edges, n_nodes)`** – Plots a directed graph representation of the decoded VRP tour using NetworkX and Matplotlib.  
  * **`evaluate_solution(x, var_map, D, Q, h, const_term, edge_map, compute_qubo_value_fn)`** – Evaluates a full solution by decoding edges, computing route cost, QUBO objective value, arc costs, and penalties, and visualizes the resulting tour. Returns a dictionary with all relevant evaluation metrics.  

  


