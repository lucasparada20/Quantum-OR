# Quantum OR

Quantum OR is a collection of research and implementation projects at the
intersection of quantum computing and operations research. The repository
contains QUBO formulations, quantum optimization experiments, and classical
algorithms for quantum hardware design problems.

## Projects

| Project | Description | Technologies |
|---|---|---|
| [QAOA for the VRP](qaoa-for-the-vrp/) | Formulates a simplified Vehicle Routing Problem as a QUBO and solves it with the Quantum Approximate Optimization Algorithm. It also includes a CPLEX implementation of the QUBO model. | Python, PennyLane, SciPy, C++, CPLEX |
| [DFS for QUBO register design](dfs-qubo-register-design/) | Uses depth-first branch-and-bound search to place atoms in a circular register so that their pairwise interactions approximate a target QUBO matrix. | C++ |

## Repository structure

```text
Quantum-OR/
├── qaoa-for-the-vrp/          # QAOA and classical QUBO model for the VRP
├── dfs-qubo-register-design/  # DFS heuristic for quantum register design
└── README.md
```

Each project has its own README with its mathematical formulation, algorithm,
requirements, and execution instructions.

## Quick start

Clone the repository:

```bash
git clone https://github.com/lucasparada20/Quantum-OR.git
cd Quantum-OR
```

### QAOA for the VRP

The Python implementation is in `qaoa-for-the-vrp/src`. From that directory,
run:

```bash
cd qaoa-for-the-vrp/src
python3 main.py
```

See the [project README](qaoa-for-the-vrp/README.md) for the QUBO formulation,
Python module descriptions, and CPLEX build instructions.

### DFS for QUBO register design

Compile and run the C++ implementation:

```bash
cd dfs-qubo-register-design
g++ -O3 -Wall -Wextra dfs_qubo.cpp -o dfs_qubo
./dfs_qubo
```

See the [project README](dfs-qubo-register-design/README.md) for the objective
function, notation, search procedure, and stopping criteria.

## Scope

The repository focuses on examples that connect familiar operations-research
models with quantum computing, including:

- QUBO modeling
- Quantum approximate optimization
- Ising Hamiltonian construction
- Classical solution methods for quantum register design
- Comparisons between classical and quantum-oriented approaches

The implementations are intended for research, experimentation, and education.
