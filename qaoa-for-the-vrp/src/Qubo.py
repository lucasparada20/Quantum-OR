import numpy as np
import pennylane as qml
from collections import defaultdict
from typing import DefaultDict, Dict, Tuple


class QUBO:
    def build_vrp_qubo(
        self,
        D: list[list[float]],
        k: int,
        A: float,
        print_matrices: bool = False
    ) -> tuple[
        DefaultDict[Tuple[int, int], float],
        DefaultDict[Tuple[Tuple[int, int], Tuple[int, int]], float],
        float,
        Dict[Tuple[int, int], int],
        np.ndarray,
        np.ndarray,
    ]:
        """
        Build QUBO from VRP data.
        Returns: qubo_lin, qubo_quad, const, var_index_map, Q, h
        """
        n = len(D)
        lin = defaultdict(float)     # {(i,j): coef}
        quad = defaultdict(float)    # {((i1,j1),(i2,j2)): coef}
        const = [0.0]
        var_index_map = {}           # {(i,j): qubit_index}
        idx_counter = 0

        # Assign qubit indices to each x_ij (i != j)
        for i in range(n):
            for j in range(n):
                if i != j:
                    var_index_map[(i, j)] = idx_counter
                    idx_counter += 1

        # --- Helpers ---
        def add_out_degree_penalty(i, target):
            vars_idx = [(i, j) for j in range(n) if j != i]
            for a in range(len(vars_idx)):
                lin[vars_idx[a]] += A
                for b in range(a + 1, len(vars_idx)):
                    v1, v2 = sorted([vars_idx[a], vars_idx[b]])
                    quad[(v1, v2)] += 2 * A
            for v in vars_idx:
                lin[v] += -2 * A * target
            const[0] += A * target * target

        def add_in_degree_penalty(i, target):
            vars_idx = [(j, i) for j in range(n) if j != i]
            for a in range(len(vars_idx)):
                lin[vars_idx[a]] += A
                for b in range(a + 1, len(vars_idx)):
                    v1, v2 = sorted([vars_idx[a], vars_idx[b]])
                    quad[(v1, v2)] += 2 * A
            for v in vars_idx:
                lin[v] += -2 * A * target
            const[0] += A * target * target

        # --- Travel cost ---
        for i in range(n):
            for j in range(n):
                if i != j:
                    lin[(i, j)] += D[i][j]

        # --- Degree penalties ---
        for i in range(1, n):
            add_out_degree_penalty(i, 1)
        for i in range(1, n):
            add_in_degree_penalty(i, 1)
        add_out_degree_penalty(0, k)
        add_in_degree_penalty(0, k)

        # --- Build dense QUBO matrices ---
        Nq = len(var_index_map)
        Q = np.zeros((Nq, Nq))
        h = np.zeros(Nq)

        # Linear → h
        for (i, j), coef in lin.items():
            q = var_index_map[(i, j)]
            h[q] += coef

        # Quadratic → Q
        for ((i1, j1), (i2, j2)), coef in quad.items():
            q1 = var_index_map[(i1, j1)]
            q2 = var_index_map[(i2, j2)]
            Q[q1, q2] += coef
            Q[q2, q1] += coef  # symmetric

        if print_matrices:
            print("Variable index map (i,j) → qubit index:")
            for v, idx in var_index_map.items():
                print(f"  {v} → {idx}")
            print("\nLinear term vector h:")
            print(h)
            print("\nQuadratic term matrix Q:")
            print(Q)
            print("\nConstant term:")
            print(const[0])

        return lin, quad, const[0], var_index_map, Q, h

    def qubo_to_hamiltonian(
        self,
        qubo_lin: Dict[Tuple[int, int], float],
        qubo_quad: Dict[Tuple[Tuple[int, int], Tuple[int, int]], float],
        const: float,
        var_index_map: Dict[Tuple[int, int], int],
    ) -> qml.Hamiltonian:
        """
        Convert QUBO (binary variables) to Ising Hamiltonian for PennyLane.
        """
        const_term = const
        z_terms = defaultdict(float)
        zz_terms = defaultdict(float)

        # Linear: c * (1 - Z)/2
        for var, coef in qubo_lin.items():
            q = var_index_map[var]
            const_term += coef / 2
            z_terms[(q,)] += -coef / 2

        # Quadratic: c * (1 - Z_i - Z_j + Z_i Z_j)/4
        for (var1, var2), coef in qubo_quad.items():
            q1 = var_index_map[var1]
            q2 = var_index_map[var2]
            const_term += coef / 4
            z_terms[(q1,)] += -coef / 4
            z_terms[(q2,)] += -coef / 4
            zz_terms[tuple(sorted((q1, q2)))] += coef / 4

        coeffs = []
        ops = []

        # Constant term
        if abs(const_term) > 1e-9:
            coeffs.append(const_term)
            ops.append(qml.Identity(0))

        # Z terms
        for (q,), c in z_terms.items():
            if abs(c) > 1e-9:
                coeffs.append(c)
                ops.append(qml.PauliZ(q))

        # ZZ terms
        for (q1, q2), c in zz_terms.items():
            if abs(c) > 1e-9:
                coeffs.append(c)
                ops.append(qml.PauliZ(q1) @ qml.PauliZ(q2))

        return qml.Hamiltonian(coeffs, ops)

    def compute_qubo_value(
        self,
        sample: np.ndarray,
        Q: np.ndarray,
        h: np.ndarray,
        constant: float
    ) -> float:
        """
        Computes the total QUBO value of a binary solution:
        x^T Q x + h^T x + constant
        """
        x = np.array(sample)
        quad_term = x @ Q @ x  # Equivalent to xᵀ Q x
        linear_term = h @ x    # Equivalent to hᵀ x
        return quad_term + linear_term + constant
