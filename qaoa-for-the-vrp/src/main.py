#!/usr/bin/env python3

import numpy as np
from collections import Counter

from Qaoa import QAOA
from Qubo import QUBO
from GraphOutils import SolutionUtils



if __name__ == "__main__":
    # ========== Instance D1 in the paper ============
    #D = [
    #    [0.0, 36.84, 5.06, 30.63],
    #    [36.84, 0.0, 24.55, 63.22],
    #    [5.06, 24.55, 0.0, 15.50],
    #    [30.63, 63.22, 15.50, 0.0]
    #]
    #k = 2
    # ========== Instance D2 in the paper ============
    D = [
        [0,      6.794,  61.653, 24.557, 47.767],
        [6.794,  0,      87.312, 47.262, 39.477],
        [61.653, 87.312, 0,      9.711,  42.887],
        [24.557, 47.262, 9.711,  0,      40.98 ],
        [47.767, 39.477, 42.887, 40.98,  0     ]
    ]
    k = 2
    A = 1000.0 
    
    # objects
    qubo = QUBO()              # build the model and Hamiltonian
    qaoa = QAOA()              # run QAOA
    utils = SolutionUtils()    # evaluate a solution

    # Build QUBO
    qubo_lin, qubo_quad, const_term, var_map, Q, h = qubo.build_vrp_qubo(D, k, A, print_matrices=True)
    H = qubo.qubo_to_hamiltonian(qubo_lin, qubo_quad, const_term, var_map)
    print("Hamiltonian:\n",H)

    # Run QAOA
    samples, min_cost, opt_params = qaoa.run(
        qubo_lin, qubo_quad, const_term, var_map, p=6
    )

    # Build edge_map once
    edge_map = {q: ij for ij, q in var_map.items()}
    
    # ===== Best solution from all samples =====
    bitstrings = ["".join(map(str, s)) for s in samples]
    best_bitstring, _ = Counter(bitstrings).most_common(1)[0]
    best_sample = [int(b) for b in best_bitstring]

    best_edges, best_route_cost = utils.decode_solution(best_sample, var_map, D)
    full_qubo_best = qubo.compute_qubo_value(best_sample, Q, h, const_term)
    arc_cost_best = utils.compute_arc_cost(best_sample, edge_map, D)
    penalty_best = full_qubo_best - arc_cost_best

    print("========= Best Solution Found =========")
    print("Bitstring       :", best_bitstring)
    print("Edges           :", best_edges)
    print(f"Route cost      : {best_route_cost:.2f}")
    print(f"Arc Cost        : {arc_cost_best:.2f}")
    print(f"Penalty Term    : {penalty_best:.2f}")
    print(f"Total QUBO Val. : {full_qubo_best:.2f}")

    utils.plot_solution_graph(best_edges, len(D))
    
    print("========= Fixed Solution Evaluation =========")
    #x1 = np.array([1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0])
    x2 = np.array([1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0]) #Best solution for D2
    print("len(bitstring):", len(x2))
    print("len(inv_map):", len(var_map))
    print("inv_map keys:", sorted(var_map.keys()))    
    results = utils.evaluate_solution(
        x2,
        var_map,
        D,
        Q,
        h,
        const_term,
        edge_map,
        compute_qubo_value_fn=qubo.compute_qubo_value
    )
    print(results)    
