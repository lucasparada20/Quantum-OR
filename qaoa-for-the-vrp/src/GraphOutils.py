import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from typing import Dict, List, Tuple, Optional, Any


class SolutionUtils:
    def decode_solution(
        self,
        bitstring: List[int],
        var_map: Dict[Tuple[int, int], int],
        D: Optional[List[List[float]]] = None
    ) -> Tuple[List[Tuple[int, int]], Optional[float]]:
        """
        Decode a QAOA bitstring into edges and optionally compute route cost.

        Args:
            bitstring: list/array of 0/1 from QAOA sample.
            var_map: mapping {(i,j): qubit_index}.
            D: optional distance matrix for computing route cost.

        Returns:
            edges: list of (i, j) edges selected.
            cost: total cost if D is provided, else None.
        """
        inv_map = {q: ij for ij, q in var_map.items()}

        edges: List[Tuple[int, int]] = [
            inv_map[q] for q, bit in enumerate(bitstring) if bit == 1
        ]

        cost: Optional[float] = None
        if D is not None:
            cost = sum(D[i][j] for (i, j) in edges)

        return edges, cost

    def compute_arc_cost(
        self,
        sample: List[int],
        edge_map: Dict[int, Tuple[int, int]],
        D: List[List[float]]
    ) -> float:
        """
        Compute total travel cost for arcs selected in a QAOA bitstring.

        Args:
            sample: list/array of 0/1, output of QAOA sample.
            edge_map: dict {qubit_index: (i, j)}.
            D: distance/cost matrix.

        Returns:
            total_cost: sum of selected arc costs.
        """
        total_cost = 0.0
        for q, bit in enumerate(sample):
            if bit == 1:
                i, j = edge_map[q]
                total_cost += D[i][j]
        return total_cost

    def plot_solution_graph(
        self,
        edges: List[Tuple[int, int]],
        n_nodes: int
    ) -> None:
        """
        Plot directed graph of a decoded VRP tour.

        Args:
            edges: list of directed edges (i, j).
            n_nodes: number of nodes in the graph.
        """
        G = nx.DiGraph()
        G.add_nodes_from(range(n_nodes))
        G.add_edges_from(edges)

        pos = nx.spring_layout(G, seed=42)

        plt.figure(figsize=(6, 5))
        nx.draw(
            G, pos,
            with_labels=True,
            node_color='lightblue',
            edge_color='red',
            node_size=700,
            arrowsize=20,
            width=2,
            font_weight='bold'
        )
        plt.title("Decoded VRP Tour")
        plt.axis('off')
        plt.show()

    def evaluate_solution(
        self,
        x: np.ndarray,
        var_map: Dict[Tuple[int, int], int],
        D: List[List[float]],
        Q: np.ndarray,
        h: np.ndarray,
        const_term: float,
        edge_map: Dict[int, Tuple[int, int]],
        compute_qubo_value_fn=None
    ) -> Dict[str, Any]:
        """
        Evaluate a given QUBO solution bitstring.
        """
        if compute_qubo_value_fn is None:
            raise ValueError("A compute_qubo_value function must be provided")

        # Decode into edges and cost
        edges_x, route_cost_x = self.decode_solution(x, var_map, D)

        # Compute QUBO-related costs
        full_qubo_x = compute_qubo_value_fn(x, Q, h, const_term)
        arc_cost_x = self.compute_arc_cost(x, edge_map, D)
        penalty_x = full_qubo_x - arc_cost_x

        results = {
            "bitstring": "".join(map(str, x)),
            "edges": edges_x,
            "route_cost": route_cost_x,
            "arc_cost": arc_cost_x,
            "penalty_term": penalty_x,
            "total_qubo_val": full_qubo_x,
        }

        # Plot
        self.plot_solution_graph(edges_x, len(D))

        return results
