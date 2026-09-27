import numpy as np
import pennylane as qml
from scipy.optimize import minimize
from Qubo import QUBO


class QAOA:
    def qaoa_layer(self, gamma: float, beta: float, H) -> None:
        """One QAOA layer."""
        qml.templates.ApproxTimeEvolution(H, gamma, 1)
        for i in range(len(H.wires)):
            qml.RX(2 * beta, wires=i)

    def qaoa_circuit(self, params: np.ndarray, H, n_wires: int) -> None:
        """Full QAOA circuit with p layers."""
        p = len(params) // 2
        gammas, betas = params[:p], params[p:]
        for i in range(n_wires):
            qml.Hadamard(wires=i)  # Initial |+> state
        for gamma, beta in zip(gammas, betas):
            self.qaoa_layer(gamma, beta, H)

    def run(
        self,
        qubo_lin,
        qubo_quad,
        const_term: float,
        var_map,
        p: int = 1,
        n_shots: int = 10000
    ):
        """Run QAOA optimization and sampling."""
        H = QUBO().qubo_to_hamiltonian(qubo_lin, qubo_quad, const_term, var_map)
        n_wires = len(var_map)
        dev = qml.device("default.qubit", wires=n_wires, shots=n_shots)

        @qml.qnode(dev)
        def qnode(params):
            self.qaoa_circuit(params, H, n_wires)
            return qml.sample()

        @qml.qnode(dev)
        def raw_energy_qnode(params):
            self.qaoa_circuit(params, H, n_wires)
            return qml.expval(H)

        def energy_qnode(params):
            return raw_energy_qnode(params)

        # Verbose cost fn
        iteration_counter = {'count': 0}
        def cost_fn_verbose(params):
            energy = energy_qnode(params)
            iteration_counter['count'] += 1
            print(f"[Iter {iteration_counter['count']:03d}] Cost: {energy:.6f}, Params: {np.round(params, 4)}")
            return energy

        best_val, best_params = float("inf"), None
        for seed in range(1):  # number of restarts
            init_params = np.random.uniform(0, np.pi, size=(2 * p,))
            print(f"\n--- Restart {seed+1} ---")
            res = minimize(cost_fn_verbose, init_params, method="COBYLA", options={"maxiter": 5000})
            if res.fun < best_val:
                best_val, best_params = res.fun, res.x

        print("\n========= Optimization Summary =========")
        print("Optimal QAOA parameters:", best_params)
        print("Minimum cost:", best_val)
        print("Constant:", const_term)

        samples = qnode(best_params)
        return samples, best_val, best_params
