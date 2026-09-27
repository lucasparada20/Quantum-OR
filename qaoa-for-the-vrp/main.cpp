#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>

ILOSTLBEGIN

// Add (sum_{j != i} x[i][j] - target)^2 to 'qubo' with weight A
static void add_out_degree_penalty(
    int i, int target, double A, int n,
    const IloArray<IloNumVarArray>& x, IloEnv& env, IloExpr& qubo)
{
    std::vector<IloNumVar> vars;
    vars.reserve(n-1);
    for (int j = 0; j < n; ++j) if (j != i) vars.push_back(x[i][j]);

    IloExpr s2(env);
    for (size_t a = 0; a < vars.size(); ++a) s2 += vars[a];
    for (size_t a = 0; a + 1 < vars.size(); ++a)
        for (size_t b = a + 1; b < vars.size(); ++b)
            s2 += 2.0 * (vars[a] * vars[b]);

    IloExpr s(env);
    for (size_t a = 0; a < vars.size(); ++a) s += vars[a];

    qubo += A * (s2 - 2.0 * target * s + double(target * target));

    s.end();
    s2.end();
}

// Add (sum_{j != i} x[j][i] - target)^2 to 'qubo' with weight A
static void add_in_degree_penalty(
    int i, int target, double A, int n,
    const IloArray<IloNumVarArray>& x, IloEnv& env, IloExpr& qubo)
{
    std::vector<IloNumVar> vars;
    vars.reserve(n-1);
    for (int j = 0; j < n; ++j) if (j != i) vars.push_back(x[j][i]);

    IloExpr s2(env);
    for (size_t a = 0; a < vars.size(); ++a) s2 += vars[a];
    for (size_t a = 0; a + 1 < vars.size(); ++a)
        for (size_t b = a + 1; b < vars.size(); ++b)
            s2 += 2.0 * (vars[a] * vars[b]);

    IloExpr s(env);
    for (size_t a = 0; a < vars.size(); ++a) s += vars[a];

    qubo += A * (s2 - 2.0 * target * s + double(target * target));

    s.end();
    s2.end();
}

int main() {
    IloEnv env;

    // ========= Data for instance D1 ==================================
    /*const int k = 2;
    const int n = 4;                 // nodes 0..3 ; depot = 0

    std::vector<std::vector<double> > D(n, std::vector<double>(n, 0.0));
    D[0][0] = 0.0;  D[0][1] = 36.84; D[0][2] = 5.06;  D[0][3] = 30.63;
    D[1][0] = 36.84;D[1][1] = 0.0;  D[1][2] = 24.55; D[1][3] = 63.22;
    D[2][0] = 5.06; D[2][1] = 24.55;D[2][2] = 0.0;  D[2][3] = 15.50;
    D[3][0] = 30.63;D[3][1] = 63.22;D[3][2] = 15.50; D[3][3] = 0.0;
	*/
	// ========= Data for instance D2 ==================================
	const int k = 2;
	const int n = 5; // nodes 0..4 ; depot = 0

	std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));

	D[0][0] = 0.0;    D[0][1] = 6.794;  D[0][2] = 61.653; D[0][3] = 24.557; D[0][4] = 47.767;
	D[1][0] = 6.794;  D[1][1] = 0.0;    D[1][2] = 87.312; D[1][3] = 47.262; D[1][4] = 39.477;
	D[2][0] = 61.653; D[2][1] = 87.312; D[2][2] = 0.0;    D[2][3] = 9.711;  D[2][4] = 42.887;
	D[3][0] = 24.557; D[3][1] = 47.262; D[3][2] = 9.711;  D[3][3] = 0.0;    D[3][4] = 40.980;
	D[4][0] = 47.767; D[4][1] = 39.477; D[4][2] = 42.887; D[4][3] = 40.980; D[4][4] = 0.0;


    const double A = 1000.0;         // penalty weight

    // ----- Model -----
    IloModel model(env);

    // x[i][j] for i != j (binary in [0,1])
    IloArray<IloNumVarArray> x(env, n);
    for (int i = 0; i < n; ++i) {
        x[i] = IloNumVarArray(env, n);
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                x[i][j] = IloNumVar(env, 0.0, 0.0); // fixed 0
            } else {
                std::ostringstream nm; nm << "x_" << i << "_" << j;
                x[i][j] = IloNumVar(env, 0.0, 1.0, nm.str().c_str());
            }
        }
    }

    // Linear travel cost: sum_{i != j} w_ij x_ij
    IloExpr travel(env);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j) travel += D[i][j] * x[i][j];

    // QUBO objective = travel + A * (degree penalties)
    IloExpr qubo(env);
    qubo += travel;

    // (2) sum_{j in source[i]} x_{ij} = 1 for i = 1..n-1  (out-degree)
    for (int i = 1; i <= n - 1; ++i)
        add_out_degree_penalty(i, 1, A, n, x, env, qubo);

    // (3) sum_{j in target[i]} x_{ji} = 1 for i = 1..n-1  (in-degree)
    for (int i = 1; i <= n - 1; ++i)
        add_in_degree_penalty(i, 1, A, n, x, env, qubo);

    // (4) sum_{j in source[0]} x_{0j} = k  (depot out-degree)
    add_out_degree_penalty(0, k, A, n, x, env, qubo);

    // (5) sum_{j in target[0]} x_{j0} = k  (depot in-degree)
    add_in_degree_penalty(0, k, A, n, x, env, qubo);

    IloObjective obj = IloMinimize(env, qubo);
    model.add(obj);

    // ----- Solve -----
    IloCplex cplex(model);
	cplex.exportModel("qubo.lp");
	cplex.setParam(IloCplex::Param::OptimalityTarget, 3); //nonconvex MIQP
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Param::Threads, 1);

    if (!cplex.solve()) {
        std::cout << "No solution\n";
        qubo.end();
        travel.end();
        obj.end();
        env.end();
        return 1;
    }

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Status: " << cplex.getStatus() << "\n";
    std::cout << "QUBO objective: " << std::setprecision(2) << cplex.getObjValue() << "\n";

    // Print chosen arcs x_ij = 1
    std::cout << "Arcs with x=1:\n";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && cplex.getValue(x[i][j]) > 0.5)
                std::cout << "  " << i << " -> " << j << "\n";

    // Simple route trace (no SECs: subtours possible)
    std::vector<std::vector<int> > succ(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && cplex.getValue(x[i][j]) > 0.5)
                succ[i].push_back(j);

    for (int r = 0; r < k; ++r) {
        int u = 0, steps = 0;
        std::cout << "Route " << (r+1) << ": " << u;
        while (!succ[u].empty() && steps <= 2*n) {
            int v = succ[u].front();
            succ[u].erase(succ[u].begin());
            std::cout << " -> " << v;
            u = v; ++steps;
            if (u == 0) break;
        }
        std::cout << "\n";
    }

    env.end();
    return 0;
}
