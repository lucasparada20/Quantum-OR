#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>

using std::size_t;

#undef INFINITY
double INFINITY = 10E8;
double EPSILON = 1E-6;

class DFS {
public:
    struct Point {
        double x;
        double y;
    };

    DFS(std::vector<std::vector<double> >& qubo, double c6, double registerRadius, double minimumSeparation, double gridStep, int searchLimit)
        : q_(qubo), c6_(c6), radius_(registerRadius), minimumSeparation_(minimumSeparation), searchLimit_(searchLimit), visitedNodes_(0), bestError_(INFINITY), startTime_(clock()), timeLimit_(60.0) 
        {
            int n = (int)q_.size();
            positions_.resize(n);
            bestPositions_.resize(n); 
            placed_.assign(n, false); //Boolean flag that determines if atom is placed

            // Create the grid inside the circular register.
            for (double x = -radius_; x <= radius_ + EPSILON; x += gridStep) {
                for (double y = -radius_; y <= radius_ + EPSILON; y += gridStep) {
                    if (x * x + y * y <= radius_ * radius_ + EPSILON) 
                    {
                        grid_.push_back(Point{x, y});
                    }
                }
            }
            std::cout << "Grid size:" << grid_.size() << std::endl;

            // Container to sort the atoms
            std::vector<PairRank> rankedPairs;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PairRank rank;
                    rank.firstAtom = i;
                    rank.secondAtom = j;
                    if (q_[i][j] >= 0) {
                        rank.strength = q_[i][j];
                    } else {
                        rank.strength = -q_[i][j];
                    }
                    rankedPairs.push_back(rank);
                }
            }
            
            // Sort pairs by decreasing |Q_ij|.
            std::sort(rankedPairs.begin(), rankedPairs.end(), PairRankComparison());
            
            //Add atoms
            std::vector<bool> added(n, false);
            for (size_t i = 0; i < rankedPairs.size(); i++) {
                int first = rankedPairs[i].firstAtom;
                int second = rankedPairs[i].secondAtom;

                if (!added[first]) {
                    atomOrder_.push_back(first);
                    added[first] = true;
                }
                if (!added[second]) {
                    atomOrder_.push_back(second);
                    added[second] = true;
                }
            }
            //Add remaining atoms, if any ...    
            for (int i = 0; i < n; i++) {
                if (!added[i]) {
                    atomOrder_.push_back(i);
                }
            }

            printf("Atom sorting: ");
            for (int i = 0; i < n; i++) {
                std::cout << atomOrder_[i] << " ";
            }
            std::cout << std::endl;
        }

    bool solve() {
        search(0, 0.0);
        return bestError_ < INFINITY;
    }

    void printSolution() {
        if (bestError_ >= INFINITY) {
            std::cout << "No feasible placement was found.\n";
            return;
        }

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Best error: " << bestError_ << "\n";
        for (size_t i = 0; i < bestPositions_.size(); i++) {
            std::cout << "Atom " << i
                      << ": (" << std::setprecision(1) << bestPositions_[i].x
                      << ", " << bestPositions_[i].y << ")\n";
        }
        std::cout << "DFS nodes visited: " << visitedNodes_ << "\n";
        std::cout << "Elapsed time: " << elapsedTime() << " seconds\n";
    }

private:
    struct Candidate {
        Point point;
        double addedError;
    };

    struct PairRank {
        int firstAtom;
        int secondAtom;
        double strength;
    };

    struct PairRankComparison {
        bool operator()(PairRank a, PairRank b) {
            return a.strength > b.strength;
        }
    };

    struct CandidateComparison {
        bool operator()(Candidate a, Candidate b) {
            return a.addedError < b.addedError;
        }
    };

    std::vector<std::vector<double> > q_;
    double c6_;
    double radius_;
    double minimumSeparation_;
    int searchLimit_;
    int visitedNodes_;
    double bestError_;
    clock_t startTime_;
    double timeLimit_;
    std::vector<Point> grid_;
    std::vector<int> atomOrder_;
    std::vector<Point> positions_;
    std::vector<Point> bestPositions_;
    std::vector<bool> placed_;

    double elapsedTime() {
        return (double)(clock() - startTime_) / CLOCKS_PER_SEC;
    }

    double distance(Point& a, Point& b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return sqrt(dx * dx + dy * dy);
    }

    //Checks the minimum distance constraint of the register design problem
    bool isFeasible(Point& candidate) {
        for (size_t i = 0; i < placed_.size(); i++) {
            if (placed_[i] && distance(candidate, positions_[i]) < minimumSeparation_ - EPSILON) {
                return false;
            }
        }
        return true;
    }
    
    //Computes the delta energy for inserting the atom in the candidate point
    double incrementalError(int atom, Point& candidate) {
        double error = 0.0;
        for (size_t other = 0; other < placed_.size(); other++) {
            if (!placed_[other]) {
                continue;
            }

            double atomDistance = distance(candidate, positions_[other]);
            double difference = c6_ / atomDistance - q_[atom][other];
            error += difference * difference;
        }
        return error;
    }

    //The implementation of the DP algorithm
    void search(int depth, double error) {
        if (visitedNodes_ >= searchLimit_ || elapsedTime() >= timeLimit_) {
            std::cout << "Search limit reached: visitedNodes=" << visitedNodes_
                      << "/" << searchLimit_
                      << ", elapsed=" << elapsedTime()
                      << "/" << timeLimit_ << '\n';
            return;
        }
        visitedNodes_++;

        //We finished!
        if (depth == (int)atomOrder_.size()) 
        {
            bestError_ = error;
            bestPositions_ = positions_;
            return;
        }

        int atom = atomOrder_[depth];

        //Printing some output
        if (visitedNodes_ % 50 == 0) {
            std::cout << "Visited nodes: " << visitedNodes_
                      << ", current atom: " << atom
                      << ", grid size: " << grid_.size()
                      << ", error: " << error
                      << ", best error: " << bestError_
                      << ", elapsed time: " << elapsedTime() << " seconds\n";
        }
        
        //Generates ALL candidate positions for the atom
        //This is very expensive even for small register sizes ....
        //In the Latex, this is C_k
        std::vector<Candidate> candidates;
        for (size_t i = 0; i < grid_.size(); i++) {
            if (isFeasible(grid_[i])) 
            {
                Candidate candidate;
                candidate.point = grid_[i];
                candidate.addedError = incrementalError(atom, grid_[i]);
                candidates.push_back(candidate);
            }
        }

        std::sort(candidates.begin(), candidates.end(), CandidateComparison());

        for (size_t i = 0; i < candidates.size(); i++) {
            //Check if you are running forever ...
            if (visitedNodes_ >= searchLimit_ || elapsedTime() >= timeLimit_)
                break;

            double newError = error + candidates[i].addedError;
            if (newError < bestError_) {
                positions_[atom] = candidates[i].point;
                placed_[atom] = true;
                search(depth + 1, newError); //Very expensive DP call due to regeneration of candidates ... but such is life.
                placed_[atom] = false;
            }
        }
    }
};

int main() {
    // Hardcoded symmetric QUBO matrix from the techical assessment. Diagonal entries are not used by DFS...
    std::vector<std::vector<double> > qubo = {
        {-17.0, 10.0, 10.0, 10.0,  0.0, 20.0},
        { 10.0,-18.0, 10.0, 10.0, 10.0, 20.0},
        { 10.0, 10.0,-29.0, 10.0, 20.0, 20.0},
        { 10.0, 10.0, 10.0,-19.0, 10.0, 10.0},
        {  0.0, 10.0, 20.0, 10.0,-17.0, 10.0},
        { 20.0, 20.0, 20.0, 10.0, 10.0,-28.0}
    };

    double c6 = 1000.0; 
    double registerRadius = 5.0;
    double minimumSeparation = 1.0;
    double gridStep = 1.0; 
    int searchLimit = 1E6;

    std::cout << "DFS solver " 
        << "nbAtoms:" << qubo.size()
        << " regRadius:" << registerRadius
        << " gridStep:" << gridStep
        << " QUBO: " << std::endl;
        
    for(size_t i=0;i<qubo.size();i++)
    {
        for(size_t j=0;j<qubo[i].size();j++)
            std::cout << std::fixed << std::setprecision(1)  << qubo[i][j] << " ";
        std::cout << std::endl;
    }

    //Creates an object that generates the input data structures
    DFS solver(qubo, c6, registerRadius, minimumSeparation, gridStep, searchLimit);
    
    //The DP algorithm
    solver.solve();
    solver.printSolution();
    return 0;
}
