#ifndef STEEPESTDESCENT_H
#define STEEPESTDESCENT_H

#include <armadillo>
#include <src/Minimizer/Minimizer.h>

using namespace arma;
using namespace std;

class SteepestDescent : public Minimizer
{
public:
    SteepestDescent(
            int myRank,
            int numprocs,
            int nParameters,
            Solver *solver);
    vec runMinimizer(vec &guess, int &nCycles, int iterMax_ = 100);
    void setCout(bool couts_) { couts = couts_; }
    void setPrintPath(bool printPath_) { printPath = printPath_; }
private:
    bool couts;
    bool printPath;
};

#endif // STEEPESTDESCENT_H
