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
            int &myRank,
            int &numprocs,
            int &nParameters,
            Solver *solver);
    vec runMinimizer(vec &guess, int &nCycles);
};

#endif // STEEPESTDESCENT_H
