#ifndef FILIP_H
#define FILIP_H

#include <armadillo>
#include <src/Minimizer/Minimizer.h>

using namespace arma;
using namespace std;

class Filip : public Minimizer
{
public:
    Filip(
            int &myRank,
            int &numprocs,
            int &nParameters,
            Solver *solver);
    vec runMinimizer(vec &guess, int &nCycles);
};

#endif // FILIP_H
