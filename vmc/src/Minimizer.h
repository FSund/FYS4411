#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <armadillo>
#include "Solver/Solver.h"
#include "Solver/SolverMCBF.h"

using namespace std;
using namespace arma;

class Minimizer
{
public:
    Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters, double *guess=0);
    double* runMinimizer();

protected:
    int myRank, numprocs;
    int nParameters;
    double *parameters;

    Solver* solver;
};

#endif // MINIMIZER_H
