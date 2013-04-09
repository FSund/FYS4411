#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include <armadillo>
#include "Solver/Solver.h"
#include "Solver/SolverMCBF.h"
#include "Solver/SolverMCIS.h"

using namespace std;
using namespace arma;

class Minimizer
{
public:
    Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters, vec &guess);
    vec runMinimizer();
protected:
    int myRank, numprocs;
    int nParameters;
    vec parameters;

    Solver* solver;
private:
    void bruteforce();

    mat energyGradientNumerical();
    vec steepestDescent();
};

#endif // MINIMIZER_H
