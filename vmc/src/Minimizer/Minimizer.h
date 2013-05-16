#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "Solver/Solver.h"
#include "Solver/SolverMCBF.h"
#include "Solver/SolverMCIS.h"

using namespace std;
using namespace arma;

class Minimizer
{
public:
    Minimizer(
            int &myRank,
            int &numprocs,
            int &nParameters,
            Solver *solver);
//    virtual vec runMinimizer(const vec &guess, const int &nCycles) = 0;

protected:
    void printToFile(const vec &parameters, const Solver *solver);

    int myRank, numprocs;
    int nParameters;
    ofstream ofile;

    Solver *solver;
};

#endif // MINIMIZER_H
