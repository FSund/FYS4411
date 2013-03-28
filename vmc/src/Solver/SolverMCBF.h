#ifndef SOLVERMCBF_H
#define SOLVERMCBF_H

#include "../lib.h"
#include "../Wavefunction.h"
#include "Solver.h"
#include <mpi.h>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class SolverMCBF : public Solver
{
public:
    SolverMCBF(int &myRank, int &numprocs, int &nParticles, int &charge);
    double runMonteCarloIntegration(const int &nCycles_);
protected:
    double stepLength;
};
#endif // SOLVERMCBF_H
