#ifndef SOLVERMCIS_H
#define SOLVERMCIS_H

#include "Solver.h"
#include "../lib.h"
#include "../Wavefunction.h"
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class SolverMCIS : public Solver
{
public:
    SolverMCIS(int &myRank, int &numprocs, int &nParticles, int &charge);
    double runMonteCarloIntegration(const int &nCycles_);
    double gaussianDeviate(long *seed);
private:
    mat qForceOld;
    mat qForceNew;
    double D, dt, Ddt;
    double omegaRatio;
};

#endif // SOLVERMCIS_H
