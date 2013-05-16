#ifndef SOLVERMCIS_H
#define SOLVERMCIS_H

#include <armadillo>
#include <iostream>
#include <src/lib.h>
#include <src/Solver/Solver.h>
#include <src/Wavefunction.h>

using namespace std;
using namespace arma;

class SolverMCIS : public Solver
{
public:
    SolverMCIS();
    SolverMCIS(int &myRank, int &numprocs, int &nParticles, int &charge, string &orbitalType);
    virtual void runCycle();
    double testSolver(const int &nCycles_);
private:
    mat qForceOld;
    mat qForceNew;
    double D, dt, Ddt;
    double omegaRatio;
};

#endif // SOLVERMCIS_H
