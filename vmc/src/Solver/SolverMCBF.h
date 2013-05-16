#ifndef SOLVERMCBF_H
#define SOLVERMCBF_H

#include <armadillo>
#include <iostream>
#include <src/lib.h>
#include <src/Solver/Solver.h>
#include <src/Wavefunction.h>

using namespace std;
using namespace arma;

class SolverMCBF : public Solver
{
public:
    SolverMCBF();
    SolverMCBF(int &myRank, int &numprocs, int &nParticles, int &charge, string &orbitalType);
    virtual void runCycle();
protected:
    double stepLength;
};
#endif // SOLVERMCBF_H
