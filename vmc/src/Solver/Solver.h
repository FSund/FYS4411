#ifndef SOLVER_H
#define SOLVER_H

#include "../lib.h"
#include "../Wavefunction.h"
#include <mpi.h>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class Solver
{
public:
    Solver(int &myRank, int &numprocs, int &nParticles, int &charge);
    ~Solver();

    virtual double runMonteCarloIntegration(const int &nCycles_) = 0;

    void setAlpha(const double &alpha);
    void setBeta(const double &beta);
    void setParameters(const vec &parameters);
protected:
    double gaussianDeviate(long *seed);

    Wavefunction *wf;

    int nDimensions;
    int nParticles;
    int charge;

    mat rOld;
    mat rNew;

    long idum;

    int nCycles;
    int nAccepted;
    int nRejected;

    int myRank;
    int numprocs;
    int local_nCycles;

    double ratio;
};

#endif // SOLVER_H
