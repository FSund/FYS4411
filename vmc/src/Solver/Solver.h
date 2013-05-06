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

    double getEnergy() { return energy; }
    double getVariance() { return variance; }
    double getAcceptanceRate() { return acceptanceRate; }

    void setClosedform(const bool &closedForm_) { closedForm = closedForm_; }
protected:
    double gaussianDeviate(long *seed);
    void finalize();

    Wavefunction *wf;

    int nDimensions;
    int nParticles;
    int charge;

    mat rOld;
    mat rNew;

    long idum;

    int nCycles;
    int nAccepted;

    int myRank;
    int numprocs;
    int local_nCycles;

    double ratio;
    double acceptanceRate;
    double energy;
    double energySquared;
    double energySum;
    double energySquaredSum;
    double variance;
    double deltaE;

    bool closedForm;
};

#endif // SOLVER_H
