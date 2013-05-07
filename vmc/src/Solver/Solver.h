#ifndef SOLVER_H
#define SOLVER_H

#include <mpi.h>
#include <armadillo>
#include <iostream>
#include <src/lib.h>
#include <src/Wavefunction.h>
#include <src/Datalogger.h>

using namespace std;
using namespace arma;

class Solver
{
public:
    Solver(int &myRank, int &numprocs, int &nParticles, int &charge);
    ~Solver();

    double runMonteCarloIntegration(const int &nCycles_);
    virtual void runCycle() = 0;

    void setAlpha(const double &alpha);
    void setBeta(const double &beta);
    void setParameters(const vec &parameters);

    double getEnergy() { return energy; }
    double getVariance() { return variance; }
    double getAcceptanceRate() { return acceptanceRate; }

    void setClosedform(const bool &closedForm_) { closedForm = closedForm_; }
    void setBlocking(const bool &blocking_) { blocking = blocking_; }
protected:
    double gaussianDeviate(long *seed);
    void finalize();

    Wavefunction *wf;
    Datalogger *logger;

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
    int nThermalize;

    double ratio;
    double acceptanceRate;
    double energy;
    double energySquared;
    double energySum;
    double energySquaredSum;
    double variance;
    double deltaE;

    bool closedForm;
    bool blocking;
};

#endif // SOLVER_H
