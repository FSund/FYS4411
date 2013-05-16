#ifndef SOLVER_H
#define SOLVER_H

#include <mpi.h>
#include <armadillo>
#include <iostream>
#include <src/lib.h>
#include <src/Wavefunction.h>
#include <src/Datalogger.h>
#include <src/Localenergy/LocalEnergy.h>
#include <src/Localenergy/SingleAtomLocalEnergy.h>
#include <src/Localenergy/DiatomicLocalEnergy.h>

using namespace std;
using namespace arma;

class Solver
{
public:
    Solver();
    Solver(const int &myRank,
           const int &numprocs,
           const int &nParticles,
           const int &charge,
           const string &orbitalType);
    ~Solver();

    double runMonteCarloIntegration(const int &nCycles_);
    virtual void runCycle() = 0;

    void setAlpha(const double &alpha);
    void setBeta(const double &beta);
    void setR(const double &dist);
    void setParameters(const vec &parameters);

    const double& getEnergy() const { return energy; }
    const double& getVariance() const { return variance; }
    const double& getAcceptanceRate() const { return acceptanceRate; }
    const vec& getVariationalGradient() const { return variationalGradient; }

    void setClosedform(const bool &closedForm_) { closedForm = closedForm_; }
    void setBlocking(const bool &blocking_) { blocking = blocking_; }
    void setMinimizing(const bool &minimizing_) { minimizing = minimizing_; }
protected:
    double gaussianDeviate(long *seed);
    void finalize();

    Wavefunction *wf;
    Datalogger *logger;
    LocalEnergy *localEnergy;

    int nDimensions;
    int nParticles;
    int charge;
    int nParameters;

    mat rOld;
    mat rNew;
    vec tempVariationalGradient;
    vec variationalGradient;
    vec variationalGradientSum;
    vec variationalGradientESum;
//    double* variationalGradient;
//    double* variationalGradientSum;
//    double* variationalGradientESum;

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
    double r12Sum, r12;

    bool closedForm;
    bool blocking;
    bool minimizing;
};

#endif // SOLVER_H
