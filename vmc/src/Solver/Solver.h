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
#include <src/OneBodyDensity.h>

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

    void setBlocking(const bool &blocking_) { blocking = blocking_; }
    void setMinimizing(const bool &minimizing_) { minimizing = minimizing_; }
    void setOnebody(const bool &onebody_) { onebody = onebody_; }
    void setThermalizationSteps(const int &steps) { nThermalize = steps; }
    void setClosedform(const bool &closedForm) {
        wf->setClosedForm(closedForm);
        localEnergy->setClosedForm(closedForm);
    }
    void setUseJastrow(const bool &useJastrow) {
        wf->setUseJastrow(useJastrow);
        localEnergy->setUseJastrow(useJastrow);
    }
protected:
    double gaussianDeviate(long *seed);
    void finalize();

    Wavefunction *wf;
    Datalogger *logger;
    LocalEnergy *localEnergy;
    OneBodyLogger *onebodylogger;

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

    bool blocking;
    bool minimizing;
    bool onebody;
};

#endif // SOLVER_H
