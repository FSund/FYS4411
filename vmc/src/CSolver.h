#ifndef CSOLVER_H
#define CSOLVER_H

#include "lib.h"
#include "CHelium.h"
#include "CBeryllium.h"
#include "CWavefunction.h"
#include <mpi.h>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

const double pi = atan(1)*4;

class Solver
{
public:
    Solver(int myRank, int numprocs, int nParticles, int charge);
    ~Solver();

    double runMonteCarloIntegration(const int &nCycles_);

//    void setParameters(const double &stepLength_);
//    void setParameters(const double &stepLength_,
//            const bool &importanceSampling_,
//            const bool &closedForm_);

    void setAlpha(const double &alpha_);
    void setBeta(const double &beta_);
protected:
    Wavefunction *wf;

    mat rOld;
    mat rNew;
    mat qForceOld;
    mat qForceNew;
    double wfOld;
    double wfNew;

    double gaussianDeviate(long *seed);

    int nDimensions;
    int nParticles;
    int charge;
    double alpha;
    double beta;
    double stepLength;

//    double h;
//    double h2;
    double D;
    double Ddt;

    long idum;

    int nCycles;
    int nAccepted;
    int nRejected;

    int myRank;
    int numprocs;
    int local_nCycles;
public:
    bool importanceSampling;
    bool closedForm;
};

#endif // CSOLVER_H
