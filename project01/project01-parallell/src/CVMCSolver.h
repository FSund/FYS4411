#ifndef CVMCSOLVER_H
#define CVMCSOLVER_H

#include "lib.h"
#include <mpi.h>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

const double pi = atan(1)*4;

class VMCSolver
{
public:
    VMCSolver(int my_rank, int numprocs, const int &charge, const int &nParticles);
    ~VMCSolver();

//    double runMonteCarloIntegration(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_, const bool closedform);
    double runMonteCarloIntegration(const int &nCycles, const bool &closedForm);
//    double runMonteCarloIntegrationImportanceSampling(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_, const bool closedform, const double dt);
    double runMonteCarloIntegrationImportanceSampling(const int &nCycles, const double &dt, const bool &closedform);
protected:
    mat rOld;
    mat rNew;

    virtual double wavefunction(const mat &r)=0;
//    virtual double localEnergy(const mat &r)=0;
    virtual double localEnergyClosedForm(const mat &r)=0;

    double localEnergy(const mat &r);
    mat quantumForce(const mat &r, const double &wf);
//    double waveFunction2(const mat &r);
//    double dr(const mat &r, const int ii, const int jj);

    double gaussianDeviate(long *seed);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    int nCycles;
    int nAccepted;
    int nRejected;

    int my_rank;
    int numprocs;
    int local_nCycles;
};

#endif // CVMCSOLVER_H
