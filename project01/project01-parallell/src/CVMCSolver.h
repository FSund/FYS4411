#ifndef CVMCSOLVER_H
#define CVMCSOLVER_H
#include <armadillo>
#include "mainapplication.h"

using namespace std;
using namespace arma;

const double pi = atan(1)*4;

class VMCSolver
{
public:
    VMCSolver();
    double runMonteCarloIntegration(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_, const bool closedform);
    double runMonteCarloIntegrationImportanceSampling(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_, const bool closedform, const double dt);
    void setParameters(const double &newStepLength, const double &newAlpha, const double &newBeta);
    ~VMCSolver();
    friend class MainApplication;
protected:
    mat rOld;
    mat rNew;

    double waveFunction(const mat &r);
//    double waveFunction2(const mat &r);
    double localEnergy(const mat &r);
    double exactLocalEnergy(const mat &r);
    double dr(const mat &r, const int ii, const int jj);
    mat quantumForce(const mat &r, const double &wf);

    double gaussianDeviate(long *seed);

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;
    int nAccepted;
    int nRejected;

    int my_rank;
    int numprocs;
    int local_nCycles;
};

#endif // CVMCSOLVER_H
