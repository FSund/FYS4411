#ifndef CVMCSOLVER_H
#define CVMCSOLVER_H
#include <armadillo>

using namespace std;
using namespace arma;

class VMCSolver
{
public:
    VMCSolver();
    double runMonteCarloIntegration(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_);
    void setParameters(const double &newStepLength, const double &newAlpha, const double &newBeta);
private:
    mat rOld;
    mat rNew;

    double waveFunction(const mat &r);
    double waveFunction2(const mat &r);
    double localEnergy(const mat &r);
    double exactLocalEnergy(const mat &r);
    double dr(const mat &r);

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
};

#endif // CVMCSOLVER_H
