#ifndef JASTROW_H
#define JASTROW_H

#include <armadillo>

using namespace std;
using namespace arma;

class Jastrow
{
public:
    Jastrow(const int &nParticles);
    void initialize(const mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setBeta(const double &newBeta);
    double getRatio();
    void acceptMove();
    void rejectMove();

    double wavefunction();
    double wavefunction(const mat &r);
    rowvec localGradient(const int &k);
    double localLaplacian(const int &k);
    double betaGradient(const int &k);

    // for debugging/tests
    rowvec localGradientNumerical(const mat &r, const int &k, const double &h);
    double localLaplacianNumerical(const mat &r, const int &k, const double &h);

    // versions using default argument r = rNew
    rowvec localGradientNumerical(const int &k, const double &h)
        { return localGradientNumerical(rNew, k, h); }
    double localLaplacianNumerical(const int &k, const double &h)
        { return localLaplacianNumerical(rNew, k, h); }
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    mat rijOld, rijNew;
    mat fijOld, fijNew;
    mat a;
    double beta;
    int currentParticle;

    // numerical gradient and laplacian
    double wfCurrent, wfMinus, wfPlus;
    mat rPlus, rMinus;
    rowvec dwavefunction;
    double ddwavefunction;

    // closed form gradient and laplacian
    double rkij, rkj, rki, brki, brkj, rij;
    rowvec grad;
    double lapl;

    void calculate_rij();
    void update_rij();
    void calculate_fij();
    void update_fij();
};

#endif // JASTROW_H
