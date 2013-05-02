#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
#include "Orbitals.h"

using namespace std;
using namespace arma;

class Slater
{
public:
    Slater(const int &nParticles);
    ~Slater();

    void initialize(const mat &r);
    void setAlpha(const double &newAlpha);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    double getRatio();
    void acceptMove();
    void rejectMove();

    double wavefunction(const mat &r);
    rowvec localGradient(const mat &r, const int &i);
    rowvec localGradientNumerical(const mat &r, const int &k, const double &h);
    double localLaplacian(const mat &r, const int &i);
    double localLaplacianNumerical(const mat &r, const int &k, const double &h);

    // versions using default argument r = rNew
    double wavefunction()
        { return wavefunction(rNew); }
    rowvec localGradient(const int &i)
        { return localGradient(rNew, i); }
    rowvec localGradientNumerical(const int &k, const double &h)
        { return localGradientNumerical(rNew, k, h); }
    double localLaplacian(const int &i)
        { return localLaplacian(rNew, i); }
    double localLaplacianNumerical(const int &k, const double &h)
        { return localLaplacianNumerical(rNew, k, h); }
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    int currentParticle;
    double ratio, ratioUP, ratioDOWN;
    bool UP;

    // for numerical gradient and laplacian
    double wfCurrent, wfMinus, wfPlus;
    double dfactor;
    mat rPlus, rMinus;
    rowvec dwavefunction;
    double ddwavefunction;

    // closed form gradient and laplacian
    rowvec grad;
    double lapl;

    int N;
    mat slaterUPold, slaterDOWNold;
    mat slaterUPnew, slaterDOWNnew;
    mat slaterUPinvOld, slaterDOWNinvOld;
    mat slaterUPinvNew, slaterDOWNinvNew;

    void updateSlater();
    void updateInverse();

    Orbitals* orbitals;
public:
    // debug stuff
    mat gradient(const double &h)
        { return gradient(rNew, h); }
    mat gradient(const mat &r, const double &h);

    mat getUPinvOld() { return slaterUPinvOld; }
    mat getDOWNinvOld() { return slaterDOWNinvOld; }
    mat getUPold() { return slaterUPold; }
    mat getDOWNold() { return slaterDOWNold; }

    mat getUPinvNew() { return slaterUPinvNew; }
    mat getDOWNinvNew() { return slaterDOWNinvNew; }
    mat getUPnew() { return slaterUPnew; }
    mat getDOWNnew() { return slaterDOWNnew; }
};

#endif // SLATER_H
