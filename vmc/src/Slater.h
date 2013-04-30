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

    void initalize(const mat &r);
    void setAlpha(const double &newAlpha);
    void updatePositionAndCurrentParticle(mat &r, int &k);

    double wavefunction() { return wavefunction(rNew); }
    double wavefunction(const mat &r);

    double getRatio();
    rowvec localGradient(const int &i);
    double localLaplacian(const int &i);

    void acceptMove();
    void rejectMove();
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew, rPlus, rMinus;
    int currentParticle;
    double ratio, ratioUP, ratioDOWN;

//    mat localGradientNumerical(const double &h);
//    double localLaplacianNumerical(const double &h);
//    double wfCurrent, wfMinus, wfPlus;
//    mat dwavefunction;
//    vec ddwavefunction;
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
    mat getUPinv() { return slaterUPinvOld; }
    mat getDOWNinv() { return slaterDOWNinvOld; }
};

#endif // SLATER_H
