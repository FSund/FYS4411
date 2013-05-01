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

    double wavefunction() { return wavefunction(rNew); }
    double wavefunction(const mat &r);

    double getRatio();
    rowvec localGradient(const int &i);
    rowvec localGradient(const mat &r, const int &i);
    rowvec localGradientNumerical(const int &k, const double &h);
    rowvec localGradientNumerical(const mat &r, const int &k, const double &h);
    double localLaplacian(const int &i);
    double localLaplacian(const mat &r, const int &i);
    double localLaplacianNumerical(const int &k, const double &h);
    double localLaplacianNumerical(const mat &r, const int &k, const double &h);

    void acceptMove();
    void rejectMove();

    // debug stuff
    mat gradient(const double &h);
    mat gradient(const mat &r, const double &h);
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    int currentParticle;
    double ratio, ratioUP, ratioDOWN;

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
