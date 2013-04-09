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

    double wavefunction(const mat &r);
    double getRatio();
    mat gradient();
    double laplacian();

    void acceptMove();
    void rejectMove();

//    double hydrogenWF(const int &i, const vec3 &rvec);
//    double phi1s(const vec3 &rvec);
//    double phi2s(const vec3 &rvec);
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    int currentParticle;
    double ratio, ratioUP, ratioDOWN;

    int N;
    mat slaterUPold, slaterDOWNold;
    mat slaterUPnew, slaterDOWNnew;
    mat slaterUPinvOld, slaterDOWNinvOld;
    mat slaterUPinvNew, slaterDOWNinvNew;

//    void update();
    void updateSlater();
    void updateInverse();

    Orbitals* orbitals;
};

#endif // SLATER_H
