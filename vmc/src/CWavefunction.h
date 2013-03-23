#ifndef CWAVEFUNCTION_H
#define CWAVEFUNCTION_H

#include <armadillo>
#include "Jastrow.h"
#include "Slater.h"

using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction(const int &nParticles, const double &charge);

    double getRatio();

    void initialize(mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setAlpha(const double &newAlpha);
    void setBeta(const double &newBeta);

    void acceptMove();
    void rejectMove();

    double localEnergyNumerical();
    double laplaceNumerical();
    mat gradientNumerical();
protected:
    double wavefunction(const mat &r);
    double electronNucleusPotential();
    double electronElectronPotential();

    int nParticles, nDimensions;
    double charge;
    double h, h2;

    int currentParticle;

//    double wfOld, wfNew;

    mat rNew, rOld;
    mat rijNew, rijOld;
    mat fijNew, fijOld;

private:
    double wfMinus, wfPlus, wfCurrent;
    mat dwavefunction;
    mat rPlus, rMinus;

    Jastrow* jastrow;
    Slater* slater;
};

#endif // CWAVEFUNCTION_H
