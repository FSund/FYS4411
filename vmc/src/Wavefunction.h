#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include "Jastrow.h"
#include "Slater.h"

using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction(const int &nParticles, const double &charge);
    ~Wavefunction();

    double getRatio();

    void initialize(mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setAlpha(const double &newAlpha);
    void setBeta(const double &newBeta);
    void setParameters(const vec &parameters);

    void acceptMove();
    void rejectMove();

    double localEnergyNumerical();
    mat localGradientNumerical();
    double localLaplacianNumerical();

    mat localGradient();
    double localLaplacian();
    //    double localEnergyClosedForm(const mat &r) const;
protected:
    double wavefunction();
    double wavefunction(const mat &r);
    double electronNucleusPotential();
    double electronElectronPotential();

    int nParticles, nDimensions;
    double charge;
    double h, h2;

    int currentParticle;

    mat rNew, rOld;
    mat rijNew, rijOld;
    mat fijNew, fijOld;

private:
    double wfMinus, wfPlus, wfCurrent;
    double dfactor;
    double ddwavefunction;
    mat dwavefunction;
    mat rPlus, rMinus;
    mat grad;
    double lapl;

    Jastrow* jastrow;
    Slater* slater;
};

#endif // WAVEFUNCTION_H
