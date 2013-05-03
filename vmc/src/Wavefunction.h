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

    double localEnergy();
    double localEnergyNumerical();
    mat localGradient();
//    mat localGradient(const mat &r);
    mat localGradientNumerical(const mat &r);
    double localLaplacian();
//    double localLaplacian(const mat &r);
    double localLaplacianNumerical(const mat &r);

    //    double localEnergyClosedForm(const mat &r) const;

    double wavefunction();
    double wavefunction(const mat &r);

    // default arguments
    mat localGradientNumerical() { return localGradientNumerical(rNew); }
    double localLaplacianNumerical() { return localLaplacianNumerical(rNew); }
protected:
    double electronNucleusPotential();
    double electronElectronPotential();

    int nParticles, nDimensions;
    double charge;
    double h, h2;

    int currentParticle;
private:
    mat rNew, rOld;
    mat rijNew, rijOld;
    mat fijNew, fijOld;
    double ratio;

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
