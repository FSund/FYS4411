#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include <src/Jastrow.h>
#include <src/Slater.h>

using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction();
    Wavefunction(
            const int &nParticles,
            const double &charge,
            const string &orbitalType);
    ~Wavefunction();

    double getRatio();

    void initialize(mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setAlpha(const double &newAlpha);
    void setBeta(const double &newBeta);
    void setR(const double &dist);
    void setParameters(const vec &parameters);

    void acceptMove();
    void rejectMove();

    double localEnergy();
    double localEnergyNumerical();
    mat localGradient();
    mat localGradientNumerical(const mat &r);
    double localLaplacian();
    double localLaplacianNumerical(const mat &r);
    vec variationalDerivatives();

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

    // numerical gradient and laplacian
    double wfMinus, wfPlus, wfCurrent;
    double dfactor;
    double ddwavefunction;
    mat dwavefunction;
    mat rPlus, rMinus;
    // closed form gradient and laplacian
    mat grad;
    double lapl;
    vec varGrad;

    Jastrow* jastrow;
    Slater* slater;

    bool useJastrow;
};

#endif // WAVEFUNCTION_H
