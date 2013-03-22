#ifndef CWAVEFUNCTION_H
#define CWAVEFUNCTION_H

#include <armadillo>

using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction(const int &nParticles, const double &charge);

    double wavefunction(const mat &r);
    double getRatio();

    void initialize(mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setAlpha(const double &alpha_);
    void setBeta(const double &beta_);

    void acceptMove();
    void rejectMove();

    double localEnergyNumerical();
    double laplaceNumerical();
    mat gradientNumerical();
    double electronNucleusPotential();
    double electronElectronPotential();

    double alpha;
    double beta;
protected:
    void calculate_rij();
    void update_rij();
    void calculate_fij();
    void update_fij();
    double slaterRatio();
    double phiSD();
    double phiSD(const mat &r);
    double jastrowWF() const;
    double jastrowRatio();

    double hydrogenWF(const int &i, const vec3 &rvec);
    double phi1s(const vec3 &rvec);
    double phi2s(const vec3 &rvec);

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
};

#endif // CWAVEFUNCTION_H
