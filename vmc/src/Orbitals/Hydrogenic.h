#ifndef HYDROGENIC_H
#define HYDROGENIC_H

#include <armadillo>
#include <src/Orbitals/Orbitals.h>

using namespace std;
using namespace arma;

class Hydrogenic : public Orbitals
{
public:
    Hydrogenic();
    virtual void setAlpha(const double &newAlpha);
    virtual void setR(const double &dist);
    virtual double wavefunction(const rowvec &rvec, const int &qNum);
    virtual rowvec gradient(const rowvec &rvec, const int &qNum);
    virtual double laplacian(const rowvec &rvec, const int &qNum);
    virtual double alphaGradient(const rowvec &rvec, const int &qNum);

private:
//    double phi1s(const rowvec &rvec);
//    double phi2s(const rowvec &rvec);
//    double phi2p(const rowvec &rvec, const int &k);
//    rowvec dphi1s(const rowvec &rvec);
//    rowvec dphi2s(const rowvec &rvec);
//    rowvec dphi2p(const rowvec &rvec, const int &k);
    double ddphi1s(const rowvec &rvec);
    double ddphi2s(const rowvec &rvec);
    double ddphi2p(const rowvec &rvec, const int &k);
    double alpha;

//private:
//    int nDimensions;
//    double wfCurrent;
//    double alpha;
//    double r, arg;
//    rowvec grad;
//    double lapl;
//    vec dphi;
};

#endif // HYDROGENIC_H
