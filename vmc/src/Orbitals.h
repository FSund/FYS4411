#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace std;
using namespace arma;

class Orbitals
{
public:
    Orbitals();
    void setAlpha(const double &newAlpha);
    double wavefunction(const rowvec &rvec, const int &qNum);
    rowvec gradient(const rowvec &rvec, const int &qNum);
    double laplacian(const rowvec &rvec, const int &qNum);
    double alphaGradient(const rowvec &rvec, const int &qNum);
protected:
//    double phi1s(const rowvec &rvec);
//    double phi2s(const rowvec &rvec);
//    double phi2p(const rowvec &rvec, const int &k);
//    rowvec dphi1s(const rowvec &rvec);
//    rowvec dphi2s(const rowvec &rvec);
//    rowvec dphi2p(const rowvec &rvec, const int &k);
    double ddphi1s(const rowvec &rvec);
    double ddphi2s(const rowvec &rvec);
    double ddphi2p(const rowvec &rvec, const int &k);
private:
    int nDimensions;
    double wfCurrent;
    double alpha;
    double r, arg;
    rowvec grad;
    double lapl;
    vec dphi;
};

#endif // ORBITALS_H
