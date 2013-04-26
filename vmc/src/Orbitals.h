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
    double wavefunction(const vec &rvec, const int &qNum);
protected:
    double phi1s(const vec &rvec);
    double phi2s(const vec &rvec);
private:
    int nDimensions;
    double wfCurrent;
    double alpha;
    double r, arg;
};

#endif // ORBITALS_H
