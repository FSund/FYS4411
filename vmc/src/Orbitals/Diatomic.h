#ifndef DIATOMIC_H
#define DIATOMIC_H

#include <armadillo>
#include <src/Orbitals/Orbitals.h>
#include <src/Orbitals/Hydrogenic.h>

using namespace std;
using namespace arma;

class Diatomic : public Orbitals
{
public:
    Diatomic();
    ~Diatomic();
    virtual void setAlpha(const double &newAlpha);
    virtual void setR(const double &dist);
    virtual double wavefunction(const rowvec &rvec, const int &qNum);
    virtual rowvec gradient(const rowvec &rvec, const int &qNum);
    virtual double laplacian(const rowvec &rvec, const int &qNum);
    virtual double alphaGradient(const rowvec &rvec, const int &qNum);

private:
    Hydrogenic *hydrogenic;
    double alpha;
    rowvec R;
};

#endif // DIATOMIC_H
