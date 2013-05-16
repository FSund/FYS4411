#ifndef LOCALENERGY_H
#define LOCALENERGY_H

#include <armadillo>
#include <src/Wavefunction.h>

using namespace std;
using namespace arma;

class LocalEnergy
{
public:
    LocalEnergy();
    LocalEnergy(const int &nParticles, const int &nDimensions, const int &charge);
    virtual ~LocalEnergy();
    virtual void setR(const double &dist) = 0;
    virtual double evaluate(const mat &r, Wavefunction *wf) = 0;
protected:
    int nParticles;
    int nDimensions;
    double charge;
//    rowvec R;
    bool useJastrow;
    bool numerical;
};

#endif // LOCALENERGY_H
