#ifndef DIATOMICLOCALENERGY_H
#define DIATOMICLOCALENERGY_H

#include <armadillo>
#include <src/Localenergy/LocalEnergy.h>

using namespace std;
using namespace arma;

class DiatomicLocalEnergy : public LocalEnergy
{
public:
    DiatomicLocalEnergy();
    DiatomicLocalEnergy(
            const int &nParticles,
            const int &nDimensions,
            const double &charge);
    virtual void setR(const double &dist);
    virtual double evaluate(const mat &r, Wavefunction *wf);
private:
    double electronElectronPotential(const mat &r);
    rowvec R;
};

#endif // DIATOMICLOCALENERGY_H
