#ifndef SINGLEATOMLOCALENERGY_H
#define SINGLEATOMLOCALENERGY_H

#include <armadillo>
#include <src/Localenergy/LocalEnergy.h>

using namespace std;
using namespace arma;

class SingleAtomLocalEnergy : public LocalEnergy
{
public:
    SingleAtomLocalEnergy();
    SingleAtomLocalEnergy(
            const int &nParticles,
            const int &nDimensions,
            const double &charge);
    virtual void setR(const double &dist);
    virtual double evaluate(const mat &r, Wavefunction *wf);
private:
    double localEnergy(const mat &r, Wavefunction *wf);
    double localEnergyNumerical(const mat &r, Wavefunction *wf);
    double electronNucleusPotential(const mat &rNew);
    double electronElectronPotential(const mat &rNew);
};

#endif // SINGLEATOMLOCALENERGY_H
