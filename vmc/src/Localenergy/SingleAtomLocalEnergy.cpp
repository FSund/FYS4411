#include <src/Localenergy/SingleAtomLocalEnergy.h>

SingleAtomLocalEnergy::SingleAtomLocalEnergy()
{
    cout << "! Error: Using default constructor in SingleAtomLocalEnergy! " << endl;
    exit(1);
}

SingleAtomLocalEnergy::SingleAtomLocalEnergy(
        const int &nParticles,
        const int &nDimensions,
        const double &charge):
    LocalEnergy(nParticles, nDimensions, charge)
{
}

void SingleAtomLocalEnergy::setR(const double &dist)
{
    (void) dist;
}

double SingleAtomLocalEnergy::evaluate(const mat &r, Wavefunction *wf)
{
    if (numerical)
        return localEnergyNumerical(r, wf);
    else
        return localEnergy(r, wf);
}

double SingleAtomLocalEnergy::localEnergy(const mat &r, Wavefunction *wf)
{
    double kinetic, potential;

    kinetic = -0.5*wf->localLaplacian();
    if (useJastrow)
        potential = electronNucleusPotential(r) + electronElectronPotential(r);
    else
        potential = electronNucleusPotential(r);

    return kinetic + potential;
}

double SingleAtomLocalEnergy::localEnergyNumerical(const mat &r, Wavefunction *wf)
{
    double kinetic, potential;

    kinetic = -0.5*wf->localLaplacianNumerical(r);
    if (useJastrow)
        potential = electronNucleusPotential(r) + electronElectronPotential(r);
    else
        potential = electronNucleusPotential(r);

    return kinetic + potential;
}

double SingleAtomLocalEnergy::electronNucleusPotential(const mat &rNew)
{
    // potential energy
    double potentialEnergy = 0.0;
    double r;
    for (int i = 0; i < nParticles; i++) {
        r = 0.0;
        for (int j = 0; j < nDimensions; j++)
            r += rNew(i,j)*rNew(i,j);
        r = sqrt(r);
        potentialEnergy += 1.0/r;
    }
    potentialEnergy *= -charge;

    return potentialEnergy;
}

double SingleAtomLocalEnergy::electronElectronPotential(const mat &rNew)
{
    // contribution from electron-electron potential
    double potentialEnergy = 0.0;
    double rij = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            rij = 0;
            for (int k = 0; k < nDimensions; k++) {
                rij += (rNew(i,k) - rNew(j,k))*(rNew(i,k) - rNew(j,k));
            }
            potentialEnergy += 1.0/sqrt(rij);
        }
    }

    return potentialEnergy;
}



