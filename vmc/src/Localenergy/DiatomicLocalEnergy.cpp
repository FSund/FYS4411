#include <src/Localenergy/DiatomicLocalEnergy.h>

DiatomicLocalEnergy::DiatomicLocalEnergy()
{
    cout << "! Error: Using default constructor in DiatomicLocalEnergy! " << endl;
    exit(1);
}

DiatomicLocalEnergy::DiatomicLocalEnergy(
        const int &nParticles,
        const int &nDimensions,
        const double &charge):
    LocalEnergy(nParticles, nDimensions, charge),
    R(zeros<rowvec>(nDimensions))
{
}

void DiatomicLocalEnergy::setR(const double &dist)
{
    R(0) = dist/2.0;

//    cout << "R local energy = " << R;
}

double DiatomicLocalEnergy::evaluate(const mat &r, Wavefunction *wf)
{
    double kinetic = -0.5*wf->localLaplacian();
    double potential = 0.0;
    double rp1, rp2, r12;
    rp1 = rp2 = r12 = 0.0;

    // contribution from the electron-nucleus potential
    for (int i = 0; i < nParticles; i++)
    {
        rp1 = 0.0;
        rp2 = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rp1 += (r(i,j) + R(j))*(r(i,j) + R(j));
            rp2 += (r(i,j) - R(j))*(r(i,j) - R(j));
        }
        potential += -charge*(1.0/sqrt(rp1) + 1.0/sqrt(rp2));
    }

    ////
//    cout << "rNew = " << endl << r;
    cout << "potential energy PE = " << potential << endl;
    cout << "potential EE        = " << electronElectronPotential(r) << endl;
    cout << "potential energy PP = " << charge*charge/(2.0*R(0)) << endl;
    cout << "kinetic energy = " << kinetic << endl;
    ////


    potential += electronElectronPotential(r);
    potential += charge*charge/(2.0*abs(R(0)));

    cout << "wf = " << wf->wavefunction() << endl;

    return kinetic + potential;
}

double DiatomicLocalEnergy::electronElectronPotential(const mat &r)
{
    // contribution from electron-electron potential
    double potentialEnergy = 0.0;
    double rij = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            rij = 0;
            for (int k = 0; k < nDimensions; k++) {
                rij += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            potentialEnergy += 1.0/sqrt(rij);
        }
    }

    return potentialEnergy;
}


