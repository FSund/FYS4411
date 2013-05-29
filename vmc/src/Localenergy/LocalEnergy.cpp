#include <src/Localenergy/LocalEnergy.h>

LocalEnergy::LocalEnergy()
{
    cout << "! Error: Using default constructor in LocalEnergy ! " << endl;
    exit(1);
}

LocalEnergy::LocalEnergy(const int &nParticles, const int &nDimensions, const int &charge):
    nParticles(nParticles),
    nDimensions(nDimensions),
    charge(charge),
    closedForm(true),
    useJastrow(true)
{
}

LocalEnergy::~LocalEnergy()
{
}
