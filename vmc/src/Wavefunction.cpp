#include <src/Wavefunction.h>

Wavefunction::Wavefunction(const int &nParticles, const double &charge):
    nParticles(nParticles),
    nDimensions(3),
    charge(charge),
    h(1e-3),
    h2(1e6)
{
    dwavefunction = zeros<mat>(nParticles, nDimensions);
    rPlus = rMinus = zeros<mat>(nParticles, nDimensions);
    rNew = rOld = zeros(nParticles, nDimensions);

    jastrow = new Jastrow(nParticles);
    slater = new Slater(nParticles);

//    useJastrow = false;
    useJastrow = true;
}

Wavefunction::~Wavefunction()
{
    delete jastrow;
    delete slater;
}

void Wavefunction::initialize(mat &r)
{
    rNew = rOld = r;

    jastrow->initialize(r);
    slater->initialize(r);
}

void Wavefunction::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;

    slater->updatePositionAndCurrentParticle(r, k);
    jastrow->updatePositionAndCurrentParticle(r, k);
}

void Wavefunction::setAlpha(const double &newAlpha)
{
    slater->setAlpha(newAlpha);
}

void Wavefunction::setBeta(const double &newBeta)
{
    jastrow->setBeta(newBeta);
}

void Wavefunction::setParameters(const vec &parameters)
{
    slater->setAlpha(parameters(0));
    jastrow->setBeta(parameters(1));
}

double Wavefunction::wavefunction()
{
    if (useJastrow)
        return slater->wavefunction()*jastrow->wavefunction();
    else
        return slater->wavefunction();
}

double Wavefunction::wavefunction(const mat &r)
{
    if (useJastrow)
        return slater->wavefunction(r)*jastrow->wavefunction(r);
    else
        return slater->wavefunction(r);
}

double Wavefunction::getRatio()
{
    if (useJastrow)
        ratio = slater->getRatio()*jastrow->getRatio();
    else
        ratio = slater->getRatio();

    return ratio; // squared in slater and jastrow
}

void Wavefunction::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    jastrow->acceptMove();
    slater->acceptMove();
}

void Wavefunction::rejectMove()
{
    rNew.row(currentParticle) = rOld.row(currentParticle);

    jastrow->rejectMove();
    slater->rejectMove();
}

double Wavefunction::localEnergy()
{
    double kinetic, potential;

    kinetic = -0.5*localLaplacian();
    if (useJastrow)
        potential = electronNucleusPotential() + electronElectronPotential();
    else
        potential = electronNucleusPotential();

    return kinetic + potential;
}

double Wavefunction::localEnergyNumerical()
{
    double kinetic, potential;

    kinetic = -0.5*localLaplacianNumerical();
    if (useJastrow)
        potential = electronNucleusPotential() + electronElectronPotential();
    else
        potential = electronNucleusPotential();

    return kinetic + potential;
}

mat Wavefunction::localGradient()
{
    grad = zeros<mat>(nParticles, nDimensions);
    if (useJastrow)
        for (int i = 0; i < nParticles; i++)
            grad.row(i) = slater->localGradient(i) + jastrow->localGradient(i);
    else
        for (int i = 0; i < nParticles; i++)
            grad.row(i) = slater->localGradient(i);

    return grad;
}

mat Wavefunction::localGradientNumerical(const mat &r)
{
    // this function is correct!
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    dfactor = 1.0/(wfCurrent*2.0*h);
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            dwavefunction(i,j) = (wfPlus - wfMinus)*dfactor;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }

    return dwavefunction;
}

double Wavefunction::localLaplacian()
{
    lapl = 0.0;

    if (useJastrow)
        for (int i = 0; i < nParticles; i++)
        {
            lapl += slater->localLaplacian(i) + jastrow->localLaplacian(i)
                    + 2.0*dot(slater->localGradient(i), jastrow->localGradient(i));
        }
    else
        for (int i = 0; i < nParticles; i++)
            lapl += slater->localLaplacian(i);

    return lapl;
}

double Wavefunction::localLaplacianNumerical(const mat &r)
{
    // this function is correct!
    ddwavefunction = 0.0;
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            ddwavefunction += wfMinus + wfPlus;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }

    ddwavefunction = h2*(ddwavefunction/wfCurrent - 2.0*double(nParticles*nDimensions));
//    ddwavefunction = h2*(ddwavefunction - 2.0*double(nParticles*nDimensions)*wfCurrent); // not-local

    return ddwavefunction;
}

double Wavefunction::electronNucleusPotential()
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

double Wavefunction::electronElectronPotential()
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
            potentialEnergy += 1/sqrt(rij);
        }
    }

    return potentialEnergy;
}
