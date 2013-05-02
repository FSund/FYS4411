#include "Slater.h"

Slater::Slater(const int &nParticles):
    nParticles(nParticles),
    nDimensions(3),
    rOld(zeros<mat>(nParticles, nDimensions)),
    rNew(zeros<mat>(nParticles, nDimensions)),
    N(nParticles/2),
    slaterUPold(zeros<mat>(N, N)),
    slaterDOWNold(zeros<mat>(N, N)),
    slaterUPnew(zeros<mat>(N, N)),
    slaterDOWNnew(zeros<mat>(N, N)),
    slaterUPinvOld(zeros<mat>(N, N)),
    slaterDOWNinvOld(zeros<mat>(N, N)),
    slaterUPinvNew(zeros<mat>(N, N)),
    slaterDOWNinvNew(zeros<mat>(N, N))
{
    orbitals = new Orbitals;
}

Slater::~Slater()
{
    delete orbitals;
}

void Slater::initialize(const mat &r)
{
    rNew = rOld = r;

    for (int i = 0; i < N; i++) // loop over particles
    {
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUPold(i,j)   = orbitals->wavefunction(rOld.row(i),   j);
            slaterDOWNold(i,j) = orbitals->wavefunction(rOld.row(i+N), j);
        }
    }
    slaterUPinvOld   = inv(slaterUPold);
    slaterDOWNinvOld = inv(slaterDOWNold);

    slaterUPnew      = slaterUPold;
    slaterDOWNnew    = slaterDOWNold;
    slaterUPinvNew   = slaterUPinvOld;
    slaterDOWNinvNew = slaterDOWNinvOld;

    ratioUP = ratioDOWN = 1.0;
}

void Slater::setAlpha(const double &newAlpha)
{
    orbitals->setAlpha(newAlpha);
}

void Slater::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;
    if (k < N) UP = true;
    else UP = false;

    // only updating the inverse slater matrix _if_ we accept the move, because
    // we only need the old inverse to find the ratio
    updateSlater();
}

double Slater::getRatio()
{
    if (UP)
    {
        int i = currentParticle;
        ratioUP = 0.0;
        for (int j = 0; j < N; j++) // loop over orbitals
            ratioUP += slaterUPnew(i,j)*slaterUPinvOld(j,i);

//        ratioUP *= ratioUP; // bad because we need not squared in inverse update
        return ratioUP*ratioUP;
    }
    else
    {
        int i = currentParticle - N;
        ratioDOWN = 0.0;
        for (int j = 0; j < N; j++) // loop over orbitals
            ratioDOWN += slaterDOWNnew(i,j)*slaterDOWNinvOld(j,i);

//        ratioDOWN *= ratioDOWN;
        return ratioDOWN*ratioDOWN;
    }
}

void Slater::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    updateInverse();
    if (UP)
    {
        slaterUPold.row(currentParticle) = slaterUPnew.row(currentParticle);
        slaterUPinvOld = slaterUPinvNew; // can't use just column/row!
    }
    else
    {
        slaterDOWNold.row(currentParticle-N) = slaterDOWNnew.row(currentParticle-N);
        slaterDOWNinvOld = slaterDOWNinvNew; // can't use just column/row!
    }
}

void Slater::rejectMove()
{
    rNew.row(currentParticle) = rOld.row(currentParticle);

    // not updating the inverse, since we haven't changed it!!!
    if (UP)
        slaterUPnew.row(currentParticle) = slaterUPold.row(currentParticle);
    else
        slaterDOWNnew.row(currentParticle-N) = slaterDOWNold.row(currentParticle-N);
}

double Slater::wavefunction(const mat &r)
{
    /* Only for use in numerical derivatives */

    mat slaterUP(N, N);
    mat slaterDOWN(N, N);

    for (int i = 0; i < N; i++) // loop over particles
    {
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUP(i,j)   = orbitals->wavefunction(r.row(i), j);
            slaterDOWN(i,j) = orbitals->wavefunction(r.row(i+N), j);
        }
    }

    return det(slaterUP)*det(slaterDOWN);
}

rowvec Slater::localGradient(const mat &r, const int &i)
{
    /* Calculates the gradient for particle k */
    // NEED TO USE 1/R !!!

    grad = zeros<rowvec>(nDimensions);
    if (i < N)
        for (int j = 0; j < N; j++)
            grad += orbitals->gradient(r.row(i),j)*slaterUPinvOld(j,i)/ratioUP;
    else
        for (int j = 0; j < N; j++)
            grad += orbitals->gradient(r.row(i),j)*slaterDOWNinvOld(j,i-N)/ratioDOWN;

    return grad;
}

rowvec Slater::localGradientNumerical(const mat &r, const int &k, const double &h)
{
    /* Calculates the the gradient for particle k */

    dwavefunction = zeros<rowvec>(nDimensions);
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nDimensions; i++)
    {
        rPlus(k,i) += h;
        rMinus(k,i) -= h;
        wfMinus = wavefunction(rMinus);
        wfPlus = wavefunction(rPlus);
        dwavefunction(i) = wfPlus - wfMinus;
        rPlus(k,i) = rMinus(k,i) = r(k,i);
    }
    dwavefunction /= (2.0*wfCurrent*h);
//    dwavefunction /= (2.0*h); // not local

    return dwavefunction;
}

double Slater::localLaplacian(const mat &r, const int &i)
{
    /* Calculates the laplacian for particle k */

    lapl = 0.0;
    if (i < N)
        for (int j = 0; j < N; j++)
            lapl += orbitals->laplacian(r.row(i),j)*slaterUPinvOld(j,i); // OLD
    else
        for (int j = 0; j < N; j++)
            lapl += orbitals->laplacian(r.row(i),j)*slaterDOWNinvOld(j,i-N); // OLD

    return lapl;
}

double Slater::localLaplacianNumerical(const mat &r, const int &k, const double &h)
{
    /* Calculates the laplacian for particle k */
    ddwavefunction = 0.0;
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nDimensions; i++) {
        rPlus(k,i) += h;
        rMinus(k,i) -= h;
        wfMinus = wavefunction(rMinus);
        wfPlus = wavefunction(rPlus);
        ddwavefunction += wfMinus + wfPlus - 2.0*wfCurrent;
        rPlus(k,i) = rMinus(k,i) = r(k,i);
    }
    ddwavefunction /= (wfCurrent*h*h);

    return ddwavefunction;
}

void Slater::updateSlater()
{
    int i = currentParticle;
    if (UP)
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUPnew(i,j) = orbitals->wavefunction(rNew.row(i), j);
        }
    else
        for (int j = 0; j < N; j++) // loop over orbitals
            slaterDOWNnew(i-N,j) = orbitals->wavefunction(rNew.row(i), j);
}

void Slater::updateInverse()
{
//    getRatio();

    if (UP) // only updating the UP matrix
    {
        int i = currentParticle;

        vec SUP = zeros<vec>(N);
        // updating all but the ith column of the inverse slater matrix
        for (int j = 0; j < N; j++) // loop over columns j != i
        {
            if (j == i) continue; // not for column i
            for (int l = 0; l < N; l++) // loop over elements in column j
            {
                SUP(j) += slaterUPnew(i,l)*slaterUPinvOld(l,j);
            }
        }
        for (int j = 0; j < N; j++) // loop over columns in the inverse slater matrix
        {
            if (j == i) continue; // not updating column i
            for (int k = 0; k < N; k++) // loop over elements in column j
            {
                // maybe we should recalculate the ratio first to be sure??
                slaterUPinvNew(k,j) = slaterUPinvOld(k,j) - SUP(j)*slaterUPinvOld(k,i)/ratioUP;
            }
        }
        // updating column i of the inverse slater matrix
        for (int k = 0; k < N; k++) // loop over elements in i'th column
        {
            slaterUPinvNew(k,i) = slaterUPinvOld(k,i)/ratioUP;
        }
    }
    else // only updating the DOWN matrix
    {
        int i = currentParticle - N;

        vec SDOWN = zeros<vec>(N);
        // updating all but the ith column of the inverse slater matrix
        for (int j = 0; j < N; j++) // loop over columns j != i
        {
            if (j == i) continue; // not for column i
            for (int l = 0; l < N; l++) // loop over elements in column j
            {
                SDOWN(j) += slaterDOWNnew(i,l)*slaterDOWNinvOld(l,j);
            }
        }
        for (int j = 0; j < N; j++) // loop over columns in the inverse slater matrix
        {
            if (j == i) continue; // not updating column i
            for (int k = 0; k < N; k++) // loop over elements in column j
            {
                // maybe we should recalculate the ratio first to be sure??
                slaterDOWNinvNew(k,j) = slaterDOWNinvOld(k,j) - SDOWN(j)*slaterDOWNinvOld(k,i)/ratioDOWN;
            }
        }
        // updating column i of the inverse slater matrix
        for (int k = 0; k < N; k++) // loop over elements in i'th column
        {
            slaterDOWNinvNew(k,i) = slaterDOWNinvOld(k,i)/ratioDOWN;
        }
    }
}

// debug stuff
mat Slater::gradient(const mat &r, const double &h)
{
    rPlus = rMinus = r;
    mat temp(nParticles, nDimensions);
    wfCurrent = wavefunction(r);
    dfactor = 1.0/(wfCurrent*2.0*h);
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            temp(i,j) = (wfPlus - wfMinus)*dfactor;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }

    return temp;
}
