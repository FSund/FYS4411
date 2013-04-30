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

void Slater::initalize(const mat &r)
{
    rNew = rOld = r;

    for (int i = 0; i < N; i++) // loop over particles
    {
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUPold(i,j)   = orbitals->wavefunction(rOld.row(i), j);
            slaterDOWNold(i,j) = orbitals->wavefunction(rOld.row(i+N), j);
        }
    }
    slaterUPinvOld = inv(slaterUPold);
    slaterDOWNinvOld = inv(slaterDOWNold);

    ratioUP = ratioDOWN = 1.0;

    slaterUPnew = slaterUPold;
    slaterDOWNnew = slaterDOWNold;
    slaterUPinvNew = slaterUPinvOld;
    slaterDOWNinvNew = slaterDOWNinvOld;
}

void Slater::setAlpha(const double &newAlpha)
{
    orbitals->setAlpha(newAlpha);
}

void Slater::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;

    // only updating the inverse slater matrix _if_ we accept the move, because
    // we only need the old inverse to find the ratio
    updateSlater();
}

//double Slater::wavefunction()
//{
//    /* Only for use in the numerical derivative in the numerical local energy */

//    mat slaterUP(N, N);
//    mat slaterDOWN(N, N);

//    for (int i = 0; i < N; i++) // loop over particles
//    {
//        for (int j = 0; j < N; j++) // loop over orbitals
//        {
//            slaterUP(i,j)   = orbitals->wavefunction(rNew.row(i), j);
//            slaterDOWN(i,j) = orbitals->wavefunction(rNew.row(i+N), j);
//        }
//    }

//    return det(slaterUP)*det(slaterDOWN);
//}

double Slater::wavefunction(const mat &r)
{
    /* Only for use in the numerical derivative in the numerical local energy */

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

double Slater::getRatio()
{
    if (currentParticle < N)
    {
        int i = currentParticle;
        ratioUP = 0.0;
        for (int j = 0; j < N; j++) // loop over orbitals
            ratioUP += slaterUPnew(i,j)*slaterUPinvOld(j,i);

        return ratioUP;
    }
    else
    {
        int i = currentParticle - N;
        ratioDOWN = 0.0;
        for (int j = 0; j < N; j++) // loop over orbitals
            ratioDOWN += slaterDOWNnew(i,j)*slaterDOWNinvOld(j,i);

        return ratioDOWN;
    }
}

rowvec Slater::localGradient(const int &i)
{
    // Note: Old == New because we have accepted/rejected before getting the energy!
    grad = zeros<rowvec>(nDimensions);
    if (i < N)
        for (int j = 0; j < N; j++)
            grad += orbitals->gradient(rNew.row(i),j)*slaterUPinvNew(j,i);
    else
        for (int j = 0; j < N; j++)
            grad += orbitals->gradient(rNew.row(i),j)*slaterDOWNinvNew(j,i-N);

    return grad;
}

//mat Slater::localGradientNumerical(const double &h)
//{
//    mat gradient = zeros<mat>(nParticles, nDimensions);
//    double wfCurrent, wfMinus, wfPlus;

//    // computing the first derivative
//    rPlus = rMinus = rNew;
//    wfCurrent = wavefunction();
//    for (int i = 0; i < nParticles; i++) {
//        for (int j = 0; j < nDimensions; j++) {
//            rPlus(i,j) += h;
//            rMinus(i,j) -= h;
//            wfMinus = wavefunction(rMinus);
//            wfPlus = wavefunction(rPlus);
//            gradient(i,j) = wfPlus - wfMinus;
//            rPlus(i,j) = rMinus(i,j) = rNew(i,j);
//        }
//    }
//    gradient /= (2.0*wfCurrent*h);

//    return gradient;
//}

double Slater::localLaplacian(const int &i)
{
    lapl = 0.0;
    if (i < N)
        for (int j = 0; j < N; j++)
            lapl += orbitals->laplacian(rNew.row(i),j)*slaterUPinvNew(j,i);
    else
        for (int j = 0; j < N; j++)
            lapl += orbitals->laplacian(rNew.row(i),j)*slaterUPinvNew(j,i-N);

//    cout << "lapl Slater = " << lapl << endl;
    return lapl;
}

//double Slater::localLaplacianNumerical(const double &h)
//{
//    // computing the second derivative
//    ddwavefunction = 0.0;
//    rPlus = rMinus = rNew;
//    wfCurrent = wavefunction();

//    for (int i = 0; i < nParticles; i++) {
//        for (int j = 0; j < nDimensions; j++) {
//            rPlus(i,j) += h;
//            rMinus(i,j) -= h;
//            wfMinus = wavefunction(rMinus);
//            wfPlus = wavefunction(rPlus);
//            ddwavefunction += wfMinus + wfPlus;
//            rPlus(i,j) = rMinus(i,j) = rNew(i,j);
//        }
//    }

//    ddwavefunction = h2*(ddwavefunction/wfCurrent - 2.0*double(nParticles*nDimensions));

//    return ddwavefunction;
//}

void Slater::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    updateInverse();
    if (currentParticle < N)
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

    // not updating the inverse, since we haven't changed it
    if (currentParticle < N)
        slaterUPnew.row(currentParticle) = slaterUPold.row(currentParticle);
    else
        slaterDOWNnew.row(currentParticle-N) = slaterDOWNold.row(currentParticle-N);
}

void Slater::updateSlater()
{
    int i = currentParticle;
    if (i < N)
        for (int j = 0; j < N; j++) // loop over orbitals
            slaterUPnew(i,j) = orbitals->wavefunction(rNew.row(i), j);
    else
        for (int j = 0; j < N; j++) // loop over orbitals
            slaterDOWNnew(i-N,j) = orbitals->wavefunction(rNew.row(i), j);
}

void Slater::updateInverse()
{
    if (currentParticle < N) // only updating the UP matrix
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
