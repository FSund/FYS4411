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

void Slater::initalize(const mat &r)
{
    rNew = rOld = r;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            slaterUPold(i,j)   = orbitals->wavefunction(i, rOld.row(j));
            slaterDOWNold(i,j) = orbitals->wavefunction(i, rOld.row(j+N));
        }
    }
    slaterUPinvOld = inv(slaterUPold);
    slaterDOWNinvOld = inv(slaterDOWNold);

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

//    updateSlater();
}

double Slater::wavefunction()
{
    mat slaterUP(N, N);
    mat slaterDOWN(N, N);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            slaterUP(i,j) = orbitals->wavefunction(i, rNew.row(j));
            slaterDOWN(i,j) = orbitals->wavefunction(i, rNew.row(j+N));
        }
    }

    return det(slaterUP)*det(slaterDOWN);
}

double Slater::wavefunction(const mat &r)
{
    /* For use in the numerical derivative, in the numerical local energy,
     * and in the temporary ratio function */

    mat slaterUP(N, N);
    mat slaterDOWN(N, N);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            slaterUP(i,j) = orbitals->wavefunction(i, r.row(j));
            slaterDOWN(i,j) = orbitals->wavefunction(i, r.row(j+N));
        }
    }

    return det(slaterUP)*det(slaterDOWN);
}

double Slater::getRatio()
{
//    ratio = 0.0;
//    int k = currentParticle;
//    for (int i = 0; i < N; i++)
//    {
//        ratio += slaterUPnew(i,k)*slaterUPinvOld(k,i);
//        ratio += slaterDOWNnew(i,k)*slaterDOWNinvOld(k,i);
//    }

//    return ratio;

    return wavefunction(rNew)*wavefunction(rNew) /
            (wavefunction(rOld)*wavefunction(rOld));
}

mat Slater::gradient()
{
    mat temp(nParticles, nDimensions);
    temp.fill(1.0);
    return temp;
}

double Slater::laplacian()
{
    return 1.0;
}

void Slater::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    // do something with the inverse (and maybe the slater matrix?)
}

void Slater::rejectMove()
{
    rNew.row(currentParticle) = rOld.row(currentParticle);

    // do something with the inverse and the slater matrix
}

void Slater::updateSlater()
{
    int k = currentParticle;
    for (int i = 0; i < N; i++)
    {
        slaterUPnew(i,k)   = orbitals->wavefunction(i, rNew.row(k));
        slaterDOWNnew(i,k) = orbitals->wavefunction(i, rNew.row(k+N));
    }
}

void Slater::updateInverse()
{
    vec SUP = zeros<vec>(N);
    vec SDOWN = zeros<vec>(N);
    int k = currentParticle;
    for (int j = 0; j < N; j++)
    {
        for (int l = 0; l < N; l++)
        {
            SUP(j) += slaterUPnew(l,k)*slaterUPinvOld(l,j);
            SDOWN(j) += slaterDOWNnew(l,k)*slaterUPinvOld(l,j);
        }
    }
}
