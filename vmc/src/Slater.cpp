#include "Slater.h"

Slater::Slater(const int &nParticles):
    nParticles(nParticles),
    nDimensions(3),
    rOld(zeros<mat>(nParticles, nDimensions)),
    rNew(zeros<mat>(nParticles, nDimensions))
{
}

void Slater::initalize(const mat &r)
{
    rNew = rOld = r;
}

void Slater::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;
}

void Slater::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

double Slater::wavefunction()
{
    int size = int(nParticles/2);
    mat slaterUP(size, size);
    mat slaterDOWN(size, size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            slaterUP(i,j)   = hydrogenWF(i, rNew.row(j));
            slaterDOWN(i,j) = hydrogenWF(i, rNew.row(j+size));
        }
    }

    return det(slaterUP)*det(slaterDOWN);
}

double Slater::wavefunction(const mat &r)
{
    int size = int(nParticles/2);
    mat slaterUP(size, size);
    mat slaterDOWN(size, size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            slaterUP(i,j)   = hydrogenWF(i, r.row(j));
            slaterDOWN(i,j) = hydrogenWF(i, r.row(j+size));
        }
    }

    return det(slaterUP)*det(slaterDOWN);
}

double Slater::getRatio()
{
    return wavefunction(rNew)*wavefunction(rNew) /
            (wavefunction(rOld)*wavefunction(rOld));
}

mat Slater::gradient()
{

}

double Slater::laplacian()
{

}

void Slater::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);
}

void Slater::rejectMove()
{
    rNew.row(currentParticle) = rOld.row(currentParticle);
}

double Slater::hydrogenWF(const int &i, const vec3 &rvec)
{
    double wfCurrent = 0;

    switch (i)
    {
    case 0:
        wfCurrent = phi1s(rvec);
//        cout << "i = " << i << ", using phi1s" << endl;
        break;
    case 1 :
        wfCurrent = phi2s(rvec);
//        cout << "i = " << i << ", using phi2s" << endl;
        break;
//    case 2 :
//        wavefunction = phi2s(r);
//        break;
//    case 3 :
//        wavefunction = phi2s(r);
//        break;
    default :
        // Process for all other cases.
        cout << "!!! We don't have this wavefunction yet!" << endl;
        exit(1);
    }

    return wfCurrent;
}

double Slater::phi1s(const vec3 &rvec)
{
    double r = 0;
    for (int i = 0; i < nDimensions; i++)
    {
        r += rvec(i)*rvec(i); // maybe store rii vector in orbital class?
    }
    r = sqrt(r);

    return exp(-alpha*r);
}

double Slater::phi2s(const vec3 &rvec)
{
    double r = 0;
    for (int i = 0; i < nDimensions; i++)
    {
        r += rvec(i)*rvec(i); // maybe store rii vector in orbital class?
    }
    r = sqrt(r);

    double arg = alpha*r*0.5;

    return (1.0 - arg)*exp(-arg);
}
