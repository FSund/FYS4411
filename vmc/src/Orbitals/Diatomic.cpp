#include <src/Orbitals/Diatomic.h>

Diatomic::Diatomic():
    Orbitals(),
    R(zeros<rowvec>(nDimensions))
{
    hydrogenic = new Hydrogenic();
}

Diatomic::~Diatomic()
{
    delete hydrogenic;
}

void Diatomic::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
    hydrogenic->setAlpha(newAlpha);
}

void Diatomic::setR(const double &dist)
{
    R(0) = dist/2.0;
}

double Diatomic::wavefunction(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0)
        return hydrogenic->wavefunction(rvec + R, qNum/2)
                + hydrogenic->wavefunction(rvec - R, qNum/2);
    else
        return hydrogenic->wavefunction(rvec + R, qNum/2)
                - hydrogenic->wavefunction(rvec - R, qNum/2);
}

rowvec Diatomic::gradient(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0)
        return hydrogenic->gradient(rvec + R, qNum/2)
                + hydrogenic->gradient(rvec - R, qNum/2);
    else
        return hydrogenic->gradient(rvec + R, qNum/2)
                - hydrogenic->gradient(rvec - R, qNum/2);
}

double Diatomic::laplacian(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0)
        return hydrogenic->laplacian(rvec + R, qNum/2)
                + hydrogenic->laplacian(rvec - R, qNum/2);
    else
        return hydrogenic->laplacian(rvec + R, qNum/2)
                - hydrogenic->laplacian(rvec - R, qNum/2);
}

double Diatomic::alphaGradient(const rowvec &rvec, const int &qNum)
{
    if (qNum%2 == 0)
        return hydrogenic->alphaGradient(rvec + R, qNum/2)
                + hydrogenic->alphaGradient(rvec - R, qNum/2);
    else
        return hydrogenic->alphaGradient(rvec + R, qNum/2)
                - hydrogenic->alphaGradient(rvec - R, qNum/2);
}

