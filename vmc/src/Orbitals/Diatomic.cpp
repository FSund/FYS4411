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

//    cout << "R orbitals = " << R;
}

double Diatomic::wavefunction(const rowvec &rvec, const int &qNum)
{
    return hydrogenic->wavefunction(rvec + R, qNum)
            + hydrogenic->wavefunction(rvec - R, qNum);
}

rowvec Diatomic::gradient(const rowvec &rvec, const int &qNum)
{
    return hydrogenic->gradient(rvec + R, qNum)
            + hydrogenic->gradient(rvec - R, qNum);
}

double Diatomic::laplacian(const rowvec &rvec, const int &qNum)
{
    return hydrogenic->laplacian(rvec + R, qNum)
            + hydrogenic->laplacian(rvec - R, qNum);
}

double Diatomic::alphaGradient(const rowvec &rvec, const int &qNum)
{
    return hydrogenic->alphaGradient(rvec + R, qNum)
            + hydrogenic->alphaGradient(rvec - R, qNum);
}

