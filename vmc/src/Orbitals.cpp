#include "Orbitals.h"

Orbitals::Orbitals():
    nDimensions(3)
{
}

void Orbitals::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

double Orbitals::wavefunction(const int &qNum, const vec &rvec)
{
    switch (qNum)
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
        cout << "! We don't have this orbital yet!" << endl;
        exit(1);
    }

    return wfCurrent;
}

double Orbitals::phi1s(const vec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    return exp(-alpha*r);
}

double Orbitals::phi2s(const vec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);
    arg = alpha*r*0.5;

    return (1.0 - arg)*exp(-arg);
}
