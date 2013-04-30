#include "Orbitals.h"

Orbitals::Orbitals():
    nDimensions(3)
{
}

void Orbitals::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

double Orbitals::wavefunction(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0 :
        wfCurrent = phi1s(rvec);
        break;
    case 1 :
        wfCurrent = phi2s(rvec);
        break;
    case 2:
        wfCurrent = phi2p(rvec, 0);
        break;
    case 3:
        wfCurrent = phi2p(rvec, 1);
        break;
    case 4:
        wfCurrent = phi2p(rvec, 2);
        break;
    default :
        // Process for all other cases.
        cout << "! We don't have this orbital yet!" << endl;
        exit(1);
    }

    return wfCurrent;
}

double Orbitals::phi1s(const rowvec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    return exp(-alpha*r);
}

double Orbitals::phi2s(const rowvec &rvec)
{
    r = 0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);
    arg = alpha*r*0.5;

    return (1.0 - arg)*exp(-arg);
}

double Orbitals::phi2p(const rowvec &rvec, const int &k)
{
    r = 0.0;
    for (int i = 0; i<nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec(k)*exp(-0.5*alpha*r);
}

rowvec Orbitals::gradient(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0:
        grad = dphi1s(rvec);
        break;
    case 1:
        grad = dphi2s(rvec);
        break;
    case 2:
        grad = dphi2p(rvec, 0);
        break;
    case 3:
        grad = dphi2p(rvec, 1);
        break;
    case 4:
        grad = dphi2p(rvec, 2);
        break;
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }

    return grad;
}

rowvec Orbitals::dphi1s(const rowvec &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return (-alpha/r*exp(-alpha*r))*rvec;
}

rowvec Orbitals::dphi2s(const rowvec &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec*(alpha*r - 4.0)*exp(-alpha*r/2.0)/(4.0*r);
}

rowvec Orbitals::dphi2p(const rowvec &rvec, const int &k)
{
    vec dphi;
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    double arg = exp(-alpha*r/2.0);
    for (int i = 0; i < nDimensions; i++)
    {
        dphi = -alpha*alpha/(2.0*r)*arg*rvec(i)*rvec(k);
    }
    dphi(k) += alpha*arg;

    return dphi;
}

double Orbitals::laplacian(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0:
        lapl = ddphi1s(rvec);
        break;
    case 1:
        lapl = ddphi2s(rvec);
        break;
    case 2:
        lapl = ddphi2p(rvec, 0);
        break;
    case 3:
        lapl = ddphi2p(rvec, 1);
        break;
    case 4:
        lapl = ddphi2p(rvec, 2);
        break;
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }

    return lapl;
}

double Orbitals::ddphi1s(const rowvec &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    return alpha*(alpha*r - 2.0)*exp(-alpha*r)/r;
}

double Orbitals::ddphi2s(const rowvec &rvec)
{
    double r2 = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r2 += rvec(i)*rvec(i);
    }
    r = sqrt(r2);

    return -alpha*(alpha*alpha*r2 - 10*alpha*r + 16)*exp(-alpha*r/2)/(8*r);
}

double Orbitals::ddphi2p(const rowvec &rvec, const int &k)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return pow(alpha, 2.0)*rvec(k)*(alpha*r - 8.0)*exp(-alpha*r/2.0)/(4.0*r);
}
