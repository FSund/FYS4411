#include <src/Orbitals.h>

Orbitals::Orbitals():
    nDimensions(3),
    dphi(vec(nDimensions))
{
}

void Orbitals::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

double Orbitals::wavefunction(const rowvec &rvec, const int &qNum)
{
    r = 0.0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    switch (qNum)
    {
    case 0 :
        return exp(-alpha*r);
    case 1 :
        arg = alpha*r*0.5;
        return (1.0 - arg)*exp(-arg);
    case 2:
        return rvec(0)*exp(-0.5*alpha*r);
    case 3:
        return rvec(1)*exp(-0.5*alpha*r);
    case 4:
        return rvec(2)*exp(-0.5*alpha*r);
    default :
        // Process for all other cases.
        cout << "! We don't have this orbital yet!" << endl;
        exit(1);
    }
}

rowvec Orbitals::gradient(const rowvec &rvec, const int &qNum)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    switch (qNum)
    {
    case 0:
        return (-alpha/r*exp(-alpha*r))*rvec;
    case 1:
        return alpha*rvec*(alpha*r - 4.0)*exp(-alpha*r/2.0)/(4.0*r);
    case 2:
        dphi = -alpha*rvec*rvec(0);
        dphi(0) += 2.0*r;
        dphi *= exp(-alpha*r/2.0)/(2.0*r);
        return dphi;
    case 3:
        dphi = -alpha*rvec*rvec(1);
        dphi(1) += 2.0*r;
        dphi *= exp(-alpha*r/2.0)/(2.0*r);
        return dphi;
    case 4:
        dphi = -alpha*rvec*rvec(2);
        dphi(2) += 2.0*r;
        dphi *= exp(-alpha*r/2.0)/(2.0*r);
        return dphi;
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }
}

double Orbitals::laplacian(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0:
        return ddphi1s(rvec);
    case 1:
        return ddphi2s(rvec);
    case 2:
        return ddphi2p(rvec, 0);
    case 3:
        return ddphi2p(rvec, 1);
    case 4:
        return ddphi2p(rvec, 2);
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }
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

    return alpha*rvec(k)*(alpha*r - 8.0)*exp(-alpha*r/2.0)/(4.0*r);
}
