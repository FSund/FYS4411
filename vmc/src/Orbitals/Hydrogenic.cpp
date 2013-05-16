#include <src/Orbitals/Hydrogenic.h>

Hydrogenic::Hydrogenic():
    Orbitals()
{
}

void Hydrogenic::setAlpha(const double &newAlpha)
{
    alpha = newAlpha;
}

void Hydrogenic::setR(const double &dist)
{
    (void) dist;
}

double Hydrogenic::wavefunction(const rowvec &rvec, const int &qNum)
{
    r = 0.0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    switch (qNum)
    {
    case 0 : // 1s
        return exp(-alpha*r);
    case 1 : // 2s
        arg = alpha*r*0.5;
        return (1.0 - arg)*exp(-arg);
    case 2 : // 2p
        return rvec(0)*exp(-0.5*alpha*r);
    case 3 : // 2p
        return rvec(1)*exp(-0.5*alpha*r);
    case 4 : // 2p
        return rvec(2)*exp(-0.5*alpha*r);
    default :
        // Process for all other cases.
        cout << "! We don't have this orbital yet!" << endl;
        exit(1);
    }

//	 return
//		 qNum == 0 ? exp(-alpha*r):
//		 qNum == 1 ? (1.0 - alpha*r*0.5)*exp(-alpha*r*0.5):
//		 qNum == 2 ? rvec(0)*exp(-0.5*alpha*r):
//		 qNum == 3 ? rvec(0)*exp(-0.5*alpha*r):
//		 qNum == 4 ? rvec(0)*exp(-0.5*alpha*r):
//			 something default
}

rowvec Hydrogenic::gradient(const rowvec &rvec, const int &qNum)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    switch (qNum)
    {
    case 0 :
        return (-alpha/r*exp(-alpha*r))*rvec;
    case 1 :
        return alpha*rvec*(alpha*r - 4.0)*exp(-alpha*r*0.5)/(4.0*r);
    case 2 :
        dphi = -alpha*rvec*rvec(0);
        dphi(0) += 2.0*r;
        dphi *= exp(-alpha*r*0.5)/(2.0*r);
        return dphi;
    case 3 :
        dphi = -alpha*rvec*rvec(1);
        dphi(1) += 2.0*r;
        dphi *= exp(-alpha*r*0.5)/(2.0*r);
        return dphi;
    case 4 :
        dphi = -alpha*rvec*rvec(2);
        dphi(2) += 2.0*r;
        dphi *= exp(-alpha*r*0.5)/(2.0*r);
        return dphi;
    default :
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }
}

double Hydrogenic::laplacian(const rowvec &rvec, const int &qNum)
{
    switch (qNum)
    {
    case 0 :
        return ddphi1s(rvec);
    case 1 :
        return ddphi2s(rvec);
    case 2 :
        return ddphi2p(rvec, 0);
    case 3 :
        return ddphi2p(rvec, 1);
    case 4 :
        return ddphi2p(rvec, 2);
    default:
        cout << "Please implement more hydrogen wavefunctions" << endl;
        exit(1);
    }
}

double Hydrogenic::ddphi1s(const rowvec &rvec)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);
    return alpha*(alpha*r - 2.0)*exp(-alpha*r)/r;
}

double Hydrogenic::ddphi2s(const rowvec &rvec)
{
    double r2 = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r2 += rvec(i)*rvec(i);
    }
    r = sqrt(r2);

    return -alpha*(alpha*alpha*r2 - 10.0*alpha*r + 16.0)*exp(-alpha*r*0.5)/(8.0*r);
}

double Hydrogenic::ddphi2p(const rowvec &rvec, const int &k)
{
    r = 0.0;
    for(int i = 0; i < nDimensions; i++){
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    return alpha*rvec(k)*(alpha*r - 8.0)*exp(-alpha*r*0.5)/(4.0*r);
}

double Hydrogenic::alphaGradient(const rowvec &rvec, const int &qNum)
{
    r = 0.0;
    for (int i = 0; i < nDimensions; i++)
        r += rvec(i)*rvec(i);
    r = sqrt(r);

    switch (qNum)
    {
    case 0 :
        return -r*exp(-alpha*r);
    case 1 :
        return r*(alpha*r - 4.0)*exp(-alpha*r*0.5)*0.25;
    case 2 :
        return -rvec(0)*r*exp(-alpha*r*0.5)*0.5;
    case 3 :
        return -rvec(1)*r*exp(-alpha*r*0.5)*0.5;
    case 4 :
        return -rvec(2)*r*exp(-alpha*r*0.5)*0.5;
    default :
        // Process for all other cases.
        cout << "! We don't have this orbital yet!" << endl;
        exit(1);
    }
}
