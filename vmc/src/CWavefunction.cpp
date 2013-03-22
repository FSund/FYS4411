#include "CWavefunction.h"

Wavefunction::Wavefunction(const int &nParticles, const double &charge):
    nParticles(nParticles),
    nDimensions(3),
    charge(charge),
    h(1e-3),
    h2(1e6)
{
//    alpha = 1.8;
//    beta = 0.36;
//    alpha = 3.9;
//    beta = 4.0;

    dwavefunction = zeros<mat>(nParticles, nDimensions);
    rPlus = rMinus = zeros<mat>(nParticles, nDimensions);
    rNew = rOld = zeros(nParticles, nDimensions);
    rijOld = rijNew = fijNew = fijOld = zeros(nParticles, nParticles);
}

void Wavefunction::initialize(mat &r)
{
    rNew = rOld = r;
    calculate_rij();
    calculate_fij();
}

void Wavefunction::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;

    update_rij();
    update_fij();
}

void Wavefunction::setAlpha(const double &alpha_)
{
    alpha = alpha_;
}

void Wavefunction::setBeta(const double &beta_)
{
    beta = beta_;
}

double Wavefunction::wavefunction(const mat &r)
{
//    return jastrowWF()*phiSD();

    return phiSD(r)*jastrowWF();
}

double Wavefunction::getRatio()
{
//    double phiNew, phiOld;
//    phiNew = phiSD(rNew);
//    phiOld = phiSD(rOld);

//    cout << "rij" << endl << rijNew;
//    cout << "jastrowRatio = " << jastrowRatio() << endl;
//    cout << "slaterRatio = " << slaterRatio() << endl;

    return slaterRatio()*jastrowRatio();
}

void Wavefunction::calculate_rij()
{
    rijNew.zeros();
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            for (int k = 0; k < nDimensions; k++)
            {
                rijNew(i,j) += (rNew(i,k) - rNew(j,k))*(rNew(i,k) - rNew(j,k));
            }
            rijNew(i,j) = sqrt(rijNew(i,j));
        }
    }
    rijOld = rijNew;
}

void Wavefunction::update_rij()
{
    int k = currentParticle;
    for (int i = 0; i < k; i++)
    {
        rijNew(i,k) = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rijNew(i,k) += (rNew(i,j) - rNew(k,j))*(rNew(i,j) - rNew(k,j));
        }
        rijNew(i,k) = sqrt(rijNew(i,k));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rijNew(k,i) = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rijNew(k,i) += (rNew(i,j) - rNew(k,j))*(rNew(i,j) - rNew(k,j));
        }
        rijNew(k,i) = sqrt(rijNew(k,i));
    }
}

void Wavefunction::calculate_fij()
{
    double a;
    double r12;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            r12 = rijNew(i,j);
            a = ((i+j)%2 == 0) ? 0.25 : 0.5;
            fijNew(i,j) = a*r12/(1.0 + beta*r12);
        }
    }
}

void Wavefunction::update_fij()
{
//    calculate_fij();

    int k = currentParticle;
    double r12, a;
    for (int i = 0; i < k; i++)
    {
        r12 = rijNew(i,k);
        a = ((i+k)%2 == 0) ? 0.25 : 0.5;
        fijNew(i,k) = a*r12/(1.0 + beta*r12);
    }
    for (int i = k+1; i < nParticles; i++)
    {
        r12 = rijNew(k,i);
        a = ((k+i)%2 == 0) ? 0.25 : 0.5;
        fijNew(k,i) = a*r12/(1.0 + beta*r12);
    }
}

double Wavefunction::slaterRatio()
{
    return phiSD(rNew)*phiSD(rNew)/(phiSD(rOld)*phiSD(rOld));
}

double Wavefunction::phiSD()
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

double Wavefunction::phiSD(const mat &r)
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


//    double wf;
//    wf = (hydrogenWF(0, r.row(0))*hydrogenWF(1, r.row(1)) - hydrogenWF(0, r.row(1))*hydrogenWF(1, r.row(0)))*
//        (hydrogenWF(0, r.row(2))*hydrogenWF(1, r.row(3)) - hydrogenWF(0, r.row(3))*hydrogenWF(1, r.row(2)));

    //    return wf;
}

double Wavefunction::jastrowWF() const
{
    double arg = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            arg += fijNew(i,j);
        }
    }

    return exp(arg);
}

double Wavefunction::jastrowRatio()
{
//    return 1.0;

    double dU = 0.0;
    int k = currentParticle;
    for (int i = 0; i < k; i++)
    {
        dU += fijNew(i, k) - fijOld(i, k);
    }
    for (int i = k+1; i < nParticles; i++)
    {
        dU += fijNew(k, i) - fijOld(k, i);
    }

    return exp(dU);
}

void Wavefunction::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    fijOld.row(currentParticle) = fijNew.row(currentParticle);
    fijOld.col(currentParticle) = fijNew.col(currentParticle);
    //    fijOld = fijNew;
    rijOld.row(currentParticle) = rijNew.row(currentParticle);
    rijOld.col(currentParticle) = rijNew.col(currentParticle);
//    rijOld = rijNew;
}

void Wavefunction::rejectMove()
{
    rNew.row(currentParticle) = rOld.row(currentParticle);

    fijNew.row(currentParticle) = fijOld.row(currentParticle);
    fijNew.col(currentParticle) = fijOld.col(currentParticle);
//    fijNew = fijOld;
    rijNew.row(currentParticle) = rijOld.row(currentParticle);
    rijNew.col(currentParticle) = rijOld.col(currentParticle);
//    rijOld = rijNew;
}

double Wavefunction::localEnergyNumerical()
{
    double kinetic = -0.5*laplaceNumerical();
    double potential = electronNucleusPotential() + electronElectronPotential();

    return kinetic + potential;
}

double Wavefunction::laplaceNumerical()
{
    // computing the second derivative
    double ddwavefunction = 0.0;
    rPlus = rMinus = rNew;
    double wf = wavefunction(rNew);

    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            ddwavefunction += wfMinus + wfPlus;
            rPlus(i,j) = rMinus(i,j) = rNew(i,j);
        }
    }
    ddwavefunction = h2*(ddwavefunction/wf - 2.0*double(nParticles*nDimensions));

    return ddwavefunction;
}

mat Wavefunction::gradientNumerical()
{
    // computing the first derivative
    rPlus = rMinus = rNew;
    wfCurrent = wavefunction(rNew);
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            dwavefunction(i,j) = wfPlus - wfMinus;
            rPlus(i,j) = rMinus(i,j) = rNew(i,j);
        }
    }
    dwavefunction /= (2.0*wfCurrent*h);

    return dwavefunction;
}

double Wavefunction::electronNucleusPotential()
{
    // potential energy
    double potentialEnergy = 0.0;
    for (int i = 0; i < nParticles; i++) {
        potentialEnergy += 1.0/norm(rNew.row(i), 2);
    }
    potentialEnergy *= -charge;

    return potentialEnergy;
}

double Wavefunction::electronElectronPotential()
{
    // contribution from electron-electron potential
    double potentialEnergy = 0.0;
    double rij = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            rij = 0;
            for (int k = 0; k < nDimensions; k++) {
                rij += (rNew(i,k) - rNew(j,k))*(rNew(i,k) - rNew(j,k));
            }
            // alternativ: rij = norm(r.row(i) - r.row(j), 2);
            potentialEnergy += 1/sqrt(rij);
        }
    }

    return potentialEnergy;
}

double Wavefunction::hydrogenWF(const int &i, const vec3 &rvec)
{
//    double wfCurrent = 0;

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

double Wavefunction::phi1s(const vec3 &rvec)
{
    double r = 0;
    for (int i = 0; i < nDimensions; i++)
    {
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

//    cout << "phi1s = " << exp(-alpha*r) << endl;
    return exp(-alpha*r);
}

double Wavefunction::phi2s(const vec3 &rvec)
{
    double r = 0;
    for (int i = 0; i < nDimensions; i++)
    {
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    double arg = alpha*r*0.5;

//    cout << "phi2s = " << (1.0 + arg)*exp(arg) << endl;
    return (1.0 - arg)*exp(-arg);
}

//double Helium::localEnergyClosedForm(const mat &r) const
//{
//    double EL1, EL2, dr;
//    double rInverseSum = 0, rSum = 0, rijSum = 0;
//    vec drvec(nDimensions);

//    for (int i = 0; i < nParticles; i++) {
//        drvec = r.row(i);
//        dr = 0.0;
//        for (int k = 0; k < nDimensions; k++)
//            dr += drvec(k)*drvec(k);
//        dr = sqrt(dr);
//        rSum += dr;
//        rInverseSum += 1.0/dr;
//        for (int j = i + 1; j < nParticles; j++) {
//            drvec -= r.row(j);
//            dr = 0.0;
//            for (int k = 0; k < nDimensions; k++)
//                dr += drvec(k)*drvec(k);
//            dr = sqrt(dr);
//            rijSum += dr;
//        }
//    }

//    vec r1vec(r.row(0).t());
//    vec r2vec(r.row(1).t());

//    double r1 = norm(r1vec, 2);
//    double r2 = norm(r2vec, 2);

//    double betafactor = 1.0/(1.0 + beta*rijSum);
//    double betafactor2 = betafactor*betafactor;
//    double rfactor = 1.0 - dot(r1vec, r2vec)/r1/r2;

//    ////
////    cout << "betafactor = " << betafactor << endl;
////    cout << "rfactor = " << rfactor << endl;
////    cout << "1 = " << alpha*rSum*rfactor/rijSum << endl;
////    cout << "2 = " << -1.0/(2.0*betafactor2) << endl;
////    cout << "3 = " << -2.0/rijSum << endl;
////    cout << "4 = " << 2.0*beta/betafactor << endl;
//    ////

//    EL1 = (alpha - charge)*rInverseSum + 1.0/rijSum - alpha*alpha;
//    EL2 = EL1 + (0.5*betafactor2)*
//            (
//                alpha*rSum*rfactor/rijSum - 0.5*betafactor2 - 2.0/rijSum
//                + 2.0*beta*betafactor
//            );

//    ////
////    cout << "braces = " << test << endl;
////    cout << "EL1 = " << EL1 << endl;
////    cout << "EL2 = " << EL2 << endl;
////    cout << endl;
//    ////
//    return EL2;
//}
