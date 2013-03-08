#include "CHelium.h"

//this my trial wave function class
//CHelium::CHelium(int myrank, int numprocs)
//{
//    int nParticles = 2;
//    int charge = 2;
//    VMCSolver(myrank, numprocs, charge, nParticles)
//}

CHelium::CHelium(int &my_rank, int &numprocs):
    VMCSolver(my_rank, numprocs, 2, 2)
{
}

//double CHelium::localEnergy(const mat &r)
//{
//    // kinetic energy
//    mat rPlus(nParticles, nDimensions); // mat rPlus(r)
//    mat rMinus(nParticles, nDimensions);
//    double waveFunctionMinus, waveFunctionPlus;
//    double waveFunctionCurrent = wavefunction(r);
//    double kineticEnergy = 0;
//    // computing the second derivative
//    for (int i = 0; i < nParticles; i++) {
//        for (int j = 0; j < nDimensions; j++) {
//            rPlus(i,j) += h;
//            rMinus(i,j) -= h;
//            waveFunctionMinus = wavefunction(rMinus);
//            waveFunctionPlus = wavefunction(rPlus);
//            //kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
//            kineticEnergy -= waveFunctionMinus + waveFunctionPlus;
//            rPlus(i,j) = rMinus(i,j) = r(i,j);
//        }
//    }
//    //kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;
//    kineticEnergy = 0.5*h2*(kineticEnergy/waveFunctionCurrent
//                            + 2.0*double(nParticles*nDimensions));

//    // potential energy
//    double potentialEnergy = 0;
//    for (int i = 0; i < nParticles; i++) {
//        // r(i) = norm(r.row(i), 2);
//        potentialEnergy += 1.0/norm(r.row(i), 2);
//    }
//    potentialEnergy *= -charge;

//    // contribution from electron-electron potential
//    double r12 = 0;
//    for (int i = 0; i < nParticles; i++) {
//        for (int j = i + 1; j < nParticles; j++) {
//            r12 = 0;
//            for (int k = 0; k < nDimensions; k++) {
//                r12 += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
//            }
//            potentialEnergy += 1/sqrt(r12);
//        }
//    }

//    return kineticEnergy + potentialEnergy;
//}

double CHelium::localEnergyClosedForm(const mat &r)
{
    double EL1, EL2, dr;
    double rInverseSum = 0, rSum = 0, rijSum = 0;
    vec drvec(nDimensions);

    for (int i = 0; i < nParticles; i++) {
        drvec = r.row(i);
        dr = 0.0;
        for (int k = 0; k < nDimensions; k++)
            dr += drvec(k)*drvec(k);
        dr = sqrt(dr);
        rSum += dr;
        rInverseSum += 1.0/dr;
        for (int j = i + 1; j < nParticles; j++) {
            drvec -= r.row(j);
            dr = 0.0;
            for (int k = 0; k < nDimensions; k++)
                dr += drvec(k)*drvec(k);
            dr = sqrt(dr);
            rijSum += dr;
        }
    }

    vec r1vec(r.row(0).t());
    vec r2vec(r.row(1).t());

    double r1 = norm(r1vec, 2);
    double r2 = norm(r2vec, 2);

    double betafactor = 1.0/(1.0 + beta*rijSum);
    double betafactor2 = betafactor*betafactor;
    double rfactor = 1.0 - dot(r1vec, r2vec)/r1/r2;

    ////
//    cout << "betafactor = " << betafactor << endl;
//    cout << "rfactor = " << rfactor << endl;
//    cout << "1 = " << alpha*rSum*rfactor/rijSum << endl;
//    cout << "2 = " << -1.0/(2.0*betafactor2) << endl;
//    cout << "3 = " << -2.0/rijSum << endl;
//    cout << "4 = " << 2.0*beta/betafactor << endl;
    ////

    EL1 = (alpha - charge)*rInverseSum + 1.0/rijSum - alpha*alpha;
    EL2 = EL1 + (0.5*betafactor2)*
            (
                alpha*rSum*rfactor/rijSum - 0.5*betafactor2 - 2.0/rijSum
                + 2.0*beta*betafactor
            );

    ////
//    cout << "braces = " << test << endl;
//    cout << "EL1 = " << EL1 << endl;
//    cout << "EL2 = " << EL2 << endl;
//    cout << endl;
    ////
    return EL2;
}

double CHelium::wavefunction(const mat &R) const
{
    double r = 0.0;
    double rSum = 0.0;

    // r1 + r2 + ...
    for (int i = 0; i < nParticles; i++) {
        r = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            r += R(i,j)*R(i,j);
        }
        rSum += sqrt(r);
    }

    double r12Sum = 0.0;
    // r1*r2*...
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            r = 0.0;
            for (int k = 0; k < nDimensions; k++)
            {
                r += (R(i,k) - R(j,k))*(R(i,k) - R(j,k));
            }
            r12Sum += sqrt(r);
        }
    }

    return exp(-alpha*rSum + r12Sum/(2.0*(1.0 + beta*r12Sum)));
}

double CHelium::wavefunction(const data &s) const
{
    return 0.0;
}

double CHelium::slaterRatio()
{
    return waveFunctionNew*waveFunctionNew/(waveFunctionOld*waveFunctionOld);
}

double CHelium::jastrowRatio(const int &k)
{
    return 1.0;
}

double CHelium::jastrowRatio(const data &s, const int &k) const
{
    return 1.0;
}

//void CHelium::setParameters(
//        const double &alpha_,
//        const double &beta_,
//        const double &stepLength_)
//{
//    alpha = alpha_;
//    beta = beta_;
//    stepLength = stepLength_;
//}

//void CHelium::setParameters(
//        const double &alpha_,
//        const double &beta_,
//        const double &stepLength_,
//        const double &h_,
//        const double &h2_)
//{
//    alpha = alpha_;
//    beta = beta_;
//    stepLength = stepLength_;
//    h = h_;
//    h2 = h2_;
//}
