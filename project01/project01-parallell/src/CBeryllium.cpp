#include "CBeryllium.h"

//this my trial wave function class
//CBeryllium::CBeryllium(int myrank, int numprocs)
//{
//    int nParticles = 2;
//    int charge = 2;
//    VMCSolver(myrank, numprocs, charge, nParticles)
//}

CBeryllium::CBeryllium(int my_rank_, int numprocs_):
    VMCSolver(my_rank_, numprocs_, 4, 4)
{
//    hydrogenWF[0] = &CBeryllium::phi1s;
//    hydrogenWF[1] = &CBeryllium::phi2s;
}

double CBeryllium::localEnergyClosedForm(const mat &r)
{
    cout << "!! Closed form isn't implemented for Beryllium !!" << endl;
    return 1.0;
}

double CBeryllium::wavefunction(const mat &r)
{
    return phiSD(r);
}

double CBeryllium::wavefunction(const mat &r, const mat &fij)
{
    return phiSD(r)*jastrowWF(fij);
}

double CBeryllium::slaterRatio()
{
    return waveFunctionNew*waveFunctionNew/(waveFunctionOld*waveFunctionOld);
}

double CBeryllium::jastrowRatio(const int &k)
{
    return 1.0;

    double dU = 0.0;
    for (int i = 0; i < k; i++)
        dU += calculate_fij_element(rijNew, i, k) - calculate_fij_element(rijOld, i, k);
    for (int i = k+1; i < nParticles; i++)
        dU += calculate_fij_element(rijNew, k, i) - calculate_fij_element(rijOld, k, i);
    return exp(dU);
}

double CBeryllium::jastrowWF(const mat &fij)
{
    double arg = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            arg += fij(i,j);
        }
    }
    return exp(arg);
}

//double CBeryllium::phiSD(const mat &r)
//{
//    double dist[nParticles];
//    for (int i = 0; i < nParticles; i++)
//    {
//        dist[i] = 0.0;
//        for (int j = 0; j < nDimensions; j++)
//        {
//            dist[i] += r(i,j)*r(i,j);
//        }
//        dist[i] = sqrt(dist[i]);
//    }
//    double answer = 1.0;
//    for (int i = 0; i < 4; i+=2)
//        answer *= (phi1s(dist[i])*phi2s(dist[i+1]) - phi1s(dist[i+1])*phi2s(dist[i]));

//    return answer;
//}

double CBeryllium::phiSD(const mat &r)
{
    // Wrong version (returns 0)
//    double answerDet;
//    mat slater(nParticles, nParticles);

//    for (int i = 0; i < nParticles; i++)
//    {
//        for (int j = 0; j < nParticles; j++)
//        {
//            slater(i,j) = hydrogenWF(i, r.row(j));
//        }
//    }
//    answerDet = det(slater, true);

    // Old version
//    double answer = 1.0;
//    for (int i = 0; i < 4; i+=2)
//        answer *= (phi1s(r.row(i))*phi2s(r.row(i+1)) - phi1s(r.row(i+1))*phi2s(r.row(i)));
//    return answer;

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

//    cout << "sdUP " << endl << slaterUP << "sdDOWN" << endl << slaterDOWN << endl;

    return det(slaterUP)*det(slaterDOWN);
}

//double CBeryllium::phi1s(const vec3 &rvec)
//{
//    double r = 0;
//    for (int i = 0; i < nDimensions; i++)
//    {
//        r += rvec(i)*rvec(i);
//    }
//    r = sqrt(r);

//    return exp(-alpha*r);
//}

//double CBeryllium::phi2s(const vec3 &rvec)
//{
//    double r = 0;
//    for (int i = 0; i < nDimensions; i++)
//    {
//        r += rvec(i)*rvec(i);
//    }
//    r = sqrt(r);

//    double arg = -alpha*r*0.5;
//    return (1.0 + arg)*exp(arg);
//}

double CBeryllium::phi1sf(const double &r) const
{
    return exp(-alpha*r);
}

double CBeryllium::phi2sf(const double &r) const
{
    double arg = -alpha*r*0.5;
    return (1.0 + arg)*exp(arg);
}

//void CBeryllium::setParameters(
//        const double &alpha_,
//        const double &beta_,
//        const double &stepLength_)
//{
//    alpha = alpha_;
//    beta = beta_;
//    stepLength = stepLength_;
//}

//void CBeryllium::setParameters(
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
