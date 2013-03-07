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

double CBeryllium::phiSD(const mat &r)
{
    double dist[nParticles];
    for (int i = 0; i < nParticles; i++)
    {
        dist[i] = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            dist[i] += r(i,j)*r(i,j);
        }
        dist[i] = sqrt(dist[i]);
    }
    double answer = 1.0;
    for (int i = 0; i < 4; i+=2)
        answer *= (phi1s(dist[i])*phi2s(dist[i+1]) - phi1s(dist[i+1])*phi2s(dist[i]));

    return answer;
}

double CBeryllium::phi1s(double r)
{
    return exp(-alpha*r);
}

double CBeryllium::phi2s(double r)
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
