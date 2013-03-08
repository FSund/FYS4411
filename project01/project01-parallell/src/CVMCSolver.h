#ifndef CVMCSOLVER_H
#define CVMCSOLVER_H

#include "lib.h"
#include <mpi.h>
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

const double pi = atan(1)*4;

class VMCSolver
{
public:
    VMCSolver(int my_rank, int numprocs, const int &charge, const int &nParticles);
    ~VMCSolver();

//    double runMonteCarloIntegration(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_, const bool closedform);
    double runMonteCarloIntegration(const int &nCycles);
//    double runMonteCarloIntegrationImportanceSampling(const int &newNCycles, const double &stepLength_, const double &alpha_, const double &beta_, const bool closedform, const double dt);
//    double runMonteCarloIntegrationImportanceSampling(const int &nCycles, const double &dt, const bool &closedform);

    void runcycle_importanceSampling(const int &i);
    void runcycle(const int &i);

    void setParameters(
            const double &alpha_,
            const double &beta_,
            const double &stepLength_);
    void setParameters(
            const double &alpha_,
            const double &beta_,
            const double &stepLength_,
            const double &h_,
            const double &h2_,
            const bool &importanceSampling_,
            const bool &closedForm_);
<<<<<<< HEAD

    double hydrogenWF(const int &i, const vec3 &rvec);
    double phi1s(const vec3 &rvec);
    double phi2s(const vec3 &rvec);
=======
>>>>>>> 11f443a3f559fada5f12f377790cb7e18191ad0c
protected:
    mat rOld;
    mat rNew;
    mat rijOld;
    mat rijNew;
    mat qForceOld;
    mat qForceNew;
<<<<<<< HEAD
    mat fijOld;
    mat fijNew;
    double waveFunctionOld;
    double waveFunctionNew;

    virtual double wavefunction(const mat&) = 0;
    virtual double wavefunction(const mat&, const mat&) = 0;
//    virtual double localEnergy(const mat &r)=0;
    virtual double localEnergyClosedForm(const mat&) = 0;
=======
    double waveFunctionOld;
    double waveFunctionNew;

    virtual double wavefunction(const mat&)=0;
    virtual double wavefunction(const mat&, const mat&)=0;
//    virtual double localEnergy(const mat &r)=0;
    virtual double localEnergyClosedForm(const mat&)=0;
>>>>>>> 11f443a3f559fada5f12f377790cb7e18191ad0c

    double localEnergy(const mat &r);
    mat quantumForce(const mat &r, const double &wf);
//    double waveFunction2(const mat &r);
//    double dr(const mat &r, const int ii, const int jj);

    void calculate_rij(const mat &r, mat &rij);
    void update_rij(const mat &r, mat &rij, const int j);
<<<<<<< HEAD
    void calculate_fij(const mat &rij, mat &fij);
    double calculate_fij_element(const mat &rij, const int &i, const int &j);
    void update_slater();

    virtual double slaterRatio() = 0;
    virtual double jastrowRatio(const int &k) = 0;

=======
    double jastrowRatio(const int k);
    double fij(const mat rij, const int i, const int j);
    void update_slater();

>>>>>>> 11f443a3f559fada5f12f377790cb7e18191ad0c
    double gaussianDeviate(long *seed);

    int nDimensions;
    int charge;
    double alpha;
    double beta;
    double stepLength;
    int nParticles;

    double h;
    double h2;
    double D;
    double Ddt;

    long idum;

    int nCycles;
    int nAccepted;
    int nRejected;

    int my_rank;
    int numprocs;
    int local_nCycles;

public:
    bool importanceSampling;
    bool closedForm;
};

#endif // CVMCSOLVER_H
