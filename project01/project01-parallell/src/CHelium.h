#ifndef CHELIUM_H
#define CHELIUM_H

#include "CVMCSolver.h"
#include "mainapplication.h"

/*It is convenient to make classes of trial wave functions, both many-body
 *wavefunctions and single-particle wave functions and the quantum numbers
 *involved, such as spin, orbital momentum and principal quantum numbers.*/

class CHelium : public VMCSolver
{
public:
    friend class MainApplication;

    CHelium();
    CHelium(int &my_rank, int &numprocs);

//    virtual double localEnergy(const mat &r);
    virtual double localEnergyClosedForm(const mat &r);
    virtual double wavefunction(const mat &R);

    void setParameters(
            const double &alpha_,
            const double &beta_,
            const double &stepLength_);
    void setParameters(
            const double &alpha_,
            const double &beta_,
            const double &stepLength_,
            const double &h_,
            const double &h2_);
protected:
    double alpha;
    double beta;
//    int charge;
//    int nParticles;
};

#endif // CHELIUM_H
