#ifndef CBERYLLIUM_H
#define CBERYLLIUM_H

#include "CVMCSolver.h"
#include "mainapplication.h"

/*It is convenient to make classes of trial wave functions, both many-body
 *wavefunctions and single-particle wave functions and the quantum numbers
 *involved, such as spin, orbital momentum and principal quantum numbers.*/

class CBeryllium : public VMCSolver
{
public:
    friend class MainApplication;

    CBeryllium();
    CBeryllium(int my_rank_, int numprocs_);

    virtual double localEnergyClosedForm(const mat &r);
    virtual double wavefunction(const mat &r);

//    double phi1s(vec3 r);
//    double phi2s(vec3 r);
    double phi1s(double dr);
    double phi2s(double dr);

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

#endif // CBERYLLIUM_H
