#ifndef CBERYLLIUM_H
#define CBERYLLIUM_H

#define CALL_MEMBER_FN(object, ptrToMember)  ((object).*(ptrToMember))

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
    virtual double wavefunction(const mat &r) const;
    virtual double wavefunction(const data &s) const;
    double jastrowWF(const mat &fij);
    double jastrowWF(const data &s) const;
    double phiSD(const mat &r) const;
    double phiSD(const data &s) const;
    virtual double slaterRatio();
    virtual double jastrowRatio(const int &k);
    virtual double jastrowRatio(const data &s, const int &k) const;


//    double (*hydrogenWF[2]) (const vec3 &position);
//    double (*hydrogenWF[2]) (const vec3 &position);
//    double (CBeryllium::*hydrogenWF[2]) (const vec3 &position);

//    double phi1s(const vec3 &r);
//    double phi2s(const vec3 &r);
    double phi1sf(const double &r) const;
    double phi2sf(const double &r) const;
//    double phi2pf(const double &r) const;

//    void setParameters(
//            const double &alpha_,
//            const double &beta_,
//            const double &stepLength_);
//    void setParameters(
//            const double &alpha_,
//            const double &beta_,
//            const double &stepLength_,
//            const double &h_,
//            const double &h2_);
protected:
//    double alpha;
//    double beta;
//    int charge;
//    int nParticles;
};

#endif // CBERYLLIUM_H
