#ifndef LOCALENERGY_H
#define LOCALENERGY_H

#include <armadillo>
#include <src/Wavefunction.h>

using namespace std;
using namespace arma;

class LocalEnergy
{
public:
    LocalEnergy();
    LocalEnergy(const int &nParticles, const int &nDimensions, const int &charge);
    virtual ~LocalEnergy();
    virtual void setR(const double &dist) = 0;
    virtual double evaluate(const mat &r, Wavefunction *wf) = 0;
    void setClosedForm(const bool &closedForm_) { closedForm = closedForm_; }
    void setUseJastrow(const bool &useJastrow_) { useJastrow = useJastrow_; }
protected:
    int nParticles;
    int nDimensions;
    double charge;
//    rowvec R;
    bool closedForm;
    bool useJastrow;
};

#endif // LOCALENERGY_H
