#ifndef JASTROW_H
#define JASTROW_H

#include <armadillo>

using namespace std;
using namespace arma;

class Jastrow
{
public:
    Jastrow(const int &nParticles);
    void initialize(const mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setBeta(const double &newBeta);

    double wavefunction();
    double wavefunction(const mat &r);
    double getRatio();
    rowvec localGradient(const int &k);
    rowvec localGradient(const mat &r, const int &k);
    double localLaplacian(const int &k);
    double localLaplacian(const mat &r, const int &k);


    void acceptMove();
    void rejectMove();

    // for debugging/tests
    // these functions return the full gradient and laplacian, not just for one
    // particle like the closed for functions
    rowvec localGradientNumerical(const int &k, const double &h);
    rowvec localGradientNumerical(const mat &r, const int &k, const double &h);
    double localLaplacianNumerical(const int &k, const double &h);
    double localLaplacianNumerical(const mat &r, const int &k, const double &h);
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    mat rijOld, rijNew; // rij only needed by Jastrow -- rii needed by H-like-orbitals
    mat fijOld, fijNew; // fij only needed by Jastrow
    mat a;
    double beta;
    int currentParticle;

    // numerical gradient and laplacian
    double wfCurrent, wfMinus, wfPlus;
    mat rPlus, rMinus;
    rowvec dwavefunction;
    double ddwavefunction;

    // closed form gradient and laplacian
    double rij;
    rowvec grad;
    double lapl;

    void calculate_rij();
    void update_rij();
    void calculate_fij();
    void update_fij();
};

#endif // JASTROW_H
