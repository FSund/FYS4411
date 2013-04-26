#ifndef JASTROW_H
#define JASTROW_H

#include <armadillo>

using namespace std;
using namespace arma;

class Jastrow
{
public:
    Jastrow(const int &nParticles);
    void initalize(const mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setBeta(const double &newBeta);

//    double wavefunction();
    double wavefunction(const mat &r);
    double getRatio();
    mat gradient();
    double laplacian();

    void acceptMove();
    void rejectMove();
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    mat rijOld, rijNew; // rij only needed by Jastrow -- rii needed by H-like-orbitals
    mat fijOld, fijNew; // fij only needed by Jastrow
    mat a;
    double beta;
    int currentParticle;

    void calculate_rij();
    void update_rij();
    void calculate_fij();
    void update_fij();
};

#endif // JASTROW_H
