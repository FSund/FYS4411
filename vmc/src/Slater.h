#ifndef SLATER_H
#define SLATER_H

#include <armadillo>

using namespace std;
using namespace arma;

class Slater
{
public:
    Slater(const int &nParticles);
    void initalize(const mat &r);
    void updatePositionAndCurrentParticle(mat &r, int &k);
    void setAlpha(const double &newAlpha);

    double wavefunction();
    double wavefunction(const mat &r);
    double getRatio();
    mat gradient();
    double laplacian();

    void acceptMove();
    void rejectMove();

    double hydrogenWF(const int &i, const vec3 &rvec);
    double phi1s(const vec3 &rvec);
    double phi2s(const vec3 &rvec);
private:
    int nParticles;
    int nDimensions;
    mat rOld, rNew;
    double alpha;
    int currentParticle;
};

#endif // SLATER_H
