#ifndef MOLECULEMINIMIZER_H
#define MOLECULEMINIMIZER_H

#include <armadillo>
#include <src/Minimizer/Minimizer.h>
#include <src/Minimizer/SteepestDescent.h>

using namespace arma;
using namespace std;

class MoleculeMinimizer : public Minimizer
{
public:
    MoleculeMinimizer(
            int &myRank,
            int &numprocs,
            int &nParameters,
            Solver *solver);
    ~MoleculeMinimizer();
    void runMinimizer(int &nCycles, vec &guess, double minR, double maxR, double &dR, int iterMax_);
private:
    SteepestDescent *minimizer;
};

#endif // MOLECULEMINIMIZER_H
