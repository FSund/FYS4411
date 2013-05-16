#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

#include <armadillo>
#include <src/Minimizer/Minimizer.h>

using namespace std;
using namespace arma;

class BruteForce : public Minimizer
{
public:
    BruteForce(
            int &myRank,
            int &numprocs,
            int &nParameters,
            Solver *solver);
    void runMinimizer(
            mat &minMax,
            vec &stepLengths,
            int &nCycles,
            const string &filename="bruteforce_minimization.dat");
};

#endif // BRUTEFORCE_H
