#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "Solver/Solver.h"
#include "Solver/SolverMCBF.h"
#include "Solver/SolverMCIS.h"

using namespace std;
using namespace arma;

class Minimizer
{
public:
    Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters);
    vec runMinimizer(const vec &guess, const int &nCycles);
protected:
    int myRank, numprocs;
    int nParameters;
    vec parameters;

    Solver* solver;
private:
    void bruteforce(const int &nCycles);

    vec energyGradientNumerical(const vec &param);
    double energyGradientNumerical(const vec &param, const int &k);
//    vec energyGradientNumerical(const vec &param);
    mat energyHessianNumerical();
    vec steepestDescent(const vec &param);
    vec newtonsMethod(const vec &param);
//    vec CGM();
};

#endif // MINIMIZER_H
