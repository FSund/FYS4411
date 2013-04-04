#ifndef VMCAPP_H
#define VMCAPP_H
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "Solver/Solver.h"
#include "Solver/SolverMCBF.h"
#include "Solver/SolverMCIS.h"
#include "Minimizer.h"

using namespace std;
using namespace arma;

class VMCApp
{
public:
    VMCApp(int &myRank, int &numprocs);
    void runApplication();
    void minimize();
private:
//    void variational_paramenters();
//    void minimum();
//    void optimal_steplength();
//    void steplength_secant();
//    void closedformBenchmark();
//    void importanceSampling();

//    void beryllium_variational_parameters();

    int myRank, numprocs;
};

#endif // VMCAPP_H