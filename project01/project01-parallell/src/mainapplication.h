#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "CVMCSolver.h"
#include "CHelium.h"
#include "CBeryllium.h"
using namespace std;
using namespace arma;

class MainApplication
{
public:
    MainApplication(int &my_rank_, int &numprocs_);
    void runApplication();
private:
    void variational_paramenters();
    void minimum();
    void optimal_steplength();
    void steplength_secant();
    void closedformBenchmark();
    void importanceSampling();

    void beryllium_variational_parameters();

    int my_rank, numprocs;
};

#endif // MAINAPPLICATION_H
