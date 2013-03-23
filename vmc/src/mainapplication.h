#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "CSolver.h"
using namespace std;
using namespace arma;

class MainApplication
{
public:
    MainApplication(int &myRank, int &numprocs);
    void runApplication();
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

#endif // MAINAPPLICATION_H
