#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "CVMCSolver.h"

using namespace std;
using namespace arma;

class MainApplication
{
public:
    MainApplication();
    void runApplication();
private:
    void variational_paramenters();
    void minimum();
    void optimal_steplength();
    void steplength_secant();
    void closedformBenchmark();
    void importanceSampling();

};

#endif // MAINAPPLICATION_H
