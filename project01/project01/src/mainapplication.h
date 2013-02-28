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
    void test();
    void minimum();
};

#endif // MAINAPPLICATION_H
