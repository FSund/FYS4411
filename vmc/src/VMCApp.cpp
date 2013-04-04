#include "VMCApp.h"

VMCApp::VMCApp(int &myRank, int &numprocs):
    myRank(myRank),
    numprocs(numprocs)
{
}

void VMCApp::runApplication()
{
//    double alpha = 1.8;
//    double beta = 0.36;
//    int nParticles = 2;
//    int charge = 2;

    double alpha = 3.9;
    double beta = 4.0;
    int nParticles = 4;
    int charge = 4;

//    SolverMCBF solver(myRank, numprocs, nParticles, charge);
//    solver.setAlpha(alpha);
//    solver.setBeta(beta);

//    double energy = solver.runMonteCarloIntegration(1e5);

//    if (myRank == 0) cout << "Energy = " << energy << endl;



    SolverMCIS solver(myRank, numprocs, nParticles, charge);
    solver.setAlpha(alpha);
    solver.setBeta(beta);

    double energy = solver.runMonteCarloIntegration(1e5);

    if (myRank == 0) cout << "Energy = " << energy << endl;
}

void VMCApp::minimize()
{
    double* guess = new double[2];
    guess[0] = 1.2; guess[1] = 4.5;
    int nParameters = 2;
    int nParticles = 4;
    int charge = 4;

    Minimizer m(myRank, numprocs, nParticles, charge, nParameters, guess);
    m.runMinimizer();
}
