#include "mainapplication.h"

MainApplication::MainApplication(int &myRank, int &numprocs):
    myRank(myRank),
    numprocs(numprocs)
{
}

void MainApplication::runApplication()
{
//    double alpha = 1.8;
//    double beta = 0.36;
//    int nParticles = 2;
//    int charge = 2;

    double alpha = 3.9;
    double beta = 4.0;
    int nParticles = 4;
    int charge = 4;

    Solver solver(myRank, numprocs, nParticles, charge);
    solver.setAlpha(alpha);
    solver.setBeta(beta);


    double energy = solver.runMonteCarloIntegration(1e5);

    if (myRank == 0) cout << "Energy = " << energy << endl;
}
