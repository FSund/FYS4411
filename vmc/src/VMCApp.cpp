#include "VMCApp.h"

VMCApp::VMCApp(int &myRank, int &numprocs):
    myRank(myRank),
    numprocs(numprocs)
{
}

void VMCApp::runApplication()
{
    double alpha, beta;
    int nParticles, charge;

    // helium
    alpha = 2.0; // 1.8
    beta = 0.36;
    nParticles = 2;
    charge = nParticles;

    // beryllium
    alpha = 4.0; // 3.8
    beta = 0.1;
    nParticles = 4;
    charge = nParticles;

    // neon
    alpha = 10.6; // 10.6
    beta = 0.1;
    nParticles = 10;
    charge = nParticles;

    SolverMCIS solver(myRank, numprocs, nParticles, charge);
//    SolverMCBF solver(myRank, numprocs, nParticles, charge);
    solver.setAlpha(alpha);
    solver.setBeta(beta);

    int nCycles = 10;
    double energy = solver.runMonteCarloIntegration(nCycles);

    if (myRank == 0)
    {
        cout << "Energy = " << energy << endl;
        cout << "Variance = " << solver.getVariance() << endl;
        cout << "Acceptance ratio = " << solver.getAcceptanceRate() << endl;
    }
}

void VMCApp::minimize()
{
    vec guess(2);
    vec minParam(2);
    guess << 3.5 << 0.5;
    int nParameters = 2;
    int nParticles = 4;
    int charge = 4;

    Minimizer m(myRank, numprocs, nParticles, charge, nParameters);
    minParam = m.runMinimizer(guess);

    cout << "minParam = " << minParam.t() << endl;
}
