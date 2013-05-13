#include <src/VMCApp.h>

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
    alpha = 10.0; // 10.6
    beta = 0.2;
    nParticles = 10;
    charge = nParticles;

    SolverMCIS solver(myRank, numprocs, nParticles, charge);
//    SolverMCBF solver(myRank, numprocs, nParticles, charge);
    solver.setAlpha(alpha);
    solver.setBeta(beta);
//    solver.setBlocking(false);

    int nCycles = 1e5;
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
    int nParameters, nParticles, charge;
    nParameters = 2;
    vec guess(nParameters);
    vec minParam(nParameters);

    // beryllium
//    guess << 3.5 << .3;
//    nParticles = 4;
//    charge = nParticles;

    // neon
    guess << 9.0 << 0.5;
    nParticles = 10;
    charge = nParticles;

    int nCycles = 4e4;
    Minimizer m(myRank, numprocs, nParticles, charge, nParameters);
    minParam = m.runMinimizer(guess, nCycles);

    cout << "minParam = " << minParam.t() << endl;
}

