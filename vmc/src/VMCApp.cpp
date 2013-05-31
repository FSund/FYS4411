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

    int nCycles = 1e7;

    // helium
    alpha = 1.8307; // 1.8
    beta = 0.4936;
    nParticles = 2;
    charge = nParticles;

//    // beryllium
//    alpha = 3.97;
//    beta = 0.094;
//    nParticles = 4;
//    charge = nParticles;

//    // neon
//    alpha = 10.26;
//    beta = 0.083;
//    nParticles = 10;
//    charge = nParticles;

    string orbitalType = "Hydrogenic";
//    SolverMCBF solver(myRank, numprocs, nParticles, charge, orbitalType);
    SolverMCIS solver(myRank, numprocs, nParticles, charge, orbitalType);

    ////
    solver.setAlpha(alpha);
    solver.setBeta(beta);
//    solver.setThermalizationSteps(1e5); // default = 1e5
//    solver.setMinimizing(false); // default = false
//    solver.setOnebody(true); // default = false
//    solver.setBlocking(true); // default = false
//    solver.setUseJastrow(false); // default = true
//    solver.setClosedform(false); // default = true
    ////

    solver.runMonteCarloIntegration(nCycles);

    if (myRank == 0)
    {
        cout << "Energy = " << solver.getEnergy() << endl;
        cout << "Variance = " << solver.getVariance() << endl;
        cout << "Acceptance ratio = " << solver.getAcceptanceRate() << endl;
        if (nParticles == 2)
            cout << "<r12> = " << solver.getr12() << endl;
    }
}

void VMCApp::runApplicationDiatomic()
{
    double alpha, beta, dist;
    int nParticles, charge;

    int nCycles = 1e4;

//    // H2
//    alpha = 1.29;
//    beta = 0.39;
//    nParticles = 2;
//    charge = 1;
//    dist = 1.4;

    // H2
    alpha = 1.29;
    beta = 0.39;
    nParticles = 2;
    charge = 1;
    dist = 2.45; // 2.45

    // Be2
    alpha = 4.0; // 10.6
    beta = 0.4;
    nParticles = 8;
    charge = 4;
    dist = 2.45;

    string orbitalType = "Diatomic";
//    SolverMCBF solver(myRank, numprocs, nParticles, charge, orbitalType);
    SolverMCIS solver(myRank, numprocs, nParticles, charge, orbitalType);

    ////
    solver.setAlpha(alpha);
    solver.setBeta(beta);
    solver.setR(dist);
//    solver.setThermalizationSteps(0); // default = 1e5
//    solver.setMinimizing(false); // default = false
//    solver.setOnebody(true); // default = false
//    solver.setBlocking(true); // default = false
//    solver.setUseJastrow(false); // default = true
//    solver.setClosedform(false); // default = true
    ////

    solver.runMonteCarloIntegration(nCycles);

    if (myRank == 0)
    {
        cout << "Energy = " << solver.getEnergy() << endl;
        cout << "Variance = " << solver.getVariance() << endl;
        cout << "Acceptance ratio = " << solver.getAcceptanceRate() << endl;
    }
}

void VMCApp::minimize()
{
    int nParameters, nParticles, charge;
    nParameters = 2;

    /* steepest descent thing */
    vec guess(nParameters);
    vec minParam(nParameters);

//    // helium
//    guess << 2.0 << 0.5;
//    nParticles = 2;
//    charge = nParticles;

//    // beryllium
//    nParticles = 4;
//    charge = nParticles;
//    guess << nParticles << 0.5;

    // neon
    nParticles = 10;
    charge = nParticles;
    guess << nParticles << 0.5;

    int nCycles = 1e6;
    string orbitalType = "Hydrogenic";
    Solver *solver = new SolverMCIS(myRank, numprocs, nParticles, charge, orbitalType);

    ////
    solver->setThermalizationSteps(1e5); // default = 1e5
    solver->setMinimizing(true);
//    solver->setOnebody(true); // default = false
//    solver->setBlocking(true); // default = false
    solver->setUseJastrow(true); // default = true
    solver->setClosedform(true); // default = true
    ////

    SteepestDescent *m = new SteepestDescent(myRank, numprocs, nParameters, solver);
    minParam = m->runMinimizer(guess, nCycles);
    if (myRank == 0) cout << "minParam = " << minParam.t() << endl;


    /* brute force*/
//    vec minMax(nParameters, 2);
//    vec stepLengths(nParameters);

    // beryllium
//    minMax << 3.6 << 4.1 << endr
//           << 0.0 << 0.3 << endr;

    // neon
//    minMax << 8.0 << 12.0 << endr
//           << 0.0 << 0.5 << endr;

//    stepLengths << 0.5 << endr << 0.05;
//    stepLengths << 0.025 << endr << 0.025;

//    BruteForce *m = new BruteForce(myRank, numprocs, nParameters, solver);
    //    m->runMinimizer(minMax, stepLengths, nCycles, filename);
}

void VMCApp::minimizeMolecules()
{
    int nParameters, nParticles, charge;
    double minR, maxR, dR;
    nParameters = 2;

    /* steepest descent thing */
    vec guess(nParameters);
    vec minParam(nParameters);

    // H2
    guess << 2.0 << 0.4;
    nParticles = 2;
    charge = 1;
    minR = 1.0;
    maxR = 1.8;
    dR = 0.005;

    // Be2
    guess << 2.0 << 0.4;
    nParticles = 8;
    charge = 4;
    minR = 0.5;
    maxR = 10;
    dR = 1.0;

    string orbitalType = "Diatomic";
    int nCycles = 1e5;
    Solver *solver = new SolverMCIS(myRank, numprocs, nParticles, charge, orbitalType);

    ////
//    solver->setThermalizationSteps(1e5); // default = 1e5
    solver->setMinimizing(true); // default = false
//    solver->setOnebody(true); // default = false
//    solver->setBlocking(true); // default = false
//    solver->setUseJastrow(false); // default = true
//    solver->setClosedform(false); // default = true
    int iterMax = 30;
    ////

    MoleculeMinimizer *m = new MoleculeMinimizer(myRank, numprocs, nParameters, solver);
    m->runMinimizer(nCycles, guess, minR, maxR, dR, iterMax);

    cout << "After runMinimizer" << endl;

    delete m;
    delete solver;
}

void VMCApp::minimizeMoleculesConstR()
{
    int nParameters, nParticles, charge;
    double R;
    nParameters = 2;

    /* steepest descent thing */
    vec guess(nParameters);
    vec minParam(nParameters);

    // H2
    guess << 2.0 << 0.4;
    nParticles = 2;
    charge = 1;
    R = 1.4;

    // Be2
    guess << 4 << 0.5;
    nParticles = 8;
    charge = 4;
    R = 2.45;

    string orbitalType = "Diatomic";
    int nCycles = 1e5;
    Solver *solver = new SolverMCIS(myRank, numprocs, nParticles, charge, orbitalType);

    ////
    solver->setR(R);
//    solver->setThermalizationSteps(1e5); // default = 1e5
    solver->setMinimizing(true); // default = false
//    solver->setOnebody(true); // default = false
//    solver->setBlocking(true); // default = false
//    solver->setUseJastrow(false); // default = true
//    solver->setClosedform(false); // default = true
    ////

    SteepestDescent m(myRank, numprocs, nParameters, solver);
    m.runMinimizer(guess, nCycles);

    delete solver;
}

void VMCApp::timing()
{
    double alpha, beta;
    int nParticles, charge;

    int nCycles = 1e5;

    // helium
    alpha = 1.8307; // 1.8
    beta = 0.4936;
    nParticles = 2;
    charge = nParticles;

    // beryllium
    alpha = 3.97;
    beta = 0.094;
    nParticles = 4;
    charge = nParticles;

    // neon
    alpha = 10.26;
    beta = 0.083;
    nParticles = 10;
    charge = nParticles;

    string orbitalType = "Hydrogenic";
//    SolverMCBF solver(myRank, numprocs, nParticles, charge, orbitalType);
    SolverMCIS solver(myRank, numprocs, nParticles, charge, orbitalType);

    ////
    solver.setAlpha(alpha);
    solver.setBeta(beta);
    solver.setThermalizationSteps(0); // default = 1e5
//    solver.setMinimizing(false); // default = false
//    solver.setOnebody(true); // default = false
//    solver.setBlocking(true); // default = false
//    solver.setUseJastrow(false); // default = true
    solver.setClosedform(false); // default = true
    ////

    solver.runMonteCarloIntegration(nCycles);

    if (myRank == 0)
    {
        cout << "Energy = " << solver.getEnergy() << endl;
        cout << "Variance = " << solver.getVariance() << endl;
        cout << "Acceptance ratio = " << solver.getAcceptanceRate() << endl;
    }
}
