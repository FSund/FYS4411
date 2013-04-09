#include "Minimizer.h"

Minimizer::Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters, vec &guess):
    myRank(myRank),
    numprocs(numprocs),
    nParameters(nParameters),
    parameters(guess)
{
//    solver = new SolverMCBF(myRank, numprocs, nParticles, charge);
    solver = new SolverMCIS(myRank, numprocs, nParticles, charge);
}

vec Minimizer::runMinimizer()
{
    cout << "Minimizer::runMinimizer()" << endl;

    bruteforce();

    cout << "Exiting Minimizer::runMinimizer()" << endl;

    return zeros<vec>(3);
}

void Minimizer::bruteforce()
{
    double minAlpha = 2.78;
    double maxAlpha = 6;
    double minBeta = 1.3456;
    double maxBeta = 6;
    double step = 0.0723;
    int nStepsAlpha = round((maxAlpha - minAlpha)/step);
    int nStepsBeta = round((maxBeta - minBeta)/step);
    int nCycles = 1e5;

    double energy = 0.0;

    ofstream ofile;
    if (myRank == 0)
    {
        ofile.open("minimization.dat");
    }

    for (double beta = minBeta; beta <= minBeta + step*nStepsBeta; beta += step)
    {
        solver->setBeta(beta);
        for (double alpha = minAlpha; alpha <= minAlpha + step*nStepsAlpha; alpha += step)
        {
            if (myRank == 0) cout << "alpha = " << alpha << ", beta = " << beta << endl;

            solver->setAlpha(alpha);
            energy = solver->runMonteCarloIntegration(nCycles);

            if (myRank == 0)
            {
                ofile << alpha << " " << beta << " " << energy << endl;
                cout << "alpha = " << alpha << ", beta = " << beta << ", energy = " << energy << endl;
            }
        }
    }
}

mat Minimizer::energyGradientNumerical()
{
    double h = 1e-3;
//    double h2 = 1e6;
    int nCycles = 1e4;
    vec paramPlus(nParameters), paramMinus(nParameters);
    mat dE;
    double energyCurrent, energyPlus, energyMinus;

    paramPlus = paramMinus = parameters;

    // computing the first derivative
    energyCurrent = solver->runMonteCarloIntegration(nCycles);
    for (int i = 0; i < nParameters; i++) {
        paramPlus(i) += h;
        solver->setParameters(paramPlus);
        energyPlus = solver->runMonteCarloIntegration(nCycles);

        paramMinus(i) -= h;
        solver->setParameters(paramMinus);
        energyMinus = solver->runMonteCarloIntegration(nCycles);

        dE(i) = energyPlus - energyMinus;

        paramPlus(i) = paramMinus(i) = parameters(i);
    }
    dE /= (2.0*energyCurrent*h);

    return dE;
}

vec Minimizer::steepestDescent()
{
//    solver->setParameters(parameters);

//    double tolerance = 1e-10;
//    int iterMax = 100;
//    int i;
//    vec f, x, z;
//    mat A;
//    double c, alpha;

//    x = parameters;
//    A = x;
//    f = A*x;
//    i = 0;
//    while (i <= iterMax || norm(f,2) < tolerance )
//    {
//        z = A*f;
//        c = dot(f,f);
//        alpha = c/dot(f,z);

//        x = x - alpha*f;
//        f = A*x;

//        i++;
//    }

//    return x;
}

