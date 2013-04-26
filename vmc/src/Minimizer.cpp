#include "Minimizer.h"

Minimizer::Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters, vec &guess):
    myRank(myRank),
    numprocs(numprocs),
    nParameters(nParameters),
    parameters(guess)
{
    solver = new SolverMCBF(myRank, numprocs, nParticles, charge);
//    solver = new SolverMCIS(myRank, numprocs, nParticles, charge);
}

vec Minimizer::runMinimizer()
{
    cout << "Minimizer::runMinimizer()" << endl;

//    bruteforce();

    vec minParam = zeros<vec>(nParameters);
//    minParam = steepestDescent();

    bruteforce();

    cout << "Exiting Minimizer::runMinimizer()" << endl;

    return minParam;
}

void Minimizer::bruteforce()
{
    double minAlpha = 3.4;
    double maxAlpha = 3.4;
    double minBeta = 0.5;
    double maxBeta = 6;
    double step = 0.05;
    int nStepsAlpha = round((maxAlpha - minAlpha)/step);
    int nStepsBeta = round((maxBeta - minBeta)/step);
    int nCycles = 1e5;

    double energy = 0.0;

//    cout << Wavefunction

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
//            if (myRank == 0) cout << "alpha = " << alpha << ", beta = " << beta << endl;

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

vec &Minimizer::energyGradientNumerical()
{
    const double h = 1e-3;
//    double h2 = 1e6;
    const int nCycles = 1e4;
    vec paramPlus(nParameters), paramMinus(nParameters);
    vec dE(nParameters);
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

//vec Minimizer::energyGradientNumerical(const vec &param)
//{
//    double h = 1e-3;
////    double h2 = 1e6;
//    int nCycles = 1e4;
//    vec paramPlus(nParameters), paramMinus(nParameters);
//    vec dE;
//    double energyCurrent, energyPlus, energyMinus;

//    paramPlus = paramMinus = param;

//    // computing the first derivative
//    energyCurrent = solver->runMonteCarloIntegration(nCycles);
//    for (int i = 0; i < nParameters; i++) {
//        paramPlus(i) += h;
//        solver->setParameters(paramPlus);
//        energyPlus = solver->runMonteCarloIntegration(nCycles);

//        paramMinus(i) -= h;
//        solver->setParameters(paramMinus);
//        energyMinus = solver->runMonteCarloIntegration(nCycles);

//        dE(i) = energyPlus - energyMinus;

//        paramPlus(i) = paramMinus(i) = param(i);
//    }
//    dE /= (2.0*energyCurrent*h);

//    return dE;
//}

mat &Minimizer::energyHessianNumerical()
{
    const double h = 1e-3;
    const double h2 = 1e6;
    const int nCycles = 1e4;
    vec param(nParameters);
    mat H = zeros<mat>(nParameters, nParameters);
    double energyCurrent, energy;
    int dx, dy;

    param = parameters;

    // computing the hessian matrix
    energyCurrent = solver->runMonteCarloIntegration(nCycles);
    for (int i = 0; i < nParameters; i++) {
        for (int j = 0; j < nParameters; j++) {
            if (i == j) // diagonal elements
            {
                for (dx = -1; dx <= 1; dx+=2)
                {
                    param(i) += dx*h;
                    solver->setParameters(param);
                    energy = solver->runMonteCarloIntegration(nCycles);
                    H(i,i) += energy;
                    param(i) = parameters(i);
                }
                H(i,i) = (H(i,i) - energyCurrent)/(h2*energyCurrent*energyCurrent);
            }
            else // off-diagonal elements
            {
                for (dx = -1; dx <= 1; dx+=2)
                {
                    for (dy = -1; dy <= 1; dy+=2)
                    {
                        param(i) += h*dx;
                        param(j) += h*dy;
                        solver->setParameters(param);

                        energy = solver->runMonteCarloIntegration(nCycles);
                        H(i,j) += energy*dx*dy;

                        param = parameters;
                    }
                }
                H(i,j) /= (4.0*h2*energyCurrent*energyCurrent);
            }
        }
    }

    return H;
}

vec Minimizer::steepestDescent()
{
    solver->setParameters(parameters);

    double tolerance = 1e-10;
    int iterMax = 10;
    int i;
    vec r(nParameters);
    vec oldParam = parameters;
    mat A(nParameters, nParameters);
    double alpha;

//    A = energyHessianNumerical();
    i = 0;
//    while (i <= iterMax || norm(oldParam-parameters, 2) < tolerance )
    while (i <= iterMax)
    {
        cout << "i = " << i << endl;

        r = -energyGradientNumerical();
        A = energyHessianNumerical();
        alpha = dot(r,r)/dot(r, A*r);
        parameters = oldParam + alpha*r;

        cout << "A = " << endl << A;
        cout << "r = " << endl << r;
        cout << "alpha = " << alpha << endl;
        cout << "alpha*r = " << endl << alpha*r;
        cout << "parameters = " << endl << parameters << endl;

        oldParam = parameters;

        i++;
    }

    return parameters;
}

//vec Minimizer::CGM()
//{
//    // http://en.wikipedia.org/wiki/Energy_minimization#Nonlinear_conjugate_gradient_method

//    vec paramOld(nParameters);
//    vec F(nParameters);
//    vec FOld(nParameters);
//    vec h(nParameters);
//    vec hOld(nParameters);
//    int iterMax = 1000;
//    paramOld = parameters;

//    F = FOld = energyGradientNumerical();
//    h = hOld = zeros<vec>(nParameters);
//    double kappa = 1.0;

//    int i = 0;
//    while (i <= iterMax || norm(f,2) < tolerance )
//    {
//        hOld = h;
//        paramOld = parameters;
//        parameters = paramOld + kappa*hOld;
//        F = energyGradientNumerical();
//        gamma = (F%F)/(FOld%FOld);
//        h = F + gamma*hOld;
//        i++;
//    }
//}

//vc Minimizer::CGM2()
//{
//    // http://en.wikipedia.org/wiki/Energy_minimization#Nonlinear_conjugate_gradient_method

//    vec paramOld(nParameters);
//    vec F(nParameters);
//    vec FOld(nParameters);
//    vec h(nParameters);
//    vec hOld(nParameters);
//    int iterMax = 1000;
//    paramOld = parameters;

//    F = FOld = energyGradientNumerical();
//    h = hOld = zeros<vec>(nParameters);
//    double kappa = 1.0;

//    int i = 0;
//    while (i <= iterMax || norm(f,2) < tolerance )
//    {
//        hOld = h;
//        paramOld = parameters;
//        parameters = paramOld + kappa*hOld;
//        F = energyGradientNumerical();
//        gamma = (F%F)/(FOld%FOld);
//        h = F + gamma*hOld;
//        i++;
//    }
//}
