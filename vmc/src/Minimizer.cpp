#include "Minimizer.h"

Minimizer::Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters):
    myRank(myRank),
    numprocs(numprocs),
    nParameters(nParameters)
{
//    solver = new SolverMCBF(myRank, numprocs, nParticles, charge);
    solver = new SolverMCIS(myRank, numprocs, nParticles, charge);
}

vec Minimizer::runMinimizer(const vec &guess)
{
    cout << "Minimizer::runMinimizer()" << endl;

//    bruteforce();

    cout << "Exiting Minimizer::runMinimizer()" << endl;

//    return zeros<vec>(nParameters);
//    return steepestDescent(parameters);
    return newtonsMethod(guess);
}

void Minimizer::bruteforce()
{
    mat minmax(nParameters, 2);
    minmax << 3.0 << 4.1 << endr
           << 0.0 << 1.0 << endr;

    // benchmark: 56s -n 2 at home, 1m52s -n 1 at home
//    minmax << 3.4 << 3.4 << endr
//           << 0.1 << 0.2 << endr;

    vec step(nParameters);
    step << 0.05 << endr << 0.05;

    ofstream ofile;
    if (myRank == 0)
    {
        ofile.open("minimization.dat");
    }

    int nCycles = 1e5;
    double energy = 0.0;
    for (double alpha = minmax(0,0); alpha <= minmax(0,1); alpha += step(0))
    {
        solver->setAlpha(alpha);
        for (double beta = minmax(1,0); beta <= minmax(1,1); beta += step(1))
        {
            solver->setBeta(beta);
            energy = solver->runMonteCarloIntegration(nCycles);

            if (myRank == 0)
            {
                ofile << alpha << " " << beta << " " << energy << endl;
                cout << "alpha = " << setw(4) << alpha;
                cout << ", beta = " << setw(4) << beta;
                cout << ", energy = " << energy << endl;
            }
        }
    }
}

vec Minimizer::steepestDescent(const vec &param)
{
    solver->setParameters(param);

//    double tolerance = 1e-10;
    int iterMax = 100;
//    int i;
    vec r(nParameters);
    vec oldParam = param;
    vec newParam = oldParam;
//    mat A(nParameters, nParameters);
//    double alpha;
    double stepSize = 0.00001;

//    A = energyHessianNumerical();
//    i = 0;
//    while (i <= iterMax || norm(oldParam-parameters, 2) < tolerance )
    for (int i = 0; i < iterMax; i++)
    {
        r = -energyGradientNumerical(oldParam);
//        alpha = dot(r,r)/dot(r, A*r);
        newParam = oldParam + stepSize*r;

//        cout << "i = " << i << endl;
//        cout << "A = " << endl << A;
//        cout << "r = " << r.t();
//        cout << "alpha = " << alpha << endl;
//        cout << "alpha*r = " << endl << alpha*r;
        cout << "parameters = " << newParam.t() << endl;
        solver->setParameters(newParam);
        cout << "energy = " << solver->runMonteCarloIntegration(1e4);

        oldParam = newParam;

//        i++;
    }

    return parameters;
}

//vec Minimizer::steepestDescent()
//{
//    solver->setParameters(parameters);

////    double tolerance = 1e-10;
//    int iterMax = 10;
//    int i;
//    vec r(nParameters);
//    vec oldParam = parameters;
//    mat A(nParameters, nParameters);
//    double alpha;

////    A = energyHessianNumerical();
//    i = 0;
////    while (i <= iterMax || norm(oldParam-parameters, 2) < tolerance )
//    while (i <= iterMax)
//    {
//        cout << "i = " << i << endl;

//        r = -energyGradientNumerical();
//        A = energyHessianNumerical();
//        alpha = dot(r,r)/dot(r, A*r);
//        parameters = oldParam + alpha*r;

//        cout << "A = " << endl << A;
//        cout << "r = " << endl << r;
//        cout << "alpha = " << alpha << endl;
//        cout << "alpha*r = " << endl << alpha*r;
//        cout << "parameters = " << endl << parameters << endl;

//        oldParam = parameters;

//        i++;
//    }

//    return parameters;
//}

vec Minimizer::newtonsMethod(const vec &guess)
{
    double tolerance = 1e-10;
    int iterMax = 100;
    int nCycles = 1e4;

//    vec fpp = zeros<vec>(nParameters);
//    vec fp = zeros<vec>(nParameters);
    double fp, fpp;
    vec parampp = guess;
    vec paramp = guess*1.1;
    vec param = paramp;

    solver->setParameters(parampp);
    solver->runMonteCarloIntegration(nCycles);
    fpp = solver->getVariance();

    cout << "fpp = " << fpp << endl;

//    while () ???
    for (int i = 0; i < iterMax; i++)
    {
        for (int j = 0; j < nParameters; j++)
        {
            solver->setParameters(paramp);
            solver->runMonteCarloIntegration(nCycles);
            fp = solver->getVariance();

            cout << "fp = " << fp << endl;

            param(j) = paramp(j) - fp*(paramp(j) - parampp(j))/(fp - fpp);

            parampp(j) = paramp(j);
            paramp(j) = param(j);
            fpp = fp;

            cout << "parameters = " << param.t();
        }

//        cout << "parameters = " << param.t();
    }

    return param;
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

vec Minimizer::energyGradientNumerical(const vec &param)
{
    const double h = 1e-3;
//    double h2 = 1e6;
    const int nCycles = 1e4;
    vec paramPlus(nParameters), paramMinus(nParameters);
    vec dE(nParameters);
    double energyPlus, energyMinus;
    double factor = 2.0*h;

    paramPlus = paramMinus = param;
    // computing the first derivative
    for (int i = 0; i < nParameters; i++) {
        paramPlus(i) += h;
        solver->setParameters(paramPlus);
        energyPlus = solver->runMonteCarloIntegration(nCycles);

        paramMinus(i) -= h;
        solver->setParameters(paramMinus);
        energyMinus = solver->runMonteCarloIntegration(nCycles);

        dE(i) = (energyPlus - energyMinus)/factor;

        paramPlus(i) = paramMinus(i) = param(i);
    }

    return dE;
}

double Minimizer::energyGradientNumerical(const vec &param, const int &k)
{
    const double h = 1e-3;
//    double h2 = 1e6;
    const int nCycles = 1e4;
    vec paramPlus(nParameters), paramMinus(nParameters);
    double dE;
    double energyPlus, energyMinus;
    double factor = 2.0*h;

    paramPlus = paramMinus = param;
    // computing the first derivative
    paramPlus(k) += h;
    solver->setParameters(paramPlus);
    energyPlus = solver->runMonteCarloIntegration(nCycles);

    paramMinus(k) -= h;
    solver->setParameters(paramMinus);
    energyMinus = solver->runMonteCarloIntegration(nCycles);

    dE = (energyPlus - energyMinus)/factor;

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

mat Minimizer::energyHessianNumerical()
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
