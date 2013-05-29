#include <src/Minimizer/Filip.h>

Filip::Filip(
        int &myRank,
        int &numprocs,
        int &nParameters,
        Solver *solver):
    Minimizer(myRank, numprocs, nParameters, solver)
{
}

vec Filip::runMinimizer(vec &guess, int &nCycles)
{
    /* stochastic gradient descent */
    cout << "Starting Filip minimization" << endl;

    int iterMax = 100;
    double gamma0 = 1.0;
    double k = 0.75;
    vec gamma(nParameters);
    gamma.fill(gamma0);

    vec r(nParameters);
    vec rOld(nParameters);
    vec oldParam = guess;
    vec newParam = oldParam;
    double energy;
    double oldEnergy;
    double minEnergy = 0.0;
    double variance;
    vec minParameters(nParameters,1);

    solver->setParameters(oldParam);
    solver->runMonteCarloIntegration(nCycles);
    oldEnergy = solver->getEnergy();
    rOld = -solver->getVariationalGradient();
    printToFile(oldParam, solver);

    if (myRank == 0)
    {
        cout << "old param  = " << oldParam.t();
        cout << "old energy = " << oldEnergy << endl << endl;
    }

    for (int i = 0; i < iterMax; i++)
    {
        r = -solver->getVariationalGradient();
        for (int j = 0; j < nParameters; j++)
        {
            if ((rOld(j) < 0 && r(j) > 0) || (rOld(j) > 0 && r(j) < 0))
                gamma(j) *= k;
        }
        r = r%gamma;
        newParam = oldParam + r;

        for (int ii = 0; ii < nParameters; ii++)
        {
            if (abs(newParam(ii) - oldParam(ii))/oldParam(ii) > 0.3)
                newParam(ii) = oldParam(ii) + 0.3*oldParam(ii)*((newParam(ii) - oldParam(ii)) > 0 ? 1 : -1);
            if (newParam(ii) < 0.0)
                newParam(ii) = oldParam(ii)*0.5;
        }

        if (myRank == 0)
        {
            cout << "i = " << i << ", gamma = " << gamma.t();
            cout << "r                   = " << r.t();
            cout << "new param           = " << newParam.t();
        }

        solver->setParameters(newParam);
        solver->runMonteCarloIntegration(nCycles);
        energy = solver->getEnergy();
        variance = solver->getVariance();
        if (myRank == 0)
            printToFile(newParam, solver);

        if (energy < minEnergy)
        {
            minEnergy = energy;
            minParameters = newParam;
            if (myRank == 0)
            {
                cout << "new minEnergy       = " << minEnergy << endl;
                cout << "new min parameters  = " << minParameters.t();
            }
        }

        oldParam = newParam;
        oldEnergy = energy;

        if (myRank == 0)
        {
            cout << "energy              = " << solver->getEnergy() << endl;
            cout << "variance            = " << solver->getVariance() << endl;
            cout << "variance/E          = " << solver->getVariance()/solver->getEnergy() << endl;
            cout << endl;
        }

        rOld = r;
    }

    return oldParam;
}
