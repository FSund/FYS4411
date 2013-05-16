#include <src/Minimizer/SteepestDescent.h>

SteepestDescent::SteepestDescent(
        int &myRank,
        int &numprocs,
        int &nParameters,
        Solver *solver):
    Minimizer(myRank, numprocs, nParameters, solver)
{
}

vec SteepestDescent::runMinimizer(vec &guess, int &nCycles)
{
    /* stochastic gradient descent */

//    double tolerance = 1e-10;
    int iterMax = 100;
//    int nCycles = 1e6;
    int n = 10;

    vec r(nParameters);
    vec oldParam = guess;
    vec newParam = oldParam;
    double energy;
    double oldEnergy;
    double minEnergy = 0.0;
    vec minParameters(nParameters,1);
    vec energyVec(1);
    mat parameterMat(nParameters, 1);

//    double minEnergyDev;
    double meanE = 0.0;
    double stdE = 100;

    double gamma0 = 0.5;
    double k = 0.75;
    double gamma;

    solver->setParameters(oldParam);
    solver->runMonteCarloIntegration(nCycles);
    oldEnergy = solver->getEnergy();
    printToFile(oldParam, solver);

    if (myRank == 0)
    {
        cout << "old param  = " << oldParam.t();
        cout << "old energy = " << oldEnergy << endl << endl;
    }

    for (int i = 0; i < iterMax; i++)
    {
        gamma = gamma0*pow((i+1),-k);
        r = -solver->getVariationalGradient();

        if (myRank == 0)
        {
            cout << "i = " << i << ", gamma = " << gamma << endl;
        }

        r = r*gamma;

        if (myRank == 0)
            cout << "r                   = " << r.t();

        newParam = oldParam + r;

        for (int ii = 0; ii < nParameters; ii++)
        {
            if (abs(newParam(ii) - oldParam(ii))/oldParam(ii) > 0.3)
                newParam(ii) = oldParam(ii) + 0.3*oldParam(ii)*((newParam(ii) - oldParam(ii)) > 0 ? 1 : -1);
            if (newParam(ii) < 0.0)
                newParam(ii) = -newParam(ii)*0.5;
        }

        solver->setParameters(newParam);
        solver->runMonteCarloIntegration(nCycles);
        energy = solver->getEnergy();
        printToFile(newParam, solver);

        parameterMat.resize(nParameters, i+1);
        parameterMat.col(i) = newParam;
        energyVec.resize(i+1, 1);
        energyVec(i) = energy;

        if (energy < minEnergy)
        {
            minEnergy = energy;
            minParameters = newParam;
            if (myRank == 0)
            {
//                cout << endl;
                cout << "new minEnergy       = " << minEnergy << endl;
                cout << "new min parameters  = " << minParameters.t();
//                cout << endl;
            }
        }

        oldParam = newParam;
        oldEnergy = energy;

        meanE = mean(energyVec.rows(i<n ? 0 : (i-n), i));
        stdE = stddev(energyVec.rows(i<n ? 0 : (i-n), i));
//        minEnergyDev = abs(mean/minEnergy - 1.0);

        if (myRank == 0)
        {
            cout << "new param           = " << newParam.t();
            cout << "energy              = " << energy << endl;
            cout << "variance            = " << solver->getVariance() << endl;
            cout << "test " << abs((meanE - minEnergy)/meanE) << " < " << gamma*stdE/abs(meanE) << endl;
            cout << endl;
        }

        if (i > n && abs((meanE - minEnergy)/meanE) < gamma*stdE/abs(meanE))
        {
            if (myRank == 0)
            {
                cout << "------------- Exiting minimization ----------------" << endl;
                cout << "minimum energy      = " << minEnergy << endl;
                cout << "min parameters      = " << minParameters.t();
                cout << "return parameters   = " << mean(parameterMat.cols(i-n,i), 1).t();
                cout << endl;
            }
//            return mean(parameterMat.cols(i-n,i), 1);
        }
    }

    return oldParam;
}

