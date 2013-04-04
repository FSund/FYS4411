#include "Minimizer.h"

Minimizer::Minimizer(int &myRank, int &numprocs, int &nParticles, int &charge, int &nParameters, double* guess):
    myRank(myRank),
    numprocs(numprocs),
    nParameters(nParameters)
{
//    cout << "myRank = " << myRank << ", guess = " << guess << endl;

    solver = new SolverMCBF(myRank, numprocs, nParticles, charge);

    parameters = new double[nParameters];
    for (int i = 0; i < nParameters; i++)
        *parameters = guess[i];
}

double* Minimizer::runMinimizer()
{
    cout << "Minimizer::runMinimizer()" << endl;

    //

    return parameters;
}
