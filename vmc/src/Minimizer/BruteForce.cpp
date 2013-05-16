#include <src/Minimizer/BruteForce.h>

void BruteForce::runMinimizer(mat &minMax,
        vec &stepLengths,
        int &nCycles,
        const string &filename)
{
    vec gradient(nParameters);
    double energy, variance;
    ofstream ofile;
    if (myRank == 0)
    {
        ofile.open(filename.c_str());
    }
    for (double alpha = minMax(0,0); alpha <= minMax(0,1); alpha += stepLengths(0))
    {
        solver->setAlpha(alpha);
        for (double beta = minMax(1,0); beta <= minMax(1,1); beta += stepLengths(1))
        {
            solver->setBeta(beta);
            solver->runMonteCarloIntegration(nCycles);
            energy = solver->getEnergy();
            variance = solver->getVariance();
            gradient = solver->getVariationalGradient();
            if (myRank == 0)
            {
                ofile << alpha << " " << beta << " ";
                ofile << energy << " " << variance;
                for (int i = 0; i < nParameters; i++)
                    ofile << " " << gradient(i);
                ofile << endl;

                cout << "alpha = " << setw(6) << alpha;
                cout << ", beta = " << setw(6) << beta;
                cout << ", energy = " << setw(8) << energy;
                cout << ", variance = " << variance;
                cout << endl;
            }
        }
    }
}




