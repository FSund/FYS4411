#include <src/Minimizer/MoleculeMinimizer.h>

MoleculeMinimizer::MoleculeMinimizer(
        int &myRank,
        int &numprocs,
        int &nParameters,
        Solver *solver) :
    Minimizer(myRank, numprocs, nParameters, solver),
    minimizer(new SteepestDescent(myRank, numprocs, nParameters, solver))
{
}

MoleculeMinimizer::~MoleculeMinimizer()
{
    cout << "before delete" << endl;
    delete minimizer;
}

void MoleculeMinimizer::runMinimizer(int &nCycles, vec &guess, double minR, double maxR, double &dR)
{
    int nR = int((maxR - minR)/dR);
    vec parameters = guess;
    double R;
    int iterMax = 50;

    ofstream ofile;
    ofile.open("BF_moleculeMinimizer.dat");
    minimizer->setCout(false);
    minimizer->setPrintPath(false);

    int minimizerCycles = nCycles/10;

    for (int i = 0; i <= nR; i++)
    {
        R = minR + i*dR;
        solver->setR(R);
        parameters = minimizer->runMinimizer(parameters, minimizerCycles, iterMax);
        solver->runMonteCarloIntegration(nCycles);

        if (myRank == 0)
        {
            for (int j = 0; j < nParameters; j++)
                ofile << parameters(j) << " ";
            ofile << R << " ";
            ofile << solver->getEnergy();
            ofile << endl;
        }
    }

    cout << "before exit MoleculeMinimizer" << endl;
}
