#include <src/Minimizer/Minimizer.h>

Minimizer::Minimizer(
        int &myRank,
        int &numprocs,
        int &nParameters,
        Solver *solver):
    myRank(myRank),
    numprocs(numprocs),
    nParameters(nParameters),
//    nParticles(nParticles),
//    charge(charge),
    solver(solver)
{
    ofile.open("minimizingPath.dat");
    solver->setMinimizing(true);
}

void Minimizer::printToFile(const vec &parameters, const Solver *solver)
{
    for (uint i = 0; i < parameters.n_rows; i++)
        ofile << parameters(i) << " ";
    ofile << solver->getEnergy() << " ";
    ofile << solver->getVariance() << " ";
    ofile << endl;
}
