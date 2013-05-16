#include <src/Minimizer/Minimizer.h>

Minimizer::Minimizer(
        int &myRank,
        int &numprocs,
        int &nParameters,
        Solver *solver):
    myRank(myRank),
    numprocs(numprocs),
//    nParticles(nParticles),
//    charge(charge),
    solver(solver)
{
    ofile.open("minimizingPath.dat");
}

void Minimizer::printToFile(const vec &parameters, const Solver *solver)
{
    for (int i = 0; i < parameters.n_rows; i++)
        ofile << parameters(i) << " ";
    ofile << solver->getEnergy() << " ";
    ofile << solver->getVariance() << " ";
    ofile << endl;
}
