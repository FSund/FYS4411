#include <src/Solver/Solver.h>

Solver::Solver(int &myRank,
        int &numprocs,
        int &nParticles,
        int &charge):
    nDimensions(3),
    nParticles(nParticles),
    charge(charge),
    rOld(zeros<mat>(nParticles, nDimensions)),
    rNew(zeros<mat>(nParticles, nDimensions)),
    nAccepted(0),
    myRank(myRank),
    numprocs(numprocs),
    nThermalize(1e3),
//    closedForm(false)
    closedForm(true),
//    blocking(false)
    blocking(true)
{
    idum = -1 - myRank;
    logger = new Datalogger(myRank, numprocs);

    if (nParticles == 2)
    {
        wf = new Wavefunction(nParticles, charge);
        cout << "Helium" << endl;
    }
    else if (nParticles == 4)
    {
        wf = new Wavefunction(nParticles, charge);
        cout << "Beryllium" << endl;
    }
    else if (nParticles == 10)
    {
        wf = new Wavefunction(nParticles, charge);
        cout << "Neon" << endl;
    }
    else
    {
        cout << "! Unknown element/number of particles, exiting." << endl;
        exit(1);
    }
}

Solver::~Solver()
{
    delete wf;
}

void Solver::setAlpha(const double &alpha)
{
    wf->setAlpha(alpha);
}

void Solver::setBeta(const double &beta)
{
    wf->setBeta(beta);
}

void Solver::setParameters(const vec &parameters)
{
    wf->setParameters(parameters);
}

double Solver::gaussianDeviate(long *seed)
{
    double R, randomNormal;
    // Box-Muller transform
    R = sqrt(-2.0*log(ran2(seed)));
    randomNormal = R*cos(2.0*pi*ran2(seed));
    return randomNormal;
}

double Solver::runMonteCarloIntegration(const int &nCycles_)
{
    energySum = 0;
    energySquaredSum = 0;
    nAccepted = 0;

    double dt = 1e-3;
    // initial trial positions
    for(int i = 0; i < nParticles; i++)
    {
        for(int j = 0; j < nDimensions; j++)
        {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(dt);
        }
    }
    rNew = rOld;
    wf->initialize(rOld);
    nCycles = nCycles_/numprocs;

    if (blocking)
        logger->initialize(nCycles, nParticles, "test2");

    runCycle();
    finalize();

    if (blocking)
        logger->finish();

    return energy;
}

void Solver::finalize()
{
    energy = energySum/(nCycles*nParticles);
    energySquared = energySquaredSum/(nCycles*nParticles);
    double totalEnergy = 0.0;
    double totalEnergySquared = 0.0;
    int totalNAccepted = 0;

    MPI_Allreduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&energySquared, &totalEnergySquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalNAccepted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    energy = totalEnergy/numprocs;
    energySquared = totalEnergySquared/numprocs;
    acceptanceRate = double(totalNAccepted)/double(nCycles*numprocs*nParticles);
    variance = energySquared - energy*energy;
}
