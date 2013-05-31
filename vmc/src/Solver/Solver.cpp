#include <src/Solver/Solver.h>

Solver::Solver()
{
    cout << "! Error: Using default constructor in Solver! " << endl;
    exit(1);
}

Solver::Solver(
        const int &myRank,
        const int &numprocs,
        const int &nParticles,
        const int &charge,
        const string &orbitalType
        ):
    nDimensions(3),
    nParticles(nParticles),
    charge(charge),
    nParameters(2),
    rOld(mat(nParticles, nDimensions)),
    rNew(mat(nParticles, nDimensions)),
    tempVariationalGradient(vec(nParameters)),
    variationalGradient(vec(nParameters)),
    variationalGradientSum(vec(nParameters)),
    variationalGradientESum(vec(nParameters)),
    nAccepted(0),
    myRank(myRank),
    numprocs(numprocs),
    nThermalize(1e5),
    blocking(false),
    minimizing(false),
    onebody(false)
{
    idum = -1 - myRank;
    logger = new Datalogger(myRank, numprocs);
    onebodylogger = new OneBodyLogger(myRank, numprocs, nParticles);

    if (orbitalType == "Hydrogenic")
    {
        localEnergy = new SingleAtomLocalEnergy(nParticles, nDimensions, charge);
        if (nParticles == 2)
        {
            wf = new Wavefunction(nParticles, charge, orbitalType);
            cout << "Helium" << endl;
        }
        else if (nParticles == 4)
        {
            wf = new Wavefunction(nParticles, charge, orbitalType);
            cout << "Beryllium" << endl;
        }
        else if (nParticles == 10)
        {
            wf = new Wavefunction(nParticles, charge, orbitalType);
            cout << "Neon" << endl;
        }
        else
        {
            cout << "! Unknown element/number of particles, exiting." << endl;
            exit(1);
        }
    }
    if (orbitalType == "Diatomic")
    {
        localEnergy = new DiatomicLocalEnergy(nParticles, nDimensions, charge);
        if (nParticles == 2)
        {
            wf = new Wavefunction(nParticles, charge, orbitalType);
            cout << "H2" << endl;
        }
        else if (nParticles == 8)
        {
            wf = new Wavefunction(nParticles, charge, orbitalType);
            cout << "Be2" << endl;
        }
        else
        {
            cout << "! Unknown element/number of particles, exiting." << endl;
            exit(1);
        }
    }
}

Solver::~Solver()
{
    delete wf;
    delete localEnergy;
    delete logger;
    delete onebodylogger;
}

void Solver::setAlpha(const double &alpha)
{
    wf->setAlpha(alpha);
}

void Solver::setBeta(const double &beta)
{
    wf->setBeta(beta);
}

void Solver::setR(const double &dist)
{
    localEnergy->setR(dist);
    wf->setR(dist);
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
    r12Sum = 0.0;
    variationalGradientSum = zeros<vec>(nParameters);
    variationalGradientESum = zeros<vec>(nParameters);
    nAccepted = 0;

    double dt = 5e-3;
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
        logger->initialize(nCycles, nParticles);
    if (onebody)
        onebodylogger->initialize(nCycles);

    runCycle();

    if (blocking)
        logger->finish();
    if (onebody)
        onebodylogger->finish();

    finalize();

    return energy;
}

void Solver::finalize()
{
    energy = energySum/nCycles;
    energySquared = energySquaredSum/nCycles;
    r12 = r12Sum/nCycles;
    double totalEnergy = 0.0;
    double totalEnergySquared = 0.0;
    double totalr12 = 0.0;
    int totalNAccepted = 0;

    MPI_Allreduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&energySquared, &totalEnergySquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&r12, &totalr12, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalNAccepted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    energy = totalEnergy/numprocs;
    energySquared = totalEnergySquared/numprocs;
    r12 = totalr12/numprocs;
    acceptanceRate = double(totalNAccepted)/double(nCycles*numprocs*nParticles);
    variance = (energySquared - energy*energy)/double(nCycles*numprocs);
//    variance = energySquared - energy*energy;

    if (minimizing)
    {
        vec totalVariationalGradient = zeros<vec>(nParameters);
        vec totalVariationalGradientE = zeros<vec>(nParameters);
        for (int i = 0; i < nParameters; i++)
        {
            MPI_Allreduce(&variationalGradientSum(i), &totalVariationalGradient(i), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&variationalGradientESum(i), &totalVariationalGradientE(i), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

        totalVariationalGradient /= (nCycles*numprocs);
        totalVariationalGradientE /= (nCycles*numprocs);
        variationalGradient = 2.0*(totalVariationalGradientE - energy*totalVariationalGradient);
    }
}
