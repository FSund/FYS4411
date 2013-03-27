#include "CSolver.h"

Solver::Solver(
        int myRank,
        int numprocs,
        int nParticles,
        int charge):
    nDimensions(3),
    nParticles(nParticles),
    charge(charge),
    stepLength(1.0),
//    h(1e-3),
//    h2(1e6),
    nAccepted(0),
    nRejected(0),
    myRank(myRank),
    numprocs(numprocs)
{
    idum = -1 - myRank;

    rNew = mat(nParticles, nDimensions);
//    oldS.qForce = mat(nParticles, nDimensions);
    wfNew = 0.0;

    rOld = mat(nParticles, nDimensions);
//    newS.qForce = mat(nParticles, nDimensions);
    wfOld = 0.0;

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
    else
    {
        cout << "! Unknown element or number of particles, exiting." << endl;
        exit(1);
    }
}

Solver::~Solver()
{
}

//void Solver::setParameters(
////        const double &alpha_,
////        const double &beta_,
//        const double &stepLength_)
//{
////    alpha = alpha_;
////    beta = beta_;
//    stepLength = stepLength_;
//}

//void Solver::setParameters(
////        const double &alpha_,
////        const double &beta_,
//        const double &stepLength_,
////        const double &h_,
////        const double &h2_,
//        const bool &importanceSampling_,
//        const bool &closedForm_)
//{
////    alpha = alpha_;
////    beta = beta_;
//    stepLength = stepLength_;
////    h = h_;
////    h2 = h2_;
//    importanceSampling = importanceSampling_;
//    closedForm = closedForm_;
//}

void Solver::setAlpha(const double &alpha_)
{
    alpha = alpha_;
    wf->setAlpha(alpha);
}

void Solver::setBeta(const double &beta_)
{
    beta = beta_;
    wf->setBeta(beta);
}

double Solver::runMonteCarloIntegration(const int &nCycles_)
{
    nCycles = nCycles_;
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    wfOld = 0;
    wfNew = 0;

    double nAccepted = 0;

    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE;

    double ratio;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    wf->initialize(rOld);

    nCycles = nCycles/numprocs;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }
            wf->updatePositionAndCurrentParticle(rNew, i);
            ratio = wf->getRatio();
            ratio *= ratio; // should be squared!

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= ratio) {
                rOld.row(i) = rNew.row(i);
                wf->acceptMove();
                nAccepted++;
            } else {
                rNew.row(i) = rOld.row(i);
                wf->rejectMove();
            }
            // update energies
            deltaE = wf->localEnergyNumerical();

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles*nParticles);
    double energySquared = energySquaredSum/(nCycles*nParticles);
    double totalEnergy = 0.0;
    double totalEnergySquared = 0.0;

    MPI_Reduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&energySquared, &totalEnergySquared, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    energy = totalEnergy/numprocs;
    energySquared = totalEnergySquared/numprocs;

//    if (myRank == 0) cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;

    return energy;
}

double Solver::gaussianDeviate(long *seed)
{
    double R, randomNormal;
    // Box-Muller transform
//    randomUniform << ran2(seed) << ran2(seed);
//    R = sqrt(-2*log(randomUniform(0)));
//    randomNormal(0) = R*cos(2*pi*randomUniform(1));
//    randomNormal(1) = R*sin(2*pi*randomUniform(1))

    R = sqrt(-2.0*log(ran2(seed)));
    randomNormal = R*cos(2.0*pi*ran2(seed));
    return randomNormal;
}
