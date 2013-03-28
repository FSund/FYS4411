#include "SolverMCIS.h"

SolverMCIS::SolverMCIS(int &myRank,
        int &numprocs,
        int &nParticles,
        int &charge):
    Solver(myRank, numprocs, nParticles, charge),
    qForceOld(zeros<mat>(nParticles, nDimensions)),
    qForceNew(zeros<mat>(nParticles, nDimensions)),
    D(0.5),
    dt(1e-3),
    Ddt(D*dt)
{
}

double SolverMCIS::runMonteCarloIntegration(const int &nCycles_)
{
    nCycles = nCycles_;

    double nAccepted = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE;
    double ratio;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(dt);
        }
    }
    rNew = rOld;
    wf->initialize(rOld);

    qForceOld = 2.0*wf->gradientNumerical();

    nCycles = nCycles/numprocs;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // loop over particles
        for(int i = 0; i < nParticles; i++) {
            // New position to test
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(dt)
                        + qForceOld(i,j)*Ddt;
            }
            wf->updatePositionAndCurrentParticle(rNew, i);
            ratio = wf->getRatio();
            ratio *= ratio; // should be squared!

            qForceNew = 2.0*wf->gradientNumerical();
            // Green's function ratio
            omegaRatio = 0.0;
            for (int j = 0; j < nDimensions; j++)
            {
                omegaRatio += (qForceOld(i,j) + qForceNew(i,j))*
                        (0.5*Ddt*(qForceOld(i,j) - qForceNew(i,j)) + rOld(i,j) - rNew(i,j));
            }
            omegaRatio = exp(0.5*omegaRatio);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= omegaRatio*ratio) {
                rOld.row(i) = rNew.row(i);
                qForceOld.row(i) = qForceNew.row(i);
                wf->acceptMove();

                nAccepted++;
            } else {
                rNew.row(i) = rOld.row(i);
                qForceNew.row(i) = qForceOld.row(i);
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

double SolverMCIS::gaussianDeviate(long *seed)
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
