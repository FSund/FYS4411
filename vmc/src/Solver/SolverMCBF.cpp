#include <src/Solver/SolverMCBF.h>

SolverMCBF::SolverMCBF(int &myRank,
        int &numprocs,
        int &nParticles,
        int &charge):
    Solver(myRank, numprocs, nParticles, charge),
    stepLength(1.0)
{
}

void SolverMCBF::runCycle()
{
    // loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles+nThermalize; cycle++)
    {
        // New position to test
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            wf->updatePositionAndCurrentParticle(rNew, i);
            ratio = wf->getRatio(); // squared in Wavefunction

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= ratio)
            {
                rOld.row(i) = rNew.row(i);
                wf->acceptMove();
                nAccepted++;
            }
            else
            {
                rNew.row(i) = rOld.row(i);
                wf->rejectMove();
            }

            if (cycle >= nThermalize)
            {
                // update energies
                if (closedForm)
                    deltaE = wf->localEnergy();
                else
                    deltaE = wf->localEnergyNumerical();

                if (blocking)
                    logger->log(deltaE);

                energySum += deltaE;
                energySquaredSum += deltaE*deltaE;
            }
        }
    }
}


