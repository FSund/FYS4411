#include <src/Solver/SolverMCBF.h>

SolverMCBF::SolverMCBF()
{
    cout << "! Error: Using default constructor in SolverMCBF! " << endl;
    exit(1);
}

SolverMCBF::SolverMCBF(int &myRank,
        int &numprocs,
        int &nParticles,
        int &charge,
        string &orbitalType):
    Solver(myRank, numprocs, nParticles, charge, orbitalType),
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
                if (cycle > nThermalize) nAccepted++;
            }
            else
            {
                rNew.row(i) = rOld.row(i);
                wf->rejectMove();
            }
        }
        if (cycle >= nThermalize)
        {
            deltaE = localEnergy->evaluate(rNew, wf);

            if (blocking)
                logger->log(deltaE);

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
            if (minimizing)
            {
                tempVariationalGradient = wf->variationalDerivatives();
                variationalGradientSum += tempVariationalGradient;
                variationalGradientESum += tempVariationalGradient*deltaE;
            }
        }
    }
}


