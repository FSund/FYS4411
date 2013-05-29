#include <src/Solver/SolverMCIS.h>

SolverMCIS::SolverMCIS()
{
    cout << "! Error: Using default constructor in SolverMCIS! " << endl;
    exit(1);
}

SolverMCIS::SolverMCIS(int &myRank,
        int &numprocs,
        int &nParticles,
        int &charge,
        string &orbitalType):
    Solver(myRank, numprocs, nParticles, charge, orbitalType),
    qForceOld(zeros<mat>(nParticles, nDimensions)),
    qForceNew(zeros<mat>(nParticles, nDimensions)),
    D(0.5),
    dt(5e-3),
    Ddt(D*dt)
{
}

void SolverMCIS::runCycle()
{
    if (closedForm)
        qForceOld = 2.0*wf->localGradient(); // rOld == rNew in Wavefunction
    else
        qForceOld = 2.0*wf->localGradientNumerical(rOld);

    // loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles+nThermalize; cycle++)
    {
        // loop over particles
        for (int i = 0; i < nParticles; i++)
        {
            // New position to test
            for (int j = 0; j < nDimensions; j++)
            {
                rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(dt)
                        + qForceOld(i,j)*Ddt;
            }
            wf->updatePositionAndCurrentParticle(rNew, i);
            ratio = wf->getRatio(); // squared in Wavefunction

            if (closedForm)
                qForceNew = 2.0*wf->localGradient();
            else
                qForceNew = 2.0*wf->localGradientNumerical();

            // Green's function ratio
            omegaRatio = 0.0;
            for (int k = 0; k < nParticles; k++)
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    omegaRatio += (qForceOld(k,j) + qForceNew(k,j))*
                            (0.5*Ddt*(qForceOld(k,j) - qForceNew(k,j)) + rOld(k,j) - rNew(k,j));
                }
            }
            omegaRatio = exp(0.5*omegaRatio);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= omegaRatio*ratio)
            {
                rOld.row(i) = rNew.row(i);
                qForceOld = qForceNew;
                wf->acceptMove();
                if (cycle > nThermalize) nAccepted++;
            }
            else
            {
                rNew.row(i) = rOld.row(i);
                qForceNew = qForceOld;
                wf->rejectMove();
            }
        }
        if (cycle >= nThermalize)
        {
            // update energies
//            if (closedForm)
//                deltaE = localEnergy->evaluate(rNew, wf);
//            else
                deltaE = localEnergy->evaluate(rNew, wf);

//                cout << "deltaE = " << deltaE << endl;
//                cout << endl;

            if (blocking)
                logger->log(deltaE);

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
            if (nParticles == 2)
            {
                r12 = 0.0;
                for (int ii = 0; ii < nDimensions; ii++)
                    r12 += (rOld(0,ii) - rOld(1,ii))*(rOld(0,ii) - rOld(1,ii));
                r12Sum += sqrt(r12);
            }
            if (minimizing)
            {
                tempVariationalGradient = wf->variationalDerivatives();
                variationalGradientSum += tempVariationalGradient;
                variationalGradientESum += tempVariationalGradient*deltaE;
            }
            if (onebody)
            {
                onebodylogger->log(rOld);
            }
        }
    }
}
