#include "SolverMCBF.h"

SolverMCBF::SolverMCBF(int &myRank,
        int &numprocs,
        int &nParticles,
        int &charge):
    Solver(myRank, numprocs, nParticles, charge),
    stepLength(1.0)
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

double SolverMCBF::runMonteCarloIntegration(const int &nCycles_)
{
//    nCycles = nCycles_;

    energySum = 0;
    energySquaredSum = 0;
    nAccepted = 0;

//    // initial trial positions
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = 0; j < nDimensions; j++) {
//            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
//        }
//    }

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

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++)
    {

        // New position to test
        for(int i = 0; i < nParticles; i++)
        {
            for(int j = 0; j < nDimensions; j++)
            {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            wf->updatePositionAndCurrentParticle(rNew, i);
            ratio = wf->getRatio(); // squared in Wavefunction

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= ratio)
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
            // update energies
            if (closedForm)
                deltaE = wf->localEnergy();
            else
                deltaE = wf->localEnergyNumerical();

//            wf->localEnergy();
//            cout << "deltaE   NUM = " << wf->localEnergyNumerical() << endl;
//            cout << "deltaE    CF = " << wf->localEnergy() << endl;
//            cout << endl;

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }

    finalize();

    return energy;
}


