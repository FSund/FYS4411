#include "CVMCSolver.h"
#include "lib.h"
#include "mpi.h"
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(0.5*charge),
    beta(0.5*charge),
    nCycles(1000000),
    nAccepted(0),
    nRejected(0)
{
    rOld(nParticles, nDimensions);
    rNew(nParticles, nDimensions);
    int argc = 1;
    char** argv;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
}

double VMCSolver::runMonteCarloIntegration(
        const int &nCycles_,
        const double &stepLength_,
        const double &alpha_,
        const double &beta_,
        const bool closedform)
{
    rOld.zeros(nParticles, nDimensions);
    rNew.zeros(nParticles, nDimensions);

    nCycles = nCycles_;
    alpha = alpha_;
    beta = beta_;
    stepLength = stepLength_;
    nAccepted = 0;
    nRejected = 0;

    double r12 = 0;
    double waveFunctionOld = 0, waveFunctionNew = 0;
    double energySum = 0, energySquaredSum = 0;
    double exactEnergySum = 0, exactEnergySquaredSum = 0;
    double deltaE = 0;

    // MPI stuff
    local_nCycles = int(nCycles/numprocs);
    int my_start = my_rank*local_nCycles;
    int my_end = my_start + local_nCycles;
    idum -= my_rank;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // loop over Monte Carlo cycles
    for (int cycle = my_start; cycle < my_end; cycle++) {
    //for (int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = waveFunction(rOld);

        for (int i = 0; i < nParticles; i++) {
            // New position to test
            for (int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = waveFunction(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                rOld.row(i) = rNew.row(i); // update position
                waveFunctionOld = waveFunctionNew;
                nAccepted++;
            } else {
                rNew = rOld; // reset position, throw away the test-position
                nRejected++;
            }

            // update energies
            if (closedform)
                deltaE = exactLocalEnergy(rNew);
            else
                deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        // calculating the optimal value of the distance between the particles
//        r12Sum += dr(rNew, 0, 1);
        r12 += norm(rNew.row(0), 2);
    }

//    cout << "My rank = " << my_rank << ", nAccepted = " << nAccepted << endl;
//    cout << "My rank = " << my_rank << ", nRejected = " << nRejected << endl;


    double energy = energySum/(local_nCycles*nParticles);
    double energySquared = energySquaredSum/(local_nCycles*nParticles);
    r12 = r12/local_nCycles;

    double total_energy = 0.0;
    double total_energySquared = 0.0;
    double total_r12 = 0.0;
    int totalnAcc, totalnRej;
    MPI_Reduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&energySquared, &total_energySquared, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&r12, &total_r12, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&r12, &total_r12, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalnAcc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nRejected, &totalnRej, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (my_rank == 0) {
//        cout << "Energy            = " << total_energy/numprocs << endl;
//        cout << "Energy squared    = " << total_energySquared/numprocs << endl;
//        cout << "r12 average       = " << total_r12/numprocs << endl;
//        cout << "nAccepted   = " << nAccepted << endl;
//        cout << "nRejected   = " << nRejected << endl;
//        cout << "Total tests = " << nAccepted + nRejected << endl;
//        cout << "nAcc = " << totalnAcc << endl;
//        cout << "nRej = " << totalnRej << endl;
//        cout << endl;
    }

    return energy;
}

double VMCSolver::runMonteCarloIntegrationImportanceSampling(
        const int &nCycles_,
        const double &stepLength_,
        const double &alpha_,
        const double &beta_,
        const bool closedform,
        const double timestep)
{
    rOld.zeros(nParticles, nDimensions);
    rNew.zeros(nParticles, nDimensions);

    nCycles = nCycles_;
    alpha = alpha_;
    beta = beta_;
    stepLength = stepLength_;
    nAccepted = 0;
    nRejected = 0;

    double r12 = 0;
    double waveFunctionOld = 0, waveFunctionNew = 0;
    double energySum = 0, energySquaredSum = 0;
    double exactEnergySum = 0, exactEnergySquaredSum = 0;
    double deltaE = 0;
    double greensfunction;
    mat qforceOld(nParticles, nDimensions);
    mat qforceNew(nParticles, nDimensions);

    // MPI stuff
    local_nCycles = int(nCycles/numprocs);
    int my_start = my_rank*local_nCycles;
    int my_end = my_start + local_nCycles;
    idum -= my_rank;

    // initial trial positions
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            // rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(timestep);
        }
    }

    rNew = rOld;
    waveFunctionOld = waveFunction(rOld);
    qforceOld = quantumForce(rOld, waveFunctionOld);

    // loop over Monte Carlo cycles
    for (int cycle = my_start; cycle < my_end; cycle++) {
    //for (int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = waveFunction(rOld);

        for (int i = 0; i < nParticles; i++) {
            // New position to test
            for (int j = 0; j < nDimensions; j++) {
                // rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
                rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(timestep)
                        + qforceOld(i,j)*timestep;
            }

            // Recalculate the value of the wave function
            waveFunctionNew = waveFunction(rNew);
            qforceNew = quantumForce(rNew, waveFunctionNew);

            // log of the ratio of the greens function
            greensfunction = 0.0;
            for (int j = 0; j < nDimensions; j++) {
                greensfunction += 0.5*(qforceOld(i,j) + qforceNew(i,j))
                        *(timestep*0.5*(qforceOld(i,j) - qforceNew(i,j))
                          - rNew(i,j) + rOld(i,j));
            }
//            greensfunction += -(rNew - rOld - qforceOld)

            // Check for step acceptance (if yes, update position, if no, reset position)
//            if (ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
            if (ran2(&idum) <= greensfunction*waveFunctionNew*waveFunctionNew/(waveFunctionOld*waveFunctionOld)) {
                rOld.row(i) = rNew.row(i); // update position
                qforceOld.row(i) = qforceNew.row(i);
                waveFunctionOld = waveFunctionNew;
                nAccepted++;
            } else {
                rNew = rOld; // reset position, throw away the test-position
                nRejected++;
            }

            // update energies
            if (closedform)
                deltaE = exactLocalEnergy(rNew);
            else
                deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        // calculating the optimal value of the distance between the particles
//        r12Sum += dr(rNew, 0, 1);
        r12 += norm(rNew.row(0), 2);
    }

//    cout << "My rank = " << my_rank << ", nAccepted = " << nAccepted << endl;
//    cout << "My rank = " << my_rank << ", nRejected = " << nRejected << endl;


    double energy = energySum/(local_nCycles*nParticles);
    double energySquared = energySquaredSum/(local_nCycles*nParticles);
    r12 = r12/local_nCycles;

    double total_energy = 0.0;
    double total_energySquared = 0.0;
    double total_r12 = 0.0;
    int totalnAcc, totalnRej;
    MPI_Reduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&energySquared, &total_energySquared, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&r12, &total_r12, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&r12, &total_r12, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Allreduce(&nAccepted, &totalnAcc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nRejected, &totalnRej, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (my_rank == 0) {
//        cout << "Energy            = " << total_energy/numprocs << endl;
//        cout << "Energy squared    = " << total_energySquared/numprocs << endl;
//        cout << "r12 average       = " << total_r12/numprocs << endl;
//        cout << "nAccepted   = " << nAccepted << endl;
//        cout << "nRejected   = " << nRejected << endl;
//        cout << "Total tests = " << nAccepted + nRejected << endl;
//        cout << "nAcc = " << totalnAcc << endl;
//        cout << "nRej = " << totalnRej << endl;
//        cout << endl;
    }

    return energy;
}

void VMCSolver::setParameters(
        const double &newStepLength,
        const double &newAlpha,
        const double &newBeta
)
{
    stepLength = newStepLength;
    alpha = newAlpha*charge;
    beta = newBeta*charge;
}

VMCSolver::~VMCSolver()
{
    MPI_Finalize();
}

double VMCSolver::localEnergy(const mat &r)
{
    // kinetic energy
    mat rPlus(r), rMinus(r);
    double waveFunctionMinus, waveFunctionPlus;
    double waveFunctionCurrent = waveFunction(r);
    double kineticEnergy = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus);
            waveFunctionPlus = waveFunction(rPlus);
            //kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            kineticEnergy -= waveFunctionMinus + waveFunctionPlus;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5*h2*(kineticEnergy/waveFunctionCurrent
                            + 2.0*double(nParticles*nDimensions));
    //kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // potential energy
    double potentialEnergy = 0;
    for (int i = 0; i < nParticles; i++) {
        // r(i) = norm(r.row(i), 2);
        potentialEnergy += 1.0/norm(r.row(i), 2);
    }
    potentialEnergy *= -charge;

    // contribution from electron-electron potential
    double r12 = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for (int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            potentialEnergy += 1/sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

double VMCSolver::exactLocalEnergy(const mat &r)
{
    double EL1, EL2, dr;
    double rInverseSum = 0, rSum = 0, rijSum = 0;
    vec drvec(nDimensions);

    for (int i = 0; i < nParticles; i++) {
        drvec = r.row(i);
        dr = 0.0;
        for (int k = 0; k < nDimensions; k++)
            dr += drvec(k)*drvec(k);
        dr = sqrt(dr);
        rSum += dr;
        rInverseSum += 1.0/dr;
        for (int j = i + 1; j < nParticles; j++) {
            drvec -= r.row(j);
            dr = 0.0;
            for (int k = 0; k < nDimensions; k++)
                dr += drvec(k)*drvec(k);
            dr = sqrt(dr);
            rijSum += dr;
        }
    }

    vec r1vec(r.row(0).t());
    vec r2vec(r.row(1).t());

    double r1 = norm(r1vec, 2);
    double r2 = norm(r2vec, 2);

    double betafactor = 1.0/(1.0 + beta*rijSum);
    double betafactor2 = betafactor*betafactor;
    double rfactor = 1.0 - dot(r1vec, r2vec)/r1/r2;

    ////
//    cout << "betafactor = " << betafactor << endl;
//    cout << "rfactor = " << rfactor << endl;
//    cout << "1 = " << alpha*rSum*rfactor/rijSum << endl;
//    cout << "2 = " << -1.0/(2.0*betafactor2) << endl;
//    cout << "3 = " << -2.0/rijSum << endl;
//    cout << "4 = " << 2.0*beta/betafactor << endl;
    ////

    EL1 = (alpha - charge)*rInverseSum + 1.0/rijSum - alpha*alpha;
    EL2 = EL1 + (0.5*betafactor2)*
            (
                alpha*rSum*rfactor/rijSum - 0.5*betafactor2 - 2.0/rijSum
                + 2.0*beta*betafactor
            );

    ////
//    cout << "braces = " << test << endl;
//    cout << "EL1 = " << EL1 << endl;
//    cout << "EL2 = " << EL2 << endl;
//    cout << endl;
    ////
    return EL2;
}

double VMCSolver::dr(const mat &r, const int ii, const int jj)
{
//    double dr = 0.0;
//    vec3 drvec;
//    drvec = r.row(ii) - r.row(jj);
//    norm()
//    for (int i = 0; i < nDimensions; i++)
//    {
//        dr += drvec(i)*drvec(i);
//    }
//    dr = sqrt(dr);
    return norm((r.row(ii) - r.row(jj)), 2);
}

mat VMCSolver::quantumForce(const mat &r, const double &wf) {
    double wfMinus, wfPlus;
    mat rPlus(nParticles, nDimensions);
    mat rMinus(nParticles, nDimensions);
    mat qforce(nParticles, nDimensions);

    // compute the first derivative
    rPlus = rMinus = r;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = waveFunction(rPlus);
            wfPlus = waveFunction(rMinus);
            qforce(i,j) = (wfPlus - wfMinus)/(wf*h);
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }
    return qforce;
}

double VMCSolver::gaussianDeviate(long *seed)
{
    return ran2(seed);
}

double VMCSolver::waveFunction2(const mat &r)
{
    double argument = 0;
    for (int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }

    return exp(-argument * alpha);
}

double VMCSolver::waveFunction(const mat &R)
{
    double r = 0.0;
    double rSum = 0.0;

    // r1 + r2 + ...
    for (int i = 0; i < nParticles; i++) {
        r = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            r += R(i,j)*R(i,j);
        }
        rSum += sqrt(r);
    }

    double r12Sum = 0.0;
    // r1*r2*...
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            r = 0.0;
            for (int k = 0; k < nDimensions; k++)
            {
                r += (R(i,k) - R(j,k))*(R(i,k) - R(j,k));
            }
            r12Sum += sqrt(r);
        }
    }

    return exp(-alpha*rSum + r12Sum/(2.0*(1.0 + beta*r12Sum)));
}
