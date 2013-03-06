#include "CVMCSolver.h"

VMCSolver::VMCSolver(
        int my_rank,
        int numprocs,
        const int &charge,
        const int &nParticles):
    nDimensions(3),
    charge(charge),
    stepLength(1.0),
    nParticles(nParticles),
    h(0.001),
    h2(1000000),
    idum(-1),
    nAccepted(0),
    nRejected(0),
    my_rank(my_rank),
    numprocs(numprocs)
{
    rOld = mat(nParticles, nDimensions);
    rNew = mat(nParticles, nDimensions);
}

VMCSolver::~VMCSolver()
{
}

//double VMCSolver::runMonteCarloIntegration(
//        const int &nCycles_,
//        const double &stepLength_,
//        const double &alpha_,
//        const double &beta_,
//        const bool closedform)
double VMCSolver::runMonteCarloIntegration(const int &nCycles, const bool &closedForm)
{
    rOld = zeros(nParticles, nDimensions);
    rNew = zeros(nParticles, nDimensions);

    nAccepted = 0;
    nRejected = 0;

    double r12 = 0;
    double waveFunctionOld = 0, waveFunctionNew = 0;
    double energySum = 0, energySquaredSum = 0;
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

        // Store the current value of the wave function
        waveFunctionOld = wavefunction(rOld);

        for (int i = 0; i < nParticles; i++) {
            // New position to test
            for (int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = wavefunction(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                rOld.row(i) = rNew.row(i); // update position
                waveFunctionOld = waveFunctionNew;
                nAccepted++;
            } else {
                rNew.row(i) = rOld.row(i); // reset position, throw away the test-position
                nRejected++;
            }

            // update energies
            if (closedForm)
                deltaE = localEnergyClosedForm(rNew);
            else
                deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        // calculating the optimal value of the distance between the particles
        r12 += norm(rNew.row(0), 2);
    }

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
//        cout << "nReject    MPI_Finalize();ed   = " << nRejected << endl;
//        cout << "Total tests = " << nAccepted + nRejected << endl;
//        cout << "nAcc = " << totalnAcc << endl;
//        cout << "nRej = " << totalnRej << endl;
//        cout << endl;
    }

    return total_energy/numprocs;
}

//double VMCSolver::runMonteCarloIntegrationImportanceSampling(
//        const int &nCycles_,
//        const double &stepLength_,
//        const double &alpha_,
//        const double &beta_,
//        const bool closedform,
//        const double dt)
double VMCSolver::runMonteCarloIntegrationImportanceSampling(
        const int &nCycles,
        const double &dt,
        const bool &closedform)
{
    rOld.zeros(nParticles, nDimensions);
    rNew.zeros(nParticles, nDimensions);

    nAccepted = 0;
    nRejected = 0;

    double r12 = 0;
    double waveFunctionOld = 0, waveFunctionNew = 0;
    double energySum = 0, energySquaredSum = 0;
    double deltaE = 0;
    double omegaRatio;
    h = dt;
    h2 = 1.0/(h*h);
    double D = 0.5;
    double Ddt = D*h;
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
            rOld(i,j) = gaussianDeviate(&idum)*sqrt(2.0*Ddt);
//            rOld(i,j) = 2.0;
        }
    }

    rNew = rOld;
    waveFunctionOld = wavefunction(rOld);
    qforceOld = quantumForce(rOld, waveFunctionOld);

    // loop over Monte Carlo cycles
    for (int cycle = my_start; cycle < my_end; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = wavefunction(rOld);

        for (int i = 0; i < nParticles; i++) {
            // New position to test
            for (int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(2.0*Ddt)
                        + qforceOld(i,j)*Ddt;
            }

            // Recalculate the value of the wave function
            waveFunctionNew = wavefunction(rNew);
            qforceNew = quantumForce(rNew, waveFunctionNew);

            // log of the ratio of the greens function
            // from slides
            omegaRatio = 0.0;
            for (int j = 0; j < nDimensions; j++) {
                omegaRatio += (qforceOld(i,j) + qforceNew(i,j))
                        * (2.0*(rOld(i,j) - rNew(i,j)) + Ddt*(qforceOld(i,j) - qforceNew(i,j)));
            }
            omegaRatio /= 4.0;
            omegaRatio = exp(omegaRatio);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= omegaRatio*waveFunctionNew*waveFunctionNew/(waveFunctionOld*waveFunctionOld)) {
                rOld.row(i) = rNew.row(i); // update position
                qforceOld.row(i) = qforceNew.row(i);
                waveFunctionOld = waveFunctionNew;
                nAccepted++;
            } else {
                rNew.row(i) = rOld.row(i); // reset position, throw away the test-position
                qforceNew.row(i) = qforceOld.row(i);
                nRejected++;
            }

            // update energies
            if (closedform)
                deltaE = localEnergyClosedForm(rNew);
            else
                deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        // calculating the optimal value of the distance between the particles
        r12 += norm(rNew.row(0), 2);
    }

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

    return total_energy/numprocs;
}

mat VMCSolver::quantumForce(const mat &r, const double &wf) {
    double wfMinus, wfPlus;
    mat rPlus(nParticles, nDimensions); // mat rPlus(r)
    mat rMinus(nParticles, nDimensions);
    mat qforce(nParticles, nDimensions);

    // computing the first derivative
    rPlus = rMinus = r;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            qforce(i,j) = wfPlus - wfMinus;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }
    qforce = qforce/(wf*h);
    return qforce;
}

double VMCSolver::gaussianDeviate(long *seed)
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

double VMCSolver::localEnergy(const mat &r)
{
    // kinetic energy
    mat rPlus(r);
    mat rMinus(r);
    double waveFunctionMinus, waveFunctionPlus;
    double waveFunctionCurrent = wavefunction(r);
    double kineticEnergy = 0;
    // computing the second derivative
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = wavefunction(rMinus);
            waveFunctionPlus = wavefunction(rPlus);
            //kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            kineticEnergy -= waveFunctionMinus + waveFunctionPlus;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }
    //kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;
    kineticEnergy = 0.5*h2*(kineticEnergy/waveFunctionCurrent
                            + 2.0*double(nParticles*nDimensions));

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

// In CWavefunction !!!
//double VMCSolver::localEnergyClosedForm(const mat &r)
//{
//    double EL1, EL2, dr;
//    double rInverseSum = 0, rSum = 0, rijSum = 0;
//    vec drvec(nDimensions);

//    for (int i = 0; i < nParticles; i++) {
//        drvec = r.row(i);
//        dr = 0.0;
//        for (int k = 0; k < nDimensions; k++)
//            dr += drvec(k)*drvec(k);
//        dr = sqrt(dr);
//        rSum += dr;
//        rInverseSum += 1.0/dr;
//        for (int j = i + 1; j < nParticles; j++) {
//            drvec -= r.row(j);
//            dr = 0.0;
//            for (int k = 0; k < nDimensions; k++)
//                dr += drvec(k)*drvec(k);
//            dr = sqrt(dr);
//            rijSum += dr;
//        }
//    }

//    vec r1vec(r.row(0).t());
//    vec r2vec(r.row(1).t());

//    double r1 = norm(r1vec, 2);
//    double r2 = norm(r2vec, 2);

//    double betafactor = 1.0/(1.0 + beta*rijSum);
//    double betafactor2 = betafactor*betafactor;
//    double rfactor = 1.0 - dot(r1vec, r2vec)/r1/r2;

//    ////
////    cout << "betafactor = " << betafactor << endl;
////    cout << "rfactor = " << rfactor << endl;
////    cout << "1 = " << alpha*rSum*rfactor/rijSum << endl;
////    cout << "2 = " << -1.0/(2.0*betafactor2) << endl;
////    cout << "3 = " << -2.0/rijSum << endl;
////    cout << "4 = " << 2.0*beta/betafactor << endl;
//    ////

//    EL1 = (alpha - charge)*rInverseSum + 1.0/rijSum - alpha*alpha;
//    EL2 = EL1 + (0.5*betafactor2)*
//            (
//                alpha*rSum*rfactor/rijSum - 0.5*betafactor2 - 2.0/rijSum
//                + 2.0*beta*betafactor
//            );

//    ////
////    cout << "braces = " << test << endl;
////    cout << "EL1 = " << EL1 << endl;
////    cout << "EL2 = " << EL2 << endl;
////    cout << endl;
//    ////
//    return EL2;
//}

//double VMCSolver::dr(const mat &r, const int ii, const int jj)
//{
////    double dr = 0.0;
////    vec3 drvec;
////    drvec = r.row(ii) - r.row(jj);
////    norm()
////    for (int i = 0; i < nDimensions; i++)
////    {
////        dr += drvec(i)*drvec(i);
////    }
////    dr = sqrt(dr);
//    return norm((r.row(ii) - r.row(jj)), 2);
//}

//double VMCSolver::waveFunction2(const mat &r)
//{
//    double argument = 0;
//    for (int i = 0; i < nParticles; i++) {
//        double rSingleParticle = 0;
//        for (int j = 0; j < nDimensions; j++) {
//            rSingleParticle += r(i,j) * r(i,j);
//        }
//        argument += sqrt(rSingleParticle);
//    }

//    return exp(-argument * alpha);
//}


// in CWavefunction !!!
//double VMCSolver::wavefunction(const mat &R)
//{
//    double r = 0.0;
//    double rSum = 0.0;

//    // r1 + r2 + ...
//    for (int i = 0; i < nParticles; i++) {
//        r = 0.0;
//        for (int j = 0; j < nDimensions; j++)
//        {
//            r += R(i,j)*R(i,j);
//        }
//        rSum += sqrt(r);
//    }

//    double r12Sum = 0.0;
//    // r1*r2*...
//    for (int i = 0; i < nParticles; i++) {
//        for (int j = i + 1; j < nParticles; j++) {
//            r = 0.0;
//            for (int k = 0; k < nDimensions; k++)
//            {
//                r += (R(i,k) - R(j,k))*(R(i,k) - R(j,k));
//            }
//            r12Sum += sqrt(r);
//        }
//    }

//    return exp(-alpha*rSum + r12Sum/(2.0*(1.0 + beta*r12Sum)));
//}
