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
    h(1e-3),
    h2(1e6),
    idum(-1),
    nAccepted(0),
    nRejected(0),
    my_rank(my_rank),
    numprocs(numprocs)
{
    oldS.r = mat(nParticles, nDimensions);
    oldS.rij = mat(nParticles, nParticles);
    oldS.qForce = mat(nParticles, nDimensions);
    oldS.waveFunction = 0.0;

    newS.r = mat(nParticles, nDimensions);
    newS.rij = mat(nParticles, nParticles);
    newS.qForce = mat(nParticles, nDimensions);
    newS.waveFunction = 0.0;

    rOld = mat(nParticles, nDimensions);
    rNew = mat(nParticles, nDimensions);
    rijOld = mat(nParticles, nParticles);
    rijNew = mat(nParticles, nParticles);
}

VMCSolver::~VMCSolver()
{
}

void VMCSolver::setParameters(
        const double &alpha_,
        const double &beta_,
        const double &stepLength_)
{
    alpha = alpha_;
    beta = beta_;
    stepLength = stepLength_;
}

void VMCSolver::setParameters(
        const double &alpha_,
        const double &beta_,
        const double &stepLength_,
        const double &h_,
        const double &h2_,
        const bool &importanceSampling_,
        const bool &closedForm_)
{
    alpha = alpha_;
    beta = beta_;
    stepLength = stepLength_;
    h = h_;
    h2 = h2_;
    importanceSampling = importanceSampling_;
    closedForm = closedForm_;
}

double VMCSolver::hydrogenWF(const int &i, const vec3 &rvec) const
{
    double wavefunction;

    switch (i)
    {
    case 0:
        wavefunction = phi1s(rvec);
//        cout << "i = " << i << ", using phi1s" << endl;
        break;
    case 1 :
        wavefunction = phi2s(rvec);
//        cout << "i = " << i << ", using phi2s" << endl;
        break;
//    case 2 :
//        wavefunction = phi2s(r);
//        break;
//    case 3 :
//        wavefunction = phi2s(r);
//        break;
    default :
        // Process for all other cases.
        cout << "!!! We don't have this wavefunction yet!" << endl;
        exit(1);
    }

    return wavefunction;
}

double VMCSolver::phi1s(const vec3 &rvec) const
{
    double r = 0;
    for (int i = 0; i < nDimensions; i++)
    {
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

//    cout << "phi1s = " << exp(-alpha*r) << endl;
    return exp(-alpha*r);
}

double VMCSolver::phi2s(const vec3 &rvec) const
{
    double r = 0;
    for (int i = 0; i < nDimensions; i++)
    {
        r += rvec(i)*rvec(i);
    }
    r = sqrt(r);

    double arg = -alpha*r*0.5;

//    cout << "phi2s = " << (1.0 + arg)*exp(arg) << endl;
    return (1.0 + arg)*exp(arg);
}

double VMCSolver::runMonteCarloIntegration(const int &nCycles)
{
//    rOld.zeros(nParticles, nDimensions);
//    rNew.zeros(nParticles, nDimensions);
//    waveFunctionOld = 0;
//    waveFunctionNew = 0;

    oldS.r.zeros(nParticles, nDimensions);
    newS.r.zeros(nParticles, nDimensions);
    oldS.waveFunction = 0;
    newS.waveFunction = 0;

    nAccepted = 0;
    nRejected = 0;

    double r12 = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE = 0;
//    double omegaRatio;
//    h = dt;
//    h2 = 1.0/(h*h);
    D = 0.5;
    Ddt = D*h;
//    mat qForceOld(nParticles, nDimensions);
//    mat qForceNew(nParticles, nDimensions);

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
            if (importanceSampling)
            {
//                rOld(i,j) = gaussianDeviate(&idum)*sqrt(2.0*Ddt);
                oldS.r(i,j) = gaussianDeviate(&idum)*sqrt(2.0*Ddt);
            }
            else
            {
//                rOld(i,j) = stepLength*ran2(&idum) - 0.5;
                oldS.r(i,j) = stepLength*ran2(&idum) - 0.5;
            }
        }
    }

//    rNew = rOld;
//    waveFunctionOld = wavefunction(rOld);

    newS.r = oldS.r;
    oldS.waveFunction = wavefunction(oldS.r);
//    calculate_rij(rOld, rijOld);
//    fijOld = calculate_fij_element()
    if (importanceSampling)
    {
//        qForceOld = quantumForce(rOld, waveFunctionOld);
//        qForceNew = qForceOld;

        oldS.qForce = quantumForce(oldS.r, oldS.waveFunction);
        newS.qForce = oldS.qForce;
    }

    // loop over Monte Carlo cycles
    for (int cycle = my_start; cycle < my_end; cycle++) {

        // Store the current value of the wave function
//        waveFunctionOld = wavefunction(rOld);
        oldS.waveFunction = wavefunction(oldS.r);

        for (int i = 0; i < nParticles; i++) {

            if (importanceSampling)
                runcycle_importanceSampling(i);
            else
                runcycle(i);

            // update energies
            if (closedForm)
            {
//                deltaE = localEnergyClosedForm(rNew);
                deltaE = localEnergyClosedForm(newS.r);
            }
            else
            {
//                deltaE = localEnergy(rNew);
                deltaE = localEnergy(newS.r);
            }
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        // calculating the optimal value of the distance between the particles
//        r12 += norm(rNew.row(0), 2);
        r12 += norm(newS.r.row(0), 2);
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

void VMCSolver::runcycle(const int &i)
{
    // New position to test
    for (int j = 0; j < nDimensions; j++) {
//        rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
        newS.r(i,j) = oldS.r(i,j) + stepLength*(ran2(&idum) - 0.5);
    }
//    update_rij(rNew, rijNew, i);
//    update_rij(newS.r, newS.rij, i);
    update_rij(newS, i);

    // Recalculate the value of the wave function
//    waveFunctionNew = wavefunction(rNew);
    newS.waveFunction = wavefunction(newS);

    // Check for step acceptance (if yes, update position, if no, reset position)
//    if (ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
    if (ran2(&idum) <= slaterRatio()*jastrowRatio(i)) {
//        rOld.row(i) = rNew.row(i); // update position SHOULD BE JUST ROW!!!
//        waveFunctionOld = waveFunctionNew;
//        rijOld = rijNew;

        oldS.r.row(i) = newS.r.row(i);
        oldS.waveFunction = newS.waveFunction;
        oldS.rij = newS.rij;

        nAccepted++;
    } else {
//        rNew.row(i) = rOld.row(i); // reset position, throw away the test-position
//        rijNew = rijOld;

        newS.r.row(i) = oldS.r.row(i);
        newS.rij = oldS.rij;

        nRejected++;
    }
}

void VMCSolver::runcycle_importanceSampling(const int &i)
{
    double omegaRatio = 0.0;

    // New position to test
    for (int j = 0; j < nDimensions; j++) {
//        rNew(i,j) = rOld(i,j) + gaussianDeviate(&idum)*sqrt(2.0*Ddt)
//                + qForceOld(i,j)*Ddt;
        newS.r(i,j) = oldS.r(i,j) + gaussianDeviate(&idum)*sqrt(2.0*Ddt)
                + oldS.qForce(i,j)*Ddt;
    }

    // Recalculate the value of the wave function
//    waveFunctionNew = wavefunction(rNew);
//    qForceNew = quantumForce(rNew, waveFunctionNew);

    newS.waveFunction = wavefunction(newS.r);
    newS.qForce = quantumForce(newS.r, newS.waveFunction);

    // log of the ratio of the greens function
    // from slides
    omegaRatio = 0.0;
    for (int j = 0; j < nDimensions; j++) {
//        omegaRatio += (qForceOld(i,j) + qForceNew(i,j))
//                * (2.0*(rOld(i,j) - rNew(i,j)) + Ddt*(qForceOld(i,j) - qForceNew(i,j)));
        omegaRatio += (oldS.qForce(i,j) + newS.qForce(i,j))
                * (2.0*(oldS.r(i,j) - newS.r(i,j)) + Ddt*(oldS.qForce(i,j) - newS.qForce(i,j)));
    }
    omegaRatio /= 4.0;
    omegaRatio = exp(omegaRatio);

    // Check for step acceptance (if yes, update position, if no, reset position)
    if (ran2(&idum) <= omegaRatio*slaterRatio()*jastrowRatio(i)) {
//        rOld.row(i) = rNew.row(i); // update position
//        qForceOld.row(i) = qForceNew.row(i); // SHOULD BE JUST THE ROW !!!
//        waveFunctionOld = waveFunctionNew;

        oldS.r.row(i) = newS.r.row(i); // update position
        oldS.qForce.row(i) = newS.qForce.row(i); // SHOULD BE JUST THE ROW !!!
        oldS.waveFunction = newS.waveFunction;

        nAccepted++;
    } else {
//        rNew.row(i) = rOld.row(i);
//        qForceNew.row(i) = qForceOld.row(i); // SHOULD BE JUST THE ROW !!!

        newS.r.row(i) = oldS.r.row(i);
        newS.qForce.row(i) = oldS.qForce.row(i); // SHOULD BE JUST THE ROW !!!

        nRejected++;
    }
}

const mat VMCSolver::quantumForce(const mat &r, const double &wf) const
{
    double wfMinus, wfPlus;
    mat rPlus(r);
    mat rMinus(r);
    mat qforce_return(nParticles, nDimensions);

    // computing the first derivative
    rPlus = rMinus = r;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            qforce_return(i,j) = wfPlus - wfMinus;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }
    qforce_return = qforce_return/(wf*h);

    return qforce_return;
}

void VMCSolver::calculate_rij(const mat &r, mat &rij)
{
    rij.zeros();
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            for (int k = 0; k < nDimensions; k++)
            {
                rij(i,j) += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            rij(i,j) = sqrt(rij(i,j));
        }
    }
}

void VMCSolver::update_rij(const mat &r, mat &rij, const int j)
{
    for (int i = 0; i < j; i++)
    {
        rij(i,j) = 0.0;
        for (int k = 0; k < nDimensions; k++)
        {
            rij(i,j) += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
        }
        rij(i,j) = sqrt(rij(i,j));
    }
    for (int i = j+1; i < nParticles; i++)
    {
        rij(i,j) = 0.0;
        for (int k = 0; k < nDimensions; k++)
        {
            rij(j,i) += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
        }
        rij(j,i) = sqrt(rij(j,i));
    }
}

void VMCSolver::update_rij(VMCSolver::data &s, const int &j) const
{
    for (int i = 0; i < j; i++)
    {
        s.rij(i,j) = 0.0;
        for (int k = 0; k < nDimensions; k++)
        {
            s.rij(i,j) += (s.r(i,k) - s.r(j,k))*(s.r(i,k) - s.r(j,k));
        }
        s.rij(i,j) = sqrt(s.rij(i,j));
    }
    for (int i = j+1; i < nParticles; i++)
    {
        s.rij(i,j) = 0.0;
        for (int k = 0; k < nDimensions; k++)
        {
            s.rij(j,i) += (s.r(i,k) - s.r(j,k))*(s.r(i,k) - s.r(j,k));
        }
        s.rij(j,i) = sqrt(s.rij(j,i));
    }
}

void VMCSolver::calculate_fij(const mat &rij, mat &fij)
{
//    double a = ((i+j)%2 == 0) ? 0.25 : 0.5;
//    double r12 = rij(i,j);

    double a;
    double r12;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            r12 = rij(i,j);
            a = ((i+j)%2 == 0) ? 0.25 : 0.5;
            fij(i,j) = a*r12/(1.0 + beta*r12);
        }
    }

//    return a*r12/(1.0 + beta*r12);
}

double VMCSolver::calculate_fij_element(const mat &rij, const int &i, const int &j) const
{
    double a = ((i+j)%2 == 0) ? 0.25 : 0.5;
    double r12 = rij(i,j);

    return a*r12/(1.0 + beta*r12);
}

void VMCSolver::update_slater()
{
//    for (int i = 0; i < nParticles; i++)
//    {
//        for (int j = 0; j < nParticles; j++)
//        {
//            S(j) = Dnew
//        }
//    }
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
