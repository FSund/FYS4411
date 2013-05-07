#include <src/Slater.h>

Slater::Slater(const int &nParticles):
    nParticles(nParticles),
    nDimensions(3),
    rOld(mat(nParticles, nDimensions)),
    rNew(mat(nParticles, nDimensions)),
    N(nParticles/2),
    slaterUPold(mat(N, N)),
    slaterDOWNold(mat(N, N)),
    slaterUPnew(mat(N, N)),
    slaterDOWNnew(mat(N, N)),
    slaterUPinvOld(mat(N, N)),
    slaterDOWNinvOld(mat(N, N)),
    slaterUPinvNew(mat(N, N)),
    slaterDOWNinvNew(mat(N, N)),
    SUP(vec(N)),
    SDOWN(vec(N)),
    dwavefunction(rowvec(nDimensions)),
    grad(rowvec(nDimensions))
{
    orbitals = new Orbitals;
}

Slater::~Slater()
{
    delete orbitals;
}

void Slater::initialize(const mat &r)
{
    rNew = rOld = r;

    for (int i = 0; i < N; i++) // loop over particles
    {
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUPold(i,j)   = orbitals->wavefunction(rOld.row(i),   j);
            slaterDOWNold(i,j) = orbitals->wavefunction(rOld.row(i+N), j);
        }
    }
    slaterUPinvOld   = inv(slaterUPold);
    slaterDOWNinvOld = inv(slaterDOWNold);

    slaterUPnew      = slaterUPold;
    slaterDOWNnew    = slaterDOWNold;
    slaterUPinvNew   = slaterUPinvOld;
    slaterDOWNinvNew = slaterDOWNinvOld;

    ratioUP = ratioDOWN = 1.0;
}

void Slater::setAlpha(const double &newAlpha)
{
    orbitals->setAlpha(newAlpha);
}

void Slater::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;
    if (k < N) UP = true;
    else UP = false;

    // updating both the determinant and the inverse, since importance sampling
    // accepts ~95% of the moves anyway, and we need the inverse in our
    // new closed form gradient
    updateSlater();
    updateInverse();
}

double Slater::getRatio()
{
    if (UP)
    {
        int i = currentParticle;
        ratioUP = 0.0;
        for (int j = 0; j < N; j++) // loop over orbitals
            ratioUP += slaterUPnew(i,j)*slaterUPinvOld(j,i);

        return ratioUP*ratioUP;
    }
    else
    {
        int i = currentParticle - N;
        ratioDOWN = 0.0;
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            ratioDOWN += slaterDOWNnew(i,j)*slaterDOWNinvOld(j,i);
        }

        return ratioDOWN*ratioDOWN;
    }
}

void Slater::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    if (UP)
    {
        slaterUPold.row(currentParticle) = slaterUPnew.row(currentParticle);
        slaterUPinvOld = slaterUPinvNew; // can't use just column/row!
    }
    else
    {
        slaterDOWNold.row(currentParticle-N) = slaterDOWNnew.row(currentParticle-N);
        slaterDOWNinvOld = slaterDOWNinvNew; // can't use just column/row!
    }
}

void Slater::rejectMove()
{
//    rNew.row(currentParticle) = rOld.row(currentParticle);
    rNew = rOld;

    // now updating the inverse as well, since we update it for every new set
    // of positions with the new closed form gradient
    if (UP)
    {
        slaterUPnew.row(currentParticle) = slaterUPold.row(currentParticle);
        slaterUPinvNew = slaterUPinvOld;
    }
    else
    {
        slaterDOWNnew.row(currentParticle-N) = slaterDOWNold.row(currentParticle-N);
        slaterDOWNinvNew = slaterDOWNinvOld;
    }
}

double Slater::wavefunction()
{
    return det(slaterUPnew)*det(slaterDOWNnew);
}

double Slater::wavefunction(const mat &r)
{
    /* Only for use in numerical derivatives */

    mat slaterUP(N, N);
    mat slaterDOWN(N, N);

    for (int i = 0; i < N; i++) // loop over particles
    {
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUP(i,j)   = orbitals->wavefunction(r.row(i), j);
            slaterDOWN(i,j) = orbitals->wavefunction(r.row(i+N), j);
        }
    }

    return det(slaterUP)*det(slaterDOWN);
}

rowvec Slater::localGradient(const int &i)
{
    /* Calculates the gradient for particle k */

    grad = zeros<vec>(nDimensions);
    if (i < N)
        for (int j = 0; j < N; j++)
//            grad += orbitals->gradient(r.row(i),j)*slaterUPinvOld(j,i)/ratioUP;
            grad += orbitals->gradient(rNew.row(i),j)*slaterUPinvNew(j,i);
    else
        for (int j = 0; j < N; j++)
//            grad += orbitals->gradient(r.row(i),j)*slaterDOWNinvOld(j,i-N)/ratioDOWN;
            grad += orbitals->gradient(rNew.row(i),j)*slaterDOWNinvNew(j,i-N);

    return grad;
}

rowvec Slater::localGradientNumerical(const mat &r, const int &k, const double &h)
{
    /* Calculates the the gradient for particle k */

    dwavefunction = zeros<vec>(nDimensions);
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nDimensions; i++)
    {
        rPlus(k,i) += h;
        rMinus(k,i) -= h;
        wfMinus = wavefunction(rMinus);
        wfPlus = wavefunction(rPlus);
        dwavefunction(i) = wfPlus - wfMinus;
        rPlus(k,i) = rMinus(k,i) = r(k,i);
    }
    dwavefunction /= (2.0*wfCurrent*h);
//    dwavefunction /= (2.0*h); // not local

    return dwavefunction;
}

double Slater::localLaplacian(const int &i)
{
    /* Calculates the laplacian for particle i */

    lapl = 0.0;
    if (i < N)
    {
        for (int j = 0; j < N; j++)
        {
            lapl += orbitals->laplacian(rNew.row(i),j)*slaterUPinvNew(j,i);
        }
    }
    else
    {
        for (int j = 0; j < N; j++)
        {
            lapl += orbitals->laplacian(rNew.row(i),j)*slaterDOWNinvNew(j,i-N);
        }
    }

    return lapl;
}

double Slater::localLaplacianNumerical(const mat &r, const int &k, const double &h)
{
    /* Calculates the laplacian for particle k */
    ddwavefunction = 0.0;
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nDimensions; i++) {
        rPlus(k,i) += h;
        rMinus(k,i) -= h;
        wfMinus = wavefunction(rMinus);
        wfPlus = wavefunction(rPlus);
        ddwavefunction += wfMinus + wfPlus - 2.0*wfCurrent;
        rPlus(k,i) = rMinus(k,i) = r(k,i);
    }
    ddwavefunction /= (wfCurrent*h*h);

    return ddwavefunction;
}

void Slater::updateSlater()
{
    int i = currentParticle;
    if (UP)
        for (int j = 0; j < N; j++) // loop over orbitals
            slaterUPnew(i,j) = orbitals->wavefunction(rNew.row(i), j);
    else
        for (int j = 0; j < N; j++) // loop over orbitals
            slaterDOWNnew(i-N,j) = orbitals->wavefunction(rNew.row(i), j);
}

void Slater::updateInverse()
{
    ratio = getRatio();

    if (UP) // only updating the UP matrix
    {
        int i = currentParticle;

        SUP.zeros();
        // updating all but the ith column of the inverse slater matrix
        for (int j = 0; j < N; j++) // loop over columns j != i
        {
            if (j == i) continue; // not for column i
            for (int l = 0; l < N; l++) // loop over elements in column j
            {
                SUP(j) += slaterUPnew(i,l)*slaterUPinvOld(l,j);
            }
        }
        for (int j = 0; j < N; j++) // loop over columns in the inverse slater matrix
        {
            if (j == i) continue; // not updating column i
            for (int k = 0; k < N; k++) // loop over elements in column j
            {
                // maybe we should recalculate the ratio first to be sure??
                slaterUPinvNew(k,j) = slaterUPinvOld(k,j) - SUP(j)*slaterUPinvOld(k,i)/ratioUP;
            }
        }
        // updating column i of the inverse slater matrix
        for (int k = 0; k < N; k++) // loop over elements in i'th column
        {
            slaterUPinvNew(k,i) = slaterUPinvOld(k,i)/ratioUP;
        }
    }
    else // only updating the DOWN matrix
    {
        int i = currentParticle - N;

        SDOWN.zeros();
        // updating all but the ith column of the inverse slater matrix
        for (int j = 0; j < N; j++) // loop over columns j != i
        {
            if (j == i) continue; // not for column i
            for (int l = 0; l < N; l++) // loop over elements in column j
            {
                SDOWN(j) += slaterDOWNnew(i,l)*slaterDOWNinvOld(l,j);
            }
        }
        for (int j = 0; j < N; j++) // loop over columns in the inverse slater matrix
        {
            if (j == i) continue; // not updating column i
            for (int k = 0; k < N; k++) // loop over elements in column j
            {
                // maybe we should recalculate the ratio first to be sure??
                slaterDOWNinvNew(k,j) = slaterDOWNinvOld(k,j) - SDOWN(j)*slaterDOWNinvOld(k,i)/ratioDOWN;
            }
        }
        // updating column i of the inverse slater matrix
        for (int k = 0; k < N; k++) // loop over elements in i'th column
        {
            slaterDOWNinvNew(k,i) = slaterDOWNinvOld(k,i)/ratioDOWN;
        }
    }
}

// debug stuff
mat Slater::gradient(const mat &r, const double &h)
{
    rPlus = rMinus = r;
    mat temp(nParticles, nDimensions);
    wfCurrent = wavefunction(r);
    dfactor = 1.0/(wfCurrent*2.0*h);
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            wfMinus = wavefunction(rMinus);
            wfPlus = wavefunction(rPlus);
            temp(i,j) = (wfPlus - wfMinus)*dfactor;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }

    return temp;
}
