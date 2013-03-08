#include "mainapplication.h"

MainApplication::MainApplication(int &my_rank_, int &numprocs_):
    my_rank(my_rank_),
    numprocs(numprocs_)
{
}

void MainApplication::runApplication()
{
    // cout << "MainApplication::runApplication()" << endl;
//    variational_paramenters();
//    optimal_steplength();
//    steplength_secant();
//    closedformBenchmark();
//    importanceSampling();

    beryllium_variational_parameters();
}

void MainApplication::variational_paramenters()
{
    mat energySurface;
    int nCycles, varsteps;
    double steplength, alpha, beta, alphamax, betamax, energy, alphamin, betamin;

    FILE* file;
    file = fopen("./energySurface.dat","w");
    FILE* xyfile;
    xyfile = fopen("./alpha_beta.dat","w");

    nCycles = 1e6;
    bool closedform = 0;

    alphamin = 1.8; // 1.8
    alphamax = 1.8;
    betamin = 0.36; // 0.36
    betamax = 0.36;

    varsteps = 1;
    steplength = 1.0;

    energySurface = zeros(varsteps, varsteps);

    CHelium solver(numprocs, my_rank);

    for (int nAlpha = 0; nAlpha < varsteps; nAlpha++)
    {
        alpha = alphamin + nAlpha*(alphamax-alphamin)/varsteps;
        if (my_rank == 0) cout << "Alpha = " << alpha << endl;
        for (int nBeta = 0; nBeta < varsteps; nBeta++)
        {
            beta = betamin + nBeta*(betamax-betamin)/varsteps;
            if (my_rank == 0) cout << "Beta = " << beta << endl;

            solver.setParameters(alpha, beta, steplength);
            // energy = solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta, closedform);
            solver.closedForm = closedform;
            energy = solver.runMonteCarloIntegration(nCycles);

            energySurface(nAlpha, nBeta) = energy;
            fprintf(file, "%f ", energy);
            fflush(0); // flush files (write now instead of waiting for program to finish)
        }
        fprintf(file, "\n");
    }

    for (int i = 0; i < varsteps; i++)
    {
        alpha = alphamin + i*(alphamax-alphamin)/varsteps;
        fprintf(xyfile, "%f ", alpha);
    }
    fprintf(xyfile, "\n");
    for (int i = 0; i < varsteps; i++)
    {
        beta = betamin + i*(betamax-betamin)/varsteps;
        fprintf(xyfile, "%f ", beta);
    }
    fprintf(xyfile, "\n");
}

//void MainApplication::minimum()
//{
//    int nCycles;
//    double alpha, beta, steplength;

//    nCycles = 1000000;
//    alpha = 1.0;
//    beta = 1.0;
//    steplength = 1.0;

//    CHelium solver(numprocs, my_rank);
//    solver.setParameters(alpha, beta, steplength);
//    energy = solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta);
//}


void MainApplication::optimal_steplength()
{
    int nCycles, varsteps;
    double steplength, alpha, beta, steplengthMin, steplengthMax, dr;

    FILE* file;
    file = fopen("./drPlot.dat","w");

    nCycles = 1e6;
    bool closedform = 0;

    alpha = 1.8;
    beta = 0.36;

    varsteps = 100;
    steplengthMin = 1.4;
    steplengthMax = 1.6;
    dr = (steplengthMax - steplengthMin)/varsteps;

    CHelium solver(numprocs, my_rank);

    for (double nSteplength = 0; nSteplength < varsteps; nSteplength++)
    {
        steplength = steplengthMin + nSteplength*dr;
        if (my_rank == 0) cout << "Steplength        = " << steplength << endl;

        solver.setParameters(alpha, beta, steplength);
        // solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta, closedform);
        solver.closedForm = closedform;
        solver.runMonteCarloIntegration(nCycles);

        fprintf(file, "%f %d %d \n", steplength, solver.nAccepted, solver.nRejected);
        fflush(0); // flush files (write now instead of waiting for program to finish)
    }
}

void MainApplication::steplength_secant()
{
    int nCycles;
    double steplength, alpha, beta, acceptanceRate, wanted_acceptanceRate,
            steplength_p, steplength_pp, fpp, fp, tolerance;
    CHelium solver(numprocs, my_rank);

//    FILE* file;
//    file = fopen("./drPlot.dat","w");

    nCycles = 1e6;
    bool closedform = 0;

    alpha = 1.8;
    beta = 0.36;

//    maxSteps = 100;
    wanted_acceptanceRate = 0.50;
    tolerance = 1e-4;

    steplength_pp = 1;
    steplength_p = 1.1;

    solver.setParameters(alpha, beta, steplength_pp);
    //solver.runMonteCarloIntegration(nCycles, steplength_pp, alpha, beta, closedform);
    solver.closedForm = closedform;
    solver.runMonteCarloIntegration(nCycles);

    acceptanceRate = double(solver.nAccepted)/double(solver.nAccepted + solver.nRejected);
    fpp = acceptanceRate - wanted_acceptanceRate;

//    if (solver.my_rank == 0)
//    {
//        cout << "acceptanceRate   = " << acceptanceRate << endl;
//        cout << "steplength_p     = " << steplength_p << endl;
//        cout << "steplength_pp    = " << steplength_pp << endl;
//        cout << "nAccepted      = " << solver.nAccepted << endl;
//        cout << "nRejected      = " << solver.nRejected << endl;
//        cout << "Total tests    = " << solver.nAccepted + solver.nRejected << endl;
//        cout << endl;
//    }

    while (abs(fpp) > tolerance)
    {
        solver.setParameters(alpha, beta, steplength_p);
        //solver.runMonteCarloIntegration(nCycles, steplength_p, alpha, beta, closedform);
        solver.closedForm = closedform;
        solver.runMonteCarloIntegration(nCycles);

        acceptanceRate = double(solver.nAccepted)/double(solver.nAccepted + solver.nRejected);
        fp = acceptanceRate - wanted_acceptanceRate;

        steplength = steplength_p - fp*(steplength_p - steplength_pp)/(fp - fpp);

//        if (solver.my_rank == 0)
//        {
//            cout << "acceptancerate = " << acceptanceRate << endl;
//            cout << "fp             = " << fp << endl;
//            cout << "fpp            = " << fpp << endl;
//            cout << "steplength_p   = " << steplength_p << endl;
//            cout << "steplength_pp  = " << steplength_pp << endl;
//            cout << "convtest       = " << abs(fpp/wanted_acceptanceRate) << endl;
//        }

        steplength_pp = steplength_p;
        steplength_p = steplength;
        fpp = fp;
    }
    cout << "Found best steplength as " << steplength
         << " within " << tolerance
         << " of the acceptance rate " << wanted_acceptanceRate << "%"
         << endl;

//    fprintf(file, "%f %d %d \n", steplength, solver.nAccepted, solver.nRejected);
    //    fflush(0); // flush files (write now instead of waiting for program to finish)
}

void MainApplication::closedformBenchmark()
{
    int nCycles;
    double steplength, alpha, beta, energy;
    CHelium solver(numprocs, my_rank);

    bool closedform = 0;

    nCycles = 1e6;
    alpha = 1.8;
    beta = 0.36;
    steplength = 1.485;

    closedform = 1;
    solver.setParameters(alpha, beta, steplength);
    // energy = solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta, closedform);
    solver.closedForm = closedform;
    energy = solver.runMonteCarloIntegration(nCycles);

    if (solver.my_rank == 0) cout << "Closed form energy = " << energy << endl;

//    closedform = 0;
//    energy = solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta, closedform);
//    if (solver.my_rank == 0) cout << "Numerical energy = " << energy << endl;


}

void MainApplication::importanceSampling()
{
    int nCycles, varsteps;
    double steplength, alpha, beta, energy, timestep, acceptanceRate;
    CHelium solver(numprocs, my_rank);

    bool closedform = 1;

    nCycles = 1e6;
    alpha = 1.8;
    beta = 0.36;
    steplength = 1.485;

    FILE* file;
    file = fopen("./dtPlot.dat","w");

    timestep = 1;
    varsteps = 100;
    for (double nTimestep = 0; nTimestep <= varsteps; nTimestep++) {
        timestep /= 10;

        solver.setParameters(alpha, beta, steplength);
        solver.closedForm = closedform;
        // NEED TO SET h = timestep SOMEHOW!!
        energy = solver.runMonteCarloIntegration(nCycles);
        acceptanceRate = double(solver.nAccepted)/double(solver.nAccepted + solver.nRejected);

        if (solver.my_rank == 0) {
            cout << "Timestep        = " << timestep << endl;
            cout << "Energy          = " << energy << endl;
            cout << "Acceptance rate = " << acceptanceRate << endl;
        }

        fprintf(file, "%f %f %f \n", timestep, energy, acceptanceRate);
        fflush(0); // flush files (write now instead of waiting for program to finish)
    }
}

void MainApplication::beryllium_variational_parameters()
{
//    mat energySurface;
    int nCycles, varsteps;
    double steplength, dt, alpha, beta, alphamax, betamax, energy, alphamin, betamin;

//    FILE* file;
//    file = fopen("./energySurface.dat","w");
//    FILE* xyfile;
//    xyfile = fopen("./alpha_beta.dat","w");

    nCycles = 1e6;

    alphamin = 3.5; // 1.8, 3.9
    alphamax = 0.0;
    betamin = 3.5; // 0.36, 4.0
    betamax = 0.0;

    varsteps = 1;
    steplength = 1.0;

//    energySurface = zeros(varsteps, varsteps);

    CBeryllium solver(my_rank, numprocs);
    solver.importanceSampling = 0;
    solver.closedForm = 0;

    for (int nAlpha = 0; nAlpha < varsteps; nAlpha++)
    {
        alpha = alphamin + nAlpha*(alphamax-alphamin)/varsteps;
        if (my_rank == 0) cout << "Alpha = " << alpha << endl;
        for (int nBeta = 0; nBeta < varsteps; nBeta++)
        {
            beta = betamin + nBeta*(betamax-betamin)/varsteps;
            if (my_rank == 0) cout << "Beta = " << beta << endl;

            solver.setParameters(alpha, beta, steplength);
            energy = solver.runMonteCarloIntegration(nCycles);
            if (my_rank == 0) cout << "Energy = " << energy << endl;

//            energySurface(nAlpha, nBeta) = energy;
//            fprintf(file, "%f ", energy);
//            fflush(0); // flush files (write now instead of waiting for program to finish)
        }
//        fprintf(file, "\n");
    }

//    for (int i = 0; i < varsteps; i++)
//    {
//        alpha = alphamin + i*(alphamax-alphamin)/varsteps;
//        fprintf(xyfile, "%f ", alpha);
//    }
//    fprintf(xyfile, "\n");
//    for (int i = 0; i < varsteps; i++)
//    {
//        beta = betamin + i*(betamax-betamin)/varsteps;
//        fprintf(xyfile, "%f ", beta);
//    }
//    fprintf(xyfile, "\n");
}
