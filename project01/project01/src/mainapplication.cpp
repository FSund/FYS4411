#include "mainapplication.h"

MainApplication::MainApplication()
{
}

void MainApplication::runApplication()
{
    // cout << "MainApplication::runApplication()" << endl;
    test();


}

void MainApplication::test()
{
    mat energySurface;
    int nCycles, varsteps;
    double steplength, alpha, beta, alphamax, betamax, energy, alphamin, betamin;

    FILE* file;
    file = fopen("./energySurface.dat","w");
    FILE* xyfile;
    xyfile = fopen("./alpha_beta.dat","w");

    nCycles = 1000000;

    alphamin = 1.6; // 1.925
    alphamax = 2;
    betamin = 0.3; // 1.3
    betamax = 0.6;

    varsteps = 20;
    steplength = 1.0;

    energySurface = zeros(varsteps, varsteps);

    VMCSolver solver;

    for (double nAlpha = 0; nAlpha < varsteps; nAlpha++)
    {
        alpha = alphamin + nAlpha*(alphamax-alphamin)/varsteps;
        cout << "Alpha = " << alpha << endl;
        for (double nBeta = 0; nBeta < varsteps; nBeta++)
        {
            beta = betamin + nBeta*(betamax-betamin)/varsteps;
            cout << "Beta = " << beta << endl;

            energy = solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta);
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

void MainApplication::minimum()
{
    int nCycles;
    double alpha, beta, steplength;

    nCycles = 1000000;
    alpha = 1.0;
    beta = 1.0;
    steplength = 1.0;

    VMCSolver solver;

//    energy = solver.runMonteCarloIntegration(nCycles, steplength, alpha, beta);

}
