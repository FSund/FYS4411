//#include <iostream>
#include <mpi.h>
#include "src/mainapplication.h"

using namespace std;

int main(int argc, char *argv[])
{
    int my_rank, numprocs;
    // MPI initialization
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    MainApplication m(my_rank, numprocs);
    m.runApplication();


    MPI_Finalize();
    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
