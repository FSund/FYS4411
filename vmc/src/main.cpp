#include <mpi.h>
#include <src/VMCApp.h>
#include <src/Minimizer/Minimizer.h>

using namespace std;

int main(int argc, char *argv[])
{
    int myRank, numprocs;
    // MPI initialization
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    VMCApp m(myRank, numprocs);
//    m.runApplication();
//    m.diatomic();
    m.minimize();

    MPI_Finalize();
    return 0;
}

// to get QtCreator to run/debug programs correctly:
// $ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
