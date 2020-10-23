#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BUF 20
#define TAG 42

int main(int argc, char* argv[])
{
    srand(time(0));
    int rank, size, n, r;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0)
    {
      n = rand();
      MPI_Send(n, 4, MPI_INT, 1, TAG, MPI_COMM_WORLD);
    }
    else if(rank<7)
    {
        n= rand();
        MPI_Send(n, 4, MPI_INT, rank+1, TAG, MPI_COMM_WORLD);
        MPI_Recv(r, 4, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
        printf("\nHello world\tRank : %d\tNumber received : %d\n",rank,r);
    }

    MPI_Finalize();
    return 0;
}
