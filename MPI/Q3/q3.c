#include <stdio.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	MPI_Init(NULL,NULL);
	int rank,world_size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	char msg[30];
	if(rank == 0)
	{
		MPI_Status status;
		for(int i=1;i<world_size;++i)
		{
			MPI_Recv(msg,30,MPI_CHAR,i,0,MPI_COMM_WORLD,&status);
			printf("%s\n",msg);
		}
	}
	else
	{
		sprintf(msg,"Hello world from %d",rank);
		MPI_Send(msg,30,MPI_CHAR,0,0,MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}