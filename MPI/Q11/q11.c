#include <stdio.h>
#include <mpi.h>
#define N 4
int arr[N][N];
int array_of_displacements[N];
int array_of_blocklengths[N];

//Functiont to display a matrix
void print_matrix(int a[][N])
{
	for(int i=0;i<N;++i)
	{
		for(int j=0;j<N;++j)
			printf("%d ",a[i][j] );
		printf("\n");
	}
}

int main(int argc, char const *argv[])
{
	int rank,world_size,namelen;
	char name[30];

	MPI_Init(NULL,NULL);

	MPI_Comm_size(MPI_COMM_WORLD,&world_size);

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Get_processor_name(name,&namelen);

	//Calculating
	for(int i = 0 ; i < N ;++i)
	{
		array_of_displacements[i] = N*i + i;
		array_of_blocklengths[i] = N - i;
	}	

	MPI_Datatype upper_triangle;
	MPI_Type_indexed(N,array_of_blocklengths,array_of_displacements,MPI_INT,&upper_triangle);
	MPI_Type_commit(&upper_triangle);
	if(rank == 0)
	{
		for(int i = 0; i < N; ++i)
			for(int j = 0; j < N; ++j)
				arr[i][j] = i*N + j;

		printf("Sent Matrix\n");
		print_matrix(arr);
		MPI_Send(arr,1,upper_triangle,1,0,MPI_COMM_WORLD);
	}
	else if(rank == 1)
	{
		
		for(int i = 0; i < N; ++i)
			for(int j = 0; j < N; ++j)
				arr[i][j] = 0;
		
		MPI_Status status;

		MPI_Recv(arr,1,upper_triangle,0,0,MPI_COMM_WORLD,&status);

		printf("Received Matrix\n");
		print_matrix(arr);

	}
	MPI_Finalize();

	return 0;
}