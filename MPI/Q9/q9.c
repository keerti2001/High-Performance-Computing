#include <stdio.h>
#include <mpi.h>

struct dd
{
	char c;
	int iA[2];
	float fA[4];
};

int count = 3;
int array_of_block_length[3]={1,2,4};
MPI_Aint array_of_displacements[3];
MPI_Datatype array_of_types[3] = {MPI_CHAR,MPI_INT,MPI_FLOAT};

struct dd get_filled_struct(char key)
{

	struct dd temp;
	temp.c = key;
	for(int i=0;i<4;++i)
	{
		if(i<2)
		{
			temp.iA[i] = key + i;
		}
		temp.fA[i] = key + i;
	}

	return temp;
}
int main(int argc, char const *argv[])
{
	int rank,world_size,namelen;
	char name[30];

	MPI_Init(NULL,NULL);

	MPI_Comm_size(MPI_COMM_WORLD,&world_size);

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Get_processor_name(name,&namelen);

	struct dd temp;

	MPI_Get_address(&temp.c,&array_of_displacements[0]);
	MPI_Get_address(temp.iA,&array_of_displacements[1]);
	MPI_Get_address(temp.fA,&array_of_displacements[2]);

	array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
	array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
	array_of_displacements[0] = 0;

	
	MPI_Datatype newtype;

	MPI_Type_create_struct(count,array_of_block_length,array_of_displacements,array_of_types,&newtype);
	MPI_Type_commit(&newtype);


	//Collective communication of the filled structure
	if(rank == 0)
	{
		temp = get_filled_struct('1');
	}

	MPI_Bcast(&temp,1,newtype,0,MPI_COMM_WORLD);
	printf("Collective Communication Rank : %d Structure : %c - %d %d - %f %f %f %f\n",rank,temp.c,temp.iA[0],temp.iA[1],temp.fA[0],temp.fA[1],temp.fA[2],temp.fA[3]);

	//Point-to-point communication of the filled structure
	if(rank == 0)
	{
		temp = get_filled_struct('a');
		for(int i=1;i<world_size;++i)
		{
			MPI_Send(&temp,1,newtype,i,0,MPI_COMM_WORLD);
		}
	}
	else
	{	
		MPI_Status status;
		MPI_Recv(&temp,1,newtype,0,0,MPI_COMM_WORLD,&status);
		printf("Point-to-point communication Rank : %d Structure : %c - %d %d - %f %f %f %f\n",rank,temp.c,temp.iA[0],temp.iA[1],temp.fA[0],temp.fA[1],temp.fA[2],temp.fA[3]);
	}

	MPI_Finalize();

	return 0;
}