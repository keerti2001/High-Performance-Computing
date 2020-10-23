#include <mpi.h>
#include <stdio.h>

#define N 1000000
long long numbers[N];

int main(){

    MPI_Init(NULL,NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    long long partial_sum[N]={0};

    for(long long i=0;i<N;i++)
      numbers[i] = i+1;

    for(long long i=rank;i<N;i+=world_size){
        partial_sum[rank] += numbers[i];
    }

    
    long long send_sum,recv_sum;
    MPI_Status status;
    int mid_point = world_size;

    MPI_Barrier(MPI_COMM_WORLD);
    do{

      	if(mid_point%2 != 0)
      	{
            if(rank==world_size-1)
            {
            	send_sum = partial_sum[rank];
                MPI_Send(&send_sum, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
            }

            if(rank==0)
            {
                MPI_Recv(&recv_sum, 1, MPI_LONG_LONG, rank+mid_point-1, 0, MPI_COMM_WORLD, &status);
                partial_sum[rank] = partial_sum[rank] + recv_sum;
                mid_point--;
                continue;
            }
      	}
      	mid_point = mid_point/2;

      	if(rank>=mid_point)
      	{
      		send_sum = partial_sum[rank];
        	MPI_Send(&send_sum, 1, MPI_LONG_LONG, rank-mid_point, 0, MPI_COMM_WORLD); 
      	}
      	else
      	{
        	MPI_Recv(&recv_sum, 1, MPI_LONG_LONG, rank+mid_point, 0, MPI_COMM_WORLD, &status);
        	partial_sum[rank] = partial_sum[rank] + recv_sum;
      	}

  	}while(mid_point>1);

  	MPI_Finalize();
  	if(rank==0)
      printf("Sum: %lld\n",partial_sum[rank]);

}
