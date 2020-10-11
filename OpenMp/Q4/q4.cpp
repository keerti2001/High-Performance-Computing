
#include <iostream>
#include <cstdlib>
#include "omp.h"
#include <vector>
#include <chrono>
using namespace std;
using namespace std::chrono;


#define TMAX   27       //Max no of threads
#define ITERATIONS 4    //No of iterations to negate random factors

#define M 250
#define N 250
#define P 250

int A[M][N], B[N][P], C[M][P];

void operation(int num_t)
{//Creates the threads and waits for them to complete execution
   omp_set_num_threads(num_t);
   #pragma omp parallel
   {
      #pragma omp for
      for(int i=0; i<M; i++)
      {
        for(int j=0; j<P; j++)
        {
            for(int k=0; k<N; k++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
      }
   }
   
   return;
}

int calcAvgTime(int n)
{
      //clock function
      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      
      // operation(); Pass the number of threads and wait for execution
      operation(n);

      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>( t2 - t1 ).count();

   return (int)duration;
}

void populateMatrices()
{
   for (int i = 0; i < M; ++i)
   {
      for (int j = 0; j < N; ++j)
      {
         A[i][j] = rand() % 5;
      }
   }

   for (int i = 0; i < N; ++i)
   {
      for (int j = 0; j < P; ++j)
      {
         B[i][j] = rand() % 5;
      }
   }
}

int main () {
   populateMatrices();
   int avg_times[TMAX+1];

   for (int i = 1; i <= TMAX; ++i)
   {//run once prior to negate any cache issues that effect exec time
      avg_times[i] = calcAvgTime(i);
   }
   for (int i = 1; i <= TMAX; ++i)
   {
      avg_times[i] = 0;
   }
   
   //run over a number of iterations
   for (int x = 0; x < ITERATIONS; ++x)
   {
      for (int i = 1; i <= TMAX; ++i)
      {
         avg_times[i] += calcAvgTime(i);
      }
   }

   for (int i = 1; i <= TMAX; ++i)
   {
      avg_times[i] /= ITERATIONS;
   }
   cout<<"Num_threads\tExec Time (micro-s)"<<endl;
   for (int i = 1; i <= TMAX; ++i)
   {
      cout<<i<<"          \t"<<avg_times[i]<<endl;
   }

   return 0;
}
