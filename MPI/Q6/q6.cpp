#include<mpi.h>
#include<bits/stdc++.h>

using namespace std;

int main(int argc, char* argv[]) {
    int n = 1e3, num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    vector<vector<double>> A(n, vector<double>(n)), original(n, vector<double>(n));
    double TOL = 1e8, div;
    if(rank == 0) {
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++) original[i][j] = 1.0 * (rand() / 1e7);
        A = original;
        div = n / num_procs;
    }
    double nettime = 0;
    bool done = false;
    while(!done) {
        double diff = 0;
        original = A;
        MPI_Barrier(MPI_COMM_WORLD);
        double t = MPI_Wtime();
        for(int i = rank*div; i < (rank+1)*div; i++) {
            for(int j = 0; j < n; j++) {
                double sum = original[i][j];
                if(i-1 >= 0) sum += original[i-1][j];
                if(i+1 < n) sum += original[i+1][j];
                if(j+1 < n) sum += original[i][j+1];
                if(j-1 >= 0) sum += original[i][j-1];
                sum *= 0.2;
                A[i][j] = sum;
                diff += abs(original[i][j] - sum);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        nettime += MPI_Wtime() - t;
        if(diff < TOL) done = true;
    }
    if(rank == 0)
        cout << "MultiProc Time: " << 1000*nettime << " milliseconds.\n";
    MPI_Finalize();
    
}