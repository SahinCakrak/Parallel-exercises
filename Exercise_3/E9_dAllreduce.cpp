#include <mpi.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <numeric>
int main(int argc, char **argv){

    

    MPI_Init(&argc, &argv);
    int root = 0;
    int mpi_size, mpi_rank, dim;
    int pdim = static_cast<int>(dim/(mpi_size));
    double random_number;
    double global_sum;
    double sumprod;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process

    
    

    for(int i=0; i<mpi_size;i++){
        std::vector<double> rfield;
        rfield.reserve(pdim);
        for(int j=0; j<pdim;j++){
            std::srand(static_cast<unsigned>(std::time(nullptr)) + j);
    double random_number = static_cast<double>(std::rand()) / RAND_MAX;
            double sumprod = std::inner_product(rfield.begin(), rfield.end(),rfield.begin(),0);
        }


    }

    MPI_Allreduce(&sumprod, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        std::cout << "Sum of ranks: " << global_sum << std::endl;
    }
   
    MPI_Finalize();
    
    return 0;
}