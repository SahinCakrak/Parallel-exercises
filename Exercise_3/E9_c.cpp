#include <mpi.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <cstdlib>
#include <ctime>

int main(int argc, char **argv){

    

    MPI_Init(&argc, &argv);
    int root = 0;
    int mpi_size, mpi_rank;
    double random_number;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process

    
    

    std::vector<double> rfield;
    if(mpi_rank == root){
        rfield.resize(mpi_size);
        std::srand(static_cast<unsigned>(std::time(nullptr)));
        for(int i=0; i<mpi_size;i++){
            random_number = static_cast<double>(std::rand()) / RAND_MAX;
            rfield[i] = random_number;
        }
    }

    double received_number;
    MPI_Scatter(rfield.data(), 1, MPI_DOUBLE, &received_number, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    
    
    for (int i = 0; i < mpi_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (mpi_rank == i) {
            std::cout << "Process " << mpi_rank << " received number: " << received_number << std::endl;
        }
    }


   
    MPI_Finalize();
    
    return 0;
}