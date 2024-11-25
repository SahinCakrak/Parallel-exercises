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
    
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process

    std::srand(static_cast<unsigned>(std::time(nullptr)) + mpi_rank);
    double random_number = static_cast<double>(std::rand()) / RAND_MAX;
    std::cout<<"Process "<<mpi_rank<<" generated number: "<<random_number<<std::endl;

    std::vector<double> rfield;


    if(mpi_rank == root){
        rfield.resize(mpi_size);
    }
    MPI_Gather(&random_number, 1, MPI_DOUBLE, rfield.data(), 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if(mpi_rank == root){
        std::cout<<"The random numbers at root process are:"<<std::endl;
        for (int i = 0; i < mpi_size; i++) {
            std::cout<<"Process "<<i<< ": "<<rfield[i]<<std::endl;
        }
    }
    MPI_Finalize();
    
    return 0;
}