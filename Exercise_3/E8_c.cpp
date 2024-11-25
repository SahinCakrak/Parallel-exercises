#include <mpi.h>
#include <iostream>
#include <string.h>


int main(int argc, char **argv){

    

    MPI_Init(&argc, &argv);
    char message[] = "Hello world";
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process
    MPI_Bcast(message, strlen(message)+1, MPI_CHAR, 0, MPI_COMM_WORLD);

    std::cout<<"The 0 rank has send the message '"<<message<<"' to the rank "<<mpi_rank<<std::endl;

    MPI_Finalize();
    
    return 0;
}