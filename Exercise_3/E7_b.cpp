#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process

    std::cout << "Hello World" << std::endl;

    MPI_Finalize();
    return 0;
}
