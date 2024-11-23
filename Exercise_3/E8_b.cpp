#include <mpi.h>
#include <iostream>
#include <string.h>


int main(int argc, char **argv){

    

    MPI_Init(&argc, &argv);
    const char message[] = "Hello world";
    int mpi_size, mpi_rank, i;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process
    if(mpi_rank == 0){
        for(int i=1; i<mpi_size; i++)
        MPI_Send(&message, strlen(message)+1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        std::cout<<"Rank 0 send the message: "<<message<<"to rank "<<i<<std::endl;

    }
    else{
        char received[20];
        MPI_Recv(&received, 20, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout<<"Rank "<<mpi_rank<<" received the message: "<<received<<std::endl;

    }


    MPI_Finalize();
    
    return 0;
}