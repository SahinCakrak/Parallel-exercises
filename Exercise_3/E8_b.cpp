#include <mpi.h>
#include <iostream>
#include <string.h>


int main(){

    char[] message = "Hello world";

    MPI_Init(&argc, &argv);

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process
    if(mpi_rank == 0){
        for(int i=1; i<=mpi_size; i++)
        MPI_Send(&message, strlen(message)+1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        printf("Rank 0 sent the message %d to rank %d\n", message, i);

    }
    else{
        char[20] received;
        MPI_Recv(received, 20, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %d received the message: '%d' from rank 0\n", mpi_rank, received);

    }


    MPI_Finalize();
    
    return 0;
}