#include <mpi.h>
#include <iostream>
#include <cstring>


int main(int argc, char **argv){

    char Hello[] = "Hello World from process";

    MPI_Init(&argc, &argv);
    char send_buff[300];
    char recv_buff[40];
    char rank_count[40];
    int mpi_size, mpi_rank;
    strncpy(send_buff, Hello, sizeof(send_buff));
    std::cout<<send_buff<<" "<<mpi_rank<<std::endl;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process
    
    for(int i=0; i<mpi_size; i++){
        int send_to = (mpi_rank + 1) % mpi_size;            
        int recv_from = (mpi_rank - 1 + mpi_size) % mpi_size;

        MPI_Send(send_buff, sizeof(send_buff)+1 , MPI_CHAR, send_to, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_buff, sizeof(recv_buff)+1, MPI_CHAR, recv_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        strcat(send_buff, recv_buff);
        std::cout<<"Process "<<mpi_rank<<" of "<<mpi_size<<" received the message: "<<send_buff<<std::endl;
    }

    MPI_Finalize();
    
    return 0;
}