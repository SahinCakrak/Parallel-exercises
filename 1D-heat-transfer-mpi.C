#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <mpi.h>
#include <cmath>

/*
void saveVectorToCSV(const std::vector<double>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return;
    }
    
    for (size_t i = 0; i < data.mpi_size(); ++i) {
        file << data[i];
        if (i != data.mpi_size() - 1) {
            file << ",";  // Separate values with a comma
        }
    }
    file << "\n";  // End with a newline
    file.close();
}

void output_T_mpi(
    int step,
    std::vector<double> T,
    int nx,
    double dt,
    int mpi_rank,
    int mpi_size)
{
    double time = step * dt;

    if (mpi_rank == 0){
    printf("The current time of simulation is: %f \n", time );
    std::cout<<"The current temperatures are the following: \n";
    }

    // get number of Temperatures (nx) per process (local_nx)
    int local_nx = nx / mpi_size;
    std::vector<double> local_T(local_nx);

    // Scatter
    MPI_Scatter(T.data(), local_nx, MPI_DOUBLE, local_T.data(), local_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Each process prints it part of the temp-array
    std::cout << "Process " << mpi_rank << " temperatures: ";
    for (double temp : local_T) {
        std::cout << temp << " ";
    }
    std::cout << std::endl;

    // Gather all parts back to the root process
    MPI_Gather(local_T.data(), local_nx, MPI_DOUBLE, T.data(), local_nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Root process prints the complete temperature array
    if (mpi_rank == 0) {
        std::cout << "Gathered temperatures: ";
        for (double temp : T) {
            std::cout << temp << " ";
        }
        std::cout << std::endl;
    } 
}
*/

//------------------------------------------------------------------------------
// Function: output_T_mpi
// Gathers all local temperature arrays from each rank to root, then prints
// them in the correct order, mimicking the original "output_T" format.
//------------------------------------------------------------------------------
void output_T_mpi(
    int   step,            // current time step (or something akin to time)
    double dt,             // time-step size
    const std::vector<double> &local_T, // local temperature array (with ghost cells)
    int   local_nx,        // "interior" number of grid points per process
    int   nxg,             // local array size = local_nx + 2
    int   nx,              // total number of points in global domain
    int   mpi_rank,
    int   mpi_size)
{
    // The simulation time
    double time = step * dt;

    // We'll gather only the interior data from each rank
    // i.e., we skip ghost cells local_T[0] and local_T[nxg-1].
    std::vector<double> sendbuf(local_nx);
    for(int i = 0; i < local_nx; i++) {
        sendbuf[i] = local_T[i+1]; 
    }

    // The root process receives data from all ranks
    std::vector<double> recvbuf;
    if (mpi_rank == 0) {
        recvbuf.resize(nx);
    }

    // Each rank has local_nx entries (the interior), so we gather them
    MPI_Gather(
        sendbuf.data(),      // send buffer
        local_nx,            // number of elements per rank to send
        MPI_DOUBLE, 
        recvbuf.data(),      // recv buffer (only valid at root)
        local_nx,            // number of elements per rank to receive
        MPI_DOUBLE,
        0,                   // root
        MPI_COMM_WORLD
    );

    // Now rank 0 can output the results
    if (mpi_rank == 0) {
        std::cout << "The current time of simulation is: " << time << std::endl;
        std::cout << "The current temperatures are the following: " << std::endl;
        for(int i = 0; i < nx; i++) {
            std::cout << recvbuf[i] << " ";
        }
        std::cout << std::endl;
    }
}


void copy_T(int nx, std::vector<double> T, std::vector<double> Ttemp)
// helpful for copying temperature in a temporary temperature 
{
    for(int i=0; i<nx; i++){
        Ttemp[i]=T[i];
    }
}

double get_2nd_der(int i, double dx, std::vector<double> T)
// for the second spatial derivative, necessary for temp calculation
{
    double d2Tdx;
    d2Tdx = (T[i+1]-2*T[i]+T[i-1])/(dx*dx);

    return d2Tdx;
}



//------------------------------------------------------------------------------
// Function: communicate_ghost_cells_mpi
// Exchanges boundary data among neighboring ranks to fill ghost cells.
// If mpi_size == 1, do nothing.
//------------------------------------------------------------------------------
void communicate_ghost_cells_mpi(
    std::vector<double> &T, // local temperature array (with ghost cells)
    int nxg,                // local array size = local_nx + 2
    int mpi_rank,
    int mpi_size)
{
    if(mpi_size == 1) {
        // No neighbor processes to communicate with
        return;
    }

    // Left neighbor = rank - 1, right neighbor = rank + 1
    int left_neighbor  = mpi_rank - 1;
    int right_neighbor = mpi_rank + 1;

    // Buffers to send
    double send_left  = T[1];         // leftmost interior cell
    double send_right = T[nxg - 2];   // rightmost interior cell

    // Buffers to receive
    double recv_left  = 0.0;          // to fill T[0]
    double recv_right = 0.0;          // to fill T[nxg-1]

    MPI_Request reqs[4];
    int count = 0;

    // Exchange with left neighbor (if it exists)
    if(left_neighbor >= 0) {
        // Send my left interior cell to neighbor's right ghost cell
        MPI_Isend(&send_left, 1, MPI_DOUBLE, left_neighbor, 0, MPI_COMM_WORLD, &reqs[count++]);
        // Receive left ghost cell from neighbor's right interior cell
        MPI_Irecv(&recv_left, 1, MPI_DOUBLE, left_neighbor, 1, MPI_COMM_WORLD, &reqs[count++]);
    }

    // Exchange with right neighbor (if it exists)
    if(right_neighbor < mpi_size) {
        // Send my right interior cell to neighbor's left ghost cell
        MPI_Isend(&send_right, 1, MPI_DOUBLE, right_neighbor, 1, MPI_COMM_WORLD, &reqs[count++]);
        // Receive right ghost cell from neighbor's left interior cell
        MPI_Irecv(&recv_right, 1, MPI_DOUBLE, right_neighbor, 0, MPI_COMM_WORLD, &reqs[count++]);
    }

    // Wait for all communications to complete
    MPI_Waitall(count, reqs, MPI_STATUSES_IGNORE);

    // Update ghost cells
    if(left_neighbor >= 0) {
        T[0] = recv_left;
    }
    if(right_neighbor < mpi_size) {
        T[nxg - 1] = recv_right;
    }
}


int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Declaration
    int s, nx, nt;
    double lambda, T0, T1, T2, dt, dx, cfl, dfl, tmax;

    //----------------------------------------------------------------------------
    // 1) Let only mpi_rank 0 read in the parameters, then broadcast to all other processes
    //----------------------------------------------------------------------------
    if (mpi_rank==0){
        s = 1000;
        nx = 1000;
        nt = 10000;
        lambda = 2;
        T0 = 0.7;
        T1 = 0;
        T2 = 1;
        tmax = 1000;
    }

    // Broadcast to everyone
    MPI_Bcast(&lambda, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T1,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T2,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T0,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nx,     1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&nt,     1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&s,      1, MPI_INT, 0, MPI_COMM_WORLD);


    //----------------------------------------------------------------------------
    // 2) Domain decomposition 
    //----------------------------------------------------------------------------
    if(nx % mpi_size != 0) {
        if(mpi_rank == 0) {
            std::cerr << "Error: nx must be divisible by mpi_size!" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int local_nx = nx / mpi_size; // interior points per process
    int nxg      = local_nx + 2;  // local array + 2 ghost cells
    dx = s / (double) nx; 
    dt = tmax / (double)nt; 
    cfl = dt / dx;
    dfl = lambda * dt / (dx * dx);

    if(mpi_rank == 0) {
            std::cout << "Number of processes = " << mpi_size << std::endl;
            std::cout << "the Courant-Friedrichs-Lewy-number is: " << cfl << std::endl;
            std::cout << "the diffusive cfl number is: " << dfl << std::endl;
    }

    //----------------------------------------------------------------------------
    // 3) Allocate local temperature arrays
    //    Each process has a local portion of the domain plus 2 ghost cells.
    //----------------------------------------------------------------------------
    std::vector<double> T(nxg, T0);
    std::vector<double> Ttemp(nxg, T0);

    //----------------------------------------------------------------------------
    // 4) Impose global boundary conditions
    //----------------------------------------------------------------------------
    if(mpi_rank == 0) {
        T[1] = T1;  // left boundary
    }
    if(mpi_rank == mpi_size - 1) {
        T[local_nx] = T2; // right boundary
    }
    
    //----------------------------------------------------------------------------
    // 5) Output the initial state
    //----------------------------------------------------------------------------
    if(mpi_rank == 0) {
        std::cout << "Initial Temperature Distribution:" << std::endl;
    }
    output_T_mpi(0, dt, T, local_nx, nxg, nx, mpi_rank, mpi_size);

    //----------------------------------------------------------------------------
    // 6) Main time-stepping loop
    //----------------------------------------------------------------------------
    
    double start_time = MPI_Wtime(); // start timing

    for(int step = 1; step <= nt; step++) {

        // 6a) Communicate ghost cells among neighbors
        communicate_ghost_cells_mpi(T, nxg, mpi_rank, mpi_size);

        // 6b) Make a copy of T into Ttemp
        std::copy(T.begin(), T.end(), Ttemp.begin());


        // 6c) Update the local domain interior points [1..local_nx]
        //     Ghost cells are at indices 0 and nxg-1 => do not update them!
        //     Also apply boundary conditions if rank=0 or rank=mpi_size-1 as needed.
        for(int i = 1; i <= local_nx; i++) {
            // If i=1 and rank=0 => left boundary => skip if you want Dirichlet T1
            // If i=local_nx and rank=mpi_size-1 => right boundary => skip if you want Dirichlet T2
            // Otherwise, do the update.
        
            // For simplicity, we can do:
            if( (i == 1 && mpi_rank == 0) || (i == local_nx && mpi_rank == mpi_size-1) ) {
                // do nothing, keep boundary
            } else {
                T[i] = Ttemp[i] + lambda*dt * get_2nd_der(i, dx, Ttemp);
            }
        }

        // 6d) (Optional) output every k steps, e.g., if step % 1000 == 0
        // if(step % 1000 == 0) {
        //    output_T_mpi(step, dt, T, local_nx, nxg, nx, mpi_rank, mpi_size);
        // }
    }

    double end_time = MPI_Wtime(); // end timing

    // Final distribution
    if(mpi_rank == 0) {
        std::cout << "\nFinal Temperature Distribution:" << std::endl;
    }
    output_T_mpi(nt, dt, T, local_nx, nxg, nx, mpi_rank, mpi_size);

    //----------------------------------------------------------------------------
    // 7) Print timing info
    //----------------------------------------------------------------------------
    double duration = end_time - start_time;
    if(mpi_rank == 0) {
        std::cout << "\nThe duration of the program was " 
                  << duration << " seconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
