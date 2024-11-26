## OpenMP-Syntax
### Set parallel execution block
**Code block**: Complete code block is parallelized:
```
#pragma omp parallel
{
    // Code executed by all threads
}

```

For-Loop
```
#pragma omp parallel for
for (int i = 0; i < N; i++) {
    // Parallelized loop body
}

```

Sections:
```
#pragma omp parallel sections
{
    #pragma omp section
    {
        // Task 1
    }
    #pragma omp section
    {
        // Task 2
    }
}

```

### Synchronize

**Barrier**: Wait for all threads to finish, before further execution
```
#pragma omp barrier
```

**Critical Section**: Define a section that needs to be computed together
```
#pragma omp critical
{
    // Critical section
}

```

**Atomic**: Prevents Race Conditions Only one thread at a time can execute the operation on a shared variable Other threads must wait until the operation is complete. Faster than critical sections
```
#pragma omp atomic
sum += value;

```

**Reduction**: Predifined Operations where communication between threads are necessary without causing Race Conditions. Operators: `+`, `-`, `*`, `min`, `max`
```
#pragma omp parallel for reduction(+:sum)
for (int i = 1; i <= 10; i++) {
	sum += i; // Safe operation
 }
```


## Exercise 5 - OpenMP - HelloWorld
```c++
#include <iostream>
#include <omp.h>
  
int main() {
    // Parallel region
    #pragma omp parallel
    {
        #pragma omp critical
        {
        std::cout << "Hello world from thread " << omp_get_thread_num() << std::endl;
        }
        
        // Synchronize to ensure all threads finish printing before moving forward
        #pragma omp barrier
        
        // Only the master thread prints the total number of threads
        #pragma omp single
        {
            int total_threads = omp_get_num_threads();
            std::cout << "Total number of threads: " << total_threads << std::endl;
        }
    }  

    return 0;
}
```

## Exercise 6 - OpenMP - Matrix Product:
Code:
```c++
#include <iostream>
#include <vector>
#include <omp.h>  

#define N 4 // Size of the matrices and vector  

void print_matrix(const std::vector<std::vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}  

void print_vector(const std::vector<int>& vec) {
    for (const auto& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}  

void matrix_multiply(const std::vector<std::vector<int>>& A,
                     const std::vector<std::vector<int>>& B,
                     std::vector<std::vector<int>>& C) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
 
void matrix_vector_multiply(const std::vector<std::vector<int>>& A,
                            const std::vector<int>& v,
                            std::vector<int>& result) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        result[i] = 0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
}

int main() {
    // -------- Initialize -------------------------------------
    std::vector<std::vector<int>> A = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };

    std::vector<std::vector<int>> B = {
        {16, 15, 14, 13},
        {12, 11, 10, 9},
        {8, 7, 6, 5},
        {4, 3, 2, 1}
    };
    std::vector<int> v = {1, 2, 3, 4};
    std::vector<std::vector<int>> C(N, std::vector<int>(N, 0)); // Resultant matrix for AB
    std::vector<int> result(N, 0); // Resultant vector for Av

    // ------ Compute ------------------------------------
    matrix_multiply(A, B, C);  

    matrix_vector_multiply(A, v, result);  

    // -------- Print results ----------------------------
    std::cout << "Matrix A:" << std::endl;
    print_matrix(A);
    std::cout << "Matrix B:" << std::endl;
    print_matrix(B);
    std::cout << "Matrix A * B:" << std::endl;
    print_matrix(C);
    std::cout << "Matrix A * Vector v:" << std::endl;
    print_vector(result);


    return 0;
}
```


## Exercise 7 - MPI:
### General:
Compile the program:
```
mpicxx hello_world.cpp -o hello_world
```
### Part (a):
- **`MPI_Init()`**:
    
    - **Purpose**: Initializes the MPI execution environment.
    - **Parameters**: Takes two arguments: `int *argc` and `char ***argv` (these can often just be passed as `NULL` if not needed).
- **`MPI_Comm_size()`**:
    
    - **Purpose**: Determines the total number of processes in the communicator.
    - **Parameters**: Takes a communicator (`MPI_Comm`) and an integer pointer to store the size of the communicator.
- **`MPI_Comm_rank()`**:
    
    - **Purpose**: Determines the rank of the calling process within the communicator.
    - **Parameters**: Takes a communicator (`MPI_Comm`) and an integer pointer to store the rank.
- **`MPI_Finalize()`**:
    
    - **Purpose**: Terminates the MPI environment and cleans up resources.
    - **Parameters**: None.
### Part (b):
Code:
```c++
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
```

Run:
```
mpirun -np 1 ./hello_world
mpirun -np 2 ./hello_world
mpirun -np 4 ./hello_world
mpirun -np 8 ./hello_world
```

Result:
"Hello World" is as often printed as the number of used processes

### Part (c)
Code:
```c++
#include <mpi.h>
#include <iostream>  

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);  

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process  

    std::cout << "Process " << mpi_rank << " of " << mpi_size << " says Hello World" << std::endl;  

    MPI_Finalize();
    return 0;
}
```


Result:
The order of the printed processes is not deterministic, because each process in an OpenMPI program runs independently and may execute at a different pace due to variations in CPU scheduling and other system-level factors. Also modern computing environments are dynamic, with other processes or system tasks competing for resources. This can unpredictably affect the runtime behavior of OpenMPI processes.

### Part (d)
Code:
```c++
#include <mpi.h>
#include <iostream>  

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);  

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get rank of current process  

if (mpi_rank == 0) {
        std::cout << "Process " << mpi_rank << " of " << mpi_size << " says Hello World" << std::endl;

    }  

    MPI_Finalize();
    return 0;
}
```
