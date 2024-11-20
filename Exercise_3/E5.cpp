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
