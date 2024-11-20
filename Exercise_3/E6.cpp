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

    // Compute matrix-vector product
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
