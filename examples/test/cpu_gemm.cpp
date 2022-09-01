// LAPACK test code
//compile with: g++ main.cpp -llapack -lblas -o testprog

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <cstring>
using namespace std;
#include "cblas.h"
#include "cblas_f77.h"

void print_matrix(const float *A, int nr_rows_A, int nr_cols_A) {

    for(int i = 0; i < nr_rows_A; ++i){
        for(int j = 0; j < nr_cols_A; ++j){
            std::cout << A[j * nr_rows_A + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int main(){
	int rowB = 150, colB = 100, rowA = 100, colA = 100, rowC = 150, colC = 100;
	// int rowB = 3, colB = 3, rowA = 3, colA = 3, rowC = 3, colC = 3;
	float* A = (float*)malloc(sizeof(float) * rowA * colA);
	float* B = (float*)malloc(sizeof(float) * rowB * colB);
    float* C = (float*)malloc(sizeof(float) * rowC * colC);

	for(int i = 0;i < rowA;i++){
		for(int j = 0;j<colA;j++){
			A[i + j * rowA] = (float)(rand() / double(RAND_MAX));
		}
	}
	for(int i = 0;i < rowB; i++){
		for(int j = 0; j < colB; j++){
			B[i + j * rowB] = (float)(rand() / double(RAND_MAX));
		}
	}
    print_matrix(A, rowA, colA);
    print_matrix(B, rowB, colB);
    float alpha = 1.0;
    float beta = 0;
    int start = clock();
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowB, colA, colB, &alpha, B, rowB, A, colB, &beta, C, rowB);
    int cost = clock() - start;
    print_matrix(C, rowC, colC);

    free(A);
    free(B);
    free(C);
	return 0;
}