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
	int rowB = 150, colB = 100, rowA = 100, colA = 100;
	float* A = (float*)malloc(sizeof(float) * rowA * colA);
	float* B = (float*)malloc(sizeof(float) * rowB * colB);
	memset(A, 0, sizeof(float) * rowA * colA);
	for(int i = 0;i < rowA;i++){
		for(int j = 0;j<i;j++){
			A[i + j * rowA] = (float)(rand() / double(RAND_MAX));
		}
		A[i + i * rowA] = 1;
	}
	for(int i = 0;i < rowB; i++){
		for(int j = 0; j < colB; j++){
			B[i + j * rowB] = (float)(rand() / double(RAND_MAX));
		}
	}
	// print_matrix(A, rowA, colA);
	// print_matrix(B, rowB, colB);
	int start = clock();
	cblas_strsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, rowB, colB, 1.0, A, colB, B, rowB);
	int cost = clock() - start;
	// print_matrix(B, rowB, colB);
	printf("\ntime cost : %f\n", cost * 1.0 / CLOCKS_PER_SEC );

	free(A);
	free(B);
	return 0;
}