#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cublas_v2.h>
#include <curand.h>
#include <cstdlib>
#include <ctime>
void destoryHandle(cublasHandle_t& handle);
void initializeHandle(cublasHandle_t& handle);
void gpu_blas_mmul(cublasHandle_t& handle, const float *A, const float *B, float *C, const int m, const int k, const int n);
void print_matrix(const float *A, int nr_rows_A, int nr_cols_A);
const int CLOCK_RATE = 1410000;
int main(int argc, char const *argv[])
{
    cublasHandle_t handle;
    initializeHandle(handle);
    int rowB = 3, colB = 3, rowA = 3, colA = 3, rowC = 3, colC = 3;
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
	int start = clock();
	gpu_blas_mmul(handle, B, A, C, rowB, colB, colA);
	int cost = clock() - start;
	print_matrix(C, rowC, colC);
	printf("\ntime cost : %f\n", cost * 1.0 / CLOCK_RATE);
    destoryHandle(handle);
    return 0;
}
