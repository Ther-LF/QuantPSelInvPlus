#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cublas_v2.h>
#include <curand.h>
#include <cstdlib>
#include <ctime>


extern "C" void strsm_( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb );
// static cublasHandle_t handle;
// void init_handle(){
// 	cublasCreate(&handle);
// }
// void del_handle(){
// 	cublasDestroy(handle);
// }
//Print matrix A(nr_rows_A, nr_cols_A) storage in column-major format
const int CLOCK_RATE = 1410000;
void print_matrix(const float *A, int nr_rows_A, int nr_cols_A) {

    for(int i = 0; i < nr_rows_A; ++i){
        for(int j = 0; j < nr_cols_A; ++j){
            std::cout << A[j * nr_rows_A + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
void gpu_blas_mmul(const float *A, const float *B, float *C, const int m, const int k, const int n) {
	int lda=m,ldb=k,ldc=m;
    const float alf = 1;
    const float bet = 0;
    const float *alpha = &alf;
    const float *beta = &bet;
	// Allocate 3 arrays on GPU
	float *d_A, *d_B, *d_C;
	cudaMalloc(&d_A,m * k * sizeof(float));
	cudaMalloc(&d_B,k * n * sizeof(float));
	cudaMalloc(&d_C,m * n * sizeof(float));

    // Copy the data to device
    cudaMemcpy(d_A, A, m * k * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, n * k * sizeof(float), cudaMemcpyHostToDevice);

    // Print the matrix
    std::cout << "A =" << std::endl;
	print_matrix(A, m, k);
	std::cout << "B =" << std::endl;
	print_matrix(B, k, n);



	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasCreate(&handle);
	// cublasOperation_t opA = CUBLAS_OP_N;
	// cublasOperation_t opB = CUBLAS_OP_N;
	// if(transA == 'T'){
	// 	opA = CUBLAS_OP_T;
	// }
	// if(transB == 'T'){
	// 	opB = CUBLAS_OP_T;
	// }
	
	// Do the actual multiplication
	cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, d_A, lda, d_B, ldb, beta, d_C, ldc);
	// Destroy the handle
	cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(C,d_C,m * n * sizeof(float),cudaMemcpyDeviceToHost);
	std::cout << "C =" << std::endl;
	print_matrix(C, m, n);

	//Free GPU memory
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);	
}

void gpu_blas_trsm( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb ){

	// settings in gpu
	cublasSideMode_t cuSide;
	cublasFillMode_t cuUplo;
	cublasOperation_t cuTrans;
	cublasDiagType_t cuUnit;

	int rowA, colA, rowB = m, colB = n;
	// modify the settings
	if(side == 'L'){
		rowA = lda;
		colA = m;
		cuSide = CUBLAS_SIDE_LEFT;
	}else if(side == 'R'){
		rowA = lda;
		colA = n;
		cuSide = CUBLAS_SIDE_RIGHT;
	}


	if(uplo == 'L'){
		cuUplo = CUBLAS_FILL_MODE_LOWER;
	}else if(uplo == 'U'){
		cuUplo = CUBLAS_FILL_MODE_UPPER;
	}else{ //这里不知道FULL是什么字符就用了else了，后面会补上的
		cuUplo = CUBLAS_FILL_MODE_FULL;
	}

	if(trans == 'T'){
		cuTrans = CUBLAS_OP_T;
	}else if(trans == 'N'){
		cuTrans = CUBLAS_OP_N;
	}

	if(unit == 'U'){
		cuUnit = CUBLAS_DIAG_UNIT;
	}else{ //这里不知道NON_UNIT是什么字符就用了else了，后面会补上的
		cuUnit = CUBLAS_DIAG_NON_UNIT;
	}

	// Allocate 3 arrays on GPU
	float *d_A, *d_B;
	cudaMalloc(&d_A,rowA * colA * sizeof(float));
	cudaMalloc(&d_B,rowB * colB * sizeof(float));

	// Copy the data to device
    cudaMemcpy(d_A, A, rowA * colA * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, rowB * colB * sizeof(float), cudaMemcpyHostToDevice);

	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasCreate(&handle);

	// Do the actual trsm
	cublasStrsm(handle, cuSide, cuUplo, cuTrans, cuUnit, m, n, &alpha, d_A, lda, d_B, ldb);
	// Destroy the handle
	cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(B, d_B, m * n * sizeof(float), cudaMemcpyDeviceToHost);

	//Free GPU memory
	cudaFree(d_A);
	cudaFree(d_B);
}


void printDeviceProp(const cudaDeviceProp &prop)
{
    printf("Device Name : %s.\n", prop.name);
    printf("totalGlobalMem : %lu.\n", prop.totalGlobalMem);
    printf("sharedMemPerBlock : %lu.\n", prop.sharedMemPerBlock);
    printf("regsPerBlock : %d.\n", prop.regsPerBlock);
    printf("warpSize : %d.\n", prop.warpSize);
    printf("memPitch : %lu.\n", prop.memPitch);
    printf("maxThreadsPerBlock : %d.\n", prop.maxThreadsPerBlock);
    printf("maxThreadsDim[0 - 2] : %d %d %d.\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    printf("maxGridSize[0 - 2] : %d %d %d.\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    printf("totalConstMem : %lu.\n", prop.totalConstMem);
    printf("major.minor : %d.%d.\n", prop.major, prop.minor);
    printf("clockRate : %d.\n", prop.clockRate);
    printf("textufloatignment : %lu.\n", prop.textureAlignment);
    printf("deviceOverlap : %d.\n", prop.deviceOverlap);
    printf("multiProcessorCount : %d.\n", prop.multiProcessorCount);
}

int main(){
	int rowB = 15000, colB = 10000, rowA = 10000, colA = 10000;
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
	gpu_blas_trsm('R', 'L', 'N', 'U', rowB, colB, 1.0, A, colB, B, rowB);
	int cost = clock() - start;
	// print_matrix(B, rowB, colB);
	printf("\ntime cost : %f\n", cost * 1.0 / CLOCK_RATE);

	free(A);
	free(B);
	return 0;
}