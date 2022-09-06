#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>

void gpu_blas_smmul(cublasHandle_t& handle, char transA, char transB, int m, int n, int k, 
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc ){

	// Allocate 3 arrays on GPU
	float *d_A, *d_B, *d_C;
	cudaMalloc(&d_A,m * k * sizeof(float));
	cudaMalloc(&d_B,k * n * sizeof(float));
	cudaMalloc(&d_C,m * n * sizeof(float));

    // Copy the data to device
  cudaMemcpy(d_A, A, m * k * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, B, n * k * sizeof(float), cudaMemcpyHostToDevice);



	// Create a handle for CUBLAS
	// cublasHandle_t handle;
	// cublasCreate(&handle);
	cublasOperation_t opA = CUBLAS_OP_N;
	cublasOperation_t opB = CUBLAS_OP_N;
	if(transA == 'T'){
		opA = CUBLAS_OP_T;
	}
	if(transB == 'T'){
		opB = CUBLAS_OP_T;
	}
	
	// Do the actual multiplication
	cublasSgemm(handle, opA, opB, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc);
	// Destroy the handle
	// cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(C,d_C,m * n * sizeof(float),cudaMemcpyDeviceToHost);

	//Free GPU memory
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);	
}

void gpu_blas_strsm(cublasHandle_t& handle, char side, char uplo, char trans, char unit, int m, int n,
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
	// cublasHandle_t handle;
	// cublasCreate(&handle);

	// Do the actual trsm
	cublasStrsm(handle, cuSide, cuUplo, cuTrans, cuUnit, m, n, &alpha, d_A, lda, d_B, ldb);
	// Destroy the handle
	// cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(B, d_B, m * n * sizeof(float), cudaMemcpyDeviceToHost);

	//Free GPU memory
	cudaFree(d_A);
	cudaFree(d_B);
}

// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
void gpu_blas_dmmul(cublasHandle_t& handle, char transA, char transB, int m, int n, int k, 
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc ){

	// Allocate 3 arrays on GPU
	double *d_A, *d_B, *d_C;
	cudaMalloc(&d_A,m * k * sizeof(double));
	cudaMalloc(&d_B,k * n * sizeof(double));
	cudaMalloc(&d_C,m * n * sizeof(double));

  // Copy the data to device
  cudaMemcpy(d_A, A, m * k * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, B, k * n * sizeof(double), cudaMemcpyHostToDevice);



	// Create a handle for CUBLAS
	// cublasHandle_t handle;
	// cublasCreate(&handle);
	cublasOperation_t opA = CUBLAS_OP_N;
	cublasOperation_t opB = CUBLAS_OP_N;
	if(transA == 'T'){
		opA = CUBLAS_OP_T;
	}
	if(transB == 'T'){
		opB = CUBLAS_OP_T;
	}
	
	// Do the actual multiplication
	cublasDgemm(handle, opA, opB, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc);
	// Destroy the handle
	// cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(C,d_C,m * n * sizeof(double),cudaMemcpyDeviceToHost);

	//Free GPU memory
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);	
}

void gpu_blas_dtrsm(cublasHandle_t& handle, char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb ){

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
	double *d_A, *d_B;
	cudaMalloc(&d_A,rowA * colA * sizeof(double));
	cudaMalloc(&d_B,rowB * colB * sizeof(double));

	// Copy the data to device
    cudaMemcpy(d_A, A, rowA * colA * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, rowB * colB * sizeof(double), cudaMemcpyHostToDevice);

	// Create a handle for CUBLAS
	// cublasHandle_t handle;
	// cublasCreate(&handle);

	// Do the actual trsm
	cublasDtrsm(handle, cuSide, cuUplo, cuTrans, cuUnit, m, n, &alpha, d_A, lda, d_B, ldb);
	// Destroy the handle
	// cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(B, d_B, m * n * sizeof(double), cudaMemcpyDeviceToHost);

	//Free GPU memory
	cudaFree(d_A);
	cudaFree(d_B);
}

void print_matrix(const float *A, int nr_rows_A, int nr_cols_A) {

    for(int i = 0; i < nr_rows_A; ++i){
        for(int j = 0; j < nr_cols_A; ++j){
            std::cout << A[j * nr_rows_A + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
void print_matrix(const double *A, int nr_rows_A, int nr_cols_A) {

    for(int i = 0; i < nr_rows_A; ++i){
        for(int j = 0; j < nr_cols_A; ++j){
            std::cout << A[j * nr_rows_A + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
const int CLOCK_RATE = 1410000;
int main(int argc, char const *argv[])
{
	// int rowB = 150, colB = 100, rowA = 100, colA = 100, rowC = 150, colC = 100;
    int rowB = 3, colB = 3, rowA = 3, colA = 3, rowC = 3, colC = 3;
	double* A = (double*)malloc(sizeof(double) * rowA * colA);
	double* B = (double*)malloc(sizeof(double) * rowB * colB);
    double* C = (double*)malloc(sizeof(double) * rowC * colC);
    
	for(int i = 0;i < rowA;i++){
		for(int j = 0;j<colA;j++){
			A[i + j * rowA] = (double)(rand() / double(RAND_MAX));
		}
	}
	for(int i = 0;i < rowB; i++){
		for(int j = 0; j < colB; j++){
			B[i + j * rowB] = (double)(rand() / double(RAND_MAX));
		}
	}
	print_matrix(A, rowA, colA);
	print_matrix(B, rowB, colB);
	int start = clock();
	cublasHandle_t handle;
	cublasCreate(&handle);
    double alpha = 1;
    double beta = 0;
	gpu_blas_dmmul(handle, 'N', 'N', rowA, colC, colA, alpha, A, rowA, B, rowB, beta, C, rowC);
	cublasDestroy(handle);
	int cost = clock() - start;
	print_matrix(C, rowC, colC);
	printf("\ntime cost : %f\n", cost * 1.0 / CLOCK_RATE);

	free(A);
	free(B);
    free(C);
    return 0;
}
