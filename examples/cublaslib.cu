#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cublas_v2.h>
#include <curand.h>


//Print matrix A(nr_rows_A, nr_cols_A) storage in column-major format
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
void gpu_blas_mmul( char transA, char transB, int m, int n, int k, 
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
	cublasHandle_t handle;
	cublasCreate(&handle);
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
	cublasDestroy(handle);

	// Copy (and print) the result on host memory
	cudaMemcpy(C,d_C,m * n * sizeof(float),cudaMemcpyDeviceToHost);

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
