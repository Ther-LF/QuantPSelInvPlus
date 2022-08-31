#include <cstdlib>
#include <ctime>
void gpu_blas_mmul(const float *A, const float *B, float *C, const int m, const int k, const int n);
void init_handle();
void del_handle();
int main(){
	// Allocate 3 arrays on CPU
	int nr_rows_A, nr_cols_A, nr_rows_B, nr_cols_B, nr_rows_C, nr_cols_C;

	// for simplicity we are going to use square arrays
	nr_rows_A = nr_cols_A = nr_rows_B = nr_cols_B = nr_rows_C = nr_cols_C = 3;
	

	float *h_A = (float *)malloc(nr_rows_A * nr_cols_A * sizeof(float));
	float *h_B = (float *)malloc(nr_rows_B * nr_cols_B * sizeof(float));
	float *h_C = (float *)malloc(nr_rows_C * nr_cols_C * sizeof(float));

    for(int i = 0;i<nr_rows_A;i++){
        for(int j = 0;j<nr_cols_A;j++){
            h_A[j * nr_rows_A + i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        }
    }
    for(int i = 0;i<nr_rows_B;i++){
        for(int j = 0;j<nr_cols_B;j++){
            h_B[j * nr_rows_B + i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        }
    }
    // init_handle();
    gpu_blas_mmul(h_A, h_B, h_C, nr_rows_A, nr_cols_A, nr_cols_B);
    // del_handle();
    free(h_A);
    free(h_B);
    free(h_C);
    return 0;   
}