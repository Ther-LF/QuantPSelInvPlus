gpu_trsm : gpu_trsm.cu
	nvcc gpu_trsm.cu -lcudart -lcublas -L/usr/local/cuda/lib64 -o gpu_trsm

gpu_gemm : gpu_gemm.cu
	nvcc gpu_gemm.cu -lcudart -lcublas -L/usr/local/cuda/lib64 -o gpu_gemm

cpu_trsm : cpu_trsm.cpp
	g++ cpu_trsm.cpp -llapack -lblas -o cpu_trsm

cpu_gemm : cpu_gemm.cpp
	g++ cpu_gemm.cpp -llapack -lblas -o cpu_gemm