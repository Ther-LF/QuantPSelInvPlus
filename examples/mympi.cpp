#include "mpi.h"
#include "stdio.h"
int main(int argc, char const *argv[])
{
    MPI_Init(NULL, NULL);
    int size;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int number;
    if(rank == 0){
        number = -1;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }else{
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("receive %d\n", number);
    }
    // 释放 MPI 的一些资源
    MPI_Finalize();
    return 0;
}
