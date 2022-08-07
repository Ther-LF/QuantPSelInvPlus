#include  "ppexsi.hpp"

#include "pexsi/timer.h"

#define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#else
#define MYSCALAR Real
#endif


using namespace PEXSI;
using namespace std;
int main(int argc, char** argv) {
	MPI_Init( &argc, &argv );

	int mpirank, mpisize;
  	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );
	try{
		MPI_Comm world_comm;
		Int nprow = 1;
    	Int npcol = mpisize;
		//Create a communicator with npcol*nprow processors
    	MPI_Comm_split(MPI_COMM_WORLD, mpirank<nprow*npcol, mpirank, &world_comm);
		std::string Hfile = "H180.csc";
		DistSparseMatrix<Real> HMat;
		ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm );
		int dimension = HMat.size;
		int nnz = HMat.nnz;
		std::cout<<"Dimension:"<<dimension<<std::endl;
		std::cout<<"Nonzero value number:"<<nnz<<std::endl;

		IntNumVec colptrs = HMat.colptrLocal; //colptr
		IntNumVec rowind = HMat.rowindLocal; //rowindex
		NumVec<Real> nzval = HMat.nzvalLocal; //values
		
		NumMat<Real> matrix(3600, 3600);
		for(int col = 0;col<dimension;col++){
			int start = colptrs[col] - 1;//start index in row and values
			int end = colptrs[col + 1] - 1;	 //end index in row and values
			// cout<<start<<","<<end<<endl;
			for(int index = start; index < end; index++){
				int row = rowind[index] - 1;
				Real value = nzval[index];
				// cout<<index<<":"<<row<<","<<value<<endl;
				matrix(row, col) = value;
			}
		}
		string fileName = "H180.csv";
		ofstream myfile;
		myfile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
		for(int i = 0;i<3600;i++){
			string line = "";
			line = line + to_string(matrix(i, 0));
			for(int j = 1;j<3600;j++){
				line += ",";
				line += to_string(matrix(i, j));
			}
			line += '\n';
			myfile.write(line.c_str(), line.length());
		}
		myfile.close();
	}catch( std::exception& e )
	{
		std::cerr << "Processor " << mpirank << " caught exception with message: "
		<< e.what() << std::endl;
	}


	MPI_Finalize();
	return 0;
}


