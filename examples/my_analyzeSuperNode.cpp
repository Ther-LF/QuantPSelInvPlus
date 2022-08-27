#include  "ppexsi.hpp"

#include "pexsi/timer.h"

// #define _MYCOMPLEX_

#ifdef _MYCOMPLEX_
#define MYSCALAR Complex
#else
#define MYSCALAR Real
#endif
struct NodeInfo{
  int row;
  int col;
  double avg;
  double var;
  int avgRank;
  int varRank;
  NodeInfo(){
    avgRank = 0;
    varRank = 0;
  }
  NodeInfo(int row_, int col_, double avg_, double var_):row(row_), col(col_), avg(avg_), var(var_){
    avgRank = 0;
    varRank = 0;
  }
};

using namespace PEXSI;
using namespace std;

bool cmpAvg(struct NodeInfo x, struct NodeInfo y){
  return abs(x.avg) < abs(y.avg);
}
bool cmpVar(struct NodeInfo x, struct NodeInfo y){
  return x.var < y.var;
}
bool cmpRank(struct NodeInfo x, struct NodeInfo y){
  return (x.avgRank + x.varRank) < (y.avgRank + y.varRank);
}
int main(int argc, char *argv[])
{
    MPI_Init( &argc, &argv );


    int mpirank, mpisize;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpisize );    
    try{
        MPI_Comm world_comm;

        // *********************************************************************
        // Input parameter
        // *********************************************************************
        std::map<std::string,std::string> options;

        OptionsCreate(argc, argv, options);

        // Default processor number
        Int nprow = 1;
        Int npcol = mpisize;
        double delta = atof(options["-delta"].c_str());
        if( options.find("-r") != options.end() ){
            if( options.find("-c") != options.end() ){
                nprow= atoi(options["-r"].c_str());
                npcol= atoi(options["-c"].c_str());
                if(nprow*npcol > mpisize){
                    ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
                } 
            }
            else{
                ErrorHandling( "When using -r option, -c also needs to be provided." );
            }
        }
        else if( options.find("-c") != options.end() ){
            if( options.find("-r") != options.end() ){
                nprow= atoi(options["-r"].c_str());
                npcol= atoi(options["-c"].c_str());
                if(nprow*npcol > mpisize){
                ErrorHandling("The number of used processors cannot be higher than the total number of available processors." );
                } 
            }
            else{
                ErrorHandling( "When using -c option, -r also needs to be provided." );
            }
        }

        //这里将processor分为两部分：在nprow * npcol内和外的，并且用world_comm表示新的COMM
        MPI_Comm_split(MPI_COMM_WORLD, mpirank<nprow*npcol, mpirank, &world_comm);
    if (mpirank<nprow*npcol){
      //这里得到了新的rank和size
      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );
#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
      TAU_PROFILE_SET_CONTEXT(world_comm);
#endif
      //打开每个processor的日至文件
      stringstream  ss;
      ss << "logTest" << mpirank;
      statusOFS.open( ss.str().c_str() );

      ///#if defined(COMM_PROFILE) || defined(COMM_PROFILE_BCAST)
      ///      stringstream  ss3;
      ///      ss3 << "comm_stat" << mpirank;
      ///      commOFS.open( ss3.str().c_str());
      ///#endif


      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  ErrorHandling( "nprow == npcol is assumed in this test routine." );
      //}

      if( mpirank == 0 )
        cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

      std::string Hfile, Sfile;
      int isCSC = true;
      if( options.find("-T") != options.end() ){
        //-T表示txt形式的CSC文件 
        isCSC= ! atoi(options["-T"].c_str());
      }

      int checkAccuracy = false;
      if( options.find("-E") != options.end() ){ 
        //-E还不清楚
        //好像是检查一下LU分解的正确性
        checkAccuracy= atoi(options["-E"].c_str());
      }


      int doSelInv = 1;
      if( options.find("-Sinv") != options.end() ){ 
        //为什么可以大于等于1还不清楚
        doSelInv= atoi(options["-Sinv"].c_str());
      }


      int doFacto = 1;
      if( options.find("-F") != options.end() ){ 
        //同上
        doFacto= atoi(options["-F"].c_str());
      }

      int doSymbolic = 1;
      if( options.find("-Symb") != options.end() ){ 
        //同上
        doSymbolic= atoi(options["-Symb"].c_str());
      }

      int doConvert = 1;
      if( options.find("-C") != options.end() ){ 
        //同上
        doConvert= atoi(options["-C"].c_str());
      }

      int doDistribute = 1;
      if( options.find("-D") != options.end() ){ 
        //同上
        doDistribute= atoi(options["-D"].c_str());
      }

      Int doConstructPattern = 1;
      if( options.find("-Pattern") != options.end() ){ 
        doConstructPattern = atoi(options["-Pattern"].c_str());
      }

      Int doPreSelinv = 1;
      if( options.find("-PreSelinv") != options.end() ){ 
        //同上
        doPreSelinv = atoi(options["-PreSelinv"].c_str());
      }



      if(doSelInv){
        doConvert=1;
      }

      if(doConvert){
        doFacto=1;
      }

      if(doFacto){
        doDistribute=1;
      }

      if(doDistribute){
        doSymbolic = 1;
      }











      int doToDist = true;
      if( options.find("-ToDist") != options.end() ){ 
        //?
        doToDist= atoi(options["-ToDist"].c_str());
      }

      int doDiag = false;
      if( options.find("-Diag") != options.end() ){ 
        //?
        doDiag = atoi(options["-Diag"].c_str());
      }


      if( options.find("-H") != options.end() ){ 
        //输入H矩阵的路径
        Hfile = options["-H"];
      }
      else{
        ErrorHandling("Hfile must be provided.");
      }

      if( options.find("-S") != options.end() ){ 
        //输入S矩阵的路径
        Sfile = options["-S"];
      }
      else{
        statusOFS << "-S option is not given. " 
          << "Treat the overlap matrix as an identity matrix." 
          << std::endl << std::endl;
      }

      Int maxPipelineDepth = -1;
      if( options.find("-P") != options.end() ){ 
        //?
        maxPipelineDepth = atoi(options["-P"].c_str());
      }
      else{
        statusOFS << "-P option is not given. " 
          << "Do not limit SelInv pipelining depth." 
          << std::endl << std::endl;
      }

      Int symmetricStorage = 0;
      if( options.find("-SS") != options.end() ){ 
        //?
        symmetricStorage = atoi(options["-SS"].c_str());
      }
      else{
        statusOFS << "-SS option is not given. " 
          << "Do not use symmetric storage." 
          << std::endl << std::endl;
      }





      Int numProcSymbFact;
      if( options.find("-npsymbfact") != options.end() ){ 
        //?
        numProcSymbFact = atoi( options["-npsymbfact"].c_str() );
      }
      else{
        statusOFS << "-npsymbfact option is not given. " 
          << "Use default value (maximum number of procs)." 
          << std::endl << std::endl;
        numProcSymbFact = 0;
      }

      //rshift和shift为最后组成z这个复数的两部分
      Real rshift = 0.0, ishift = 0.0;
      if( options.find("-rshift") != options.end() ){ 
        //?
        rshift = atof(options["-rshift"].c_str());
      }
      if( options.find("-ishift") != options.end() ){ 
        //?
        ishift = atof(options["-ishift"].c_str());
      }


      std::string ColPerm;
      if( options.find("-colperm") != options.end() ){ 
        //这是用在对列进行置换的东西，用来加强LU的sparsity pattern
        ColPerm = options["-colperm"];
      }
      else{
        statusOFS << "-colperm option is not given. " 
          << "Use MMD_AT_PLUS_A." 
          << std::endl << std::endl;
        ColPerm = "MMD_AT_PLUS_A";
      }

      // *********************************************************************
      // Read input matrix
      // *********************************************************************

      // Setup grid.
      SuperLUGrid<MYSCALAR> g( world_comm, nprow, npcol );

      int      m, n;
      DistSparseMatrix<MYSCALAR>  AMat;

      DistSparseMatrix<Real> HMat;
      DistSparseMatrix<Real> SMat;
      Real timeSta, timeEnd;
      GetTime( timeSta );
      if(isCSC)
        ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm );  
      else{
        ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
        ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
      }

      if( Sfile.empty() ){
        // Set the size to be zero.  This will tell PPEXSI.Solve to treat
        // the overlap matrix as an identity matrix implicitly.
        SMat.size = 0;  
      }
      else{
        if(isCSC)
          ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
        else{
          ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, world_comm ); 
          ParaWriteDistSparseMatrix( "S.csc", SMat, world_comm ); 
        }
      }

      GetTime( timeEnd );
      LongInt nnzH = HMat.Nnz();
      if( mpirank == 0 ){
        cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
        cout << "H.size = " << HMat.size << endl;
        cout << "H.nnz  = " << nnzH  << endl;
      }

      // Get the diagonal indices for H and save it n diagIdxLocal_
      //虽然不是很懂这一步是为啥，但是还是跟着看吧
      //好像是用来在后面减去zshift的时候要减对角线的值用的
      std::vector<Int>  diagIdxLocal;
      { 
        Int numColLocal      = HMat.colptrLocal.m() - 1; //本地包含的column的数量
        Int numColLocalFirst = HMat.size / mpisize;
        Int firstCol         = mpirank * numColLocalFirst;//得到本processor的第一个列的序号

        diagIdxLocal.clear();

        for( Int j = 0; j < numColLocal; j++ ){
          Int jcol = firstCol + j + 1;
          for( Int i = HMat.colptrLocal(j)-1; 
              i < HMat.colptrLocal(j+1)-1; i++ ){
            Int irow = HMat.rowindLocal(i);
            if( irow == jcol ){
              diagIdxLocal.push_back( i );//把所有这个processor保存的矩阵在原来矩阵中为对角线的本地坐标保留了下来
            }
          }
        } // for (j)
      }


      GetTime( timeSta );

      AMat.size          = HMat.size;
      AMat.nnz           = HMat.nnz;
      AMat.nnzLocal      = HMat.nnzLocal;
      AMat.colptrLocal   = HMat.colptrLocal;
      AMat.rowindLocal   = HMat.rowindLocal;
      AMat.nzvalLocal.Resize( HMat.nnzLocal );
      AMat.comm = world_comm;

      MYSCALAR *ptr0 = AMat.nzvalLocal.Data();
      Real *ptr1 = HMat.nzvalLocal.Data();
      Real *ptr2 = SMat.nzvalLocal.Data();

#ifdef _MYCOMPLEX_
      Complex zshift = Complex(rshift, ishift);
#else
      Real zshift = Real(rshift);
#endif

      if( SMat.size != 0 ){
        // S is not an identity matrix
        for( Int i = 0; i < HMat.nnzLocal; i++ ){
          AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
        }
      }
      else{
        // S is an identity matrix
        for( Int i = 0; i < HMat.nnzLocal; i++ ){
          AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
        }

        for( Int i = 0; i < diagIdxLocal.size(); i++ ){
          AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
        }
      } // if (SMat.size != 0 )

      LongInt nnzA = AMat.Nnz();
      if( mpirank == 0 ){
        cout << "nonzero in A (DistSparseMatrix format) = " << nnzA << endl;
      }


      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;

      //目前为止就得到了A矩阵的所有内容
      // *********************************************************************
      // Symbolic factorization ：找到关于矩阵A的LU分解中的sparsity pattern，从而对supernode进行了分析
      // *********************************************************************
    

      GetTime( timeSta );
      SuperLUOptions luOpt;//一个传递Super_LU的接口
      luOpt.ColPerm = ColPerm;
      luOpt.numProcSymbFact = numProcSymbFact;


      SuperLUMatrix<MYSCALAR> luMat(g, luOpt );//设置初始化SuperLU Matrix，用来存储和Super_LU相关的数据


      luMat.DistSparseMatrixToSuperMatrixNRloc( AMat, luOpt );//将CSC格式转化为SuperLU的CSR压缩格式，结果保留在luMat的SuperMatrix中
      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;

      if(doSymbolic){   
        GetTime( timeSta );
        luMat.SymbolicFactorize();//进行symbolic factorization
        luMat.DestroyAOnly();//消除矩阵A，保留LU矩阵
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for performing the symbolic factorization is " << timeEnd - timeSta << endl;
      }

      // *********************************************************************
      // Numerical factorization only : LU分解
      // *********************************************************************

      Real timeTotalFactorizationSta, timeTotalFactorizationEnd; 


      // Important: the distribution in pzsymbfact is going to mess up the
      // A matrix.  Recompute the matrix A here.
      //调用pzsymbfact会吧SuperMatrix搞乱，这里重新计算
      luMat.DistSparseMatrixToSuperMatrixNRloc( AMat ,luOpt);

      GetTime( timeTotalFactorizationSta );

      if(doDistribute){
        GetTime( timeSta );
        luMat.Distribute();//再次并行分发数据，为numerical factorization做准备，说实话不是很懂
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl; 
      }


      if(doFacto){
        GetTime( timeSta );
        luMat.NumericalFactorize();//进行LU分解，不知道为什么要写成Numerical Factorization，结果记录在LUStruct中了
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 

        GetTime( timeTotalFactorizationEnd );
        if( mpirank == 0 )
          cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 


        // *********************************************************************
        // Test the accuracy of factorization by solve
        // *********************************************************************
        //这一个检查就暂时不看了，看不懂
        if( checkAccuracy ) {
          SuperLUMatrix<MYSCALAR> A1( g, luOpt );
          SuperLUMatrix<MYSCALAR> GA( g, luOpt );
          A1.DistSparseMatrixToSuperMatrixNRloc( AMat,luOpt );
          A1.ConvertNRlocToNC( GA );

          int n = A1.n();
          int nrhs = 5;
          NumMat<MYSCALAR> xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
          NumMat<MYSCALAR> xTrueLocal, bLocal;
          DblNumVec berr;
          UniformRandom( xTrueGlobal );

          GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );

          A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
          A1.DistributeGlobalMultiVector( bGlobal,     bLocal );

          luMat.SolveDistMultiVector( bLocal, berr );
          luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );
        }

        // *********************************************************************
        // Selected inversion
        // *********************************************************************

        if(doConvert || doSelInv>=1)
        {
          Real timeTotalSelInvSta, timeTotalSelInvEnd;

          NumVec<MYSCALAR> diag;
          PMatrix<MYSCALAR> * PMlocPtr;
          SuperNodeType * superPtr;
          GridType * g1Ptr;

          GetTime( timeTotalSelInvSta );

          g1Ptr = new GridType( world_comm, nprow, npcol );//重点来了，这里申请了一个nprow * npcol的网格，可能是supernode划分了
          GridType &g1 = *g1Ptr;

          superPtr = new SuperNodeType();
          SuperNodeType & super = *superPtr;

          GetTime( timeSta );
          luMat.SymbolicToSuperNode( super ); //根据sparsity pattern得到supernode

          PSelInvOptions selInvOpt;//一些PSelInv的设置，不懂
          selInvOpt.maxPipelineDepth = maxPipelineDepth;
          selInvOpt.symmetricStorage = symmetricStorage;

          FactorizationOptions factOpt;
          factOpt.ColPerm = ColPerm;//一些factorization的设置，不懂
          PMlocPtr = new PMatrix<MYSCALAR>( &g1, &super, &selInvOpt, &factOpt  );//PSelInv的主要数据结构
          PMatrix<MYSCALAR> & PMloc = *PMlocPtr;

          if(doConvert){
            luMat.LUstructToPMatrix( PMloc );//把LU也放进来了？这个意思吧，应该是=。=
            GetTime( timeEnd );
          }

          LongInt nnzLU = PMloc.Nnz();
          if( mpirank == 0 ){
            cout << "nonzero in L+U  (PMatrix format) = " << nnzLU << endl;
          }

          //记录本processor所拥有的所有supernode，记录其下标和属性
          vector<struct NodeInfo> nodeInfos;
          int superNum = PMloc.NumSuper();
          const GridType* grid_ = PMloc.Grid();
          IntNumVec superPtr_ = super.superPtr;
          if(mpirank == 0){
            cout<<"superNum:"<<superNum<<endl;
          }
          for(int ksup = 0;ksup<superNum;ksup++){
            if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){//如果本processor和supernode在同一列，说明本processor确实包含了部分supernode信息
              vector<LBlock<Real> >& Lcol = PMloc.L( LBj( ksup, grid_ ) );//得到所有的LBlock
              for(LBlock<Real> LB : Lcol){//对于每个L Block，记录它的属性
                int col = ksup;//记录L Block的列
                int row = LB.blockIdx;//记录L Block的行
                NumMat<Real> nzval = LB.nzval;
                //下面开始计算总和
                double sum = 0;
                for(int j = 0;j<LB.numCol;j++){
                  for(int i = 0;i<LB.numRow;i++){
                    sum += abs(nzval(i, j));
                  }
                }
                //计算平均值
                int numCol = superPtr_(col + 1) - superPtr_(col);//这个supernode具有的列数
                int numRow = superPtr_(row + 1) - superPtr_(row);//这个supernode具有的行数
                // int numCol = LB.numCol;//用nzval的col
                // int numRow = LB.numRow;//用nzval的row
                double avg = sum / (numCol * numRow);
                //计算方差
                double var = 0;
                for(int j = 0;j<LB.numCol;j++){
                  for(int i = 0;i<LB.numRow;i++){
                    var += (abs(nzval(i, j)) - avg) * (abs(nzval(i, j)) - avg);
                  }
                }
                var /= (numCol * numRow);
                //将这个L Block放入数组中用于后续发送
                struct NodeInfo tmp(row, col, avg, var);
                nodeInfos.push_back(tmp);
                // cout<<"("<<tmp.row<<","<<tmp.col<<"):"<<tmp.avg<<" "<<tmp.var<<endl;
              }
            }
          }
          //现在申请MPI_Type来发送NodeInfo类型的数组
          int blockLength[] = {1, 1, 1, 1, 1, 1};
          MPI::Datatype oldTypes[] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
          MPI::Aint addressOffsets[6];
          struct NodeInfo tmp;
          MPI_Address(&tmp.row, &addressOffsets[0]);
          MPI_Address(&tmp.col, &addressOffsets[1]);
          MPI_Address(&tmp.avg, &addressOffsets[2]);
          MPI_Address(&tmp.var, &addressOffsets[3]);
          MPI_Address(&tmp.avgRank, &addressOffsets[4]);
          MPI_Address(&tmp.varRank, &addressOffsets[5]);
          addressOffsets[5] = addressOffsets[5] - addressOffsets[0];
          addressOffsets[4] = addressOffsets[4] - addressOffsets[0];
          addressOffsets[3] = addressOffsets[3] - addressOffsets[0];
          addressOffsets[2] = addressOffsets[2] - addressOffsets[0];
          addressOffsets[1] = addressOffsets[1] - addressOffsets[0];
          addressOffsets[0] = 0;
          MPI::Datatype newType = MPI::Datatype::Create_struct(6, blockLength, addressOffsets, oldTypes);
          newType.Commit();
          //首先发送数组的大小，方便后面使用MPI_Gatherv函数
          int *recvCount = NULL;
          if(mpirank == 0){
            recvCount = (int*)malloc(sizeof(int) * mpisize);
          }
          int cnt = nodeInfos.size();
          MPI_Gather(&cnt, 1, MPI_INT, recvCount, 1, MPI_INT, 0, world_comm);
          //然后每个processor发送不同尺寸的数组
          struct NodeInfo* recvBuf = NULL;
          int* displs = NULL;
          int recvSize = 0;
          if(mpirank == 0){
            displs = (int*)malloc(sizeof(int) * mpisize);
            displs[0] = 0;
            for(int i = 1;i<mpisize;i++){
              displs[i] = displs[i - 1] + recvCount[i - 1];
            }
            recvSize = displs[mpisize - 1] + recvCount[mpisize - 1];
            recvBuf = (struct NodeInfo*)malloc(sizeof(struct NodeInfo) * recvSize);
          }
          MPI_Gatherv(nodeInfos.data(), nodeInfos.size(), newType, recvBuf, recvCount, displs, newType, 0, world_comm);
          if(mpirank == 0){
            // for(int i = 0;i<recvSize;i++){
            //   cout<<"(" <<recvBuf[i].row << "," <<recvBuf[i].col <<"):"<<recvBuf[i].avg<<" "<<recvBuf[i].var<<endl;
            // }
             //得到数据后，首先将他们按照平均值的绝对值从小到大排序
             cout<<"quant number:"<<recvSize<<endl;
             sort(recvBuf, recvBuf + recvSize, cmpAvg);
//             for(int i = 0;i<recvSize;i++){
//               recvBuf[i].avgRank = i;
//             }
//             //然后按照方差的顺序从小到达排序
//             sort(recvBuf, recvBuf + recvSize, cmpVar);
//             for(int i = 0;i<recvSize;i++){
//               recvBuf[i].varRank = i; 
//             }
             //最后综合排名排序
//             sort(recvBuf, recvBuf + recvSize, cmpRank);
             ofstream superCol;
             superCol.open("my_quant3600Info");

            // for(int i = 0;i<recvSize;i++){//保留所有supernode
            //   string pos = "(" + to_string(recvBuf[i].row) + "," + to_string(recvBuf[i].col) + ")";
            //   superCol << pos << "\t\t\t\t" << recvBuf[i].avg << "\t\t\t\t" << recvBuf[i].var << endl;
            // }

            // for(int i = 0;i<1000;i++){//保留前1000个平均值小的supernode就行了
            //   string pos = "(" + to_string(recvBuf[i].row) + "," + to_string(recvBuf[i].col) + ")";
            //   superCol << pos << "\t\t\t\t" << recvBuf[i].avg << "\t\t\t\t" << recvBuf[i].var << endl;
            // }
            // superCol << endl;
            // superCol << "------------------------------------The last 1000 supernode-----------------------------------------"<<endl;
            // superCol << endl;
            for(int i = 0; i < recvSize; i++){//保留后1000个平均值小的supernode
              string pos = "(" + to_string(recvBuf[i].row) + "," + to_string(recvBuf[i].col) + ")";
              superCol << pos << "\t\t\t\t" << recvBuf[i].avg << "\t\t\t\t" << recvBuf[i].var << endl;
            }   
            // superCol << endl;
            superCol.close();
  
 }
          //最后释放这些资源
         // superCol.close();
          newType.Free();
          free(recvBuf);
          free(recvCount);
          free(displs);
          
          delete PMlocPtr;
          delete superPtr;
          delete g1Ptr;

        }
      }
        statusOFS.close();
    
  }
    }catch( std::exception& e ){
        std::cerr << "Processor " << mpirank << " caught exception with message: "
        << e.what() << std::endl;
    }

    MPI_Finalize();

    return 0;
}
