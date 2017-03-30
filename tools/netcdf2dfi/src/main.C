/*
###################################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

/**
 * @file   main.C
 * @brief  netcdf2dfi main関数
 * @author aics
 */

#ifdef DISABLE_MPI
 #include "mpi_stubs.h"
#else
 #include "mpi.h"
#endif

#include "netcdf2dfi.h"

/** @brief netcdf2dfiメイン関数
 *  @param[in] argc プログラム引数の数
 *  @param[in] argv プログラム引数
 *  @return 終了コード(0:正常、1:エラーあり)
 */
int main( int argc, char **argv )
{
  int ret;

  // MPI初期化
  if( MPI_Init(&argc,&argv) != MPI_SUCCESS )
  {
      std::cerr << "MPI_Init error." << std::endl;
      return false;
  }
  int nrank, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  printf( "[%d] nrank = %d\n", myrank, nrank );

  // インスタンス
  NetCDF2DFI* ncPtr = NetCDF2DFI::get_instance( argc, argv, nrank, myrank );
  if( !ncPtr )
  {
    MPI_Abort(MPI_COMM_WORLD, 9999);
    return 9;
  }
  if( myrank==0 )
  {
    ncPtr->PrintInputParam();
  }

  // ヘッダーレコードの読み込み
  ret = ncPtr->ReadHeader();
  if( ret != 0 )
  {
    MPI_Abort(MPI_COMM_WORLD, 9999);
    return 9;
  }

  // データレコードの読み込み
  ret = ncPtr->ReadData();
  if( ret != 0 )
  {
    MPI_Abort(MPI_COMM_WORLD, 9999);
    return 9;
  }

  // DFIファイルの出力
  ret = ncPtr->WriteDFI();
  if( ret != 0 )
  {
    MPI_Abort(MPI_COMM_WORLD, 9999);
    return 9;
  }

  // 終了
  if( myrank == 0 )
  {
    printf("\nnormal end\n");
  }

  MPI_Finalize();
  return 0;
}
