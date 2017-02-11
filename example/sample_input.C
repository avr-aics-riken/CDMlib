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

/*
 *サンプルコードsample_output_Uni.Cおよびsample_output_NonUni.Cで出力されるフィールドデータを、読込み用配列d_vに読み込む。
 *sample_output_NonUni.Cで出力される不等間隔格子のデータを読込む際は、座標ファイルcoord.crdも合わせて読込む。
 */

#include "cdm_DFI.h"
int main( int argc, char **argv )
{
  //CDMのエラーコード
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;

  //MPI Initialize
  if( MPI_Init(&argc,&argv) != MPI_SUCCESS )
  {
     std::cerr << "MPI_Init error." << std::endl;
     return 0;
  }

  //引数で渡されたdfiファイル名をセット
  if( argc != 2 ) {
    //エラー、DFIファイル名が引数で渡されない
    std::cerr << "Error undefined DFI file name." << std::endl;
    return CDM::E_CDM_ERROR;
  }
  std::string dfi_fname = argv[1];

  //計算空間の定義
  int GVoxel[3] = {10, 5, 7};   ///<計算空間全体のボクセルサイズ
  int GDiv[3] = {1, 1, 1};      ///<領域分割数（並列数)
  int head[3] = {1, 1, 1};      ///<計算領域の開始位置
  int tail[3] = {10, 5, 7};     ///<計算領域の終了位置
  int gsize     = 1;            ///<計算空間の仮想セル数
  //読込み配列のサイズ
  size_t size=(GVoxel[0]+2*gsize)*(GVoxel[1]+2*gsize)*(GVoxel[2]+2*gsize);

  //読込み用インスタンスのポインタ取得
  cdm_DFI* DFI_IN = cdm_DFI::ReadInit(MPI_COMM_WORLD, ///<MPIコミュニケータ
                                      dfi_fname,      ///<dfiファイル名
                                      GVoxel,         ///<計算空間全体のボクセルサイズ
                                      GDiv,           ///<領域分割数
                                      ret);           ///<エラーコード
  //エラー処理
  if( ret != CDM::E_CDM_SUCCESS || DFI_IN == NULL ) {
    //エラーインスタンス失敗
    std::cerr << "Error Readinit." << std::endl;
    return ret;
  }

  //読込みフィールドデータ型のチェック
  if( DFI_IN->GetDataType() != CDM::E_CDM_FLOAT32 ) {
    //データの型違い
    std::cerr << "Error Datatype unmatch." << std::endl;
    return CDM::E_CDM_ERROR;
  }

  //読込みフィールドデータの変数の個数を取得
  int nvari=DFI_IN->GetNumVariables();

  //単位系の取得
  std::string Lunit;
  double Lref,Ldiff;
  bool LBset;
  ret=DFI_IN->GetUnit("Length",Lunit,Lref,Ldiff,LBset);
  if( ret==CDM::E_CDM_SUCCESS ) {
    printf("Length\n");
    printf("  Unit      : %s\n",Lunit.c_str());
    printf("  reference : %e\n",Lref);
    if( LBset ) {
      printf("  difference: %e\n",Ldiff);
    }
  }

  //読込み配列のアロケート
  float *d_v = new float[size*nvari];
  //読込み配列のゼロクリア
  memset(d_v, 0, sizeof(float)*size*nvari);
  //読込みフィールドデータのステップ番号をセット
  unsigned step = 10;

  float r_time;       ///<dfiから読込んだ時間
  unsigned i_dummy;   ///<平均化ステップ
  float f_dummy;      ///<平均時間

  //フィールドデータの読込み
  ret =  DFI_IN->ReadData(d_v,    ///<読込み先配列のポインタ
                          step,   ///<読込みフィールドデータのステップ番号
                          gsize,  ///<計算空間の仮想セル数
                          GVoxel, ///<計算空間全体のボクセルサイズ
                          GDiv,   ///<領域分割数
                          head,   ///<計算領域の開始位置
                          tail,   ///<計算領域の終了位置
                          r_time, ///<dfiから読込んだ時間
                          true,   ///<平均を読込まない
                          i_dummy,
                          f_dummy );
  //エラー処理
  if( ret != CDM::E_CDM_SUCCESS ) {
    //フィールドデータの読込み失敗
    std::cerr << "Error ReadData." << std::endl;
    delete [] d_v;
    delete DFI_IN;
    return ret;
  }

  //読込んだ配列のチェック
  /*
  int index_i, index_j, index_k, index_n;
  for (int i=0; i<size*nvari; i++) {
    index_i = i%(GVoxel[0]+2*gsize);
    index_j = (i/(GVoxel[0]+2*gsize))%(GVoxel[1]+2*gsize);
    index_k = (i/((GVoxel[0]+2*gsize)*(GVoxel[1]+2*gsize)))%(GVoxel[2]+2*gsize);
    index_n = i/size;
    std::cout << "d_v for i,j,k,n = " << index_i << " " << index_j << " " << index_k << " " << index_n << " " << d_v[i] << std::endl;
  }
  */

  //正常終了処理
  delete [] d_v;  ///<配列ポインタの削除
  delete DFI_IN;  ///<読込み用インスタンスのポインタの削除
  std::cout << "Normal End." << std::endl;

  return CDM::E_CDM_SUCCESS;
}
