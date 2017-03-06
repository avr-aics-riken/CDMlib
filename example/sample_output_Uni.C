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
 *フィールドデータおよびdfiファイルを出力するサンプルコード。格子は等間隔。
 *フィールドデータ出力API関数：WriteFieldDataFile
 *index.dfiファイル出力API関数：WriteIndexDfiFile
 */

#include "cdm_DFI.h"
int main( int argc, char **argv )
{
  //IOのエラーコード
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
  int GDiv[3]   = {1, 1, 1};    ///<領域分割数（並列数）
  int head[3]   = {1, 1, 1};    ///<計算領域の開始位置
  int tail[3]   = {10, 5, 7};   ///<計算領域の終了位置
  int gsize     = 1;            ///<計算空間の仮想セル数
  float pit[3]  = {1.0/10.0, 1.0/10.0, 1.0/10.0} ; ///<ピッチ
  float org[3]  = {-0.5, -0.5, -0.5};              ///<原点座標値
  //配列のサイズ
  size_t size=(GVoxel[0]+2*gsize)*(GVoxel[1]+2*gsize)*(GVoxel[2]+2*gsize);

  std::string path  = "./output_Uni";     ///<出力ディレクトリ
  std::string prefix= "field";  ///<ベースファイル名
  int out_gc        = 1;        ///<出力仮想セル数
  int nvari         = 4;        ///<データの変数の個数
  CDM::E_CDM_FORMAT format  = CDM::E_CDM_FMT_PLOT3D; ///<出力フォーマット
  CDM::E_CDM_DTYPE datatype = CDM::E_CDM_FLOAT32;    ///<データ型
  std::string proc_fname    = "proc.dfi";            ///<procファイル名
  std::string hostname      = "";                    ///<ホスト名
  CDM::E_CDM_ONOFF TimeSliceOnOff = CDM::E_CDM_OFF;  ///<タイムスライス出力指示

  //出力用インスタンスのポインタ取得
  cdm_DFI* DFI_OUT = cdm_DFI::WriteInit<float>(MPI_COMM_WORLD,  ///<MPIコミュニケータ
                                               dfi_fname,       ///<dfiファイル名
                                               path,            ///<出力ディレクトリ
                                               prefix,          ///<ベースファイル名
                                               format,          ///<出力フォーマット
                                               out_gc,          ///<出力仮想セル数
                                               datatype,        ///<データ型
                                               nvari,           ///<データの変数の個数
                                               proc_fname,      ///<procファイル名
                                               GVoxel,          ///<計算空間全体のボクセルサイズ
                                               pit,             ///<ピッチ
                                               org,             ///<原点座標値
                                               GDiv,            ///<領域分割数
                                               head,            ///<計算領域の開始位置
                                               tail,            ///<計算領域の終了位置
                                               hostname,        ///<ホスト名
                                               TimeSliceOnOff); ///<タイムスライス出力オプション
  //エラー処理
  if( DFI_OUT == NULL )
  {
    //エラーインスタンス失敗
    std::cerr << "Error Writeinit." << std::endl;
    return CDM::E_CDM_ERROR;
  }
  //unitの登録
  DFI_OUT->AddUnit("Length","NonDimensional",1.0);
  DFI_OUT->AddUnit("Velocity","NonDimensional",1.0);
  DFI_OUT->AddUnit("Pressure","NonDimensional",0.0,0.0,true);
  //procファイル出力
  int c_id = 3;
  int bc_id = -1;
  DFI_OUT->WriteProcDfiFile(MPI_COMM_WORLD, ///<MPIコミュニケータ
                            false,          ///<ホスト名出力指示
                            c_id,           ///<cell id
                            bc_id);         ///<境界ID

  //IBLANK用配列の設定(PLOT3D形式)
  size_t size_ib=(GVoxel[0]+2*out_gc)*(GVoxel[1]+2*out_gc)*(GVoxel[2]+2*out_gc);
  int *iblank_v;
  if (format == CDM::E_CDM_FMT_PLOT3D) {
    //配列(iblank_v)のアロケート
    iblank_v = new int[size_ib];
    //配列(iblank_v)に値をセット
    for(int i=0; i<size_ib; i++) {
      iblank_v[i] = 1;
    }
  }

  //格子ファイルの出力(PLOT3D形式)
  ret = DFI_OUT->WriteGridFile(iblank_v); ///<plot3dのxyzファイル出力(IBLANK値含む)
  //エラー処理
  if( ret != CDM::E_CDM_SUCCESS ) {
    //格子ファイルの出力失敗
    std::cerr << "Error WriteGridFile." << std::endl;
    return CDM::E_CDM_ERROR;
  }

  // IBLANK用配列のデアロケート
  delete [] iblank_v; ///<IBLANK用配列のポインタの削除

  unsigned step=10;  ///<出力ステップ番号
  float r_time=0.0; ///<出力時間
  float minmax[8];  ///<minmax
  //minmaxのゼロクリア
  for(int i=0; i<8; i++) minmax[i]=0.0;

  //変数名の登録
  DFI_OUT->setVariableName(0,"u");
  DFI_OUT->setVariableName(1,"v");
  DFI_OUT->setVariableName(2,"w");
  DFI_OUT->setVariableName(3,"P");

  //TimeSliceの登録
  DFI_OUT->AddTimeSlice(step,    ///<登録ステップ番号
                        r_time,  ///<登録時間
                        minmax,  ///<最小値，最大値
                        true,    ///<平均出力なし
                        0,       ///<平均をとったステップ数
                        0.0);    ///<平均をとった時刻

  //indexファイル出力
  ret = DFI_OUT->WriteIndexDfiFile();
  //エラー処理
  if( ret != CDM::E_CDM_SUCCESS ) {
    //indexファイルの出力失敗
    std::cerr << "Error WriteIndexDfiFile." << std::endl;
    return CDM::E_CDM_ERROR;
  }

  //配列のアロケート
  float *d_v = new float[size*nvari];

  //配列に値をセット
  for(int i=0; i<size*nvari; i++) {
    d_v[i] = i;
  }

  //フィールドデータの出力
  ret = DFI_OUT->WriteFieldDataFile(step,   ///<出力ステップ番号
                                    r_time, ///<出力時間
                                    GVoxel, ///<d_v の実ボクセル数
                                    nvari,  ///<d_v の成分数
                                    gsize,  ///<d_v の仮想セル数
                                    d_v,    ///<出力するフィールドデータポインタ
                                    true,   ///<平均出力なし
                                    0,      ///<平均をとったステップ数
                                    0.0);   ///<平均をとった時刻
  //エラー処理
  if( ret != CDM::E_CDM_SUCCESS ) {
    //フィールドデータの出力失敗
    std::cerr << "Error WriteFieldDataFile." << std::endl;
    delete [] d_v;
    delete DFI_OUT;
    return ret;
  }

  //正常終了処理
  delete [] d_v;  ///<配列ポインタの削除
  delete DFI_OUT; ///<出力インスタンスのポインタ削除
  std::cout << "Normal End." << std::endl;

  return CDM::E_CDM_SUCCESS;
}
