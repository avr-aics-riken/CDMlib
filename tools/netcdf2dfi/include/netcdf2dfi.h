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
 * @file   netcdf2dfi.h
 * @brief  netcdf2dfiヘッダー
 * @author aics
 */

#ifndef _NETCDF2DFI_H_
#define _NETCDF2DFI_H_

#include "cdm_DFI.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "TextParser.h"

class NetCDF2DFI
{
protected:

  /** @brief コンストラクタ */
  NetCDF2DFI();

  /** @brief usageの表示
   *  @param[in] prog プログラムのパス
   */
  static void Usage( const char *prog );

  /** @brief 入力パラメータファイルの読み込み
   *  @param[in] ifname 入力tpファイル名
   *  @return 終了状態(true:正常、false:エラーあり)
   */
  bool ReadParam( const char *ifname );


public:

  /** @brief デストラクタ */
  ~NetCDF2DFI();

  /** @brief NetCDF2DFIクラスのインスタンス
   *  @param[in] argc   プログラム引数の数
   *  @param[in] argv   プログラム引数
   *  @param[in] nrank  MPI並列数
   *  @param[in] myrank 自身のMPIランク番号
   *  @return インスタンスしたNetCDF2DFIポインタ
   */
  static NetCDF2DFI* get_instance( int argc, char **argv, int nrank, int myrank );

  /** @brief 入力パラメータファイルの印刷
   *  @return 終了状態(true:正常、false:エラーあり)
   */
  bool PrintInputParam();

  /** @brief ncファイル群のヘッダーレコードの読み込み
   *  @return CDMエラーコード
   */
  int ReadHeader();

  /** @brief ncファイル群のデータレコードの読み込み
   *  @return CDMエラーコード
   */
  int ReadData();

  /** @brief DFIファイルの出力
   *  @return CDMエラーコード
   */
  int WriteDFI();

  /** @brief NetCDFファイル出力(変数抽出時のみ)
   *  @param[in]    ID         分散ファイルのファイル番号
   *  @param[in]    step       ステップ番号
   *  @param[in]    time       ステップ時刻
   *  @param[in]    VoxelSize  格子数(ローカル)
   *  @param[inout] varInfoX   x座標値のvar情報(入力:入力ファイルのvar情報、出力:出力ファイルのvar情報)
   *  @param[inout] varInfoY   y座標値のvar情報(入力:入力ファイルのvar情報、出力:出力ファイルのvar情報)
   *  @param[inout] varInfoZ   z座標値のvar情報(入力:入力ファイルのvar情報、出力:出力ファイルのvar情報)
   *  @param[inout] varInfoT   時刻のvar情報(入力:入力ファイルのvar情報、出力:出力ファイルのvar情報)
   *  @param[inout] vecVarInfo 各変数のvar情報(入力:入力ファイルのvar情報、出力:出力ファイルのvar情報)
   *  @param[in]    pArray     出力する配列
   *  @return CDMエラーコード
   */
  int WriteFieldData(int ID, int step, double time, int VoxelSize[3],
                     cdm_DFI_NETCDF::stVarInfo &varInfoX, cdm_DFI_NETCDF::stVarInfo &varInfoY, cdm_DFI_NETCDF::stVarInfo &varInfoZ,
                     cdm_DFI_NETCDF::stVarInfo &varInfoT,
                     vector<cdm_DFI_NETCDF::stVarInfo> &vecVarInfo, cdm_Array *pArray);

protected:

  /** @brief 読み込んだ配列からmin/maxを取得してSlice情報を生成
   *  @param[in] step   ステップ番号
   *  @param[in] time   ステップ時刻
   *  @param[in] pArray 配列
   *  @return Slice情報
   */
  cdm_Slice
  GetSliceInfo(int step, double time, cdm_Array *pArray);

  /** @brief 配列からmin/maxを取得
   *  @param[in]  pArray 配列
   *  @param[out] Min    各成分の最小値
   *  @param[out] Max    各成分の最大値
   */
  template<class T>
  void
  GetMinMax(cdm_TypeArray<T> *pArray, int gc, vector<double> &Min, vector<double> &Max);

  /** @brief 分散ファイル毎に生成したTimeSlice情報をマージ
   *  @param[in] vecFileTimeSlice 自ランクが担当した分散ファイルごとのTimeSlice情報
   */
  void
  mergeTimeSlice(vector<cdm_TimeSlice> &vecFileTimeSlice);

  /** @brief 分散ファイル毎に生成した Domain情報, Process情報のマージ
   *  @param[in] vecFileRank 分散ファイル毎に生成したRank情報
   *  @param[in] coord1      最小原点座標値
   *  @param[in] coord2      最大座標値
   */
  void
  mergeRank(vector<cdm_Rank> &vecFileRank, double coord1[3], double coord2[3]);

  // 変数
public:

  int m_nrank;  ///< MPI並列数
  int m_myrank; ///< 自身のMPIランク番号

  int         m_num_ncfiles;          ///<netcdfファイル数
  std::string m_DirectoryPath;        ///<netcdfファイルが存在するディレクトリ
  std::string m_Prefix;               ///<netcdfファイル名のprefix
  std::string m_RankNoPrefix;         ///<netcdfファイル名文字列内のランク番号前文字列(デフォルト:_id)
  int         m_GuideCell;            ///<GuideCellの数
  int         m_num_steps;            ///<ステップ数(FieldFilenameFormatがstep_rank、rank_stepの場合に必要)
  std::string m_OutputDirectoryPath;  ///<netcdfファイルの出力先(指定されたときのみ使用)
  std::string m_index_fname;          ///<出力index.dfiファイル名
  std::string m_proc_fname;           ///<出力proc.dfiファイル名
  std::vector<std::string> m_vecVariable; ///< 取り出す変数名
  CDM::E_CDM_OUTPUT_FNAME m_FieldFilenameFormat;  ///<netcdfファイル名のファイル名命名方式

  std::string m_nameX; ///< ncファイル内のx座標名
  std::string m_nameY; ///< ncファイル内のy座標名
  std::string m_nameZ; ///< ncファイル内のz座標名
  std::string m_nameT; ///< ncファイル内のtime名

protected:
  int m_indexFsta; ///< 自ランクが担当するファイル番号の開始番号
  int m_indexFend; ///< 自ランクが担当するファイル番号の終了番号:

  set<double> m_setX; ///< 各領域位置の先頭X座標値
  set<double> m_setY; ///< 各領域位置の先頭Y座標値
  set<double> m_setZ; ///< 各領域位置の先頭Z座標値

  // dfi出力用の情報
  cdm_FileInfo      DFI_Finfo;       ///< FileInfo class
  cdm_FilePath      DFI_Fpath;       ///< FilePath class
  cdm_VisIt         DFI_VisIt;       ///< VisIt class
  cdm_Unit          DFI_Unit;        ///< Unit class
  cdm_TimeSlice     DFI_TimeSlice;   ///< TimeSlice class
  cdm_Domain*       DFI_Domain;      ///< Domain class
  cdm_MPI           DFI_MPI;         ///< MPI class
  cdm_Process       DFI_Process;     ///< Process class
};

// MinMaxを取得
template<class T>
void
NetCDF2DFI::GetMinMax(cdm_TypeArray<T> *pArray, int gc, vector<double> &Min, vector<double> &Max)
{
  // 変数の個数
  int nvari = pArray->getNvari();

  // 配列サイズ(仮想セルを含まない)
  const size_t *sz = pArray->_getArraySize();

  // 開始、終了
  size_t sta[3] = {gc, gc, gc};
  size_t end[3] = {sz[0]-gc, sz[1]-gc, sz[2]-gc};

  for( int n=0;n<nvari;n++ )
  {
    double minv = DBL_MAX;
    double maxv = -DBL_MAX;
    for( int k=sta[2];k<end[2];k++ ){
    for( int j=sta[1];j<end[1];j++ ){
    for( int i=sta[0];i<end[0];i++ ){
      double val = pArray->_val(i,j,k,n);
      minv = min(minv, val);
      maxv = max(maxv, val);
    }}}
    Min.push_back(minv);
    Max.push_back(maxv);
  }
}

#endif /* _NETCDF2DFI_H_ */
