#ifndef _CONVOUTPUT_NETCDF_H_
#define _CONVOUTPUT_NETCDF_H_

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
 * @file   convOutput_NETCDF.h
 * @brief  convOutput_NETCDF Class Header
 * @author aics
 * @date   2013/11/7
 */
#ifdef _WITH_NETCDF4_

#include "convOutput.h"

class convOutput_NETCDF : public convOutput {

public:

  /** コンストラクタ */
  convOutput_NETCDF();

  /** デストラクタ */
  ~convOutput_NETCDF();

public:


  /**
   * @brief 出力ファイルをオープンする
   * @param [in] prefix ファイル接頭文字
   * @param [in] step ステップ数
   * @param [in] id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  cdm_FILE* OutputFile_Open(
                        const std::string prefix,
                        const unsigned step,
                        const int id,
                        const bool mio);


  /**
   * @brief NetCDFファイルのheaderの書き込み
   * @param[in] step    ステップ数
   * @param[in] dim     変数の個数
   * @param[in] d_type  データ型タイプ
   * @param[in] imax    x方向ボクセルサイズ
   * @param[in] jmax    y方向ボクセルサイズ
   * @param[in] kmax    z方向ボクセルサイズ
   * @param[in] time    時間
   * @param[in] org     原点座標
   * @param[in] pit     ピッチ
   * @param[in] prefix  ファイル接頭文字
   * @param[in] pFile   出力ファイルポインタ
   */
  bool
  WriteHeaderRecord(int step,
                    int dim,
                    CDM::E_CDM_DTYPE d_type,
                    int imax,
                    int jmax,
                    int kmax,
                    double time,
                    double* org,
                    double* pit,
                    const std::string prefix,
                    cdm_FILE *pFile);


  /**
   * @brief Field Datat 出力
   * @param [in] pFile         出力ファイルポインタ
   * @param [in] src           出力データ配列ポインタ
   * @param [in] dLen          出力データサイズ
   * @param [in] d_type        データ型タイプ
   * @param [in] flag_variname VTK形式における変数名出力フラグ
   * @param [in] variname      VTK形式で出力する変数名
   */
  bool
  WriteFieldData(cdm_FILE* pFile,
                 cdm_Array* src,
                 size_t dLen,
                 CDM::E_CDM_DTYPE d_type,
                 bool flag_variname,
                 std::string variname);


protected:

  cdm_DFI_NETCDF::stVarInfo m_varInfoX; ///< xの変数情報
  cdm_DFI_NETCDF::stVarInfo m_varInfoY; ///< yの変数情報
  cdm_DFI_NETCDF::stVarInfo m_varInfoZ; ///< zの変数情報
  cdm_DFI_NETCDF::stVarInfo m_varInfoT; ///< tの変数情報

  /** データレコード出力用のvar id(変数定義順に格納)
   *  write_HeaderRecordで格納し、write_DataRecordで配列出力時に使用
   */
  vector<cdm_DFI_NETCDF::stVarInfo> m_vecVarInfo;

  /** 追記モードを取得
   *  rank & 初回ステップで無いときに追記モード
   *  @param[in] step 出力ステップ番号
   *  @return 追記モードフラグ(true:追記モード,false:時系列分割ファイルモード):
   */
  bool
  CheckAddMode(const unsigned step);
};

#endif /*_WITH_NETCDF4_*/
#endif // _CONVOUTPUT_NETCDF_H_
