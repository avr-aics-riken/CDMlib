#ifndef _CONVOUTPUT_H_
#define _CONVOUTPUT_H_

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
 * @file   convOutput.h
 * @brief  convOutput Class Header
 * @author aics
 * @date   2013/11/7
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <errno.h>

#include "cpm_ParaManager.h"
#include "cdm_DFI.h"
#include "conv_Define.h"
#include "InputParam.h"

using namespace std;

class convOutput {

public:

  InputParam* m_InputCntl;
  const cdm_TimeSlice* m_pTSlice;
  const cdm_FileInfo*  m_pFinfo;
  const cdm_Unit*  m_pUnit;

public:

  /** コンストラクタ */
  convOutput();

  /** デストラクタ */
  ~convOutput();

  /**
   * @brief InputParamのポインタをコピー
   * @param [in] InputCntl InputParamクラスポインタ
   * @return  エラーコード
   */
  bool importInputParam(InputParam* InputCntl);

  /**
   * @brief 出力クラスのインスタンス
   * @param[in] out_format 出力ファイルフォーマット
   * @return convOutputクラスポインタ
   */
  static convOutput*
  OutputInit(const CDM::E_CDM_FORMAT out_format);

  /**
   * @brief 出力ファイルをオープンする
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   ステップ数
   * @param [in] id     ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  virtual cdm_FILE*
  OutputFile_Open(
                  const std::string prefix,
                  const unsigned step,
                  const int id,
                  const bool mio=false)
  { return NULL; };

  /**
   * @brief 出力ファイルをクローズする
   * @param [in] fp ファイルポインタ
   */
  virtual void
  OutputFile_Close(cdm_FILE* fp)
  { cdm_FILE::CloseFile(fp); };

  /**
   * @brief grid 出力(plot3d用)
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   step番号
   * @param [in] myRank ランク番号
   * @param [in] dType  dfiのデータタイプ
   * @param [in] guide  ガイドセル数
   * @param [in] org    原点座標値
   * @param [in] pit    ピッチ
   * @param [in] sz     ボクセルサイズ
   */
  virtual void
  WriteGridData(std::string prefix,
                int step,
                int myRank,
                int dType,
                int guide,
                double org[3],
                double pit[3],
                int sz[3]) {};

  /**
   * @brief grid 出力(plot3d用, 不等間隔格子対応版)
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   step番号
   * @param [in] myRank ランク番号
   * @param [in] dType  dfiのデータタイプ
   * @param [in] guide  ガイドセル数
   * @param [in] out_domain  出力用のdomainインスタンス
   * @param [in] out_process 出力用のprocessインスタンス
   */
  virtual void
  WriteGridData(std::string prefix,
                int step,
                int myRank,
                int dType,
                int guide,
                cdm_Domain* out_domain,
                cdm_Process* out_proces) {};

  /**
   * @brief ファイルのheaderの書き込み
   * @param[in] step     ステップ数
   * @param[in] dim      変数の個数
   * @param[in] d_type   データ型タイプ
   * @param[in] imax     x方向ボクセルサイズ
   * @param[in] jmax     y方向ボクセルサイズ
   * @param[in] kmax     z方向ボクセルサイズ
   * @param[in] time     時間
   * @param[in] org      原点座標
   * @param[in] pit      ピッチ
   * @param[in] prefix   ファイル接頭文字
   * @param[in] pFile    ファイルポインタ
   */
  virtual bool
  WriteHeaderRecord(int step,
                    int dim,
                    CDM::E_CDM_DTYPE d_type,
                    int imax,
                    int jmax,
                    int kmax,
                    double time,
                    double* org,
                    double* pit,
                    std::string prefix,
                    cdm_FILE *pFile)
  { return true; };

  /**
   * @brief ファイルのheaderの書き込み(vtk用、不等間隔格子対応版)
   * @param[in] step        ステップ数
   * @param[in] time        時間
   * @param[in] dim         変数の個数
   * @param[in] d_type      データ型タイプ
   * @param[in] dfi_type    dfi種別
   * @param[in] out_domain  出力用のdomainインスタンス
   * @param[in] out_process 出力用のprocessインスタンス
   * @param[in] gc          出力ガイドセル数
   * @param[in] prefix      ファイル接頭文字
   * @param[in] pFile       出力ファイルポインタ
   */
  virtual bool
  WriteHeaderRecord(int step,
                    double time,
                    int dim,
                    CDM::E_CDM_DTYPE d_type,
                    CDM::E_CDM_DFITYPE dfi_type,
                    cdm_Domain* out_domain,
                    cdm_Process* out_process,
                    int gc,
                    const std::string prefix,
                    cdm_FILE *pFile)
  { return true; };

  /**
   * @brief Field Datat 出力
   * @param [in] pFile  出力ファイルポインタ
   * @param [in] src    出力データ配列ポインタ
   * @param [in] dLen   出力データサイズ
   */
  virtual bool
  WriteFieldData(cdm_FILE* pFile,
                 cdm_Array* src,
                 size_t dLen);

  /**
   * @brief Field Datat 出力 (vtk用)
   * @param [in] pFile         出力ファイルポインタ
   * @param [in] src           出力データ配列ポインタ
   * @param [in] dLen          出力データサイズ
   * @param [in] d_type        データ型タイプ
   * @param [in] flag_variname VTK形式における変数名出力フラグ
   * @param [in] variname      VTK形式で出力する変数名
   */
  virtual bool
  WriteFieldData(cdm_FILE* pFile,
                 cdm_Array* src,
                 size_t dLen,
                 CDM::E_CDM_DTYPE d_type,
                 bool flag_variname,
                 std::string variname)
  { return true; };

  /**
   * @brief マーカーの書き込み
   * @param[in] dmy   マーカー
   * @param[in] pFile ファイルポインタ
   * @param[in] out   plot3d用Fortran出力フラグ　通常はfalse
   */
  virtual bool
  WriteDataMarker(int dmy, cdm_FILE* pFile, bool out=false) { return true; };


  /**
   * @brief avs の　ヘッダーレコード出力コントロール
   * @param[in] myRank    ランクID
   * @param[in] in_dfi    dfiのポインタ配列
   */
  virtual void
  output_avs(
             int myRank,
             vector<cdm_DFI *>in_dfi){};

  /**
   * @brief avs のヘッダーレコード出力コントロール (不等間隔格子対応版)
   * @param [in] myRank      rankID
   * @param [in] in_dfi      dfiのポインター
   * @param [in] out_domain  出力用のdomainインスタンス
   * @param [in] out_process 出力用のprocessインスタンス
   * @param [in] gc          出力ガイドセル数
   */
  virtual void
  output_avs(int myRank,
             vector<cdm_DFI *>in_dfi,
             cdm_Domain* out_domain,
             cdm_Process* out_process,
             int gc){};

protected :

};

#endif //_CONVOUTPUT_H_
