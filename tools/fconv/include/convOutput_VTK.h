#ifndef _CONVOUTPUT_VTK_H_
#define _CONVOUTPUT_VTK_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   convOutput_VTK.h
 * @brief  convOutput_VTK Class Header
 * @author aics
 * @date   2013/11/12
 */

#include "convOutput.h"

class convOutput_VTK : public convOutput {

public:

  /** コンストラクタ */
  convOutput_VTK();

  /** デストラクタ */
  ~convOutput_VTK();

public:


  /**
   * @brief 出力ファイルをオープンする
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   ステップ数
   * @param [in] id     ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  cdm_FILE* OutputFile_Open(
                        const std::string prefix,
                        const unsigned step,
                        const int id,
                        const bool mio);
                  

  /**
   * @brief vtkファイルのheaderの書き込み
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
   * @brief vtkファイルのheaderの書き込み(不等間隔格子対応版)
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
  bool
  WriteHeaderRecord(int step,
                    double time,
                    int dim,
                    CDM::E_CDM_DTYPE d_type,
                    CDM::E_CDM_DFITYPE dfi_type,
                    cdm_Domain* out_domain,
                    cdm_Process* out_process,
                    int gc,
                    const std::string prefix,
                    cdm_FILE *pFile); 

  /**
   * @brief Field Datat 出力
   * @param [in] pFile  出力ファイルポインタ
   * @param [in] src    出力データ配列ポインタ
   * @param [in] dLen   出力データサイズ
   */
  bool
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
  bool
  WriteFieldData(cdm_FILE* pFile,
                 cdm_Array* src,
                 size_t dLen,
                 CDM::E_CDM_DTYPE d_type,
                 bool flag_variname,
                 std::string variname);

  /**
   * @brief マーカーの書き込み
   * @param[in] dmy    マーカー
   * @param[in] pFile  出力ファイルポインタ 
   * @param[in] out    出力フラグ
   */
  bool
  WriteDataMarker(int dmy, cdm_FILE* pFile, bool out);

  /**
   * @brief 出力ファイルをクローズする
   * @param [in] pFile ファイルポインタ
   */
  void OutputFile_Close(cdm_FILE* pFile); 

protected:

};

#endif // _CONVOUTPUT_VTK_H_ 
