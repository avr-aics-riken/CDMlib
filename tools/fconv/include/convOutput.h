#ifndef _CONVOUTPUT_H_
#define _CONVOUTPUT_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   convOutput.h
 * @brief  convOutput Class Header
 * @author aics
 * @date   2013/11/7
 */

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
  virtual FILE*
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
  OutputFile_Close(FILE* fp)
  { fclose(fp); };

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
   * @brief ファイルのheaderの書き込み
   * @param[in] step     ステップ数
   * @param[in] dim      成分数
   * @param[in] d_type   データ型タイプ
   * @param[in] imax     x方向ボクセルサイズ
   * @param[in] jmax     y方向ボクセルサイズ
   * @param[in] kmax     z方向ボクセルサイズ
   * @param[in] time     時間
   * @param[in] org      原点座標
   * @param[in] pit      ピッチ
   * @param[in] prefix   ファイル接頭文字
   * @param[in] fp       ファイルポインタ
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
                    FILE *fp) 
  { return true; };

  /**
   * @brief Field Datat 出力
   * @param [in] fp     出力ファイルポインタ
   * @param [in] src    出力データ配列ポインタ
   * @param [in] dLen   出力データサイズ
   */
  virtual bool
  WriteFieldData(FILE* fp, 
                 cdm_Array* src, 
                 size_t dLen);

  /**
   * @brief マーカーの書き込み
   * @param[in] dmy マーカー
   * @param[in] fp  ファイルポインタ
   * @param[in] out plot3d用Fortran出力フラグ　通常はfalse
   */
  virtual bool
  WriteDataMarker(int dmy, FILE* fp, bool out=false) { return true; }; 


  /**
   * @brief avs の　ヘッダーレコード出力コントロール
   * @param[in] myRank    ランクID
   * @param[in] in_dfi    dfiのポインタ配列
   */
  virtual void
  output_avs( 
             int myRank,
             vector<cdm_DFI *>in_dfi){};

protected :

};

#endif //_CONVOUTPUT_H_
