#ifndef _CONVOUTPUT_VTK_H_
#define _CONVOUTPUT_VTK_H_

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
  FILE* OutputFile_Open(
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
   * @param[in] fp      出力ファイルポインタ
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
                    FILE *fp); 

  /**
   * @brief Field Datat 出力
   * @param [in] fp     出力ファイルポインタ
   * @param [in] src    出力データ配列ポインタ
   * @param [in] dLen   出力データサイズ
   */
  bool
  WriteFieldData(FILE* fp,
                 cdm_Array* src,
                 size_t dLen);

  /**
   * @brief マーカーの書き込み
   * @param[in] dmy マーカー
   * @param[in] fp  出力ファイルポインタ 
   * @param[in] out 出力フラグ
   */
  bool
  WriteDataMarker(int dmy, FILE* fp, bool out);

  /**
   * @brief 出力ファイルをクローズする
   * @param [in] fp ファイルポインタ
   */
  void OutputFile_Close(FILE* fp); 

protected:

};

#endif // _CONVOUTPUT_VTK_H_ 
