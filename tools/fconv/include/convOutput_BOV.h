#ifndef _CONVOUTPUT_BOV_H_
#define _CONVOUTPUT_BOV_H_

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
 * @file   convOutput_BOV.h
 * @brief  convOutput_BOV Class Header
 * @author aics
 * @date   2013/11/7
 */

#include "convOutput.h"

class convOutput_BOV : public convOutput {

public:

  /** コンストラクタ */
  convOutput_BOV();

  /** デストラクタ */
  ~convOutput_BOV();

public:


  /**
   * @brief 出力ファイルをオープンする
   * @param [in] prefix ファイル接頭文字
   * @param [in] step ステップ数
   * @param [in] id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  FILE* OutputFile_Open(
                        const std::string prefix,
                        const unsigned step,
                        const int id,
                        const bool mio);

  /**
   * @brief bovファイルのheaderの書き込み
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
   * @param[in] fp       ファイルポインタ
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


protected:

};

#endif // _CONVOUTPUT_BOV_H_ 
