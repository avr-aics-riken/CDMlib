#ifndef _CONVOUTPUT_FUB_H_
#define _CONVOUTPUT_FUB_H_

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
 * @file   convOutput_FUB.h
 * @brief  convOutput_FUB Class Header
 * @author aics
 * @date   2014/04/11
 */

#include "convOutput.h"

class convOutput_FUB : public convOutput {

public:

  /** コンストラクタ */
  convOutput_FUB();

  /** デストラクタ */
  ~convOutput_FUB();

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
   * @brief sphファイルのheaderの書き込み
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
   * @brief マーカーの書き込み
   * @param[in] dmy   マーカー
   * @param[in] pFile 出力ファイルポインタ
   * @param[in] out   plot3d用
   */
  bool
  WriteDataMarker(int dmy, cdm_FILE* pFile, bool out);

protected:

};

#endif // _CONVOUTPUT_FUB_H_
