#ifndef _CONVOUTPUT_AVS_H_
#define _CONVOUTPUT_AVS_H_

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
 * @file   convOutput_AVS.h
 * @brief  convOutput_AVS Class Header
 * @author aics
 * @date   2013/11/7
 */

#include "convOutput.h"

class convOutput_AVS : public convOutput {

public:

  /** コンストラクタ */
  convOutput_AVS();

  /** デストラクタ */
  ~convOutput_AVS();

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
   * @brief Field Datat 出力
   * @param [in] fp     出力ファイルポインタ
   * @param [in] src    出力データ配列ポインタ
   * @param [in] dLen   出力データサイズ
   */
  bool WriteFieldData(FILE* fp,
                      cdm_Array* src,
                      size_t dLen); 

  /**
   * @brief avsファイルのヘッダー処理
   * @param [in] myRank   rankID
   * @param [in] in_dfi   dfiのポインター
   */
  void output_avs(
                  int myRank, 
                  vector<cdm_DFI *>in_dfi);

protected:
  /**
   * @brief avsファイルのヘッダー処理（Mx1)
   * @param [in] myRank rankID
   * @param [in] in_dfi dfiのポインター
   */
/*
  void output_avs_Mx1(
                      int myRank, 
                      vector<cdm_DFI *>in_dfi); 
*/
  /**
   * @brief avsファイルのヘッダー処理（MxM)
   * @param [in] myRank rankID
   * @param [in] in_dfi dfiのポインター
   */
  void output_avs_MxM(
                      int myRank, 
                      vector<cdm_DFI *>in_dfi); 

  /**
   * @brief avsファイルのヘッダー処理（MxN)
   * @param [in] myRank   rankID
   * @param [in] in_dfi   dfiのポインター
   * @param [in] paraMngr パラマネージャー
   * @param [in] head     headインデックス
   */
  void output_avs_MxN(
                      int myRank, 
                      vector<cdm_DFI *>in_dfi,
                      cpm_ParaManager* paraMngr,
                      int *head); 

  /**
   * @brief avs coord data ファイル出力
   * @param[in] RankID  ランクID
   * @param[in] mio     分割出力指示
   * @param[in] min_ext 座標値の最小値
   * @param[in] max_ext 座標値の最大値
   */
  void output_avs_coord(int RankID,
                        bool mio,
                        double min_ext[3],
                        double max_ext[3]);

  /**
   * @brief avsファイルヘッダー出力
   * @param[in] dfi     cdm_DFIクラスポインタ
   * @param[in] RankID  ランクID
   * @param[in] mio     分割出力指示
   * @param[in] ndim    3
   * @param[in] nspace  3
   * @param[in] dims    サイズ
   */
  void output_avs_header(cdm_DFI* dfi,
                         int RankID,
                         bool mio,
                         int ndim,
                         int nspace,
                         int dims[3]);  

};

#endif // _CONVOUTPUT_AVS_H_ 
