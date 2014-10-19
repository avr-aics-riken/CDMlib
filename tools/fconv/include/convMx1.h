#ifndef _CONVMX1_H_
#define _CONVMX1_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "conv.h"
#include "convOutput.h"

/**
 * @file   convMx1.h
 * @brief  convMx1 Class Header
 * @author aics
 * @date   2013/11/14
 */

class convMx1 : public CONV {

public:

  convOutput *ConvOut;

  typedef std::map<int,int> headT;

  //vector<stepinfo> m_stepList;   ///<並列処理用インデックスリスト
  vector<step_rank_info> m_StepRankList;   ///<並列処理用インデックスリスト

public:

  /** コンストラクタ */
  convMx1();

  /** デストラクタ */
  ~convMx1();

public:

  /**
   * @brief Mx1の実行
   * @return エラーコード
   */ 
  bool exec();

  /**
   * @brief 並列形状 nijkをnijkでコンバートして出力
   * @param[in] fp     出力ファイルポインタ
   * @param[in] inPath dfiのディレクトリパス
   * @param[in] l_step 出力step番号
   * @param[in] l_dtime 出力時刻
   * @param[in] d_type データタイプ
   * @param[in] mio    分割出力指示
   * @param[in] div    分割数
   * @param[in] sz     サイズ
   * @param[in] dfi    dfi
   * @param[in] DFI_Process cdm_Process
   * @param[in] mapHeadX 
   * @param[in] mapHeadY
   * @param[in] mapHeadZ
   * @param[out] min 最小値
   * @param[out] max 最大値
   */
  bool
  convMx1_out_nijk(FILE* fp,
                    std::string inPath,
                    int l_step,
                    double l_dtime,
                    CDM::E_CDM_DTYPE d_type,
                    bool mio,
                    int div[3],
                    int sz[3],
                    cdm_DFI* dfi,
                    cdm_Process* DFI_Process,
                    headT mapHeadX, 
                    headT mapHeadY, 
                    headT mapHeadZ,
                    double* min,
                    double* max);

  /**
   * @brief 並列形状 nijkをijknまたは ijknをijknにコンバートして出力
   * @param[in] fp     出力ファイルポインタ
   * @param[in] inPath dfiのディレクトリパス
   * @param[in] l_step 出力step番号
   * @param[in] l_dtime 出力時刻
   * @param[in] d_type データタイプ
   * @param[in] mio    分割出力指示
   * @param[in] div    分割数
   * @param[in] sz     サイズ
   * @param[in] dfi    dfi
   * @param[in] DFI_Process cdm_Process
   * @param[in] mapHeadX
   * @param[in] mapHeadY
   * @param[in] mapHeadZ
   * @param[out] min 最小値
   * @param[out] max 最大値
   */ 
  bool
  convMx1_out_ijkn(FILE* fp,
                    std::string inPath,
                    int l_step,
                    double l_dtime,
                    CDM::E_CDM_DTYPE d_type,
                    bool mio,
                    int div[3],
                    int sz[3],
                    cdm_DFI* dfi,
                    cdm_Process* DFI_Process,
                    headT mapHeadX,
                    headT mapHeadY,
                    headT mapHeadZ,
                    double* min,
                    double* max);

  /**
   * @brief 補間処理
   * @param[in]  src_old  1つ前の層
   * @param[in]  src      処理する層
   * @param[out] outArray 足しこむ配列
   * @param[in]  ivar_src 図心データの コンポーネント位置
   * @param[in]  ivar_out 格子データの コンポーネント位置
   */
  bool InterPolate(cdm_Array* src_old, cdm_Array* src, cdm_Array* outArray,
                   int ivar_src, int ivar_out); 

  /**
   * @brief 配列のゼロクリア
   * @param[out] data     配列
   * @param[in]  ivar_out コンポーネント位置
   */
  template<class T>
  void zeroClearArray(cdm_TypeArray<T>* data, int ivar_out);

  /**
   * @brief 図心データを格子点に補間
   * @param[out] O 格子点データ
   * @param[in]  S 図心データ
   * @param[in]  ivar_out 格子データの コンポーネント位置
   * @param[in]  ivar_src 図心データの コンポーネント位置
   */
  template<class T>
  bool setGridData_XY(cdm_TypeArray<T>* O,
                      cdm_TypeArray<T>* S,
                      int ivar_out,
                      int ivar_src); 

  /**
   * @brief 内部の格子点のデータを重み付けでで割る
   * @param[out] O        格子点データ
   * @param[in]  ivar_out コンポーネント位置
   */
  template<class T>
  void VolumeDataDivide8(cdm_TypeArray<T> *O, int ivar_out); 

  /**
   * @brief NIJK配列をスカラーのIJK配列にコピーコントロール
   * @param[in] src  コピー元配列
   * @param[in] ivar コピーするコンポーネント位置
   * @return IJKにコピーされて配列ポインタ
   */
  cdm_Array* nijk_to_ijk(cdm_Array* src, int ivar); 

  /**
   * @brief NIJK配列をスカラーのIJK配列にコピー 
   * @param[in] S コピー元配列
   * @param[in] O コピー先配列
   * @param[in] ivar コピーするコンポーネント位置
   */
  template<class T>
  void copyArray_nijk_ijk(cdm_TypeArray<T> *S, cdm_TypeArray<T> *O, int ivar);

};

//inline関数
#include "inline/convMx1_inline.h"

#endif // _CONVMX1_H_
