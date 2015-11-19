#ifndef _CONVOUTPUT_PLOT3D_H_
#define _CONVOUTPUT_PLOT3D_H_

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
 * @file   convOutput_PLOT3D.h
 * @brief  convOutput_PLOT3D Class Header
 * @author aics
 * @date   2013/11/7
 */

#include "convOutput.h"
#include <typeinfo>

class convOutput_PLOT3D : public convOutput {

public:

public:

  /** コンストラクタ */
  convOutput_PLOT3D();

  /** デストラクタ */
  ~convOutput_PLOT3D();

public:

  /**
   * @brief GRID ファイル出力
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   step番号
   * @param [in] myRank ランク番号
   * @param [in] dType  dfiのデータタイプ
   * @param [in] guide  ガイドセル数
   * @param [in] org    原点座標値
   * @param [in] pit    ピッチ
   * @param [in] sz     ボクセルサイズ
   */
  void WriteGridData(std::string prefix,
                     int step,
                     int myRank,
                     int dType,
                     int guide,
                     double org[3],
                     double pit[3],
                     int sz[3]); 

  /**
   * @brief GRID ファイル出力 (不等間隔格子対応版)
   * @param [in] prefix      ファイル接頭文字
   * @param [in] step        step番号
   * @param [in] myRank      ランク番号
   * @param [in] dType       dfiのデータタイプ
   * @param [in] guide       ガイドセル数
   * @param [in] out_domain  出力用のdomainインスタンス
   * @param [in] out_process 出力用のprocessインスタンス
   */
  void WriteGridData(std::string prefix,
                     int step,
                     int myRank,
                     int dType,
                     int guide,
                     cdm_Domain* out_domain,
                     cdm_Process* out_process);

  /**
   * @brief xyzファイルの出力(template 関数)
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   ステップ
   * @param [in] rank   ランク
   * @param [in] guide  ガイドセル数
   * @param [in] origin 基点座標
   * @param [in] pitch  ピッチ
   * @param [in] size   ボクセルサイズ
   * @param [in] x      x方向座標ワーク
   * @param [in] y      y方向座標ワーク
   * @param [in] z      z方向座標ワーク
   */
  template<class T1, class T2>
  void OutputPlot3D_xyz(std::string prefix,
                        int step, 
                        int rank, 
                        int guide, 
                        T1* origin, 
                        T1* pitch, 
                        int* size, 
                        T2* x, 
                        T2* y, 
                        T2* z);

  /**
   * @brief xyzファイルの出力(不等間隔格子対応版)
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   ステップ
   * @param [in] rank   ランク
   * @param [in] guide  ガイドセル数
   * @param [in] out_domain  出力用のdomainインスタンス
   * @param [in] out_process 出力用のprocessインスタンス
   */
  template<class T>
  void OutputPlot3D_xyz(std::string prefix,
                        int step,
                        int rank,
                        int guide,
                        cdm_Domain* out_domain,
                        cdm_Process* out_process);

  /**
   * @brief グリッド数の書き出し
   * @param [in] pFile 出力ファイルポインタ
   * @param [in] ngrid グリッド数
   */
  void WriteNgrid(cdm_FILE* pFile, int ngrid);

  /**
   * @brief ブロックデータの書き出し
   * @param [in] pFile 出力ファイルポインタ
   * @param [in] id i方向のサイズ
   * @param [in] jd j方向のサイズ
   * @param [in] kd k方向のサイズ
   */
  void WriteBlockData(cdm_FILE* pFile,
                      int id,
                      int jd,
                      int kd); 

  /**
   * @brief gridデータ出力
   * @param [in] pFile      出力ファイルポインタ
   * @param [in] id         i方向サイズ
   * @param [in] jd         j方向のサイズ
   * @param [in] kd         k方向のサイズ
   * @param [in] ngrid      1
   * @param [in] x          x座標値
   * @param [in] y          y座標値
   * @param [in] z          z座標値
   */ 
  template<class T>
  bool
  WriteXYZData(cdm_FILE* pFile,
               int id, 
               int jd, 
               int kd,
               int ngrid,
               T*  x,
               T*  y,
               T*  z);

  /**
   * @brief gridデータ出力 (不等間隔格子対応版)
   * @param [in] pFile   出力ファイルポインタ
   * @param [in] sz3d    データサイズ
   * @param [in] x       x座標値データポインタ
   * @param [in] y       y座標値データポインタ
   * @param [in] z       z座標値データポインタ
   * @param [in] iblank  iblankデータポインタ
   */ 
  template<class T>
  bool
  WriteXYZData(cdm_FILE* pFile,
               size_t sz3d,
               T* x,
               T* y,
               T* z,
               int* iblank);

  /**
   * @brief Formatted 出力
   * @param [in] pFile 出力ファイルポインタ
   * @param [in] id  i方向サイズ
   * @param [in] jd  j方向のサイズ
   * @param [in] kd  k方向のサイズ
   * @param [in] x   出力座標値配列
   */
  template<class T>
  void
  WriteXYZ_FORMATTED(cdm_FILE *pFile,
                     int id,
                     int jd,
                     int kd,
                     T* x); 

  /**
   * @brief Formatted 出力 (IBLANK対応版)
   * @param [in] pFile 出力ファイルポインタ
   * @param [in] sz3d データサイズ
   * @param [in] tmp  出力データポインタ
   */
  template<class T>
  void
  WriteXYZ_FORMATTED(cdm_FILE *pFile,
                     size_t sz3d,
                     T* tmp); 

  /**
   * @brief 出力ファイルをオープンする
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   ステップ数
   * @param [in] id     ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */ 
  cdm_FILE*
  OutputFile_Open(
                  const std::string prefix,
                  const unsigned step,
                  const int id,
                  const bool mio);

  /**
   * @brief funcデータファイルののheader部の書き込み 
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
   * @param [in] pFile  出力ファイルポインタ
   * @param [in] src    出力データ配列ポインタ
   * @param [in] dLen   出力データサイズ
   */
  bool
  WriteFieldData(cdm_FILE* pFile, 
                 cdm_Array* src, 
                 size_t dLen);

  /**
   * @brief Function ブロックデータの書き出し 
   * @param [in] pFile 出力ファイルポインタ
   * @param [in] id i方向のサイズ
   * @param [in] jd j方向のサイズ
   * @param [in] kd k方向のサイズ
   * @param [in] nvar 出力項目数
   */
  void 
  WriteFuncBlockData(cdm_FILE* pFile,
                     int id,
                     int jd,
                     int kd,
                     int nvar); 

  /**
   * @brief func data の出力
   * @param [in] pFile 出力ファイルポインタ
   * @param [in] p3src plot3d func データ配列ポインタ
   */
  void
  WriteFuncData(cdm_FILE* pFile, cdm_Array* p3src);

  /**
   * @brief マーカーの書き込み
   * @param[in] dmy    マーカー
   * @param[in] pFile  出力ファイルポインタ 
   * @param[in] out    Fortranマーカー出力フラグ
   */
  bool
  WriteDataMarker(int dmy, cdm_FILE* pFile, bool out);

};

//inline 関数
#include "inline/conv_plot3d_inline.h"

#endif // _CONVOUTPUT_PLOT3D_H_ 
