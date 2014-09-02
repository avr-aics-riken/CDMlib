#ifndef _CDM_DFI_BOV_H_
#define _CDM_DFI_BOV_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_BOV.h
 * @brief  cdm_DFI_BOV Class Header
 * @author aics    
 */

#include "cdm_DFI.h"

class cdm_DFI_BOV : public cdm_DFI {

public:
  /** コンストラクタ */
  cdm_DFI_BOV();

  /** 
   * @brief コンストラクタ 
   * @param [in] F_Info  FileInfo
   * @param [in] F_Path  FilePath
   * @param [in] unit    Unit
   * @param [in] domain  Domain
   * @param [in] mpi     MPI
   * @param [in] TSlice  TimeSlice
   * @param [in] process Process
   */
  cdm_DFI_BOV(const cdm_FileInfo F_Info, 
              const cdm_FilePath F_Path, 
              const cdm_Unit unit, 
              const cdm_Domain domain, 
              const cdm_MPI mpi,
              const cdm_TimeSlice TSlice, 
              const cdm_Process process)
  {
    DFI_Finfo      = F_Info;
    DFI_Fpath      = F_Path;
    DFI_Unit       = unit;
    DFI_Domain     = domain;
    DFI_MPI        = mpi;
    DFI_TimeSlice  = TSlice;
    DFI_Process    = process;
    m_bgrid_interp_flag = false;
  };

  /**　デストラクタ */
  ~cdm_DFI_BOV();

public:

protected:

  /**
   * @brief bovファイルのヘッダーレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian エンディアンチェックフラグ true:合致
   * @param[in]  step        ステップ番号
   * @param[in]  head        dfiのHeadIndex
   * @param[in]  tail        dfiのTailIndex
   * @param[in]  gc          dfiのガイドセル数
   * @param[out] voxsize     voxsize
   * @param[out] time        時刻
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  read_HeaderRecord(FILE* fp,
                    bool matchEndian,
                    unsigned step,
                    const int head[3],
                    const int tail[3],
                    int gc,
                    int voxsize[3],
                    double &time);

  /**
   * @brief フィールドデータファイルのデータレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  buf         読込み用バッファ
   * @param[in]  head        読込みバッファHeadIndex
   * @param[in]  nz          z方向のボクセルサイズ（実セル＋ガイドセル＊２）
   * @param[out] src         読み込んだデータを格納した配列のポインタ
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  read_Datarecord(FILE* fp,
                  bool matchEndian,
                  cdm_Array* buf,
                  int head[3],
                  int nz,
                  cdm_Array* &src);

  /**
   * @brief bovファイルのAverageデータレコードの読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        読込みstep番号
   * @param[out] avr_step    平均ステップ
   * @param[out] avr_time    平均タイム
   * @return errorcode
   */
  CDM::E_CDM_ERRORCODE
  read_averaged(FILE* fp,
                bool matchEndian,
                unsigned step,
                unsigned &avr_step,
                double &avr_time);

  /**
   * @brief bovヘッダファイルの出力
   * @param[in] fp     ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  write_HeaderRecord(FILE* fp,
                     const unsigned step,
                     const double time,
                     const int RankID);

  /**
   * @brief bovデータ出力 
   * @param[in]  fp     ファイルポインタ
   * @param[in]  val    データポインタ
   * @param[in]  gc     仮想セル数  
   * @param[in]  RankID ランク番号
   * @return error code 
   */
  CDM::E_CDM_ERRORCODE
  write_DataRecord(FILE* fp, 
                   cdm_Array* val, 
                   const int gc, 
                   const int RankID);

  /**
   * @brief Averageレコードの出力
   * @param[in] fp       ファイルポインタ
   * @param[in] step_avr 平均ステップ番号
   * @param[in] time_avr 平均時刻
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  write_averaged(FILE* fp,
                 const unsigned step_avr,
                 const double time_avr);

 
  /**
   * @brief ヘッダーデータファイルの出力
   * @param [in] step step番号
   * @param [in] time time
   */
  bool
  write_ascii_header(const unsigned step,
                     const double time);
  
};

#endif // _cdm_DFI_BOV_H_
