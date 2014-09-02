#ifndef _CDM_MPI_H_
#define _CDM_MPI_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_MPI.h
 * @brief  cdm_MPI Class Header
 * @author aics    
 */

/** proc.dfi ファイルの MPI */
class cdm_MPI {

public:

  int NumberOfRank;                      ///<プロセス数
  int NumberOfGroup;                     ///<グループ数

  /** コンストラクタ **/
  cdm_MPI();

  /**
   * @brief コンストラクタ
   * @param [in] _NumberOfRank  プロセス数
   * @param [in] _NumberOfGroup グループ数
   */ 
  cdm_MPI(const int _NumberOfRank, int _NumberOfGroup=0);

  /** デストラクタ **/
  ~cdm_MPI();

  /**
   * @brief read MPI(proc.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @param [in]   domain  Domain
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl, 
       const cdm_Domain domain); 

  /**
   * @brief DFIファイル:MPIを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

};

#endif // _CDM_MPI_H_
