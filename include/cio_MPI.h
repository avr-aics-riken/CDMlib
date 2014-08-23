#ifndef _CIO_MPI_H_
#define _CIO_MPI_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_MPI.h
 * @brief  cio_MPI Class Header
 * @author aics    
 */

/** proc.dfi ファイルの MPI */
class cio_MPI {

public:

  int NumberOfRank;                      ///<プロセス数
  int NumberOfGroup;                     ///<グループ数

  /** コンストラクタ **/
  cio_MPI();

  /**
   * @brief コンストラクタ
   * @param [in] _NumberOfRank  プロセス数
   * @param [in] _NumberOfGroup グループ数
   */ 
  cio_MPI(const int _NumberOfRank, int _NumberOfGroup=0);

  /** デストラクタ **/
  ~cio_MPI();

  /**
   * @brief read MPI(proc.dfi)
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [in]   domain  Domain
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Read(cio_TextParser tpCntl, 
       const cio_Domain domain); 

  /**
   * @brief DFIファイル:MPIを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

};

#endif // _CIO_MPI_H_
