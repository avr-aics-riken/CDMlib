#ifndef _CIO_FILEPATH_H_
#define _CIO_FILEPATH_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_FilePath.h
 * @brief  cio_FilePath Class Header
 * @author aics    
 */

/** index.dfi ファイルの FilePath */
class cio_FilePath {

public:

  std::string ProcDFIFile;                       ///<proc.dfi ファイル名

  /** コンストラクタ **/
  cio_FilePath();

  /** 
   * @brief コンストラクタ
   * @param [in] _ProcDFIFile proc.dfiファイル名
   */ 
  cio_FilePath(const std::string _ProcDFIFile);

  /** デストラクタ **/

  ~cio_FilePath();

  /**
   * @brief read FilePath(inde.dfi)
   * @details proc.dfiファイル名の読込み
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Read(cio_TextParser tpCntl);

  /**
   * @brief DFIファイル:Processを出力する
   * @details proc.dfiファイル名の出力
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab); 

};

#endif // _CIO_FILEPATH_H_
