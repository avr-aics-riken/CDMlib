#ifndef _CDM_FILEPATH_H_
#define _CDM_FILEPATH_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_FilePath.h
 * @brief  cdm_FilePath Class Header
 * @author aics    
 */

/** index.dfi ファイルの FilePath */
class cdm_FilePath {

public:

  std::string ProcDFIFile;                       ///<proc.dfi ファイル名

  /** コンストラクタ **/
  cdm_FilePath();

  /** 
   * @brief コンストラクタ
   * @param [in] _ProcDFIFile proc.dfiファイル名
   */ 
  cdm_FilePath(const std::string _ProcDFIFile);

  /** デストラクタ **/

  ~cdm_FilePath();

  /**
   * @brief read FilePath(inde.dfi)
   * @details proc.dfiファイル名の読込み
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  /**
   * @brief DFIファイル:Processを出力する
   * @details proc.dfiファイル名の出力
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab); 

};

#endif // _CDM_FILEPATH_H_
