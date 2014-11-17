#ifndef _CDM_VISIT_H_
#define _CDM_VISIT_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_VisIt.h
 * @brief  cdm_VisIt Class Header
 * @author aics    
 */

class cdm_VisIt {

public:
  
  std::string PlotGC;           ///< ガイドセル描画オプション

  /** コンストラクタ **/
  cdm_VisIt();

  
  /**
   * @brief コンストラクタ
   * @param [in] _PlotGC ガイドセル描画オプション
   */
  cdm_VisIt(const std::string _PlotGC);
  

  /** デストラクタ **/
  ~cdm_VisIt();

  /**
   * @brief read VisItオプションの読み込み
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  /**
   * @brief DFIファイル:VisItオプションを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab); 

};

#endif // _CDM_VISIT_H_
