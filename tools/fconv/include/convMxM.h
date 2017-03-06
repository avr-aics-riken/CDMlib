#ifndef _CONVMXM_H_
#define _CONVMXM_H_

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

#include "conv.h"
#include "convOutput.h"

/**
 * @file   convMxM.h
 * @brief  convMxM Class Header
 * @author aics
 * @date   2013/11/14
 */

class convMxM : public CONV {

public:

  //vector<stepinfo> m_stepList;   ///<並列処理用インデックスリスト
  //vector<rankinfo> m_rankList;   ///<並列処理用インデックスリスト
  vector<step_rank_info> m_StepRankList; ///<並列処理用インデックスリスト

public:

  /** コンストラクタ */
  convMxM();

  /** デストラクタ */
  ~convMxM();

public:

  /**
   * @brief MxMの実行
   * @return エラーコード
   */
  bool exec();

  /**
   * @brief
   */
  bool mxmsolv(std::string dfiname,
               cdm_DFI* dfi,
               int l_step,
               double l_time,
               int rankID,
               double* min,
               double* max);

};

#endif // _CONVMXM_H_
