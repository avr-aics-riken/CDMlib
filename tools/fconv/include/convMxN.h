#ifndef _CONVMXN_H_
#define _CONVMXN_H_

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
 * @file   convMxN.h
 * @brief  convMxN Class Header
 * @author aics
 * @date   2013/11/14
 */

class convMxN : public CONV {

public:
  int m_Gvoxel[3];
  int m_Gdiv[3];
  int m_Head[3];
  int m_Tail[3];

  vector<cdm_DFI *>m_out_dfi; //出力DFIのポインタ

public:

  /** コンストラクタ */
  convMxN();

  /** デストラクタ */
  ~convMxN();

public:

  /**
   * @brief 領域分割と出力DFIのインスタンス
   */
  void VoxelInit();

  /**
   * @brief MxNの実行
   * @return エラーコード
   */
  bool exec();

};

#endif // _CONVEXEC_H_
