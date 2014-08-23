#ifndef _CONV_DEFINE_H_
#define _CONV_DEFINE_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   conv_Define.h
 * @brief  CONV Definition Header
 * @author aics
 * @date   2013/11/7
 */


#include <stdlib.h>

#define Exit(x) \
((void)printf("exit at %s:%u\n", __FILE__, __LINE__), exit((x)))

#define message() printf("\t%s (%d):\n",__FILE__, __LINE__)
#define mark() printf("%s (%d):\n",__FILE__, __LINE__)

#define stamped_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#define stamped_fprintf fprintf(fp, "%s (%d):  ",__FILE__, __LINE__), fprintf

#define Hostonly_ if(m_paraMngr->GetMyRankID()==0)


#define LOG_OUT_   if(m_lflag)
#define LOG_OUTV_  if(m_lflagv)
#define STD_OUT_   if(m_pflag) 
#define STD_OUTV_  if(m_pflagv) 


/** コンバート形式 */
enum E_CONV_OUTPUT_CONV_TYPE
{
  E_CONV_OUTPUT_UNKNOWN = -1, ///未定義
  E_CONV_OUTPUT_Mx1 = 0,   ///M対1
  E_CONV_OUTPUT_MxN,       ///M対N
  E_CONV_OUTPUT_MxM        ///M対M
};

/** 並列処理時のファイル割振り方法 */
enum E_CONV_OUTPUT_MULTI_FILE_CAST
{
  E_CONV_OUTPUT_CAST_UNKNOWN = -1, //未定義
  E_CONV_OUTPUT_STEP         = 0,  //step基準
  E_CONV_OUTPUT_RANK               //rank基準
};


#define ON          1
#define OFF         0

#define REAL_UNKNOWN 0
#define SPH_FLOAT    1
#define SPH_DOUBLE   2

#define SPH_DATA_UNKNOWN 0
#define SPH_SCALAR       1
#define SPH_VECTOR       2

/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _F_IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (size_t)(_K+_VC-1) * (size_t)(_NI+2*_VC) * (size_t)(_NJ+2*_VC) \
+ (size_t)(_J+_VC-1) * (size_t)(_NI+2*_VC) \
+ (size_t)(_I+_VC-1) \
)

#endif // _CONV_DEFINE_H_
