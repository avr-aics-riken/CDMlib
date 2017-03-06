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

/**
 * @file   fub_Define.h
 * Fubの定義マクロ記述ヘッダーファイル
 * @date   2016/03/17
 */
#ifndef _FUB_DEFINE_H_
#define _FUB_DEFINE_H_

/** コーナー位置 */
#define _P000 0
#define _P100 1
#define _P110 2
#define _P010 3
#define _P001 4
#define _P101 5
#define _P111 6
#define _P011 7

/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (long long)(_K+_VC) * (long long)(_NI+2*_VC) * (long long)(_NJ+2*_VC) \
+ (long long)(_J+_VC) * (long long)(_NI+2*_VC) \
+ (long long)(_I+_VC) \
)

/** 4次元インデクス(i,j,k,n) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
( (long long)(_N) * (long long)(_NI+2*_VC) * (long long)(_NJ+2*_VC) * (long long)(_NK+2*_VC) \
+ _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
)

/** 3次元インデクス(i,j,k,3) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 */
#define _IDX_V3D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) (_IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC))

/** 4次元インデクス(n,i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NN 成分数
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S4DEX(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
( (long long)(_NN) * _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
+ (long long)(_N) )

/** 3次元インデクス(3,i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 */
#define _IDX_V3DEX(_N,_I,_J,_K,_NI,_NJ,_NK,_VC) (_IDX_S4DEX(_N,_I,_J,_K,3,_NI,_NJ,_NK,_VC))

enum fub_FaceFlag
{
  X_MINUS = 0  ///< -X face
, X_PLUS  = 1  ///< +X face
, Y_MINUS = 2  ///< -Y face
, Y_PLUS  = 3  ///< +Y face
, Z_MINUS = 4  ///< -Z face
, Z_PLUS  = 5  ///< +Z face
};

/** 軸方向フラグ */
enum fub_DirFlag
{
  X_DIR = 0 ///< X direction
, Y_DIR = 1 ///< Y direction
, Z_DIR = 2 ///< Z direction
};

/** 方向フラグ */
enum fub_PMFlag
{
  PLUS2MINUS = 0 ///< plus   -> minus direction
, MINUS2PLUS = 1 ///< minus  -> plus  direction
, BOTH       = 2 ///< plus  <-> minus direction
};

#endif /* _FUB_DEFINE_H_ */
