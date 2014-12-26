#ifndef _CDM_ARRAY_H_
#define _CDM_ARRAY_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "cdm_Define.h"

extern "C"
{
  void cdm_interp_ijkn_r4_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *nvari, float *src, float *dst );
  void cdm_interp_ijkn_r8_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *nvari, double *src, double *dst );
  void cdm_interp_nijk_r4_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *nvari, float *src, float *dst );
  void cdm_interp_nijk_r8_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *nvari, double *src, double *dst );
}

class cdm_Array
{
/*
 * enum
 */
public:
/*
  /// データタイプ
  enum E_CDM_DTYPE
  {
    E_CDM_DTYPE_UNKNOWN=0,
    E_CDM_INT8,
    E_CDM_INT16,
    E_CDM_INT32,
    E_CDM_INT64,
    E_CDM_UINT8,
    E_CDM_UINT16,
    E_CDM_UINT32,
    E_CDM_UINT64,
    E_CDM_FLOAT32,
    E_CDM_FLOAT64
  };

  /// 配列形状タイプ
  enum E_CDM_ARRAYSHAPE
  {
    E_CDM_ARRAYSHAPE_UNKNOWN=0,
    E_CDM_IJKN,
    E_CDM_NIJK,
  };
*/

/*
 *  メンバ関数
 */
public:
  /// デストラクタ
  virtual ~cdm_Array()
  {
  }

  /// インスタンス
  static cdm_Array*
  instanceArray( CDM::E_CDM_DTYPE dtype
               , CDM::E_CDM_ARRAYSHAPE shape
               , size_t ix
               , size_t jx
               , size_t kx
               , size_t gc
               , size_t nvari=1 );

  /// インスタンス
  static cdm_Array*
  instanceArray( CDM::E_CDM_DTYPE dtype
               , CDM::E_CDM_ARRAYSHAPE shape
               , size_t sz[3]
               , size_t gc
               , size_t nvari=1 );

  /// インスタンス
  static cdm_Array*
  instanceArray( CDM::E_CDM_DTYPE dtype
               , CDM::E_CDM_ARRAYSHAPE shape
               , int ix
               , int jx
               , int kx
               , int gc
               , int nvari=1 );

  /// インスタンス
  static cdm_Array*
  instanceArray( CDM::E_CDM_DTYPE dtype
               , CDM::E_CDM_ARRAYSHAPE shape
               , int sz[3]
               , int gc
               , int nvari=1 );

  /// インスタンス
  template<class T>
  static cdm_Array*
  instanceArray( T *data
               , CDM::E_CDM_ARRAYSHAPE shape
               , size_t ix
               , size_t jx
               , size_t kx
               , size_t gc
               , size_t nvari=1 );

  /// インスタンス
  template<class T>
  static cdm_Array*
  instanceArray( T *data
               , CDM::E_CDM_ARRAYSHAPE shape
               , size_t sz[3]
               , size_t gc
               , size_t nvari=1 );

  /// インスタンス
  template<class T>
  static cdm_Array*
  instanceArray( T *data
               , CDM::E_CDM_ARRAYSHAPE shape
               , int ix
               , int jx
               , int kx
               , int gc
               , int nvari=1 );

  /// インスタンス
  template<class T>
  static cdm_Array*
  instanceArray( T *data
               , CDM::E_CDM_ARRAYSHAPE shape
               , int sz[3]
               , int gc
               , int nvari=1 );

  /// データポインタを取得
  void* getData( bool extract=false );

  /// データタイプの取得
  CDM::E_CDM_DTYPE getDataType() const
  {
    return m_dtype;
  }

  /// データタイプ文字列の取得
  char* getDataTypeString() const
  {
    switch( m_dtype )
    {
    case CDM::E_CDM_INT8:
      return "INT8";
      break;
    case CDM::E_CDM_INT16:
      return "INT16";
      break;
    case CDM::E_CDM_INT32:
      return "INT32";
      break;
    case CDM::E_CDM_INT64:
      return "INT64";
      break;
    case CDM::E_CDM_UINT8:
      return "UINT8";
      break;
    case CDM::E_CDM_UINT16:
      return "UINT16";
      break;
    case CDM::E_CDM_UINT32:
      return "UINT32";
      break;
    case CDM::E_CDM_UINT64:
      return "UINT64";
      break;
    case CDM::E_CDM_FLOAT32:
      return "FLOAT32";
      break;
    case CDM::E_CDM_FLOAT64:
      return "FLOAT64";
      break;
    }
    return "Unknown";
  }

  /// 配列形状の取得
  CDM::E_CDM_ARRAYSHAPE getArrayShape() const
  {
    return m_shape;
  }

  /// 配列形状文字列の取得
  char* getArrayShapeString() const
  {
    switch(m_shape)
    {
    case CDM::E_CDM_IJKN:
      return "IJKN";
      break;
    case CDM::E_CDM_NIJK:
      return "NIJK";
      break;
    }
    return "Unknown";
  }

  /// ガイドセル数を取得
  size_t getGc() const
  {
    return m_gc;
  }

  /// ガイドセル数を取得(int版)
  int getGcInt() const
  {
    return m_gcI;
  }

  /// 変数の個数を取得
  size_t getNvari() const
  {
    return m_nvari;
  }

  /// 変数の個数を取得(int版)
  int getNvariInt() const
  {
    return m_nvariI;
  }

  /// 格子数を取得
  const size_t* getArraySize()
  {
    return m_sz;
  }

  /// 格子数を取得(int版)
  const int* getArraySizeInt()
  {
    return m_szI;
  }

  /// headインデクスを取得
  const int* getHeadIndex()
  {
    switch(m_shape)
    {
    case CDM::E_CDM_IJKN:
      return m_headIndex;
      break;
    case CDM::E_CDM_NIJK:
      return m_headIndex + 1;
      break;
    }
    return NULL;
  }

  /// tailインデクスを取得
  const int* getTailIndex()
  {
    switch(m_shape)
    {
    case CDM::E_CDM_IJKN:
      return m_tailIndex;
      break;
    case CDM::E_CDM_NIJK:
      return m_tailIndex + 1;
      break;
    }
    return NULL;
  }

  /// ガイドセルを含んだ格子数を取得
  const size_t* _getArraySize()
  {
    switch(m_shape)
    {
    case CDM::E_CDM_IJKN:
      return m_Sz;
      break;
    case CDM::E_CDM_NIJK:
      return m_Sz + 1;
      break;
    }
    return NULL;
  }

  /// ガイドセルを含んだ格子数を取得(int版)
  const int* _getArraySizeInt()
  {
    switch(m_shape)
    {
    case CDM::E_CDM_IJKN:
      return m_SzI;
      break;
    case CDM::E_CDM_NIJK:
      return m_SzI + 1;
      break;
    }
    return NULL;
  }

  /// 配列長を取得
  size_t getArrayLength() const
  {
    size_t nw = 1;
    for( int i=0;i<4;i++ )
    {
      nw *= m_Sz[i];
    }
    return nw;
  }

  /// head/tailをセット
  void setHeadIndex( int head[3] )
  {
    switch(m_shape)
    {
    case CDM::E_CDM_IJKN:
      m_headIndex[0] = head[0];
      m_headIndex[1] = head[1];
      m_headIndex[2] = head[2];
      m_headIndex[3] = 0;
      m_tailIndex[0] = m_headIndex[0] + m_sz[0] - 1;
      m_tailIndex[1] = m_headIndex[1] + m_sz[1] - 1;
      m_tailIndex[2] = m_headIndex[2] + m_sz[2] - 1;
      m_tailIndex[3] = 0;
      break;
    case CDM::E_CDM_NIJK:
      m_headIndex[0] = 0;
      m_headIndex[1] = head[0];
      m_headIndex[2] = head[1];
      m_headIndex[3] = head[2];
      m_tailIndex[0] = 0;
      m_tailIndex[1] = m_headIndex[1] + m_sz[0] - 1;
      m_tailIndex[2] = m_headIndex[2] + m_sz[1] - 1;
      m_tailIndex[3] = m_headIndex[3] + m_sz[2] - 1;
    }
  }

  /// 配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArray( cdm_Array *dst, bool ignoreGc=false ) = 0;

  /// 範囲指定での配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArray( int sta[3], int end[3], cdm_Array *dst ) = 0;

//FCONV 20131216.s
  /// 指定変数の配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArrayNvari( cdm_Array *dst, int vari, bool ignoreGc=false ) = 0;

  /// 指定変数の範囲指定での配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArrayNvari( int sta[3], int end[3], cdm_Array *dst, int vari ) = 0;
//FCONV 20131216.e

  /// 指定した変数の配列のみ取得し、IJK配列でdstにコピー
  virtual int copyArrayNvari_to_ijk( cdm_Array *dst, int vari, bool ignoreGc=false ) = 0;

  /// 指定した変数の配列のみ範囲指定で取得し、IJK配列でdstにコピー
  virtual int copyArrayNvari_to_ijk( int sta[3], int end[3], cdm_Array *dst, int vari ) = 0;

  /// 粗密データの補間処理を行う
  static cdm_Array* interp_coarse( cdm_Array *src, int &err, bool head0start=true );

  /// 配列サイズ分のバイナリデータを読み込み(戻り値は読み込んだ要素数)
  virtual size_t readBinary( FILE *fp , bool bMatchEndian) = 0;

  /// 配列サイズ分のバイナリデータを書き出す(戻り値は読み込んだ要素数)
  virtual size_t writeBinary( FILE *fp ) = 0;

//20131213.s
  /// 配列サイズ分のasciiデータを書き出す(戻り値は読み込んだ要素数)
  virtual size_t writeAscii( FILE *fp ) = 0;
//20131213.e

protected:
  /// デフォルトコンストラクタ
  cdm_Array()
  {
    m_dtype = CDM::E_CDM_DTYPE_UNKNOWN;
    m_shape = CDM::E_CDM_ARRAYSHAPE_UNKNOWN;
    m_sz[0] = m_sz[1] = m_sz[2] = 0;
    m_Sz[0] = m_Sz[1] = m_Sz[2] = m_Sz[3] = 0;
    m_gc = 0;
    m_gcl[0] = m_gcl[1] = m_gcl[2] = m_gcl[3] = 0;
    m_nvari = 1;
    m_headIndex[0] = m_headIndex[1] = m_headIndex[2] = m_headIndex[3] = 0;
    m_tailIndex[0] = m_tailIndex[1] = m_tailIndex[2] = m_tailIndex[3] = 0;
  }

  /// コンストラクタ
  cdm_Array( CDM::E_CDM_DTYPE dtype
            , CDM::E_CDM_ARRAYSHAPE shape
            , size_t ix
            , size_t jx
            , size_t kx
            , size_t gc
            , size_t nvari=1 )
  {
    m_sz[0] = m_szI[0] = ix;
    m_sz[1] = m_szI[1] = jx;
    m_sz[2] = m_szI[2] = kx;

    switch(shape)
    {
    case CDM::E_CDM_IJKN:
      m_Sz[0] = m_SzI[0] = ix+2*gc;
      m_Sz[1] = m_SzI[1] = jx+2*gc;
      m_Sz[2] = m_SzI[2] = kx+2*gc;
      m_Sz[3] = m_SzI[3] = nvari;
      m_gcl[0] = gc;
      m_gcl[1] = gc;
      m_gcl[2] = gc;
      m_gcl[3] = 0;
      break;
    case CDM::E_CDM_NIJK:
      m_Sz[0] = m_SzI[0] = nvari;
      m_Sz[1] = m_SzI[1] = ix+2*gc;
      m_Sz[2] = m_SzI[2] = jx+2*gc;
      m_Sz[3] = m_SzI[3] = kx+2*gc;
      m_gcl[0] = 0;
      m_gcl[1] = gc;
      m_gcl[2] = gc;
      m_gcl[3] = gc;
    }

    m_gc = m_gcI = gc;
    m_nvari = m_nvariI = nvari;
    m_dtype = dtype;
    m_shape = shape;

    int head[3]={0,0,0};
    setHeadIndex(head);
  }

/*
 * メンバ変数
 */
public:


protected:
  /// データタイプ
  CDM::E_CDM_DTYPE m_dtype;

  /// 配列形状
  CDM::E_CDM_ARRAYSHAPE m_shape;

  /// ガイドセル数
  size_t m_gc;

  /// 格子数
  size_t m_sz[3];

  /// ガイドセルを含んだ格子数
  size_t m_Sz[4];

  /// ガイドセル数(インデクス毎)
  size_t m_gcl[4];

  /// 変数の個数
  size_t m_nvari;


  /// ガイドセル数(int)
  int m_gcI;

  /// 格子数(int)
  int m_szI[3];

  /// ガイドセルを含んだ格子数(int)
  int m_SzI[4];

  /// 変数の個数(int)
  int m_nvariI;


  /// headインデックス
  int m_headIndex[4];

  /// tailインデックス
  int m_tailIndex[4];
};

#endif /* _CDM_ARRAY_H_ */

