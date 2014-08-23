#ifndef _CIO_ARRAY_H_
#define _CIO_ARRAY_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "cio_Define.h"

extern "C"
{
  void cio_interp_ijkn_r4_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *ncomp, float *src, float *dst );
  void cio_interp_ijkn_r8_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *ncomp, double *src, double *dst );
  void cio_interp_nijk_r4_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *ncomp, float *src, float *dst );
  void cio_interp_nijk_r8_(const int *szS, const int* gcS, const int *szD, const int *gcD, const int *ncomp, double *src, double *dst );
}

class cio_Array
{
/*
 * enum
 */
public:
/*
  /// データタイプ
  enum E_CIO_DTYPE
  {
    E_CIO_DTYPE_UNKNOWN=0,
    E_CIO_INT8,
    E_CIO_INT16,
    E_CIO_INT32,
    E_CIO_INT64,
    E_CIO_UINT8,
    E_CIO_UINT16,
    E_CIO_UINT32,
    E_CIO_UINT64,
    E_CIO_FLOAT32,
    E_CIO_FLOAT64
  };

  /// 配列形状タイプ
  enum E_CIO_ARRAYSHAPE
  {
    E_CIO_ARRAYSHAPE_UNKNOWN=0,
    E_CIO_IJKN,
    E_CIO_NIJK,
  };
*/

/*
 *  メンバ関数
 */
public:
  /// デストラクタ
  virtual ~cio_Array()
  {
  }

  /// インスタンス
  static cio_Array*
  instanceArray( CIO::E_CIO_DTYPE dtype
               , CIO::E_CIO_ARRAYSHAPE shape
               , size_t ix
               , size_t jx
               , size_t kx
               , size_t gc
               , size_t ncomp=1 );

  /// インスタンス
  static cio_Array*
  instanceArray( CIO::E_CIO_DTYPE dtype
               , CIO::E_CIO_ARRAYSHAPE shape
               , size_t sz[3]
               , size_t gc
               , size_t ncomp=1 );

  /// インスタンス
  static cio_Array*
  instanceArray( CIO::E_CIO_DTYPE dtype
               , CIO::E_CIO_ARRAYSHAPE shape
               , int ix
               , int jx
               , int kx
               , int gc
               , int ncomp=1 );

  /// インスタンス
  static cio_Array*
  instanceArray( CIO::E_CIO_DTYPE dtype
               , CIO::E_CIO_ARRAYSHAPE shape
               , int sz[3]
               , int gc
               , int ncomp=1 );

  /// インスタンス
  template<class T>
  static cio_Array*
  instanceArray( T *data
               , CIO::E_CIO_ARRAYSHAPE shape
               , size_t ix
               , size_t jx
               , size_t kx
               , size_t gc
               , size_t ncomp=1 );

  /// インスタンス
  template<class T>
  static cio_Array*
  instanceArray( T *data
               , CIO::E_CIO_ARRAYSHAPE shape
               , size_t sz[3]
               , size_t gc
               , size_t ncomp=1 );

  /// インスタンス
  template<class T>
  static cio_Array*
  instanceArray( T *data
               , CIO::E_CIO_ARRAYSHAPE shape
               , int ix
               , int jx
               , int kx
               , int gc
               , int ncomp=1 );

  /// インスタンス
  template<class T>
  static cio_Array*
  instanceArray( T *data
               , CIO::E_CIO_ARRAYSHAPE shape
               , int sz[3]
               , int gc
               , int ncomp=1 );

  /// データポインタを取得
  void* getData( bool extract=false );

  /// データタイプの取得
  CIO::E_CIO_DTYPE getDataType() const
  {
    return m_dtype;
  }

  /// データタイプ文字列の取得
  char* getDataTypeString() const
  {
    switch( m_dtype )
    {
    case CIO::E_CIO_INT8:
      return "INT8";
      break;
    case CIO::E_CIO_INT16:
      return "INT16";
      break;
    case CIO::E_CIO_INT32:
      return "INT32";
      break;
    case CIO::E_CIO_INT64:
      return "INT64";
      break;
    case CIO::E_CIO_UINT8:
      return "UINT8";
      break;
    case CIO::E_CIO_UINT16:
      return "UINT16";
      break;
    case CIO::E_CIO_UINT32:
      return "UINT32";
      break;
    case CIO::E_CIO_UINT64:
      return "UINT64";
      break;
    case CIO::E_CIO_FLOAT32:
      return "FLOAT32";
      break;
    case CIO::E_CIO_FLOAT64:
      return "FLOAT64";
      break;
    }
    return "Unknown";
  }

  /// 配列形状の取得
  CIO::E_CIO_ARRAYSHAPE getArrayShape() const
  {
    return m_shape;
  }

  /// 配列形状文字列の取得
  char* getArrayShapeString() const
  {
    switch(m_shape)
    {
    case CIO::E_CIO_IJKN:
      return "IJKN";
      break;
    case CIO::E_CIO_NIJK:
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

  /// 成分数を取得
  size_t getNcomp() const
  {
    return m_ncomp;
  }

  /// 成分数を取得(int版)
  int getNcompInt() const
  {
    return m_ncompI;
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
    case CIO::E_CIO_IJKN:
      return m_headIndex;
      break;
    case CIO::E_CIO_NIJK:
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
    case CIO::E_CIO_IJKN:
      return m_tailIndex;
      break;
    case CIO::E_CIO_NIJK:
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
    case CIO::E_CIO_IJKN:
      return m_Sz;
      break;
    case CIO::E_CIO_NIJK:
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
    case CIO::E_CIO_IJKN:
      return m_SzI;
      break;
    case CIO::E_CIO_NIJK:
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
    case CIO::E_CIO_IJKN:
      m_headIndex[0] = head[0];
      m_headIndex[1] = head[1];
      m_headIndex[2] = head[2];
      m_headIndex[3] = 0;
      m_tailIndex[0] = m_headIndex[0] + m_sz[0] - 1;
      m_tailIndex[1] = m_headIndex[1] + m_sz[1] - 1;
      m_tailIndex[2] = m_headIndex[2] + m_sz[2] - 1;
      m_tailIndex[3] = 0;
      break;
    case CIO::E_CIO_NIJK:
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
  virtual int copyArray( cio_Array *dst, bool ignoreGc=false ) = 0;

  /// 範囲指定での配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArray( int sta[3], int end[3], cio_Array *dst ) = 0;

//FCONV 20131216.s
  /// 指定成分の配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArrayNcomp( cio_Array *dst, int comp, bool ignoreGc=false ) = 0;

  /// 指定成分の範囲指定での配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArrayNcomp( int sta[3], int end[3], cio_Array *dst, int comp ) = 0;
//FCONV 20131216.e

  /// 粗密データの補間処理を行う
  static cio_Array* interp_coarse( cio_Array *src, int &err, bool head0start=true );

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
  cio_Array()
  {
    m_dtype = CIO::E_CIO_DTYPE_UNKNOWN;
    m_shape = CIO::E_CIO_ARRAYSHAPE_UNKNOWN;
    m_sz[0] = m_sz[1] = m_sz[2] = 0;
    m_Sz[0] = m_Sz[1] = m_Sz[2] = m_Sz[3] = 0;
    m_gc = 0;
    m_gcl[0] = m_gcl[1] = m_gcl[2] = m_gcl[3] = 0;
    m_ncomp = 1;
    m_headIndex[0] = m_headIndex[1] = m_headIndex[2] = m_headIndex[3] = 0;
    m_tailIndex[0] = m_tailIndex[1] = m_tailIndex[2] = m_tailIndex[3] = 0;
  }

  /// コンストラクタ
  cio_Array( CIO::E_CIO_DTYPE dtype
            , CIO::E_CIO_ARRAYSHAPE shape
            , size_t ix
            , size_t jx
            , size_t kx
            , size_t gc
            , size_t ncomp=1 )
  {
    m_sz[0] = m_szI[0] = ix;
    m_sz[1] = m_szI[1] = jx;
    m_sz[2] = m_szI[2] = kx;

    switch(shape)
    {
    case CIO::E_CIO_IJKN:
      m_Sz[0] = m_SzI[0] = ix+2*gc;
      m_Sz[1] = m_SzI[1] = jx+2*gc;
      m_Sz[2] = m_SzI[2] = kx+2*gc;
      m_Sz[3] = m_SzI[3] = ncomp;
      m_gcl[0] = gc;
      m_gcl[1] = gc;
      m_gcl[2] = gc;
      m_gcl[3] = 0;
      break;
    case CIO::E_CIO_NIJK:
      m_Sz[0] = m_SzI[0] = ncomp;
      m_Sz[1] = m_SzI[1] = ix+2*gc;
      m_Sz[2] = m_SzI[2] = jx+2*gc;
      m_Sz[3] = m_SzI[3] = kx+2*gc;
      m_gcl[0] = 0;
      m_gcl[1] = gc;
      m_gcl[2] = gc;
      m_gcl[3] = gc;
    }

    m_gc = m_gcI = gc;
    m_ncomp = m_ncompI = ncomp;
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
  CIO::E_CIO_DTYPE m_dtype;

  /// 配列形状
  CIO::E_CIO_ARRAYSHAPE m_shape;

  /// ガイドセル数
  size_t m_gc;

  /// 格子数
  size_t m_sz[3];

  /// ガイドセルを含んだ格子数
  size_t m_Sz[4];

  /// ガイドセル数(インデクス毎)
  size_t m_gcl[4];

  /// 成分数
  size_t m_ncomp;


  /// ガイドセル数(int)
  int m_gcI;

  /// 格子数(int)
  int m_szI[3];

  /// ガイドセルを含んだ格子数(int)
  int m_SzI[4];

  /// 成分数(int)
  int m_ncompI;


  /// headインデックス
  int m_headIndex[4];

  /// tailインデックス
  int m_tailIndex[4];
};

#endif /* _CIO_ARRAY_H_ */

