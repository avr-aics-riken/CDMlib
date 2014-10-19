#ifndef _CDM_TYPEARRAY_H_
#define _CDM_TYPEARRAY_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "cdm_Array.h"

template<class T>
class cdm_TypeArray : public cdm_Array
{
/*
 * enum
 */
public:


/*
 *  メンバ関数
 */
public:
  /// コンストラクタ
  cdm_TypeArray( CDM::E_CDM_DTYPE dtype
               , CDM::E_CDM_ARRAYSHAPE shape
               , size_t ix
               , size_t jx
               , size_t kx
               , size_t gc
               , size_t ncomp=1 )
    : cdm_Array(dtype,shape,ix,jx,kx,gc,ncomp)
  {

    m_outptr=false;

    size_t nw=1;
    for( int i=0;i<3;i++ )
    {
      nw *= (m_sz[i]+2*m_gc);
    }
    nw *= m_ncomp;
    m_data = new T[nw];
    memset(m_data,0,sizeof(T)*nw);
  }

  /// コンストラクタ
  cdm_TypeArray( T *data
               , CDM::E_CDM_DTYPE dtype
               , CDM::E_CDM_ARRAYSHAPE shape
               , size_t ix
               , size_t jx
               , size_t kx
               , size_t gc
               , size_t ncomp=1 )
    : cdm_Array(dtype,shape,ix,jx,kx,gc,ncomp)
  {

    m_outptr=true;

    m_data = data;
  }

  /// デストラクタ
  virtual ~cdm_TypeArray()
  {
    if( m_data )
    {
      if( m_data && !m_outptr )
      {
        delete [] m_data;
      }
    }
  }

  /// 実データのポインタを取得
  T* getData( bool extract=false )
  {
    T *ptr = m_data;
    if( extract )
    {
      m_data = NULL;
    }
    return ptr;
  }

  /** 参照
   *  実セルの最小インデクスを(0,0,0)とする
   *  IJKNのときval(i,j,k,n)
   *  NIJKのときval(n,i,j,k)
   */   
  const T&  val(int i, int j, int k, int l=0) const;
        T&  val(int i, int j, int k, int l=0);

  /** 参照(headインデクス考慮版)
   *  実セルの最小インデクスを(head[0],head[1],head[2])とする
   *  IJKNのときval(i,j,k,n)
   *  NIJKのときval(n,i,j,k)
   */   
  const T& hval(int i, int j, int k, int l=0) const;
        T& hval(int i, int j, int k, int l=0);

  /** 参照(ガイドセルを含む)
   *  ガイドセルを含む配列全体の最小インデクスを(0,0,0)とする
   *  IJKNのときval(i,j,k,n)
   *  NIJKのときval(n,i,j,k)
   */   
  const T& _val(size_t i, size_t j, size_t k, size_t l=0) const;
        T& _val(size_t i, size_t j, size_t k, size_t l=0);

  /// 配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArray( cdm_Array *dst, bool ignoreGc=false );

  /// 範囲指定での配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArray( int sta[3], int end[3], cdm_Array *dst );

//FCONV 20131216.s
  /// 指定成分の配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArrayNcomp( cdm_Array *dst, int comp, bool ignoreGc=false );

  /// 指定成分の範囲指定での配列コピー(自信をdstにコピー。head/tailを考慮した重複範囲をコピー)
  virtual int copyArrayNcomp( int sta[3], int end[3], cdm_Array *dst, int comp );
//FCONV 20131216.e

  /// 配列サイズ分のバイナリデータを読み込み(戻り値は読み込んだ要素数)
  virtual size_t readBinary( FILE *fp, bool bMatchEndian );

  /// 配列サイズ分のバイナリデータを書き出す(戻り値は読み込んだ要素数)
  virtual size_t writeBinary( FILE *fp );

//20131213.s
  /// 配列サイズ分のasciiデータを書き出す(戻り値は読み込んだ要素数)
  virtual size_t writeAscii( FILE *fp );
//20131213.e

protected:
  /// デフォルトコンストラクタ
  cdm_TypeArray() : cdm_Array()
  {
    m_data  = NULL;
  }

/*
 * メンバ変数
 */
public:


protected:

  /// 実データポインタタイプ
  bool m_outptr;

  /// 実データ配列
  T *m_data;

};

#if 0
#ifdef CDM_INLINE
 #undef CDM_INLINE
#endif

#ifndef CDM_NO_INLINE
 #define CDM_INLINE inline
#else
 #define CDM_INLINE
#endif

#define CDM_MEMFUN(rettype) \
        CDM_INLINE rettype 
#endif

#include "inline/cdm_Array_inline.h"

#endif /* _CDM_TYPEARRAY_H_ */

