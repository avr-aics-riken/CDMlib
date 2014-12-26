/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef _CDM_ARRAY_INLINE_H_
#define _CDM_ARRAY_INLINE_H_

#include "cdm_Array.h"
#include <typeinfo>

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

// インスタンス
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( CDM::E_CDM_DTYPE dtype
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , size_t ix
                        , size_t jx
                        , size_t kx
                        , size_t gc
                        , size_t nvari )
{
  cdm_Array *ptr = NULL;
  switch( dtype )
  {
  case CDM::E_CDM_INT8:
    ptr = new cdm_TypeArray<char>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_INT16:
    ptr = new cdm_TypeArray<short>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_INT32:
    ptr = new cdm_TypeArray<int>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_INT64:
    ptr = new cdm_TypeArray<long long>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_UINT8:
    ptr = new cdm_TypeArray<unsigned char>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_UINT16:
    ptr = new cdm_TypeArray<unsigned short>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_UINT32:
    ptr = new cdm_TypeArray<unsigned int>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_UINT64:
    ptr = new cdm_TypeArray<unsigned long long>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_FLOAT32:
    ptr = new cdm_TypeArray<float>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  case CDM::E_CDM_FLOAT64:
    ptr = new cdm_TypeArray<double>(dtype,shape,ix,jx,kx,gc,nvari);
    break;
  }

#ifdef _CDM_DEBUG
  if( ptr )
  {
    printf("dtype = %d\n",(int)dtype);
    printf("shape = %d\n",(int)shape);
    printf("ixjxkx = %d %d %d\n",(int)ix,(int)jx,(int)kx);
    printf("gc = %d\n",(int)gc);
    printf("nvari = %d\n",(int)nvari);
    size_t *m_Sz=ptr->m_Sz;
    size_t *m_gcl=ptr->m_gcl;
    printf("Sz  = %d %d %d %d\n",(int)m_Sz[0],(int)m_Sz[1],(int)m_Sz[2],(int)m_Sz[3]);
    printf("gcl = %d %d %d %d\n",(int)m_gcl[0],(int)m_gcl[1],(int)m_gcl[2],(int)m_gcl[3]);
  }
#endif

  return ptr;
}

// インスタンス
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( CDM::E_CDM_DTYPE dtype
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , size_t sz[3]
                        , size_t gc
                        , size_t nvari )
{
  return instanceArray(dtype,shape,sz[0],sz[1],sz[2],gc,nvari);
}

// インスタンス
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( CDM::E_CDM_DTYPE dtype
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , int ix
                        , int jx
                        , int kx
                        , int gc
                        , int nvari )
{
  return instanceArray(dtype,shape,size_t(ix),size_t(jx),size_t(kx),size_t(gc),size_t(nvari));
}

// インスタンス
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( CDM::E_CDM_DTYPE dtype
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , int sz[3]
                        , int gc
                        , int nvari )
{
  return instanceArray(dtype,shape,size_t(sz[0]),size_t(sz[1]),size_t(sz[2]),size_t(gc),size_t(nvari));
}

// インスタンス
template<class T>
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( T *data
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , size_t ix
                        , size_t jx
                        , size_t kx
                        , size_t gc
                        , size_t nvari )
{
  cdm_Array *ptr = NULL;
  CDM::E_CDM_DTYPE dtype = CDM::E_CDM_DTYPE_UNKNOWN; 

  if( typeid(data) == typeid(char*) )
  {
    dtype = CDM::E_CDM_INT8;
  }
  else if( typeid(data) == typeid(short*) )
  {
    dtype = CDM::E_CDM_INT16;
  }
  else if( typeid(data) == typeid(int*) )
  {
    dtype = CDM::E_CDM_INT32;
  }
  else if( typeid(data) == typeid(long long*) )
  {
    dtype = CDM::E_CDM_INT64;
  }
  else if( typeid(data) == typeid(unsigned char*) )
  {
    dtype = CDM::E_CDM_UINT8;
  }
  else if( typeid(data) == typeid(unsigned short*) )
  {
    dtype = CDM::E_CDM_UINT16;
  }
  else if( typeid(data) == typeid(unsigned int*) )
  {
    dtype = CDM::E_CDM_UINT32;
  }
  else if( typeid(data) == typeid(unsigned long long*) )
  {
    dtype = CDM::E_CDM_UINT64;
  }
  else if( typeid(data) == typeid(float*) )
  {
    dtype = CDM::E_CDM_FLOAT32;
  }
  else if( typeid(data) == typeid(double*) )
  {
    dtype = CDM::E_CDM_FLOAT64;
  }

  if( dtype != CDM::E_CDM_DTYPE_UNKNOWN )
  {
    ptr = new cdm_TypeArray<T>(data,dtype,shape,ix,jx,kx,gc,nvari);
  }

#ifdef _CDM_DEBUG
  if( ptr )
  {
    printf("dtype = %d\n",(int)dtype);
    printf("shape = %d\n",(int)shape);
    printf("ixjxkx = %d %d %d\n",(int)ix,(int)jx,(int)kx);
    printf("gc = %d\n",(int)gc);
    printf("nvari = %d\n",(int)nvari);
    size_t *m_Sz=ptr->m_Sz;
    size_t *m_gcl=ptr->m_gcl;
    printf("Sz  = %d %d %d %d\n",(int)m_Sz[0],(int)m_Sz[1],(int)m_Sz[2],(int)m_Sz[3]);
    printf("gcl = %d %d %d %d\n",(int)m_gcl[0],(int)m_gcl[1],(int)m_gcl[2],(int)m_gcl[3]);
  }
#endif

  return ptr;
}

// インスタンス
template<class T>
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( T *data
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , size_t sz[3]
                        , size_t gc
                        , size_t nvari )
{
  return instanceArray(data,shape,sz[0],sz[1],sz[2],gc,nvari);
}

// インスタンス
template<class T>
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( T *data
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , int ix
                        , int jx
                        , int kx
                        , int gc
                        , int nvari )
{
  return instanceArray(data,shape,size_t(ix),size_t(jx),size_t(kx),size_t(gc),size_t(nvari));
}

// インスタンス
template<class T>
CDM_MEMFUN(cdm_Array*)
cdm_Array::instanceArray( T *data
                        , CDM::E_CDM_ARRAYSHAPE shape
                        , int sz[3]
                        , int gc
                        , int nvari )
{
  return instanceArray(data,shape,size_t(sz[0]),size_t(sz[1]),size_t(sz[2]),size_t(gc),size_t(nvari));
}

/// データポインタを取得
CDM_MEMFUN(void*)
cdm_Array::getData( bool extract )
{
  switch( m_dtype )
  {
  case CDM::E_CDM_INT8:
    {
      cdm_TypeArray<char> *ptr = dynamic_cast<cdm_TypeArray<char>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_INT16:
    {
      cdm_TypeArray<short> *ptr = dynamic_cast<cdm_TypeArray<short>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_INT32:
    {
      cdm_TypeArray<int> *ptr = dynamic_cast<cdm_TypeArray<int>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_INT64:
    {
      cdm_TypeArray<long long> *ptr = dynamic_cast<cdm_TypeArray<long long>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_UINT8:
    {
      cdm_TypeArray<unsigned char> *ptr = dynamic_cast<cdm_TypeArray<unsigned char>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_UINT16:
    {
      cdm_TypeArray<unsigned short> *ptr = dynamic_cast<cdm_TypeArray<unsigned short>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_UINT32:
    {
      cdm_TypeArray<unsigned int> *ptr = dynamic_cast<cdm_TypeArray<unsigned int>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_UINT64:
    {
      cdm_TypeArray<unsigned long long> *ptr = dynamic_cast<cdm_TypeArray<unsigned long long>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_FLOAT32:
    {
      cdm_TypeArray<float> *ptr = dynamic_cast<cdm_TypeArray<float>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CDM::E_CDM_FLOAT64:
    {
      cdm_TypeArray<double> *ptr = dynamic_cast<cdm_TypeArray<double>*>(this);
      return ptr->getData( extract );
    }
    break;
  }

  return NULL;
}

// 参照
template<class T>
CDM_MEMFUN(T&)
cdm_TypeArray<T>::val(int i, int j, int k, int n)
{
  return m_data[ m_Sz[2] * m_Sz[1] * m_Sz[0] * size_t(n+m_gcl[3])
                         + m_Sz[1] * m_Sz[0] * size_t(k+m_gcl[2])
                                   + m_Sz[0] * size_t(j+m_gcl[1])
                                             + size_t(i+m_gcl[0]) ];
}

// 参照(const)
template<class T>
CDM_MEMFUN(const T&)
cdm_TypeArray<T>::val(int i, int j, int k, int n) const
{
  return val(i,j,k,n);
}

// 参照(headインデクス考慮版)
template<class T>
CDM_MEMFUN(T&)
cdm_TypeArray<T>::hval(int i, int j, int k, int n)
{
  return m_data[ m_Sz[2] * m_Sz[1] * m_Sz[0] * size_t(n-m_headIndex[3]+m_gcl[3])
                         + m_Sz[1] * m_Sz[0] * size_t(k-m_headIndex[2]+m_gcl[2])
                                   + m_Sz[0] * size_t(j-m_headIndex[1]+m_gcl[1])
                                             + size_t(i-m_headIndex[0]+m_gcl[0]) ];
}

// 参照(headインデクス考慮版,const)
template<class T>
CDM_MEMFUN(const T&)
cdm_TypeArray<T>::hval(int i, int j, int k, int n) const
{
  return hval(i,j,k,n);
}

// ガイドセルを含んだ参照
template<class T>
CDM_MEMFUN(T&)
cdm_TypeArray<T>::_val(size_t i, size_t j, size_t k, size_t n)
{
  return m_data[ m_Sz[2] * m_Sz[1] * m_Sz[0] * n
                         + m_Sz[1] * m_Sz[0] * k
                                   + m_Sz[0] * j
                                             + i ];
}

// ガイドセルを含んだ参照(const)
template<class T>
CDM_MEMFUN(const T&)
cdm_TypeArray<T>::_val(size_t i, size_t j, size_t k, size_t n) const
{
  return _val(i,j,k,n);
}

// 配列コピー
template<class T>
CDM_MEMFUN(int)
cdm_TypeArray<T>::copyArray( cdm_Array *dst, bool ignoreGc )
{
  cdm_TypeArray<T> *src = this;

  // コピーの範囲
  int       gcS    = src->getGcInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();
  int       gcD    = dst->getGcInt();
  const int *headD = dst->getHeadIndex();
  const int *tailD = dst->getTailIndex();
  if( ignoreGc )
  {
    gcS = gcD = 0;
  }
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headS[i]-gcS>=headD[i]-gcD) ? headS[i]-gcS : headD[i]-gcD;
    end[i] = (tailS[i]+gcS<=tailD[i]+gcD) ? tailS[i]+gcS : tailD[i]+gcD;
  }

  return copyArray(sta,end,dst);
}

// 範囲指定での配列コピー
template<class T>
CDM_MEMFUN(int)
cdm_TypeArray<T>::copyArray( int _sta[3], int _end[3], cdm_Array *dstptr )
{
  cdm_TypeArray<T> *src = this;

  //mod.s 
  cdm_TypeArray<T> *dst = dynamic_cast<cdm_TypeArray<T>*>(dstptr);
  if( !dst )
  {
    return 1;
  }

  // データタイプのチェック
  if( src->getDataType() != dst->getDataType() )
  {
    return 2;
  }
  CDM::E_CDM_DTYPE dtype = src->getDataType();

  // 配列形状
  if( src->getArrayShape() != dst->getArrayShape() )
  {
    return 3;
  }
  CDM::E_CDM_ARRAYSHAPE shape = src->getArrayShape();

  // 変数の個数
  if( src->getNvari() != src->getNvari() )
  {
    return 4;
  }
  int nvari = src->getNvari();

  // コピーの範囲
  int       gcS    = src->getGcInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();
  int       gcD    = dst->getGcInt();
  const int *headD = dst->getHeadIndex();
  const int *tailD = dst->getTailIndex();
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headS[i]-gcS>=headD[i]-gcD) ? headS[i]-gcS : headD[i]-gcD;
    end[i] = (tailS[i]+gcS<=tailD[i]+gcD) ? tailS[i]+gcS : tailD[i]+gcD;
  }
  for( int i=0;i<3;i++ )
  {
    sta[i] = (_sta[i]>=sta[i]) ? _sta[i] : sta[i];
    end[i] = (_end[i]<=end[i]) ? _end[i] : end[i];
  }

  // コピー
  if( m_shape == CDM::E_CDM_IJKN )
  {
    for( int n=0;n<nvari;n++ ){
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(i,j,k,n) = src->hval(i,j,k,n);
    }}}}
  }
  else
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
    for( int n=0;n<nvari;n++ ){
      dst->hval(n,i,j,k) = src->hval(n,i,j,k);
    }}}}
  }

  return 0;
}

//FCONV 20131216.s
//変数指定の配列コピー
template<class T>
CDM_MEMFUN(int)
cdm_TypeArray<T>::copyArrayNvari( cdm_Array *dst, int vari, bool ignoreGc )
{
  cdm_TypeArray<T> *src = this;

  // コピーの範囲
  int       gcS    = src->getGcInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();
  int       gcD    = dst->getGcInt();
  const int *headD = dst->getHeadIndex();
  const int *tailD = dst->getTailIndex();
  if( ignoreGc )
  {
    gcS = gcD = 0;
  }
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headS[i]-gcS>=headD[i]-gcD) ? headS[i]-gcS : headD[i]-gcD;
    end[i] = (tailS[i]+gcS<=tailD[i]+gcD) ? tailS[i]+gcS : tailD[i]+gcD;
  }

  return copyArrayNvari(sta,end,dst,vari);
}

//変数指定の範囲指定での配列コピー
template<class T>
CDM_MEMFUN(int)
cdm_TypeArray<T>::copyArrayNvari( int _sta[3], int _end[3], cdm_Array *dstptr, int vari )
{
  cdm_TypeArray<T> *src = this;

  cdm_TypeArray<T> *dst = dynamic_cast<cdm_TypeArray<T>*>(dstptr);
  if( !dst )
  {
    return 1;
  }

  // データタイプのチェック
  if( src->getDataType() != dst->getDataType() )
  {
    return 2;
  }
  CDM::E_CDM_DTYPE dtype = src->getDataType();

  //配列形状
  if( src->getArrayShape() != dst->getArrayShape() )
  {
    return 3;
  }
  CDM::E_CDM_ARRAYSHAPE shape = src->getArrayShape();

  //変数の個数
  if( src->getNvari() != src->getNvari() )
  {
    return 4;
  }

  //コピーの範囲
  int       gcS    = src->getGcInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();
  int       gcD    = dst->getGcInt();
  const int *headD = dst->getHeadIndex();
  const int *tailD = dst->getTailIndex();
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headS[i]-gcS>=headD[i]-gcD) ? headS[i]-gcS : headD[i]-gcD;
    end[i] = (tailS[i]+gcS<=tailD[i]+gcD) ? tailS[i]+gcS : tailD[i]+gcD;
  }
  for( int i=0;i<3;i++ )
  {
    sta[i] = (_sta[i]>=sta[i]) ? _sta[i] : sta[i];
    end[i] = (_end[i]<=end[i]) ? _end[i] : end[i];
  }

  // コピー
  if( m_shape == CDM::E_CDM_IJKN )
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(i,j,k,vari) = src->hval(i,j,k,0);
    }}}
  }
  else
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(vari,i,j,k) = src->hval(0,i,j,k);
    }}}
  }

  return 0;

}

//指定した変数の配列のみ取得し、IJK配列でコピー
template<class T>
CDM_MEMFUN(int)
cdm_TypeArray<T>::copyArrayNvari_to_ijk( cdm_Array *dst, int vari, bool ignoreGc )
{
  cdm_TypeArray<T> *src = this;

  // コピーの範囲
  int       gcS    = src->getGcInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();
  int       gcD    = dst->getGcInt();
  const int *headD = dst->getHeadIndex();
  const int *tailD = dst->getTailIndex();
  if( ignoreGc )
  {
    gcS = gcD = 0;
  }
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headS[i]-gcS>=headD[i]-gcD) ? headS[i]-gcS : headD[i]-gcD;
    end[i] = (tailS[i]+gcS<=tailD[i]+gcD) ? tailS[i]+gcS : tailD[i]+gcD;
  }

  return copyArrayNvari_to_ijk(sta,end,dst,vari);
}

//指定した変数の配列のみ範囲指定で取得し、IJK配列でコピー
template<class T>
CDM_MEMFUN(int)
cdm_TypeArray<T>::copyArrayNvari_to_ijk( int _sta[3], int _end[3], cdm_Array *dstptr, int vari )
{
  cdm_TypeArray<T> *src = this;

  cdm_TypeArray<T> *dst = dynamic_cast<cdm_TypeArray<T>*>(dstptr);
  if( !dst )
  {
    return 1;
  }

  // データタイプのチェック
  if( src->getDataType() != dst->getDataType() )
  {
    return 2;
  }
  CDM::E_CDM_DTYPE dtype = src->getDataType();

  //配列形状
  if( dst->getArrayShape() != CDM::E_CDM_IJKN )
  {
    return 3;
  }
  CDM::E_CDM_ARRAYSHAPE shape = src->getArrayShape();

  //変数の個数
  if( src->getNvari() != src->getNvari() )
  {
    return 4;
  }

  //コピーの範囲
  int       gcS    = src->getGcInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();
  int       gcD    = dst->getGcInt();
  const int *headD = dst->getHeadIndex();
  const int *tailD = dst->getTailIndex();
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headS[i]-gcS>=headD[i]-gcD) ? headS[i]-gcS : headD[i]-gcD;
    end[i] = (tailS[i]+gcS<=tailD[i]+gcD) ? tailS[i]+gcS : tailD[i]+gcD;
  }
  for( int i=0;i<3;i++ )
  {
    sta[i] = (_sta[i]>=sta[i]) ? _sta[i] : sta[i];
    end[i] = (_end[i]<=end[i]) ? _end[i] : end[i];
  }

  // vari番目の変数の値のみ取得して、dstにコピー
  if( m_shape == CDM::E_CDM_IJKN )
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(i,j,k,0) = src->hval(i,j,k,vari);
    }}}
  }
  else
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(i,j,k,0) = src->hval(vari,i,j,k);
    }}}
  }

  return 0;

}

// 粗密データの補間処理を行う
CDM_MEMFUN(cdm_Array*)
cdm_Array::interp_coarse( cdm_Array *src, int &err, bool head0start )
{
  err = 1;

  // データタイプ
  // 実数型のみ対応
  CDM::E_CDM_DTYPE dtype = src->getDataType();
  if( dtype != CDM::E_CDM_FLOAT32 && dtype != CDM::E_CDM_FLOAT64 )
  {
    err = -1;
    return NULL;
  }

  // 配列形状
  CDM::E_CDM_ARRAYSHAPE shape = src->getArrayShape();

  // 変数の個数
  int nvari = src->getNvari();

  // その他の情報の取得
  int gcS = src->getGc();
  void *ptrS = src->getData();
  const int *szS = src->getArraySizeInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();

  // 密配列のインスタンス
  int gcD = gcS*2;
  int szD[3] = {szS[0]*2, szS[1]*2, szS[2]*2};
  cdm_Array *dst = cdm_Array::instanceArray( dtype, shape, szD, gcD, nvari );
  void *ptrD = dst->getData();

  // headインデクスのセット
  int headD[3];
  for( int i=0;i<3;i++ )
  {
    headD[i] = headS[i]*2;
    if( !head0start )
    {
      headD[i] -= 1;
    }
  }
  dst->setHeadIndex( headD );

  // f90コードのコール(配列形状、実数型毎)
  if( shape == CDM::E_CDM_IJKN )
  {
    if( dtype == CDM::E_CDM_FLOAT32 )
    {
      cdm_interp_ijkn_r4_(szS,&gcS,szD,&gcD,&nvari,(float*)ptrS,(float*)ptrD);
    }
    else
    {
      cdm_interp_ijkn_r8_(szS,&gcS,szD,&gcD,&nvari,(double*)ptrS,(double*)ptrD);
    }
  }
  else
  {
    if( dtype == CDM::E_CDM_FLOAT32 )
    {
      cdm_interp_nijk_r4_(szS,&gcS,szD,&gcD,&nvari,(float*)ptrS,(float*)ptrD);
    }
    else
    {
      cdm_interp_nijk_r8_(szS,&gcS,szD,&gcD,&nvari,(double*)ptrS,(double*)ptrD);
    }
  }

  return dst;
}

// 配列サイズ分のバイナリデータを読み込み(戻り値は読み込んだ要素数)
template<class T>
size_t cdm_TypeArray<T>::readBinary( FILE *fp, bool bMatchEndian )
{
  if( !fp ) return size_t(0);
  size_t ndata = getArrayLength();
#ifdef CDM_BUFFER_MB_SIZE
  int bufferMBSize = CDM_BUFFER_MB_SIZE;
  int bufferSize = bufferMBSize * 1024 * 1024 / sizeof(T);
//  T* buffer = new T[bufferSize];
  int bufferIndex = 0;
  size_t nread = 0;
  while( bufferIndex + bufferSize < ndata ){
    nread += fread((m_data+bufferIndex),sizeof(T),bufferSize,fp);
//    nread += fread(buffer,sizeof(T),bufferSize,fp);
//    memcpy(&m_data[bufferIndex],buffer,bufferSize*sizeof(T));
    bufferIndex += bufferSize;
  }
  nread += fread((m_data+bufferIndex),sizeof(T),ndata-bufferIndex,fp);
//  nread += fread(buffer,sizeof(T),ndata-bufferIndex,fp);
//  memcpy(&m_data[bufferIndex],buffer,(ndata-bufferIndex)*sizeof(T));
//  delete[] buffer;
#else
  size_t nread = fread(m_data,sizeof(T),ndata,fp);
#endif
  if( !bMatchEndian )
  {
    size_t bsz = sizeof(T);
    if( bsz == 2 )
    {
      SBSWAPVEC(m_data,nread);
    }
    else if( bsz == 4 )
    {
      BSWAPVEC(m_data,nread);
    }
    else if( bsz == 8 )
    {
      DBSWAPVEC(m_data,nread);
    }
  }
  return nread;
}

// 配列サイズ分のバイナリデータを書き出す(戻り値は読み込んだ要素数)
template<class T>
size_t cdm_TypeArray<T>::writeBinary( FILE *fp )
{
  if( !fp ) return size_t(0);
  size_t ndata = getArrayLength();
#ifdef CDM_BUFFER_MB_SIZE
  int bufferMBSize = CDM_BUFFER_MB_SIZE;
  int bufferSize = bufferMBSize * 1024 * 1024 / sizeof(T);
//  T* buffer = new T[bufferSize];
  int bufferIndex = 0;
  size_t nwrite = 0;
  while( bufferIndex + bufferSize < ndata ){
//    memcpy(buffer,&m_data[bufferIndex],bufferSize*sizeof(T));
//    nwrite += fwrite(buffer,sizeof(T),bufferSize,fp);
    nwrite += fwrite(m_data+bufferIndex,sizeof(T),bufferSize,fp);
    bufferIndex += bufferSize;
  }
  nwrite += fwrite(m_data+bufferIndex,sizeof(T),ndata-bufferIndex,fp);
//  memcpy(buffer,&m_data[bufferIndex],(ndata-bufferIndex)*sizeof(T));
//  nwrite += fwrite(buffer,sizeof(T),ndata-bufferIndex,fp);
//  delete[] buffer;
#else
  size_t nwrite = fwrite(m_data,sizeof(T),ndata,fp);
#endif
  return nwrite;
}

// 配列サイズ分のasciiデータを書き出す(戻り値は読み込んだ要素数)
template<class T>
size_t cdm_TypeArray<T>::writeAscii( FILE *fp )
{
  if( !fp ) return size_t(0);
  //return fwrite(m_data,sizeof(T),getArrayLength(),fp);

  for(int i=0; i<getArrayLength(); i++) {
    fprintf(fp,"%e\n",(float)m_data[i]);
  }

  return getArrayLength();

}

#endif /* _CDM_ARRAY_INLINE_H_ */
