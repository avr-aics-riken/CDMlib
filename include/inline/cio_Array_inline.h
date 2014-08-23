/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef _CIO_ARRAY_INLINE_H_
#define _CIO_ARRAY_INLINE_H_

#include "cio_Array.h"
#include <typeinfo>

#ifdef CIO_INLINE
 #undef CIO_INLINE
#endif

#ifndef CIO_NO_INLINE
 #define CIO_INLINE inline
#else
 #define CIO_INLINE
#endif

#define CIO_MEMFUN(rettype) \
        CIO_INLINE rettype 

// インスタンス
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( CIO::E_CIO_DTYPE dtype
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , size_t ix
                        , size_t jx
                        , size_t kx
                        , size_t gc
                        , size_t ncomp )
{
  cio_Array *ptr = NULL;
  switch( dtype )
  {
  case CIO::E_CIO_INT8:
    ptr = new cio_TypeArray<char>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_INT16:
    ptr = new cio_TypeArray<short>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_INT32:
    ptr = new cio_TypeArray<int>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_INT64:
    ptr = new cio_TypeArray<long long>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_UINT8:
    ptr = new cio_TypeArray<unsigned char>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_UINT16:
    ptr = new cio_TypeArray<unsigned short>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_UINT32:
    ptr = new cio_TypeArray<unsigned int>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_UINT64:
    ptr = new cio_TypeArray<unsigned long long>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_FLOAT32:
    ptr = new cio_TypeArray<float>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  case CIO::E_CIO_FLOAT64:
    ptr = new cio_TypeArray<double>(dtype,shape,ix,jx,kx,gc,ncomp);
    break;
  }

#ifdef _CIO_DEBUG
  if( ptr )
  {
    printf("dtype = %d\n",(int)dtype);
    printf("shape = %d\n",(int)shape);
    printf("ixjxkx = %d %d %d\n",(int)ix,(int)jx,(int)kx);
    printf("gc = %d\n",(int)gc);
    printf("ncomp = %d\n",(int)ncomp);
    size_t *m_Sz=ptr->m_Sz;
    size_t *m_gcl=ptr->m_gcl;
    printf("Sz  = %d %d %d %d\n",(int)m_Sz[0],(int)m_Sz[1],(int)m_Sz[2],(int)m_Sz[3]);
    printf("gcl = %d %d %d %d\n",(int)m_gcl[0],(int)m_gcl[1],(int)m_gcl[2],(int)m_gcl[3]);
  }
#endif

  return ptr;
}

// インスタンス
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( CIO::E_CIO_DTYPE dtype
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , size_t sz[3]
                        , size_t gc
                        , size_t ncomp )
{
  return instanceArray(dtype,shape,sz[0],sz[1],sz[2],gc,ncomp);
}

// インスタンス
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( CIO::E_CIO_DTYPE dtype
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , int ix
                        , int jx
                        , int kx
                        , int gc
                        , int ncomp )
{
  return instanceArray(dtype,shape,size_t(ix),size_t(jx),size_t(kx),size_t(gc),size_t(ncomp));
}

// インスタンス
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( CIO::E_CIO_DTYPE dtype
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , int sz[3]
                        , int gc
                        , int ncomp )
{
  return instanceArray(dtype,shape,size_t(sz[0]),size_t(sz[1]),size_t(sz[2]),size_t(gc),size_t(ncomp));
}

// インスタンス
template<class T>
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( T *data
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , size_t ix
                        , size_t jx
                        , size_t kx
                        , size_t gc
                        , size_t ncomp )
{
  cio_Array *ptr = NULL;
  CIO::E_CIO_DTYPE dtype = CIO::E_CIO_DTYPE_UNKNOWN; 

  if( typeid(data) == typeid(char*) )
  {
    dtype = CIO::E_CIO_INT8;
  }
  else if( typeid(data) == typeid(short*) )
  {
    dtype = CIO::E_CIO_INT16;
  }
  else if( typeid(data) == typeid(int*) )
  {
    dtype = CIO::E_CIO_INT32;
  }
  else if( typeid(data) == typeid(long long*) )
  {
    dtype = CIO::E_CIO_INT64;
  }
  else if( typeid(data) == typeid(unsigned char*) )
  {
    dtype = CIO::E_CIO_UINT8;
  }
  else if( typeid(data) == typeid(unsigned short*) )
  {
    dtype = CIO::E_CIO_UINT16;
  }
  else if( typeid(data) == typeid(unsigned int*) )
  {
    dtype = CIO::E_CIO_UINT32;
  }
  else if( typeid(data) == typeid(unsigned long long*) )
  {
    dtype = CIO::E_CIO_UINT64;
  }
  else if( typeid(data) == typeid(float*) )
  {
    dtype = CIO::E_CIO_FLOAT32;
  }
  else if( typeid(data) == typeid(double*) )
  {
    dtype = CIO::E_CIO_FLOAT64;
  }

  if( dtype != CIO::E_CIO_DTYPE_UNKNOWN )
  {
    ptr = new cio_TypeArray<T>(data,dtype,shape,ix,jx,kx,gc,ncomp);
  }

#ifdef _CIO_DEBUG
  if( ptr )
  {
    printf("dtype = %d\n",(int)dtype);
    printf("shape = %d\n",(int)shape);
    printf("ixjxkx = %d %d %d\n",(int)ix,(int)jx,(int)kx);
    printf("gc = %d\n",(int)gc);
    printf("ncomp = %d\n",(int)ncomp);
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
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( T *data
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , size_t sz[3]
                        , size_t gc
                        , size_t ncomp )
{
  return instanceArray(data,shape,sz[0],sz[1],sz[2],gc,ncomp);
}

// インスタンス
template<class T>
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( T *data
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , int ix
                        , int jx
                        , int kx
                        , int gc
                        , int ncomp )
{
  return instanceArray(data,shape,size_t(ix),size_t(jx),size_t(kx),size_t(gc),size_t(ncomp));
}

// インスタンス
template<class T>
CIO_MEMFUN(cio_Array*)
cio_Array::instanceArray( T *data
                        , CIO::E_CIO_ARRAYSHAPE shape
                        , int sz[3]
                        , int gc
                        , int ncomp )
{
  return instanceArray(data,shape,size_t(sz[0]),size_t(sz[1]),size_t(sz[2]),size_t(gc),size_t(ncomp));
}

/// データポインタを取得
CIO_MEMFUN(void*)
cio_Array::getData( bool extract )
{
  switch( m_dtype )
  {
  case CIO::E_CIO_INT8:
    {
      cio_TypeArray<char> *ptr = dynamic_cast<cio_TypeArray<char>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_INT16:
    {
      cio_TypeArray<short> *ptr = dynamic_cast<cio_TypeArray<short>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_INT32:
    {
      cio_TypeArray<int> *ptr = dynamic_cast<cio_TypeArray<int>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_INT64:
    {
      cio_TypeArray<long long> *ptr = dynamic_cast<cio_TypeArray<long long>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_UINT8:
    {
      cio_TypeArray<unsigned char> *ptr = dynamic_cast<cio_TypeArray<unsigned char>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_UINT16:
    {
      cio_TypeArray<unsigned short> *ptr = dynamic_cast<cio_TypeArray<unsigned short>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_UINT32:
    {
      cio_TypeArray<unsigned int> *ptr = dynamic_cast<cio_TypeArray<unsigned int>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_UINT64:
    {
      cio_TypeArray<unsigned long long> *ptr = dynamic_cast<cio_TypeArray<unsigned long long>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_FLOAT32:
    {
      cio_TypeArray<float> *ptr = dynamic_cast<cio_TypeArray<float>*>(this);
      return ptr->getData( extract );
    }
    break;
  case CIO::E_CIO_FLOAT64:
    {
      cio_TypeArray<double> *ptr = dynamic_cast<cio_TypeArray<double>*>(this);
      return ptr->getData( extract );
    }
    break;
  }

  return NULL;
}

// 参照
template<class T>
CIO_MEMFUN(T&)
cio_TypeArray<T>::val(int i, int j, int k, int n)
{
  return m_data[ m_Sz[2] * m_Sz[1] * m_Sz[0] * size_t(n+m_gcl[3])
                         + m_Sz[1] * m_Sz[0] * size_t(k+m_gcl[2])
                                   + m_Sz[0] * size_t(j+m_gcl[1])
                                             + size_t(i+m_gcl[0]) ];
}

// 参照(const)
template<class T>
CIO_MEMFUN(const T&)
cio_TypeArray<T>::val(int i, int j, int k, int n) const
{
  return val(i,j,k,n);
}

// 参照(headインデクス考慮版)
template<class T>
CIO_MEMFUN(T&)
cio_TypeArray<T>::hval(int i, int j, int k, int n)
{
  return m_data[ m_Sz[2] * m_Sz[1] * m_Sz[0] * size_t(n-m_headIndex[3]+m_gcl[3])
                         + m_Sz[1] * m_Sz[0] * size_t(k-m_headIndex[2]+m_gcl[2])
                                   + m_Sz[0] * size_t(j-m_headIndex[1]+m_gcl[1])
                                             + size_t(i-m_headIndex[0]+m_gcl[0]) ];
}

// 参照(headインデクス考慮版,const)
template<class T>
CIO_MEMFUN(const T&)
cio_TypeArray<T>::hval(int i, int j, int k, int n) const
{
  return hval(i,j,k,n);
}

// ガイドセルを含んだ参照
template<class T>
CIO_MEMFUN(T&)
cio_TypeArray<T>::_val(size_t i, size_t j, size_t k, size_t n)
{
  return m_data[ m_Sz[2] * m_Sz[1] * m_Sz[0] * n
                         + m_Sz[1] * m_Sz[0] * k
                                   + m_Sz[0] * j
                                             + i ];
}

// ガイドセルを含んだ参照(const)
template<class T>
CIO_MEMFUN(const T&)
cio_TypeArray<T>::_val(size_t i, size_t j, size_t k, size_t n) const
{
  return _val(i,j,k,n);
}

// 配列コピー
template<class T>
CIO_MEMFUN(int)
cio_TypeArray<T>::copyArray( cio_Array *dst, bool ignoreGc )
{
  cio_TypeArray<T> *src = this;

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
CIO_MEMFUN(int)
cio_TypeArray<T>::copyArray( int _sta[3], int _end[3], cio_Array *dstptr )
{
  cio_TypeArray<T> *src = this;

  //mod.s 
  cio_TypeArray<T> *dst = dynamic_cast<cio_TypeArray<T>*>(dstptr);
  if( !dst )
  {
    return 1;
  }

  // データタイプのチェック
  if( src->getDataType() != dst->getDataType() )
  {
    return 2;
  }
  CIO::E_CIO_DTYPE dtype = src->getDataType();

  // 配列形状
  if( src->getArrayShape() != dst->getArrayShape() )
  {
    return 3;
  }
  CIO::E_CIO_ARRAYSHAPE shape = src->getArrayShape();

  // 成分数
  if( src->getNcomp() != src->getNcomp() )
  {
    return 4;
  }
  int ncomp = src->getNcomp();

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
  if( m_shape == CIO::E_CIO_IJKN )
  {
    for( int n=0;n<ncomp;n++ ){
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
    for( int n=0;n<ncomp;n++ ){
      dst->hval(n,i,j,k) = src->hval(n,i,j,k);
    }}}}
  }

  return 0;
}

//FCONV 20131216.s
//成分指定の配列コピー
template<class T>
CIO_MEMFUN(int)
cio_TypeArray<T>::copyArrayNcomp( cio_Array *dst, int comp, bool ignoreGc )
{
  cio_TypeArray<T> *src = this;

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

  return copyArrayNcomp(sta,end,dst,comp);
}

//成分指定の範囲指定での配列コピー
template<class T>
CIO_MEMFUN(int)
cio_TypeArray<T>::copyArrayNcomp( int _sta[3], int _end[3], cio_Array *dstptr, int comp )
{
  cio_TypeArray<T> *src = this;

  cio_TypeArray<T> *dst = dynamic_cast<cio_TypeArray<T>*>(dstptr);
  if( !dst )
  {
    return 1;
  }

  // データタイプのチェック
  if( src->getDataType() != dst->getDataType() )
  {
    return 2;
  }
  CIO::E_CIO_DTYPE dtype = src->getDataType();

  //配列形状
  if( src->getArrayShape() != dst->getArrayShape() )
  {
    return 3;
  }
  CIO::E_CIO_ARRAYSHAPE shape = src->getArrayShape();

  //成分数
  if( src->getNcomp() != src->getNcomp() )
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
  if( m_shape == CIO::E_CIO_IJKN )
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(i,j,k,comp) = src->hval(i,j,k,0);
    }}}
  }
  else
  {
    for( int k=sta[2];k<=end[2];k++ ){
    for( int j=sta[1];j<=end[1];j++ ){
    for( int i=sta[0];i<=end[0];i++ ){
      dst->hval(comp,i,j,k) = src->hval(0,i,j,k);
    }}}
  }

  return 0;

}

// 粗密データの補間処理を行う
CIO_MEMFUN(cio_Array*)
cio_Array::interp_coarse( cio_Array *src, int &err, bool head0start )
{
  err = 1;

  // データタイプ
  // 実数型のみ対応
  CIO::E_CIO_DTYPE dtype = src->getDataType();
  if( dtype != CIO::E_CIO_FLOAT32 && dtype != CIO::E_CIO_FLOAT64 )
  {
    err = -1;
    return NULL;
  }

  // 配列形状
  CIO::E_CIO_ARRAYSHAPE shape = src->getArrayShape();

  // 成分数
  // 成分数は1か3のみ対応
  int ncomp = src->getNcomp();
  if( ncomp != 1 && ncomp != 3 )
  {
    err = -1;
    return NULL;
  }

  // その他の情報の取得
  int gcS = src->getGc();
  void *ptrS = src->getData();
  const int *szS = src->getArraySizeInt();
  const int *headS = src->getHeadIndex();
  const int *tailS = src->getTailIndex();

  // 密配列のインスタンス
  int gcD = gcS*2;
  int szD[3] = {szS[0]*2, szS[1]*2, szS[2]*2};
  cio_Array *dst = cio_Array::instanceArray( dtype, shape, szD, gcD, ncomp );
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
  if( shape == CIO::E_CIO_IJKN )
  {
    if( dtype == CIO::E_CIO_FLOAT32 )
    {
      cio_interp_ijkn_r4_(szS,&gcS,szD,&gcD,&ncomp,(float*)ptrS,(float*)ptrD);
    }
    else
    {
      cio_interp_ijkn_r8_(szS,&gcS,szD,&gcD,&ncomp,(double*)ptrS,(double*)ptrD);
    }
  }
  else
  {
    if( dtype == CIO::E_CIO_FLOAT32 )
    {
      cio_interp_nijk_r4_(szS,&gcS,szD,&gcD,&ncomp,(float*)ptrS,(float*)ptrD);
    }
    else
    {
      cio_interp_nijk_r8_(szS,&gcS,szD,&gcD,&ncomp,(double*)ptrS,(double*)ptrD);
    }
  }

  return dst;
}

// 配列サイズ分のバイナリデータを読み込み(戻り値は読み込んだ要素数)
template<class T>
size_t cio_TypeArray<T>::readBinary( FILE *fp, bool bMatchEndian )
{
  if( !fp ) return size_t(0);
  size_t ndata = getArrayLength();
  size_t nread = fread(m_data,sizeof(T),ndata,fp);
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
size_t cio_TypeArray<T>::writeBinary( FILE *fp )
{
  if( !fp ) return size_t(0);
  return fwrite(m_data,sizeof(T),getArrayLength(),fp);
}

// 配列サイズ分のasciiデータを書き出す(戻り値は読み込んだ要素数)
template<class T>
size_t cio_TypeArray<T>::writeAscii( FILE *fp )
{
  if( !fp ) return size_t(0);
  //return fwrite(m_data,sizeof(T),getArrayLength(),fp);

  for(int i=0; i<getArrayLength(); i++) {
    fprintf(fp,"%e\n",(float)m_data[i]);
  }

  return getArrayLength();

}

#endif /* _CIO_ARRAY_INLINE_H_ */
