#ifndef _CONV_INLINE_H_
#define _CONV_INLINE_H_

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
 * @file   conv_inline.h
 * @brief  CONV クラスのinline関数ヘッダーファイル
 * @author aics
 * @date   2013/11/7
 */

#ifndef CONV_NO_INLINE
 #define CONV_INLINE inline
#else
 #define CONV_INLINE
#endif

// ##################################################################################
template<class T>
CONV_INLINE
bool CONV::copyArray(cio_TypeArray<T> *B,
                     cio_Array *&src,
                     int sta[3],
                     int end[3],
                     int n)
{

  CIO::E_CIO_DTYPE src_dtype =  src->getDataType();

  //uint8
  if( src_dtype == CIO::E_CIO_UINT8 ) {
    cio_TypeArray<unsigned char> *S = dynamic_cast<cio_TypeArray<unsigned char>*>(src);
    return copyArray(B,S,sta,end,n);
  } 
  //int8
  else if( src_dtype== CIO::E_CIO_INT8 ) {
    cio_TypeArray<char> *S = dynamic_cast<cio_TypeArray<char>*>(src);
    return copyArray(B,S,sta,end,n);
  } 
  //uint16
  else if( src_dtype== CIO::E_CIO_UINT16 ) {
    cio_TypeArray<unsigned short> *S = dynamic_cast<cio_TypeArray<unsigned short>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //int16
  else if( src_dtype== CIO::E_CIO_INT16 ) {
    cio_TypeArray<short> *S = dynamic_cast<cio_TypeArray<short>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //uint32
  else if( src_dtype== CIO::E_CIO_UINT32 ) {
    cio_TypeArray<unsigned int> *S = dynamic_cast<cio_TypeArray<unsigned int>*>(src);
    return copyArray(B,S,sta,end,n);
  } 
  //int32
  else if( src_dtype== CIO::E_CIO_INT32 ) { 
    cio_TypeArray<int> *S = dynamic_cast<cio_TypeArray<int>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //uint64
  else if( src_dtype== CIO::E_CIO_UINT64 ) {
    cio_TypeArray<unsigned long long> *S = dynamic_cast<cio_TypeArray<unsigned long long>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //int64
  else if( src_dtype== CIO::E_CIO_INT64 ) {
    cio_TypeArray<long long> *S = dynamic_cast<cio_TypeArray<long long>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //float32
  else if( src_dtype== CIO::E_CIO_FLOAT32 ) {
    cio_TypeArray<float> *S = dynamic_cast<cio_TypeArray<float>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  else if( src_dtype== CIO::E_CIO_FLOAT64 ) { 
    cio_TypeArray<double> *S = dynamic_cast<cio_TypeArray<double>*>(src);
    return copyArray(B,S,sta,end,n);
  }

  return false;

}

// ##################################################################################
// 配列のコピー
template<class T1, class T2>
CONV_INLINE
bool CONV::copyArray(cio_TypeArray<T1> *buf,
               cio_TypeArray<T2> *&src,
               int sta[3],
               int end[3],
               int n)
{

  //配列形状
  CIO::E_CIO_ARRAYSHAPE buf_shape = buf->getArrayShape();

  CIO::E_CIO_ARRAYSHAPE src_shape = src->getArrayShape();

  //入力指示領域の取得
  int IndexStart[3];
  if( m_param->Get_CropIndexStart_on() ) {
    const int *cropIndexStart = m_param->Get_CropIndexStart();
    for(int i=0; i<3; i++) IndexStart[i]=cropIndexStart[i];
  } else {
    //for(int i=0; i<3; i++) IndexStart[i]=sta[i]+1;
    for(int i=0; i<3; i++) IndexStart[i]=1;
  }

  const int* headS = src->getHeadIndex();

  //間引き数の取得
  int thin_count = m_param->Get_ThinOut();

  //IJKN&IJKN
  if( buf_shape == CIO::E_CIO_IJKN && src_shape == CIO::E_CIO_IJKN )
  {
    for( int k=sta[2];k<=end[2];k++ ){
    if( (k-(IndexStart[2]-1))%thin_count != 0 ) continue;

    for( int j=sta[1];j<=end[1];j++ ){
    if( (j-(IndexStart[1]-1))%thin_count != 0 ) continue;

    for( int i=sta[0];i<=end[0];i++ ){
      if( (i-(IndexStart[0]-1))%thin_count != 0 ) continue;

      src->hval(i/thin_count,j/thin_count,k/thin_count,n) = buf->hval(i,j,k,n);
    }}}
  }
  else if( buf_shape == CIO::E_CIO_NIJK && src_shape == CIO::E_CIO_NIJK )
  //NIJK&NIJK
  {
    for( int k=sta[2];k<=end[2];k++ ){
    if( (k-(IndexStart[2]-1))%thin_count != 0 ) continue;

    for( int j=sta[1];j<=end[1];j++ ){
    if( (j-(IndexStart[1]-1))%thin_count != 0 ) continue;

    for( int i=sta[0];i<=end[0];i++ ){
      if( (i-(IndexStart[0]-1))%thin_count != 0 ) continue;

      src->hval(n,i/thin_count,j/thin_count,k/thin_count) = buf->hval(n,i,j,k);
    }}}
  }
  else if( buf_shape == CIO::E_CIO_IJKN && src_shape == CIO::E_CIO_NIJK )
  //IJNK&NIJK
  {

    for( int k=sta[2];k<=end[2];k++ ){
    if( (k-(IndexStart[2]-1))%thin_count != 0 ) continue;
    for( int j=sta[1];j<=end[1];j++ ){
    if( (j-(IndexStart[1]-1))%thin_count != 0 ) continue;
    for( int i=sta[0];i<=end[0];i++ ){
      if( (i-(IndexStart[0]-1))%thin_count != 0 ) continue;

      src->hval(n,i/thin_count,j/thin_count,k/thin_count) = buf->hval(i,j,k,n);

    }}}
  }
  else if( buf_shape == CIO::E_CIO_NIJK && src_shape == CIO::E_CIO_IJKN )
  //NIJK&IJKN
  {
    for( int k=sta[2];k<=end[2];k++ ){
    if( (k-(IndexStart[2]-1))%thin_count != 0 ) continue;

    for( int j=sta[1];j<=end[1];j++ ){
    if( (j-(IndexStart[1]-1))%thin_count != 0 ) continue;

    for( int i=sta[0];i<=end[0];i++ ){
      if( (i-(IndexStart[0]-1))%thin_count != 0 ) continue;

      src->hval(i/thin_count,j/thin_count,k/thin_count,n) = buf->hval(n,i,j,k);
    }}}  
  }

  return true;
}

// ##################################################################################
// minmax
template<class T>
CONV_INLINE
bool CONV::calcMinMax(cio_TypeArray<T> *src,
                      double *min,
                      double *max)
{

  if( src == NULL ) return false;

  double CompVal;

  //sizeの取得
  const int *sz = src->getArraySizeInt();
  //配列形状の取得
  CIO::E_CIO_ARRAYSHAPE shape = src->getArrayShape();
  //成分数の取得
  int nComp = src->getNcomp();


  if( nComp > 1 ) {
    //nijkの処理
    if( shape == CIO::E_CIO_NIJK ) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        CompVal=(double)0.0;
        for(int n=0; n<nComp; n++) {
          if( min[n] > (double)src->val(n,i,j,k) ) min[n] = (double)src->val(n,i,j,k);
          if( max[n] < (double)src->val(n,i,j,k) ) max[n] = (double)src->val(n,i,j,k);
          CompVal = CompVal + (double)src->val(n,i,j,k)*(double)src->val(n,i,j,k); 
        }
        CompVal = sqrt(CompVal);
        if( min[nComp] > CompVal ) min[nComp]=CompVal;
        if( max[nComp] < CompVal ) max[nComp]=CompVal;
      }}}
    
    } 
    else if( shape == CIO::E_CIO_IJKN ) 
    //ijknの処理
    {  
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        CompVal=(double)0.0;
        for(int n=0; n<nComp; n++) {
          if( min[n] > (double)src->val(i,j,k,n) ) min[n] = (double)src->val(i,j,k,n);
          if( max[n] < (double)src->val(i,j,k,n) ) max[n] = (double)src->val(i,j,k,n);
          CompVal = CompVal + (double)src->val(i,j,k,n)*(double)src->val(i,j,k,n); 
        }
        CompVal = sqrt(CompVal);
        if( min[nComp] > CompVal ) min[nComp]=CompVal;
        if( max[nComp] < CompVal ) max[nComp]=CompVal;
      }}}
    } else return false;

  } else {
    //nijkの処理
    if( shape == CIO::E_CIO_NIJK ) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        if( min[0] > (double)src->val(0,i,j,k) ) min[0] = (double)src->val(0,i,j,k);
        if( max[0] < (double)src->val(0,i,j,k) ) max[0] = (double)src->val(0,i,j,k);
      }}}
    } 
    else if( shape == CIO::E_CIO_IJKN ) 
    //ijknの処理
    {  
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        if( min[0] > (double)src->val(i,j,k,0) ) min[0] = (double)src->val(i,j,k,0);
        if( max[0] < (double)src->val(i,j,k,0) ) max[0] = (double)src->val(i,j,k,0);
      }}}
    } else return false;

  }
  return true;

}
#endif // _CONV_INLINE_H
