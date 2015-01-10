#ifndef _CONV_INLINE_H_
#define _CONV_INLINE_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
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
bool CONV::copyArray(cdm_TypeArray<T> *B,
                     cdm_Array *&src,
                     int sta[3],
                     int end[3],
                     int n)
{

  CDM::E_CDM_DTYPE src_dtype =  src->getDataType();

  //uint8
  if( src_dtype == CDM::E_CDM_UINT8 ) {
    cdm_TypeArray<unsigned char> *S = dynamic_cast<cdm_TypeArray<unsigned char>*>(src);
    return copyArray(B,S,sta,end,n);
  } 
  //int8
  else if( src_dtype== CDM::E_CDM_INT8 ) {
    cdm_TypeArray<char> *S = dynamic_cast<cdm_TypeArray<char>*>(src);
    return copyArray(B,S,sta,end,n);
  } 
  //uint16
  else if( src_dtype== CDM::E_CDM_UINT16 ) {
    cdm_TypeArray<unsigned short> *S = dynamic_cast<cdm_TypeArray<unsigned short>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //int16
  else if( src_dtype== CDM::E_CDM_INT16 ) {
    cdm_TypeArray<short> *S = dynamic_cast<cdm_TypeArray<short>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //uint32
  else if( src_dtype== CDM::E_CDM_UINT32 ) {
    cdm_TypeArray<unsigned int> *S = dynamic_cast<cdm_TypeArray<unsigned int>*>(src);
    return copyArray(B,S,sta,end,n);
  } 
  //int32
  else if( src_dtype== CDM::E_CDM_INT32 ) { 
    cdm_TypeArray<int> *S = dynamic_cast<cdm_TypeArray<int>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //uint64
  else if( src_dtype== CDM::E_CDM_UINT64 ) {
    cdm_TypeArray<unsigned long long> *S = dynamic_cast<cdm_TypeArray<unsigned long long>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //int64
  else if( src_dtype== CDM::E_CDM_INT64 ) {
    cdm_TypeArray<long long> *S = dynamic_cast<cdm_TypeArray<long long>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  //float32
  else if( src_dtype== CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *S = dynamic_cast<cdm_TypeArray<float>*>(src);
    return copyArray(B,S,sta,end,n);
  }
  else if( src_dtype== CDM::E_CDM_FLOAT64 ) { 
    cdm_TypeArray<double> *S = dynamic_cast<cdm_TypeArray<double>*>(src);
    return copyArray(B,S,sta,end,n);
  }

  return false;

}

// ##################################################################################
// 配列のコピー
template<class T1, class T2>
CONV_INLINE
bool CONV::copyArray(cdm_TypeArray<T1> *buf,
               cdm_TypeArray<T2> *&src,
               int sta[3],
               int end[3],
               int n)
{

  //配列形状
  CDM::E_CDM_ARRAYSHAPE buf_shape = buf->getArrayShape();

  CDM::E_CDM_ARRAYSHAPE src_shape = src->getArrayShape();

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
  if( buf_shape == CDM::E_CDM_IJKN && src_shape == CDM::E_CDM_IJKN )
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
  else if( buf_shape == CDM::E_CDM_NIJK && src_shape == CDM::E_CDM_NIJK )
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
  else if( buf_shape == CDM::E_CDM_IJKN && src_shape == CDM::E_CDM_NIJK )
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
  else if( buf_shape == CDM::E_CDM_NIJK && src_shape == CDM::E_CDM_IJKN )
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
bool CONV::calcMinMax(cdm_TypeArray<T> *src,
                      double *min,
                      double *max)
{

  if( src == NULL ) return false;

  double VariVal;

  //sizeの取得
  const int *sz = src->getArraySizeInt();
  //配列形状の取得
  CDM::E_CDM_ARRAYSHAPE shape = src->getArrayShape();
  //変数の個数の取得
  int nVari = src->getNvari();


  if( nVari > 1 ) {
    //nijkの処理
    if( shape == CDM::E_CDM_NIJK ) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        VariVal=(double)0.0;
        for(int n=0; n<nVari; n++) {
          if( min[n] > (double)src->val(n,i,j,k) ) min[n] = (double)src->val(n,i,j,k);
          if( max[n] < (double)src->val(n,i,j,k) ) max[n] = (double)src->val(n,i,j,k);
          VariVal = VariVal + (double)src->val(n,i,j,k)*(double)src->val(n,i,j,k); 
        }
        VariVal = sqrt(VariVal);
        if( min[nVari] > VariVal ) min[nVari]=VariVal;
        if( max[nVari] < VariVal ) max[nVari]=VariVal;
      }}}
    
    } 
    else if( shape == CDM::E_CDM_IJKN ) 
    //ijknの処理
    {  
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        VariVal=(double)0.0;
        for(int n=0; n<nVari; n++) {
          if( min[n] > (double)src->val(i,j,k,n) ) min[n] = (double)src->val(i,j,k,n);
          if( max[n] < (double)src->val(i,j,k,n) ) max[n] = (double)src->val(i,j,k,n);
          VariVal = VariVal + (double)src->val(i,j,k,n)*(double)src->val(i,j,k,n); 
        }
        VariVal = sqrt(VariVal);
        if( min[nVari] > VariVal ) min[nVari]=VariVal;
        if( max[nVari] < VariVal ) max[nVari]=VariVal;
      }}}
    } else return false;

  } else {
    //nijkの処理
    if( shape == CDM::E_CDM_NIJK ) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        if( min[0] > (double)src->val(0,i,j,k) ) min[0] = (double)src->val(0,i,j,k);
        if( max[0] < (double)src->val(0,i,j,k) ) max[0] = (double)src->val(0,i,j,k);
      }}}
    } 
    else if( shape == CDM::E_CDM_IJKN ) 
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
