#ifndef _CONV_EXEC_INLINE_H_
#define _CONV_EXEC_INLINE_H_

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
 * @file   convMx1_inline.h
 * @brief  convMx1 クラスのinline関数ヘッダーファイル
 * @author aics
 * @date   2013/11/7
 */

#ifndef CONV_NO_INLINE
 #define CONV_INLINE inline
#else
 #define CONV_INLINE
#endif

/**
 * @brief 配列のゼロクリア
 * @param[out] data 配列
 * @param[in]  ivar_out コンポーネント位置
 */
template<class T>
CONV_INLINE
void convMx1::zeroClearArray(cdm_TypeArray<T>* data, int ivar_out)
{

  const int *sz = data->getArraySizeInt();
  CDM::E_CDM_ARRAYSHAPE shape = data->getArrayShape();

  if( shape == CDM::E_CDM_NIJK ) {
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      data->val(ivar_out,i,j,k) = (T)0.0;
    }}}
  } else {
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      data->val(i,j,k,ivar_out) = (T)0.0;
    }}}
  }
};

/**
 * @brief Scalarの格子点での値をセット
 * @param [out] O 格子点data
 * @param [in]  S セル中心data
 * @param[in]  ivar_out 格子データの コンポーネント位置
 * @param[in]  ivar_src 図心データの コンポーネント位置
 */
template<class T>
CONV_INLINE
bool convMx1::setGridData_XY( 
                             cdm_TypeArray<T>* O,
                             cdm_TypeArray<T>* S,
                             int ivar_out,
                             int ivar_src)
{
  if( O->getArrayShape() != S->getArrayShape() ) return false;

  //ガイドセル数の取得
  int gc = S->getGcInt();
  if( gc < 1 ) return false; ///ガイドセルがない場合は処理しない

  //S(図心）の配列サイズをセット
  const int* size = S->getArraySizeInt();
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];

  //仮想セル領域へのコピー
  if( S->getArrayShape() == CDM::E_CDM_NIJK ) {
    for(int j=0; j<jx; j++) {
      S->val(ivar_src,-1,j,0) = S->val(ivar_src,0,j,0);
      S->val(ivar_src,ix,j,0) = S->val(ivar_src,ix-1,j,0);
    }
    for(int i=-1; i<ix+1; i++) {
      S->val(ivar_src,i,-1,0) = S->val(ivar_src,i,0,0);
      S->val(ivar_src,i,jx,0) = S->val(ivar_src,i,jx-1,0);
    }
  } else {
    for(int j=0; j<jx; j++) {
      S->val(-1,j,0,ivar_src) = S->val(0,j,0,ivar_src);
      S->val(ix,j,0,ivar_src) = S->val(ix-1,j,0,ivar_src);
    }
    for(int i=-1; i<ix+1; i++) {
      S->val(i,-1,0,ivar_src) = S->val(i,0,0,ivar_src);
      S->val(i,jx,0,ivar_src) = S->val(i,jx-1,0,ivar_src);
    }
  }

  //O（格子点）の配列サイズをセット
  const int *Osz = O->getArraySizeInt();
  int id = Osz[0];
  int jd = Osz[1];
  int kd = Osz[2];

  //図心データを格子点に加える
  if( O->getArrayShape() == CDM::E_CDM_NIJK ) {
    for (int km=0; km<kx; km++) {
    for (int jm=0-gc; jm<jx+gc; jm++) {
    for (int im=0-gc; im<ix+gc; im++) {
      O->val(ivar_out, im  ,jm  ,km) = O->val(ivar_out, im  ,jm  ,km)+S->val(ivar_src,im,jm,km);
      O->val(ivar_out, im+1,jm  ,km) = O->val(ivar_out, im+1,jm  ,km)+S->val(ivar_src,im,jm,km);
      O->val(ivar_out, im  ,jm+1,km) = O->val(ivar_out, im  ,jm+1,km)+S->val(ivar_src,im,jm,km);
      O->val(ivar_out, im+1,jm+1,km) = O->val(ivar_out, im+1,jm+1,km)+S->val(ivar_src,im,jm,km);
    }}}
  } else {
    for (int km=0; km<kx; km++) {
    for (int jm=0-gc; jm<jx+gc; jm++) {
    for (int im=0-gc; im<ix+gc; im++) {
      O->val(im  ,jm  ,km ,ivar_out) = O->val(im  ,jm  ,km,ivar_out)+S->val(im,jm,km,ivar_src);
      O->val(im+1,jm  ,km ,ivar_out) = O->val(im+1,jm  ,km,ivar_out)+S->val(im,jm,km,ivar_src);
      O->val(im  ,jm+1,km ,ivar_out) = O->val(im  ,jm+1,km,ivar_out)+S->val(im,jm,km,ivar_src);
      O->val(im+1,jm+1,km ,ivar_out) = O->val(im+1,jm+1,km,ivar_out)+S->val(im,jm,km,ivar_src);
    }}}
  }
  return true;
};

// #################################################################
// 内部の格子点のデータを重み付けで割る(XY)
template<class T>
CONV_INLINE
void
convMx1::VolumeDataDivide8(cdm_TypeArray<T> *O, int n)
{
  const int* szO = O->getArraySizeInt();
  int id = szO[0];
  int jd = szO[1];
  int kd = szO[2];


  //NIJK
  if( O->getArrayShape() == CDM::E_CDM_NIJK ) {
    //I
    for(int k=0; k<kd; k++) {
    for(int j=0; j<jd; j++) {
    for(int i=0; i<id; i++) {
    //for(int n=0; n<nVari; n++) {
      O->val(n,i,j,k) = O->val(n,i,j,k)*0.125;
    //}}}}
    }}}

  //IJKN
  } else {
    //I
    //for (int n=0; n<nVari; n++){
    for (int k=0; k<kd; k++){
    for (int j=0; j<jd; j++){
    for (int i=0; i<id; i++){
      O->val(i,j,k,n) = O->val(i,j,k,n)*0.125;
    //}}}}
    }}}
  }

};

// #################################################################
// NIJK配列をスカラーのIJK配列にコピー
template<class T>
CONV_INLINE
void
convMx1::copyArray_nijk_ijk(cdm_TypeArray<T> *S, cdm_TypeArray<T> *O, int ivar)
{
  const int* sz = S->getArraySizeInt();
  int gc = O->getGcInt();
  if( S->getArrayShape() == CDM::E_CDM_NIJK ) {
    for(int k=0-gc; k<sz[2]+gc; k++) {
    for(int j=0-gc; j<sz[1]+gc; j++) {
    for(int i=0-gc; i<sz[0]+gc; i++) {
       O->val(i,j,k,0) = S->val(ivar,i,j,k);
    }}}
  }
  else {
    for(int k=0-gc; k<sz[2]+gc; k++) {
    for(int j=0-gc; j<sz[1]+gc; j++) {
    for(int i=0-gc; i<sz[0]+gc; i++) {
       O->val(i,j,k,0) = S->val(i,j,k,ivar);
    }}}
  }
};

#endif // _CONV_EXEC_INLINE_H_ 
