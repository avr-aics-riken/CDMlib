#ifndef _CDM_DFI_INLINE_H_
#define _CDM_DFI_INLINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI.h
 * @brief  cdm_DFI Class Header
 * @author aics    
 */

#ifdef CDM_INLINE
 #undef CDM_INLINE
#endif

#ifndef CDM_NO_INLINE
 #define CDM_INLINE inline
#else
 #define CDM_INLINE
#endif

// #################################################################
// フィールドデータの読込み(読み込んだデータのポインタを戻り値として
// 返す）

//template<class T, class TimeT, class TimeAvrT> 
//CDM_INLINE T*
template<class TimeT, class TimeAvrT> 
CDM_INLINE void*
cdm_DFI::ReadData(CDM::E_CDM_ERRORCODE &ret,
                  const unsigned step, 
                  const int gc, 
                  const int Gvoxel[3], 
                  const int Gdivision[3], 
                  const int head[3], 
                  const int tail[3],
                  TimeT &time,
                  const bool mode, 
                  unsigned &step_avr, 
                  TimeAvrT &time_avr)
{

   int sz[3];
   for(int i=0; i<3; i++) sz[i]=tail[i]-head[i]+1;
   cdm_Array *data = cdm_Array::instanceArray
                     ( DFI_Finfo.DataType
                     , DFI_Finfo.ArrayShape
                     , sz
                     , gc
                     , DFI_Finfo.Component);

   double d_time = (double)time;
   double d_time_avr = (double)time_avr;

//   int ret = ReadData(data, step, gc, Gvoxel, Gdivision, head, tail,
   ret = ReadData(data, step, gc, Gvoxel, Gdivision, head, tail,
                       d_time, mode, step_avr, d_time_avr);

   if( ret != CDM::E_CDM_SUCCESS ) {
     delete data;
     return NULL;
   }

//   T* ptr = (T*)data->getData(true);
   void* ptr = data->getData(true);
   delete data;
   time = d_time;
   time_avr = d_time_avr;

   return ptr;
}

// #################################################################
// フィールドデータの読込み(引数で渡された配列にデータを読込む）
template<class T, class TimeT, class TimeAvrT> 
CDM_INLINE
CDM::E_CDM_ERRORCODE cdm_DFI::ReadData(T *val,
                                       const unsigned step,
                                       const int gc,
                                       const int Gvoxel[3],
                                       const int Gdivision[3],
                                       const int head[3],
                                       const int tail[3],
                                       TimeT &time,
                                       const bool mode,
                                       unsigned &step_avr,
                                       TimeAvrT &time_avr)
{

   int sz[3];
   for(int i=0; i<3; i++) sz[i]=tail[i]-head[i]+1;

   cdm_Array *data = cdm_Array::instanceArray
                     ( val
                     , DFI_Finfo.ArrayShape
                     , sz
                     , gc
                     , DFI_Finfo.Component);

   double d_time = (double)time;
   double d_time_avr = (double)time_avr;

   CDM::E_CDM_ERRORCODE ret;
   ret = ReadData(data, step, gc, Gvoxel, Gdivision, head, tail,
                  d_time, mode, step_avr, d_time_avr);

   if( ret == CDM::E_CDM_SUCCESS ) {
     time = d_time;
     time_avr = d_time_avr;
   }

   //data->getData(true);
   delete data;

   return ret; 
}

// #################################################################
// フィールドデータの出力
template<class T, class TimeT, class TimeAvrT> 
CDM_INLINE
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteData(const unsigned step, 
                   TimeT time, 
                   const int sz[3],
                   const int nComp,
                   const int gc, 
                   T* val, 
                   T* minmax,
                   const bool avr_mode, 
                   const unsigned step_avr, 
                   TimeAvrT time_avr)
{

  cdm_Array *data = cdm_Array::instanceArray
                    ( val
                    , DFI_Finfo.ArrayShape
                    , DFI_Process.RankList[m_RankID].VoxelSize[0]
                    , DFI_Process.RankList[m_RankID].VoxelSize[1]
                    , DFI_Process.RankList[m_RankID].VoxelSize[2]
                    , gc
                    , DFI_Finfo.Component);

  double d_time = (double)time;
  double d_time_avr = (double)time_avr;
  double *d_minmax=NULL;
  if( minmax ) {
    if( DFI_Finfo.Component>1 ) {
      d_minmax = new double[DFI_Finfo.Component*2+2];
      for(int i=0; i<DFI_Finfo.Component*2+2; i++) {
        d_minmax[i] = minmax[i];
      }
    } else { 
      d_minmax = new double[2];
      d_minmax[0] = minmax[0];
      d_minmax[1] = minmax[1];
    }
  }

  CDM::E_CDM_ERRORCODE ret;
  ret = WriteData(step, gc, d_time, data, d_minmax, avr_mode, step_avr, d_time_avr);

  //val = (T*)data->getData(true);
  //data->getData(true);

  if( d_minmax ) delete [] d_minmax;

  delete data;
  return ret;
                                          
}

// #################################################################
//セル中心データを 格子点にセット
template<class T1, class T2>
CDM_INLINE
bool
cdm_DFI::setGridData(cdm_TypeArray<T1>* P,
                     cdm_TypeArray<T2>* S)
{

  if( P->getArrayShape() != S->getArrayShape() ) return false;

  //成分数をセット
  if( P->getNcompInt() != S->getNcompInt() ) return false;
  int nComp = P->getNcompInt();

  //S(セル中心）の配列サイズを取得セット
  //T2* data = S->getData();
  const int* size = S->getArraySizeInt();
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];

  //P（格子点）の配列サイズを取得セット
  T1* d    = P->getData();
  const int* Psz = P->getArraySizeInt();
  int id = Psz[0];
  int jd = Psz[1];
  int kd = Psz[2];

  //Pの配列をゼロクリア
  size_t dsize = (size_t)(id*jd*kd*nComp);
  for (size_t l=0; l<dsize; l++) d[l]=0.0;

  //S(セル中心）のデータをP(格子点)に加える
  //NIJKの処理
  if( P->getArrayShape() == CDM::E_CDM_NIJK ) {
    for (int km=0; km<kx; km++) {
    for (int jm=0; jm<jx; jm++) {
    for (int im=0; im<ix; im++) {
    for (int n=0; n<nComp; n++) {
      P->val(n, im  ,jm  ,km  ) = P->val(n, im  ,jm  ,km  )+S->val(n, im,jm,km); ///<0,0,0>
      P->val(n, im+1,jm  ,km  ) = P->val(n, im+1,jm  ,km  )+S->val(n, im,jm,km); ///<1,0,0>
      P->val(n, im+1,jm  ,km+1) = P->val(n, im+1,jm  ,km+1)+S->val(n, im,jm,km); ///<1,0,1>
      P->val(n, im  ,jm  ,km+1) = P->val(n, im  ,jm  ,km+1)+S->val(n, im,jm,km); ///<0,0,1>
      P->val(n, im  ,jm+1,km  ) = P->val(n, im  ,jm+1,km  )+S->val(n, im,jm,km); ///<0,1,0>
      P->val(n, im+1,jm+1,km  ) = P->val(n, im+1,jm+1,km  )+S->val(n, im,jm,km); ///<1,1,0>
      P->val(n, im+1,jm+1,km+1) = P->val(n, im+1,jm+1,km+1)+S->val(n, im,jm,km); ///<1,1,1>
      P->val(n, im  ,jm+1,km+1) = P->val(n, im  ,jm+1,km+1)+S->val(n, im,jm,km); ///<0,1,1>
    }}}}
  } else {
  //IJKNの処理
    for (int n=0; n<nComp; n++) {
    for (int km=0; km<kx; km++) {
    for (int jm=0; jm<jx; jm++) {
    for (int im=0; im<ix; im++) {
      P->val(im  ,jm  ,km  ,n) = P->val(im  ,jm  ,km  ,n)+S->val(im,jm,km,n); ///<0,0,0>
      P->val(im+1,jm  ,km  ,n) = P->val(im+1,jm  ,km  ,n)+S->val(im,jm,km,n); ///<1,0,0>
      P->val(im+1,jm  ,km+1,n) = P->val(im+1,jm  ,km+1,n)+S->val(im,jm,km,n); ///<1,0,1>
      P->val(im  ,jm  ,km+1,n) = P->val(im  ,jm  ,km+1,n)+S->val(im,jm,km,n); ///<0,0,1>
      P->val(im  ,jm+1,km  ,n) = P->val(im  ,jm+1,km  ,n)+S->val(im,jm,km,n); ///<0,1,0>
      P->val(im+1,jm+1,km  ,n) = P->val(im+1,jm+1,km  ,n)+S->val(im,jm,km,n); ///<1,1,0>
      P->val(im+1,jm+1,km+1,n) = P->val(im+1,jm+1,km+1,n)+S->val(im,jm,km,n); ///<1,1,1>
      P->val(im  ,jm+1,km+1,n) = P->val(im  ,jm+1,km+1,n)+S->val(im,jm,km,n); ///<0,1,1>
    }}}}
  }

  //内部の格子点のデータを重み付けで割る
  VolumeDataDivide(P);

  return true;

}

// #################################################################
//内部の格子点のデータを重み付けで割る
template<class T>
CDM_INLINE
void
cdm_DFI::VolumeDataDivide(cdm_TypeArray<T> *P)
{
  int i,j,k,n;
  const int* szP = P->getArraySizeInt();
  int id = szP[0];
  int jd = szP[1];
  int kd = szP[2];

  int ncomp = P->getNcompInt();

  //NIJK
  if( P->getArrayShape() == CDM::E_CDM_NIJK ) {

    //I
    for (k=0; k<kd;    k++){
    for (j=0; j<jd;    j++){
    for (i=1; i<id-1;  i++){
    for (n=0; n<ncomp; n++){
      P->val(n,i,j,k) = P->val(n,i,j,k)*0.5;
    }}}}

    //J
    for (k=0; k<kd;    k++){
    for (j=1; j<jd-1;  j++){
    for (i=0; i<id;    i++){
    for (n=0; n<ncomp; n++){
      P->val(n,i,j,k) = P->val(n,i,j,k)*0.5;
    }}}}

    //K
    for (k=1; k<kd-1;  k++){
    for (j=0; j<jd;    j++){
    for (i=0; i<id;    i++){
    for (n=0; n<ncomp; n++){
      P->val(n,i,j,k) = P->val(n,i,j,k)*0.5;
    }}}}

  //IJKN
  } else {

    //I
    for (n=0; n<ncomp; n++){
    for (k=0; k<kd;    k++){
    for (j=0; j<jd;    j++){
    for (i=1; i<id-1;  i++){
      P->val(i,j,k,n) = P->val(i,j,k,n)*0.5;
    }}}}

    //J
    for (n=0; n<ncomp; n++){
    for (k=0; k<kd;    k++){
    for (j=1; j<jd-1;  j++){
    for (i=0; i<id;    i++){
      P->val(i,j,k,n) = P->val(i,j,k,n)*0.5;
    }}}}

    //K
    for (n=0; n<ncomp; n++){
    for (k=1; k<kd-1;  k++){
    for (j=0; j<jd;    j++){
    for (i=0; i<id;    i++){
      P->val(i,j,k,n) = P->val(i,j,k,n)*0.5;
    }}}}

  }
}


#endif // _CDM_DFI_INLINE_H_
