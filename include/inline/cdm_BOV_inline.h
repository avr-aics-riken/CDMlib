#ifndef _CDM_BOV_INLINE_H_
#define _CDM_BOV_INLINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_BOV_inline.h
 * @brief  cdm_DFI_BOV template Header
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
// ファイルのデータレコード読込み
CDM_INLINE
CDM::E_CDM_ERRORCODE
cdm_DFI_BOV::read_Datarecord(FILE* fp,
                             bool matchEndian,
                             cdm_Array* buf,
                             int head[3],
                             int nz,
                             cdm_Array* &src)
{

  //１層ずつ読み込み
  int hzB = head[2];

  CDM::E_CDM_ARRAYSHAPE shape = buf->getArrayShape();

#ifdef CDM_BUFFER_MB_SIZE
  size_t ndata = src->getArrayLength();
  if( src->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
#else
  //NIJKの読込み
  if( shape == CDM::E_CDM_NIJK ) {
    for( int k=0; k<nz; k++ ) {
      //headインデクスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      //１層読み込み
      size_t ndata = buf->getArrayLength();
      if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;

      // コピー
      buf->copyArray(src);
    }
  }
  //IJKNの読込み 
  else if( shape == CDM::E_CDM_IJKN ) {
    for(int n=0; n<src->getNvari(); n++) {
    for(int k=0; k<nz; k++) {
      //headインデックスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      //１層読み込み
      size_t ndata = buf->getArrayLength();
      if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;

      //コピー
      buf->copyArrayNvari(src,n);
    }}
  }
#endif

  return CDM::E_CDM_SUCCESS;

}

#endif // _CDM_BOV_INLINE_H_
