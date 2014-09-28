/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_BOV.C
 * @brief  cdm_DFI_BOV Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include "cdm_DFI_BOV.h"

// #################################################################
// コンストラクタ
cdm_DFI_BOV::cdm_DFI_BOV()
{

}


// #################################################################
// デストラクタ
cdm_DFI_BOV::~cdm_DFI_BOV()
{

}

// #################################################################
// ファイルのヘッダーレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_BOV::read_HeaderRecord(FILE* fp,
                               bool matchEndian,
                               unsigned step,
                               const int head[3],
                               const int tail[3],
                               int gc,
                               int voxsize[3],
                               double &time)
{

  time=0.0;
  for(int i=0; i<DFI_TimeSlice.SliceList.size(); i++) {
     if( DFI_TimeSlice.SliceList[i].step == step ) {
       time=(double)DFI_TimeSlice.SliceList[i].time;
     }
  }

  for(int i=0; i<3; i++) voxsize[i]=tail[i]-head[i]+1+(2*gc);

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ファイルのデータレコード読込み
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

  //NIJKの読込み
  if( shape == CDM::E_CDM_NIJK ) {
    for( int k=0; k<nz; k++ ) {
      //headインデクスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      //１層読み込
      size_t ndata = buf->getArrayLength();
      if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;

      // コピー
      buf->copyArray(src);
    }
  }
  //IJKNの読込み 
  else if( shape == CDM::E_CDM_IJKN ) {
    for(int n=0; n<src->getNcomp(); n++) {
    for(int k=0; k<nz; k++) {
      //headインデックスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      //１層読み込
      size_t ndata = buf->getArrayLength();
      if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;

      //コピー
      buf->copyArrayNcomp(src,n);
    }}
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// Averaged レコードの読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_BOV::read_averaged(FILE* fp,
                           bool matchEndian,
                           unsigned step,
                           unsigned &step_avr,
                           double &time_avr)
{

  step_avr=0;
  time_avr=0.0;

  for(int i=0; i<DFI_TimeSlice.SliceList.size(); i++) {
     if( DFI_TimeSlice.SliceList[i].step == step ) {
       step_avr=(int)DFI_TimeSlice.SliceList[i].AveragedStep;
       time_avr=(double)DFI_TimeSlice.SliceList[i].AveragedTime;
     }
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ヘッダーレコード出力BOVは何も出力しない
CDM::E_CDM_ERRORCODE
cdm_DFI_BOV::write_HeaderRecord(FILE* fp,
                                const unsigned step,
                                const double time,
                                const int n)
{
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// BOVデータレコード出力
CDM::E_CDM_ERRORCODE
cdm_DFI_BOV::write_DataRecord(FILE* fp, 
                              cdm_Array* val, 
                              const int gc, 
                              const int n)
{

  CDM::E_CDM_DTYPE Dtype = (CDM::E_CDM_DTYPE)DFI_Finfo.DataType;
  int Real_size = get_cdm_Datasize(Dtype);

  int size[3];
  for(int i=0; i<3; i++ ) size[i] = (int)DFI_Process.RankList[n].VoxelSize[i]+(int)(2*gc);

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.Component > 1 ) dLen *= 3;

  unsigned int dmy = dLen * Real_size;

  if( val->writeBinary(fp) != dLen ) return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// 平均の出力BOVは何も出力しない
CDM::E_CDM_ERRORCODE
cdm_DFI_BOV::write_averaged(FILE* fp,
                            const unsigned step_avr,
                            const double time_avr)
{
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// BOV ascii header file output
bool
cdm_DFI_BOV::write_ascii_header(const unsigned step,
                                const double time)
{

  FILE* fp=NULL;

  //ファイル名生成
  bool mio=false;
  if( DFI_MPI.NumberOfRank > 1 ) mio = true;

  std::string fname,tmp;
  tmp = Generate_FileName(DFI_Finfo.Prefix,
                          m_RankID,
                          step,
                          "bov",
                          m_output_fname,
                          mio,
                          DFI_Finfo.TimeSliceDirFlag);

  if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
    fname = DFI_Finfo.DirectoryPath + "/" + tmp;
  } else {
    fname = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
  }

  //bov ヘッダーファイルオープン
  if( (fp = fopen(fname.c_str(),"w")) == NULL ) {
    printf("\tCan't open file.(%s)\n",fname.c_str());
    return false;
  }

  //TIME: 
  fprintf(fp,"Time: %e\n",time);

  //DATA_FILE:
  std::string o_fname;
  o_fname = Generate_FileName(DFI_Finfo.Prefix,
                              m_RankID,
                              step,
                              "dat",
                              m_output_fname,
                              mio,
                              DFI_Finfo.TimeSliceDirFlag);
  fprintf(fp,"DATA_FILE: %s\n",o_fname.c_str());

  //DATA_SIZE:
  fprintf(fp,"DATA_SIZE: %d %d %d\n",DFI_Process.RankList[m_RankID].VoxelSize[0],
                                     DFI_Process.RankList[m_RankID].VoxelSize[1],
                                     DFI_Process.RankList[m_RankID].VoxelSize[2]);

  //DATA_FORMAT
  std::string dType;
  if(      GetDataType() == CDM::E_CDM_INT8    ) dType=D_CDM_BYTE;
  else if( GetDataType() == CDM::E_CDM_UINT8   ) dType=D_CDM_UINT8;
  else if( GetDataType() == CDM::E_CDM_INT16   ) dType=D_CDM_INT16;
  else if( GetDataType() == CDM::E_CDM_UINT16  ) dType=D_CDM_UINT16;
  else if( GetDataType() == CDM::E_CDM_INT32   ) dType=D_CDM_INT;
  else if( GetDataType() == CDM::E_CDM_UINT32  ) dType=D_CDM_UINT32;
  else if( GetDataType() == CDM::E_CDM_INT64   ) dType=D_CDM_INT64;
  else if( GetDataType() == CDM::E_CDM_UINT64  ) dType=D_CDM_UINT64;
  else if( GetDataType() == CDM::E_CDM_FLOAT32 ) dType=D_CDM_FLOAT;
  else if( GetDataType() == CDM::E_CDM_FLOAT64 ) dType=D_CDM_DOUBLE;
  fprintf(fp,"DATA_FORMAT: %s\n",dType.c_str());

  //DATA_COMPONENT
  fprintf(fp,"DATA_COMPONENT: %d\n",DFI_Finfo.Component);

  //VARIABLE:
  fprintf(fp,"VARIABLE: %s\n",DFI_Finfo.Prefix.c_str());

  //DATA_ENDIAN
  if( DFI_Finfo.Endian == CDM::E_CDM_LITTLE ) {
    fprintf(fp,"DATA_ENDIAN: LITTLE\n");
  } else {
    fprintf(fp,"DATA_ENDIAN: BIG\n");
  }

  //CENTERING
  fprintf(fp,"CENTERING: zonal\n");

  //pchを計算
  double pch[3];
  for(int i=0; i<3; i++) {
    pch[i]=(DFI_Domain->GlobalRegion[i]/DFI_Domain->GlobalVoxel[i]);
  }

  //BRICK_ORIGN
  double org[3];
  for(int i=0; i<3; i++) org[i]=DFI_Domain->GlobalOrigin[i]+0.5*pch[i];
  if( DFI_Finfo.GuideCell>1 ) for(int i=0; i<3; i++) org[i]=org[i]-pch[i]*(double)DFI_Finfo.GuideCell;
  /*
  fprintf(fp,"BRICK_ORIGN: %e %e %e\n",DFI_Domain.GlobalOrigin[0],
                                       DFI_Domain.GlobalOrigin[1],
                                       DFI_Domain.GlobalOrigin[2]);
  */
  fprintf(fp,"BRICK_ORIGN: %e %e %e\n",org[0],org[1],org[2]);

  //BRICK_SIZE
  fprintf(fp,"BRICK_SIZE: %e %e %e\n",
          DFI_Process.RankList[m_RankID].VoxelSize[0]*pch[0],
          DFI_Process.RankList[m_RankID].VoxelSize[1]*pch[1],
          DFI_Process.RankList[m_RankID].VoxelSize[2]*pch[2]);

  //#CDM_ARRAY_SHAPE
  if( DFI_Finfo.ArrayShape == CDM::E_CDM_IJKN ) {
    fprintf(fp,"#CDM_ARRAY_SHAPE: IJKN\n");
  } else {
    fprintf(fp,"#CDM_ARRAY_SHAPE: NIJK\n");
  }

  //file close
  fclose(fp);

  return true;
}

