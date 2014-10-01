/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_PLOT3D.C
 * @brief  cdm_DFI_PLOT3D Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include "cdm_DFI_PLOT3D.h"

// #################################################################
// コンストラクタ
cdm_DFI_PLOT3D::cdm_DFI_PLOT3D()
{

}


// #################################################################
// デストラクタ
cdm_DFI_PLOT3D::~cdm_DFI_PLOT3D()
{

}

// #################################################################
// ファイルのヘッダーレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::read_HeaderRecord(FILE* fp,
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
// データレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::read_Datarecord(FILE* fp,
                                bool matchEndian,
                                cdm_Array* buf,
                                int head[3],
                                int nz,
                                cdm_Array* &src)
{
  //ngrid,nblock読込み
  int ngrid,szVal[3],ncomp;

  //ascii
  if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
    fscanf(fp,"%d\n",&ngrid);
    fscanf(fp,"%d%d%d%d\n",&szVal[0],&szVal[1],&szVal[2],&ncomp);
  //Fortran Binary
  } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
    unsigned int dmy;
    dmy = sizeof(int);
    fread(&dmy, sizeof(int), 1, fp);
    fread(&ngrid, sizeof(int), 1, fp);
    fread(&dmy, sizeof(int), 1, fp);
    dmy = sizeof(int)*4;
    fread(&dmy, sizeof(int), 1, fp);
    fread(&szVal[0], sizeof(int), 1, fp);
    fread(&szVal[1], sizeof(int), 1, fp);
    fread(&szVal[2], sizeof(int), 1, fp);
    fread(&ncomp, sizeof(int), 1, fp);
    fread(&dmy, sizeof(int), 1, fp);
  //Bunary
  } else {
    fread(&ngrid, sizeof(int), 1, fp);
    fread(&szVal[0], sizeof(int), 1, fp);
    fread(&szVal[1], sizeof(int), 1, fp);
    fread(&szVal[2], sizeof(int), 1, fp);
    fread(&ncomp, sizeof(int), 1, fp);
  }

  //フィールドデータ読込み
  CDM::E_CDM_ERRORCODE ret;

  if( buf->getDataType() == CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *dataS = dynamic_cast<cdm_TypeArray<float>*>(src);
    cdm_TypeArray<float> *dataB = dynamic_cast<cdm_TypeArray<float>*>(buf);
    ret = read_Func(fp, dataS, dataB, head, matchEndian);
  } else if( buf->getDataType() == CDM::E_CDM_FLOAT64 ) {
    cdm_TypeArray<double> *dataS = dynamic_cast<cdm_TypeArray<double>*>(src);
    cdm_TypeArray<double> *dataB = dynamic_cast<cdm_TypeArray<double>*>(buf);
    ret = read_Func(fp, dataS, dataB, head, matchEndian);
  }

  return ret;
}

// #################################################################
// Averaged レコードの読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::read_averaged(FILE* fp,
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
// ヘッダーレコード出力
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::write_HeaderRecord(FILE* fp,
                                   const unsigned step,
                                   const double time,
                                   const int n)
{
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// データレコード出力
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::write_DataRecord(FILE* fp, 
                                 cdm_Array* val, 
                                 const int gc, 
                                 const int n)
{

  //GRID データファイル出力処理
  if( m_OutputGrid == true ) {
    write_GridData();
    m_OutputGrid = false;
  }

  //フィールドデータの配列サイズ取得(ガイドセル含む)
  const int *szVal = val->_getArraySizeInt();

  //配列成分の取得
  int ncomp = val->getNcomp();

  //printf("ID : %d prefix : %s size : %d %d %d ncomp : %d\n",n,
  //        DFI_Finfo.Prefix.c_str(),szVal[0],szVal[1],szVal[2],
  //        ncomp);

  //ngrid,nblock出力
  int ngrid=1;
  //ascii
  if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
    fprintf(fp,"%5d\n",ngrid);
    fprintf(fp,"%5d%5d%5d%5d\n",szVal[0],szVal[1],szVal[2],ncomp);
  //Fortran Binary
  } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
    unsigned int dmy;
    dmy = sizeof(int);
    fwrite(&dmy, sizeof(int), 1, fp);
    fwrite(&ngrid, sizeof(int), 1, fp);
    fwrite(&dmy, sizeof(int), 1, fp);

    dmy = sizeof(int)*4;
    fwrite(&dmy, sizeof(int), 1, fp);
    fwrite(&szVal[0], sizeof(int), 1, fp);
    fwrite(&szVal[1], sizeof(int), 1, fp);
    fwrite(&szVal[2], sizeof(int), 1, fp);
    fwrite(&ncomp, sizeof(int), 1, fp);
    fwrite(&dmy, sizeof(int), 1, fp);
  //Bunary
  } else {
    fwrite(&ngrid, sizeof(int), 1, fp);
    fwrite(&szVal[0], sizeof(int), 1, fp);
    fwrite(&szVal[1], sizeof(int), 1, fp);
    fwrite(&szVal[2], sizeof(int), 1, fp);
    fwrite(&ncomp, sizeof(int), 1, fp);
  }
  //フィールドデータ出力

  if( val->getDataType() == CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *data = dynamic_cast<cdm_TypeArray<float>*>(val);
    write_Func(fp, data, szVal, ncomp);
  } else if( val->getDataType() == CDM::E_CDM_FLOAT64 ) {
    cdm_TypeArray<double> *data = dynamic_cast<cdm_TypeArray<double>*>(val);
    write_Func(fp, data, szVal, ncomp);
  }   

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// 平均の出力PLOT3Dは何も出力しない
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::write_averaged(FILE* fp,
                            const unsigned step_avr,
                            const double time_avr)
{
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// GRID データファイル出力コントロール
bool
cdm_DFI_PLOT3D::write_GridData()
{

  bool mio = false;
  if( DFI_MPI.NumberOfRank > 1 ) mio = true;

  //出力ファイル名の生成
  std::string fname,tmp;
  tmp = Generate_FileName(DFI_Finfo.Prefix,m_RankID,-1,"xyz",m_output_fname,mio,
                          DFI_Finfo.TimeSliceDirFlag);
  if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
    fname = DFI_Finfo.DirectoryPath + "/" + tmp;
  } else {
    fname = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
  }

  //GRID data file open
  FILE* fp=NULL;
  if( (fp = fopen(fname.c_str(),"w"))  == NULL ) {
    printf("\tCan't open file.(%s)\n",fname.c_str());
    return false;
  }

  //xyzを求めて出力
  int sz[3];
  for(int i=0; i<3; i++) sz[i] = DFI_Process.RankList[m_RankID].VoxelSize[i];
  if( DFI_Finfo.GuideCell>1 ) for(int i=0; i<3; i++) sz[i] = sz[i]+2*DFI_Finfo.GuideCell;

  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) {
    float pit[3],org[3];
    for(int i=0; i<3; i++) {
      pit[i]=(float)DFI_Domain.GlobalRegion[i]/(float)DFI_Domain.GlobalVoxel[i];
      org[i]=(float)DFI_Domain.GlobalOrigin[i]+pit[i]*0.5;
      if( DFI_Finfo.GuideCell>1 ) org[i]=org[i]-pit[i]*(float)DFI_Finfo.GuideCell;
    }
    //xyzを計算して出力
    write_XYZ(fp,org,pit,sz);
  }else if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) {
    double pit[3],org[3];
    for(int i=0; i<3; i++) {
      pit[i]=(double)DFI_Domain.GlobalRegion[i]/(double)DFI_Domain.GlobalVoxel[i];
      org[i]=(double)DFI_Domain.GlobalOrigin[i]+pit[i]*0.5;
      if( DFI_Finfo.GuideCell>1 ) org[i]=org[i]-pit[i]*(double)DFI_Finfo.GuideCell;
    }
    //xyzを計算して出力
    write_XYZ(fp,org,pit,sz);
  }

  //file close
  fclose(fp);

  return true;

}

