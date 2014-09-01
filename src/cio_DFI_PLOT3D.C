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

#include "cio_DFI.h"
#include "cio_DFI_PLOT3D.h"

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

  //フィールドデータの配列サイズ取得
  const int *szVal = val->getArraySizeInt();

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
// GIRD データファイル出力コントロール
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
  for(int i=0; i<3; i++) sz[i] = DFI_Process.RankList[m_RankID].VoxelSize[i]+1;

  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) {
    float pit[3],org[3];
    for(int i=0; i<3; i++) {
      pit[i]=(float)DFI_Domain.GlobalRegion[i]/(float)DFI_Domain.GlobalVoxel[i];
      org[i]=(float)DFI_Domain.GlobalOrigin[i]-pit[i]*0.5;
    }
    //xyzを計算して出力
    write_XYZ(fp,org,pit,sz);
  }else if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) {
    double pit[3],org[3];
    for(int i=0; i<3; i++) {
      pit[i]=(double)DFI_Domain.GlobalRegion[i]/(double)DFI_Domain.GlobalVoxel[i];
      org[i]=(double)DFI_Domain.GlobalOrigin[i]-pit[i]*0.5;
    }
    //xyzを計算して出力
    write_XYZ(fp,org,pit,sz);
  }

  //file close
  fclose(fp);

  return true;

}

