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
 * @file   convOutput_SPH.C
 * @brief  convOutput_SPH Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_SPH.h"

// #################################################################
// コンストラクタ
convOutput_SPH::convOutput_SPH()
{


}

// #################################################################
// デストラクタ
convOutput_SPH::~convOutput_SPH()
{


}

// #################################################################
// 出力ファイルをオープンする。
FILE* convOutput_SPH::OutputFile_Open(
                                      const std::string prefix,
                                      const unsigned step,
                                      const int id,
                                      const bool mio)
{
  FILE* fp;

  //ファイル名の生成
  std::string outfile;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir()+"/"+
            cdm_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       "sph",
                                       fnameformat,
                                       mio,
                                       CDM::E_CDM_OFF);


  //ファイルオープン
  if( (fp = fopen(outfile.c_str(), "wb")) == NULL ) {
    printf("\tCan't open file.(%s)\n",outfile.c_str());
    Exit(0);
  }

  return fp;

}

// #################################################################
//
bool convOutput_SPH::WriteHeaderRecord(
                                       int step, 
                                       int dim, 
                                       CDM::E_CDM_DTYPE out_type, 
                                       int imax, 
                                       int jmax, 
                                       int kmax,
                                       double time, 
                                       double* org, 
                                       double* pit, 
                                       std::string prefix,
                                       FILE *fp)
{
  if( !fp ) return false;
 
  unsigned int dmy;

  //出力データ種別フラグの設定
  int sv_type;
  if( dim == 1 ) sv_type = SPH_SCALAR;
  else if( dim == 3 ) sv_type = SPH_VECTOR;

  //出力データ型フラグの設定
  int d_type;
  if( out_type == CDM::E_CDM_FLOAT32 ) d_type = SPH_FLOAT;
  else if( out_type == CDM::E_CDM_FLOAT64 ) d_type = SPH_DOUBLE;

  dmy = 2 * sizeof(int);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&sv_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&d_type, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( d_type == SPH_FLOAT ) {
    dmy = 3 * sizeof(int);
  } else if( d_type == SPH_DOUBLE ) {
    dmy = 3 * sizeof(long long);
  }
  
  //dmy = 3 * sizeof(long long);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    if( fwrite(&imax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&jmax, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&kmax, sizeof(int), 1, fp) != 1 ) return false;
  } else {
    if( fwrite(&imax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&jmax, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&kmax, sizeof(long long), 1, fp) != 1 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( d_type == SPH_FLOAT ) {
    dmy = 3 * sizeof(float);
  } else {
    dmy = 3 * sizeof(double);
  }
  //dmy = 3 * sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    float tmp[3];
    tmp[0] = (float)org[0];
    tmp[1] = (float)org[1];
    tmp[2] = (float)org[2];
    if( fwrite(tmp, sizeof(float), 3, fp) != 3 ) return false;
  } else {
    if( fwrite(org, sizeof(double), 3, fp) != 3 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    float tmp[3];
    tmp[0] = (float)pit[0];
    tmp[1] = (float)pit[1];
    tmp[2] = (float)pit[2];
    if( fwrite(tmp, sizeof(float), 3, fp) != 3 ) return false;
  } else {
    if( fwrite(pit, sizeof(double), 3, fp) != 3 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  if( d_type == SPH_FLOAT ) {
    dmy = sizeof(int) + sizeof(float);
  } else {
    dmy = sizeof(long long) + sizeof(double);
  }
  //dmy = sizeof(long long) + sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( d_type == SPH_FLOAT ) {
    float ftmp = (float)time;
    if( fwrite(&step, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(&ftmp, sizeof(float), 1, fp) != 1 ) return false;
  } else {
    if( fwrite(&step, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(&time, sizeof(double), 1, fp) != 1 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
 
  return true;
}

// #################################################################
//
bool convOutput_SPH::WriteDataMarker(int dmy, FILE* fp, bool out)
{
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  return true;
}
