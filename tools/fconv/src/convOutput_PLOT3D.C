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
 * @file   convOutput_PLOT3D.C
 * @brief  convOutput_PLOT3D Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_PLOT3D.h"

// #################################################################
// コンストラクタ
convOutput_PLOT3D::convOutput_PLOT3D()
{


}

// #################################################################
// デストラクタ
convOutput_PLOT3D::~convOutput_PLOT3D()
{


}

// #################################################################
// grid 出力
void convOutput_PLOT3D::WriteGridData(std::string prefix,
                                      int step,
                                      int myRank,
                                      int dType,
                                      int guide,
                                      double org[3],
                                      double pit[3],
                                      int sz[3])
{

  //step 0 以外は出力しない
  if( step != 0 ) return; 

  int id,jd,kd;
  id=sz[0]; //+1;
  jd=sz[1]; //+1;
  kd=sz[2]; //+1;
  //格子点への補間は行わず、双対セルとして扱うため、+1は不要
  size_t maxsize = (size_t)id*(size_t)jd*(size_t)kd*3;
  size_t outsize = (size_t)id*(size_t)jd*(size_t)kd;

  int outDtype = m_InputCntl->Get_OutputDataType();
  if( outDtype == CDM::E_CDM_DTYPE_UNKNOWN ) outDtype = dType;

  if( outDtype == CDM::E_CDM_FLOAT64 ) {
    double* x = new double[maxsize];
    OutputPlot3D_xyz(prefix, step, myRank, guide, org, pit, sz, &x[0], 
                     &x[outsize] , &x[outsize*2] );
  }else if( outDtype == CDM::E_CDM_FLOAT32 ) {
    float* x = new float[maxsize];
    OutputPlot3D_xyz(prefix, step, myRank, guide, org, pit, sz, &x[0], 
                     &x[outsize] , &x[outsize*2] );
  }
}

// #################################################################
// grid 出力 (不等間隔格子対応版)
void convOutput_PLOT3D::WriteGridData(std::string prefix,
                                      int step,
                                      int myRank,
                                      int dType,
                                      int guide,
                                      cdm_Domain* out_domain,
                                      cdm_Process* out_process)
{

  //step 0 以外は出力しない
  if( step != 0 ) return; 

  int outDtype = m_InputCntl->Get_OutputDataType();
  if( outDtype == CDM::E_CDM_DTYPE_UNKNOWN ) outDtype = dType;

  if( outDtype == CDM::E_CDM_FLOAT32 ) {
    OutputPlot3D_xyz<float>(prefix, step, myRank, guide, out_domain, out_process);
  }else if( outDtype == CDM::E_CDM_FLOAT64 ) {
    OutputPlot3D_xyz<double>(prefix, step, myRank, guide, out_domain, out_process);
  }

}

// #################################################################
// グリッド数の 出力
void convOutput_PLOT3D::WriteNgrid(cdm_FILE* pFile, int ngrid)
{
  FILE *fp = pFile->m_fp;

  switch (m_InputCntl->Get_OutputFileType()) {
    case CDM::E_CDM_FILE_TYPE_FBINARY:
      unsigned int dmy;
      dmy = sizeof(int);
      WriteDataMarker(dmy,pFile,true);
      fwrite(&ngrid, sizeof(int), 1, fp);
      WriteDataMarker(dmy,pFile,true);
      break;
    case CDM::E_CDM_FILE_TYPE_ASCII:
      fprintf(fp,"%5d\n",ngrid);
      break;
    case CDM::E_CDM_FILE_TYPE_BINARY:
      fwrite(&ngrid, sizeof(int), 1, fp);
      break;      
    default:
      break;
    }
}

// #################################################################
// ブロックデータの 出力
void convOutput_PLOT3D::WriteBlockData(cdm_FILE* pFile, int id, int jd, int kd)
{
  FILE *fp = pFile->m_fp;

  switch (m_InputCntl->Get_OutputFileType()) {
    case CDM::E_CDM_FILE_TYPE_FBINARY:
      unsigned int dmy;
      dmy = sizeof(int)*3;
      WriteDataMarker(dmy,pFile,true);
      fwrite(&id, sizeof(int), 1, fp);
      fwrite(&jd, sizeof(int), 1, fp);
      fwrite(&kd, sizeof(int), 1, fp);
      WriteDataMarker(dmy,pFile,true);
      break;
    case CDM::E_CDM_FILE_TYPE_ASCII:
      fprintf(fp,"%5d%5d%5d\n",id,jd,kd);
      break;
    case CDM::E_CDM_FILE_TYPE_BINARY:
      fwrite(&id, sizeof(int), 1, fp);
      fwrite(&jd, sizeof(int), 1, fp);
      fwrite(&kd, sizeof(int), 1, fp);
      break;
   default:
      break;
   }      

}

// #################################################################
// フィールドデータ出力ファイルオープン
cdm_FILE* convOutput_PLOT3D::OutputFile_Open(
                         const std::string prefix,
                         const unsigned step,
                         const int id,
                         const bool mio)
{
  cdm_FILE* pFile;

  //ファイル名の生成
  std::string outfile;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir() +"/"+
            cdm_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       "fun",
                                       fnameformat,
                                       mio,
                                       CDM::E_CDM_OFF);

  //出力ファイルオープン
  // ascii
  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {
    if( (pFile = cdm_FILE::OpenWriteAscii(outfile, CDM::E_CDM_FMT_PLOT3D)) == NULL ) {
      printf("\tCan't open file.(%s)\n",outfile.c_str());
      Exit(0);
    }
  } else {
  //binary
    if( (pFile = cdm_FILE::OpenWriteBinary(outfile, CDM::E_CDM_FMT_PLOT3D)) == NULL ) {
      printf("\tCan't open file.(%s)\n",outfile.c_str());
      Exit(0);
    }
  }

  return pFile;

}

// #################################################################
// ヘッダーレコードを出力
bool convOutput_PLOT3D::WriteHeaderRecord(
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
                                       cdm_FILE *pFile)
{
  if( !pFile ) return false;
  if( !pFile->m_fp ) return false;

  //ngirdの初期化
  int ngrid=1;

  //ngrid の出力
  //WriteNgrid(fp,ngrid);

  //block data の出力
  
  WriteFuncBlockData(pFile,imax,jmax,kmax,dim);
  //WriteFuncBlockData(fp,imax+1,jmax+1,kmax+1,dim);
  //格子点への補間は行わず、双対セルとして扱うため、+1は不要

  return true;

}



// #################################################################
// func 出力
bool convOutput_PLOT3D::WriteFieldData(cdm_FILE* pFile, 
                                       cdm_Array* src, 
                                       size_t dLen)
{

  const int* sz = src->getArraySizeInt();

  cdm_Array *out = cdm_Array::instanceArray
                   (src->getDataType(),
                    src->getArrayShape(),
                    (int *)sz,
                    0,
                    src->getNvari());

  int ret = src->copyArray(out);

  WriteFuncData(pFile, out);
  delete out;

  return true;

}

// #################################################################
// func ブロックデータの 出力
void convOutput_PLOT3D::WriteFuncBlockData(cdm_FILE* pFile, int id, int jd, int kd, int nvar)
{
  FILE *fp = pFile->m_fp;

  switch (m_InputCntl->Get_OutputFileType()) {
    case CDM::E_CDM_FILE_TYPE_FBINARY:
      unsigned int dmy;
      dmy = sizeof(int)*4;
      WriteDataMarker(dmy,pFile,true);
      fwrite(&id, sizeof(int), 1, fp);
      fwrite(&jd, sizeof(int), 1, fp);
      fwrite(&kd, sizeof(int), 1, fp);
      fwrite(&nvar, sizeof(int), 1, fp);
      WriteDataMarker(dmy,pFile,true);
      break;
    case CDM::E_CDM_FILE_TYPE_ASCII:
      fprintf(fp,"%5d%5d%5d%5d\n",id,jd,kd,nvar);
      break;
    case CDM::E_CDM_FILE_TYPE_BINARY:
      fwrite(&id, sizeof(int), 1, fp);
      fwrite(&jd, sizeof(int), 1, fp);
      fwrite(&kd, sizeof(int), 1, fp);
      fwrite(&nvar, sizeof(int), 1, fp);
      break;
    default:
      break;
  }

}

// #################################################################
void convOutput_PLOT3D::WriteFuncData(cdm_FILE* pFile, cdm_Array* p3src)
{
  FILE *fp = pFile->m_fp;

  const int* sz = p3src->getArraySizeInt();
  size_t dLen = (size_t)sz[0]*(size_t)sz[1]*(size_t)sz[2]*p3src->getNvari();

  switch (m_InputCntl->Get_OutputFileType()) {
    case CDM::E_CDM_FILE_TYPE_FBINARY:
      p3src->writeBinary(fp);      
      break;
    case CDM::E_CDM_FILE_TYPE_ASCII:
      if( p3src->getDataType() == CDM::E_CDM_FLOAT32) {
        float *data = (float*)p3src->getData();
        for(int i=0; i<dLen; i++) fprintf(fp,"%15.6E\n",data[i]);
      } else if( p3src->getDataType() == CDM::E_CDM_FLOAT64) {
        double *data = (double*)p3src->getData();
        for(int i=0; i<dLen; i++) fprintf(fp,"%15.6E\n",data[i]);
      }
      break;
    case CDM::E_CDM_FILE_TYPE_BINARY:
      p3src->writeBinary(fp);
      break;
    default:
      break;
  }
    
}

// #################################################################
//
bool convOutput_PLOT3D::WriteDataMarker(int dmy, cdm_FILE* pFile, bool out)
{
  FILE *fp = pFile->m_fp;
  if( !out ) return true;
  if( m_InputCntl->Get_OutputFileType() != CDM::E_CDM_FILE_TYPE_FBINARY ) return true;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  return true;
}
                                         
