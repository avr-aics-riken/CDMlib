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
 * @file   convOutput_FUB.C
 * @brief  convOutput_FUB Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_FUB.h"

// #################################################################
// コンストラクタ
convOutput_FUB::convOutput_FUB()
{


}

// #################################################################
// デストラクタ
convOutput_FUB::~convOutput_FUB()
{


}

// #################################################################
// 出力ファイルをオープンする。
cdm_FILE* convOutput_FUB::OutputFile_Open(
                                      const std::string prefix,
                                      const unsigned step,
                                      const int id,
                                      const bool mio)
{
  cdm_FILE* pFile;

  std::string fmt;
  if( m_pFinfo->FileFormat == CDM::E_CDM_FMT_FUB ) fmt = "fub";
  else  fmt = "xyz";

  //ファイル名の生成
  std::string outfile;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir()+"/"+
            cdm_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       //"fub",
                                       fmt,
                                       fnameformat,
                                       mio,
                                       CDM::E_CDM_OFF);


  //ファイルオープン
  if( (pFile = cdm_FILE::OpenWriteBinary(outfile, CDM::E_CDM_FMT_FUB)) == NULL ) {
    printf("\tCan't open file.(%s)\n",outfile.c_str());
    Exit(0);
  }

  return pFile;

}

// #################################################################
//
bool convOutput_FUB::WriteHeaderRecord(
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
  FILE *fp = pFile->m_fp;
  if( !fp ) return false;
  return true;
}

// #################################################################
//
bool convOutput_FUB::WriteDataMarker(int dmy, cdm_FILE* pFile, bool out)
{
  FILE *fp = pFile->m_fp;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  return true;
}
