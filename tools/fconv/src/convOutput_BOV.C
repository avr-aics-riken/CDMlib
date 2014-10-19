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
 * @file   convOutput_BOV.C
 * @brief  convOutput_BOV Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_BOV.h"

// #################################################################
// コンストラクタ
convOutput_BOV::convOutput_BOV()
{


}

// #################################################################
// デストラクタ
convOutput_BOV::~convOutput_BOV()
{


}

// #################################################################
// 出力ファイルをオープンする。
FILE* convOutput_BOV::OutputFile_Open(
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
                                       "dat",
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
bool convOutput_BOV::WriteHeaderRecord(int step,
                                        int dim,
                                        CDM::E_CDM_DTYPE d_type,
                                        int imax,
                                        int jmax,
                                        int kmax,
                                        double time,
                                        double* org,
                                        double* pit,
                                        std::string prefix,
                                        FILE *fdumy)
{

  FILE* fp;

  //HeadFile名の生成とオープン
  std::string head_file;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  head_file = m_InputCntl->Get_OutputDir()+"/"+
            cdm_DFI::Generate_FileName(prefix,
                                       0,
                                       step,
                                       "bov",
                                       fnameformat,
                                       false,
                                       CDM::E_CDM_OFF);

  if( (fp = fopen(head_file.c_str(), "wb")) == NULL ) {
    printf("\tCan't open file.(%s)\n",head_file.c_str());
    Exit(0);
  }

  //データファイル名の生成
  std::string data_file;
  data_file = cdm_DFI::Generate_FileName(prefix,
                                       0,
                                       step,
                                       "dat",
                                       fnameformat,
                                       false,
                                       CDM::E_CDM_OFF);

  fprintf( fp, "Time: %e\n",time);
  fprintf( fp, "DATA_FILE: %s\n",data_file.c_str());
  fprintf( fp, "DATA_SIZE: %d %d %d\n",imax,jmax,kmax);
  std::string dataType;
  if     ( d_type == CDM::E_CDM_INT8    ) dataType = D_CDM_BYTE;
  else if( d_type == CDM::E_CDM_UINT8   ) dataType = D_CDM_UINT8;
  else if( d_type == CDM::E_CDM_INT16   ) dataType = D_CDM_INT16;
  else if( d_type == CDM::E_CDM_UINT16  ) dataType = D_CDM_UINT16;
  else if( d_type == CDM::E_CDM_INT32   ) dataType = D_CDM_INT;
  else if( d_type == CDM::E_CDM_UINT32  ) dataType = D_CDM_UINT32;
  else if( d_type == CDM::E_CDM_INT64   ) dataType = D_CDM_INT64;
  else if( d_type == CDM::E_CDM_UINT64  ) dataType = D_CDM_UINT64;
  else if( d_type == CDM::E_CDM_FLOAT32 ) dataType = D_CDM_FLOAT;
  else if( d_type == CDM::E_CDM_FLOAT64 ) dataType = D_CDM_DOUBLE;
  fprintf( fp, "DATA_FORMAT: %s\n",dataType.c_str());
  fprintf( fp, "DATA_COMPONENTS: %d\n",dim);
  fprintf( fp, "VARIABLE: %s\n",prefix.c_str());

  std::string endian;
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  if     ( cdumy[0] == 0x01 ) endian = "LITTLE";
  else if( cdumy[0] == 0x00 ) endian = "BIG";
  fprintf( fp, "DATA_ENDIAN: %s\n",endian.c_str());
  fprintf( fp, "CENTERING: zonal\n");
  //fprintf( fp, "BRICK_ORIGIN: %e %e %e\n",org[0]-0.5*pit[0],org[1]-0.5*pit[1],org[2]-0.5*pit[2]); 
  fprintf( fp, "BRICK_ORIGIN: %e %e %e\n",org[0],org[1],org[2] ); 
  fprintf( fp, "BRICK_SIZE: %e %e %e\n",pit[0]*imax,pit[1]*jmax,pit[2]*kmax);

  std::string arrayShape;
  if( m_InputCntl->Get_OutputArrayShape() == CDM::E_CDM_IJKN ) arrayShape="IJKN";
  else arrayShape="NIJK";
  fprintf( fp, "#CDM_ARRAY_SHAPE: %s\n",arrayShape.c_str());

  fclose(fp);
  return true;
}
