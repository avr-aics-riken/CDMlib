/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_FilePath.C
 * @brief  cio_FilePath Class
 * @author aics    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cio_FilePath::cio_FilePath()
{
  ProcDFIFile="";
}

// #################################################################
// コンストラクタ
cio_FilePath::cio_FilePath(const std::string _ProcDFIFile)
{
  ProcDFIFile=_ProcDFIFile;
}

// #################################################################
// デストラクタ
cio_FilePath::~cio_FilePath()
{

}

// #################################################################
// proc.dfi ファイル名の読込み
CIO::E_CIO_ERRORCODE
cio_FilePath::Read(cio_TextParser tpCntl) 
{

  std::string str;
  std::string label;

  //Process
  label = "/FilePath/Process";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_FILEPATH_PROCESS;
  }
  ProcDFIFile=str;

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// proc.dfi ファイル名の出力
CIO::E_CIO_ERRORCODE
cio_FilePath::Write(FILE* fp, 
                    const unsigned tab)
{

  fprintf(fp, "FilePath {\n");
  fprintf(fp, "\n");

  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "Process = \"%s\"\n",ProcDFIFile.c_str());

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CIO::E_CIO_SUCCESS;

}

