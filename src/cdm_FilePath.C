/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_FilePath.C
 * @brief  cdm_FilePath Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cdm_FilePath::cdm_FilePath()
{
  ProcDFIFile="";
}

// #################################################################
// コンストラクタ
cdm_FilePath::cdm_FilePath(const std::string _ProcDFIFile)
{
  ProcDFIFile=_ProcDFIFile;
}

// #################################################################
// デストラクタ
cdm_FilePath::~cdm_FilePath()
{

}

// #################################################################
// proc.dfi ファイル名の読込み
CDM::E_CDM_ERRORCODE
cdm_FilePath::Read(cdm_TextParser tpCntl) 
{

  std::string str;
  std::string label;

  //Process
  label = "/FilePath/Process";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_FILEPATH_PROCESS;
  }
  ProcDFIFile=str;

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// proc.dfi ファイル名の出力
CDM::E_CDM_ERRORCODE
cdm_FilePath::Write(FILE* fp, 
                    const unsigned tab)
{

  fprintf(fp, "FilePath {\n");
  fprintf(fp, "\n");

  _CDM_WRITE_TAB(fp, tab);
  fprintf(fp, "Process = \"%s\"\n",ProcDFIFile.c_str());

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}

