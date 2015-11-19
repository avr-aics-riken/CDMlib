/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_VisIt.C
 * @brief  cdm_VisIt Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K

// #################################################################
// コンストラクタ
cdm_VisIt::cdm_VisIt()
{
  PlotGC = "off";
}


// #################################################################
// コンストラクタ
cdm_VisIt::cdm_VisIt(const std::string _PlotGC)
{
  PlotGC = _PlotGC;
}


// #################################################################
// デストラクタ
cdm_VisIt::~cdm_VisIt()
{

}


// #################################################################
// オプションの読込み
CDM::E_CDM_ERRORCODE
cdm_VisIt::Read(cdm_TextParser tpCntl)
{

  std::string str;
  std::string label;

  //Process
  label = "/VisIt/PlotGuideCell";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_FILEPATH_PROCESS;
  }
  PlotGC = str;

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// proc.dfi ファイル名の出力
CDM::E_CDM_ERRORCODE
cdm_VisIt::Write(FILE* fp,
                    const unsigned tab)
{

  fprintf(fp, "VisIt {\n");
  fprintf(fp, "\n");

  _CDM_WRITE_TAB(fp, tab);
  fprintf(fp, "PlotGuideCell = \"%s\"\n", PlotGC.c_str());
  
  _CDM_WRITE_TAB(fp, tab);
  fprintf(fp, "ResultFormat  =  \"FBinary\"\n");

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}

