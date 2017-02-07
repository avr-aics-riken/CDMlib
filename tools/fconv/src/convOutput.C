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
 * @file   convOutput.C
 * @brief  convOutput Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_SPH.h"
#include "convOutput_BOV.h"
#include "convOutput_AVS.h"
#include "convOutput_VTK.h"
#include "convOutput_PLOT3D.h"
#include "convOutput_NETCDF.h"
#include "convOutput_FUB.h"

// #################################################################
// コンストラクタ
convOutput::convOutput()
{
  m_pTSlice = NULL;
  m_pFinfo = NULL;
  m_pUnit = NULL;
}

// #################################################################
// デストラクタ
convOutput::~convOutput()
{


}

// #################################################################
// InputParamのポインタをコピー
bool convOutput::importInputParam(InputParam* InputCntl)
{
  if( !InputCntl ) return false;
  m_InputCntl = InputCntl;
  return true;
}

//#################################################################
// convOutput インスタンス
convOutput*
convOutput::OutputInit(const CDM::E_CDM_FORMAT out_format)
{

  convOutput *OutConv = NULL;

  if     ( out_format == CDM::E_CDM_FMT_SPH    ) OutConv = new convOutput_SPH();
  else if( out_format == CDM::E_CDM_FMT_BOV    ) OutConv = new convOutput_BOV();
  else if( out_format == CDM::E_CDM_FMT_AVS    ) OutConv = new convOutput_AVS();
  else if( out_format == CDM::E_CDM_FMT_VTK    ) OutConv = new convOutput_VTK();
  else if( out_format == CDM::E_CDM_FMT_PLOT3D ) OutConv = new convOutput_PLOT3D();
#ifdef _WITH_NETCDF4_
  else if( out_format == CDM::E_CDM_FMT_NETCDF4 ) OutConv = new convOutput_NETCDF();
#endif
//20160411.fub.s
  else if( out_format == CDM::E_CDM_FMT_FUB ) OutConv = new convOutput_FUB();
//20160411.fub.e

  return OutConv;

}

//#################################################################
//
bool convOutput::WriteFieldData(cdm_FILE* pFile, cdm_Array* src, size_t dLen)
{
  FILE *fp = pFile->m_fp;

  if( src->writeBinary(fp) != dLen ) Exit(0);
  //if( src->writeAscii(fp) != dLen ) Exit(0);
  return true;
}

