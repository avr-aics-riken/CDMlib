/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   convOutput_NETCDF.C
 * @brief  convOutput_NETCDF Class
 * @author aics
 */
#ifdef _WITH_NETCDF4_

#include "convOutput.h"
#include "convOutput_NETCDF.h"

// #################################################################
// コンストラクタ
convOutput_NETCDF::convOutput_NETCDF()
{


}

// #################################################################
// デストラクタ
convOutput_NETCDF::~convOutput_NETCDF()
{


}

// #################################################################
// 出力ファイルをオープンする。
cdm_FILE* convOutput_NETCDF::OutputFile_Open(
                                      const std::string prefix,
                                      const unsigned step,
                                      const int id,
                                      const bool mio)
{
  cdm_FILE* pFile;

  //ファイル名の生成
  std::string outfile;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir()+"/"+
            cdm_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       D_CDM_EXT_NC,
                                       fnameformat,
                                       mio,
                                       CDM::E_CDM_OFF);

  // 追記モード
  bool addMode = CheckAddMode(step);

  //ファイルオープン
  if( (pFile = cdm_FILE::OpenWriteBinary(outfile, CDM::E_CDM_FMT_NETCDF4, addMode)) == NULL ) {
    printf("\tCan't open file.(%s)\n",outfile.c_str());
    Exit(0);
  }

  return pFile;

}

// #################################################################
// 追記モードを取得
bool
convOutput_NETCDF::CheckAddMode(const unsigned step)
{
  if( m_InputCntl->Get_OutputFilenameFormat() ==CDM::E_CDM_FNAME_RANK && step != m_pTSlice->SliceList[0].step )
  {
    return true;
  }
  return false;
}

// #################################################################
//
bool convOutput_NETCDF::WriteHeaderRecord(
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
  // varInfoのクリア
  m_varInfoX.clear();
  m_varInfoY.clear();
  m_varInfoZ.clear();
  m_varInfoT.clear();

  // 変数名
  m_varInfoX.name = "x";
  m_varInfoY.name = "y";
  m_varInfoZ.name = "z";
  m_varInfoT.name = "time";

  // セル座標値
  if( !m_varInfoX.data )
  {
    double *cod = new double[imax];
    for( int i=0;i<imax;i++ )
    {
      cod[i] = org[0] + pit[0]*(i+0.5);
    }
    m_varInfoX.data = cod;
  }
  if( !m_varInfoY.data )
  {
    double *cod = new double[jmax];
    for( int i=0;i<jmax;i++ )
    {
      cod[i] = org[1] + pit[1]*(i+0.5);
    }
    m_varInfoY.data = cod;
  }
  if( !m_varInfoZ.data )
  {
    double *cod = new double[kmax];
    for( int i=0;i<kmax;i++ )
    {
      cod[i] = org[2] + pit[2]*(i+0.5);
    }
    m_varInfoZ.data = cod;
  }

  // 出力用の情報
  int VoxelSize[3] ={imax, jmax, kmax};
  int GuideCell = 0;

//20160509,fub.s
//CDM::E_CDM_DTYPE DataType = m_InputCntl->Get_OutputDataType();
  CDM::E_CDM_DTYPE DataType = out_type;
//20160509.fub.e

  vector<string> VariableName;
  for( int i=0;i<m_pFinfo->VariableName.size();i++ )
  {
    VariableName.push_back(m_pFinfo->VariableName[i]);
  }
  cdm_Unit *pUnit, Unit;
  if( m_pUnit )
  {
    pUnit = (cdm_Unit*)m_pUnit;
  }
  else
  {
    pUnit = &Unit;
  }

  // 出力
  CDM::E_CDM_ERRORCODE ret = cdm_DFI_NETCDF::write_HeaderRecord(pFile, step, time, VoxelSize, GuideCell,
                                                                DataType, VariableName, *pUnit,
                                                                m_varInfoX, m_varInfoY, m_varInfoZ, m_varInfoT, m_vecVarInfo);

  return true;
}

// #################################################################
// data出力(成分ごと)
bool convOutput_NETCDF::WriteFieldData(cdm_FILE* pFile, 
                                       cdm_Array* src, 
                                       size_t dLen,
                                       CDM::E_CDM_DTYPE d_type,
                                       bool flag_variname,
                                       std::string variname)
{
  // ncid
  int ncid = pFile->m_ncid;

  // 成分の情報の取得
  cdm_DFI_NETCDF::stVarInfo *varInfo = NULL;
  for( int i=0;i<m_vecVarInfo.size();i++ )
  {
    if( variname == m_vecVarInfo[i].name )
    {
      varInfo = &m_vecVarInfo[i];
      break;
    }
  }

  // ncデータ型
  nc_type type = cdm_DFI_NETCDF::GetNcType( d_type );

  // 格子数
  const int *size = src->getArraySizeInt();

  // head
  const int *head = src->getHeadIndex();

  // timeのシフト量を取得
  int step = m_varInfoT.dims[0].len - 1;

  // start, count
  size_t start[4] = {step, head[2], head[1], head[0]};
  size_t count[4] = {1, size[2], size[1], size[0]};

  // 配列の出力
  int varid = varInfo->id;

  if( type == NC_BYTE )
  {
    signed char *ptr = (signed char*)src->getData();
    nc_put_vara_schar( ncid, varid, start, count, ptr );
  }
  else if( type == NC_SHORT )
  {
    short *ptr = (short*)src->getData();
    nc_put_vara_short( ncid, varid, start, count, ptr );
  }
  else if( type == NC_INT )
  {
    int *ptr = (int*)src->getData();
     nc_put_vara_int( ncid, varid, start, count, ptr );
  }
  else if( type == NC_FLOAT )
  {
    float *ptr = (float*)src->getData();
    nc_put_vara_float( ncid, varid, start, count, ptr );
  }
  else if( type == NC_DOUBLE )
  {
    double *ptr = (double*)src->getData();
    nc_put_vara_double( ncid, varid, start, count, ptr );
  }
  else if( type == NC_UBYTE )
  {
    unsigned char *ptr = (unsigned char*)src->getData();
    nc_put_vara_uchar( ncid, varid, start, count, ptr );
  }
  else if( type == NC_USHORT )
  {
    unsigned short *ptr = (unsigned short*)src->getData();
    nc_put_vara_ushort( ncid, varid, start, count, ptr );
  }
  else if( type == NC_UINT )
  {
    unsigned int *ptr = (unsigned int*)src->getData();
    nc_put_vara_uint( ncid, varid, start, count, ptr );
  }
  else if( type == NC_INT64 )
  {
    long long *ptr = (long long*)src->getData();
    nc_put_vara_longlong( ncid, varid, start, count, ptr );
  }
  else if( type == NC_UINT64 )
  {
    unsigned long long *ptr = (unsigned long long*)src->getData();
    nc_put_vara_ulonglong( ncid, varid, start, count, ptr );
  }

  return true;
}

#endif /*_WITH_NETCDF4_*/
