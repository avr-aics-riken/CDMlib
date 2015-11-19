#ifndef _CDM_FILE_H_
#define _CDM_FILE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_FILE.h
 * @brief  cdm_FILE Class Header
 * @author aics
 */

#include "cdm_Define.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string>

#ifdef _WITH_NETCDF4_
#include "netcdf.h"
#endif

/** CDM file pointer control class */
class cdm_FILE
{

public:
  FILE* m_fp;                 ///< file poiter
  int   m_ncid;               ///< netcdf file ID
  CDM::E_CDM_FORMAT m_format; ///< file type
  bool m_addMode;             ///< 追記モードフラグ(true:追記モード)
  string m_fname;             ///< ファイル名

protected :
  /** コンストラクタ */
  cdm_FILE()
  {
    m_fp = NULL;
    m_ncid = -1;
    m_format = CDM::E_CDM_FMT_UNKNOWN;
    m_addMode = false;
    m_fname = "";
  }
  
  /** デストラクタ */
  virtual ~cdm_FILE()
  {
    if( m_format == CDM::E_CDM_FMT_NETCDF4 )
    {
#ifdef _WITH_NETCDF4_
      nc_close(m_ncid);
#endif
    }
    else
    {
      fclose(m_fp);
    }
  }

public:
  /**
   * @brief 読み込みファイルをオープン
   * @param [in]  fname  ファイル名
   * @param [in]  format ファイルタイプ
   * @return インスタンスされたクラスのポインタ
   */
  static cdm_FILE*
  OpenReadBinary(const std::string fname,
                 CDM::E_CDM_FORMAT format)
  {
    int ncid = -1;
    FILE *fp = NULL;

    if( format == CDM::E_CDM_FMT_NETCDF4 )
    {
#ifdef _WITH_NETCDF4_
      int iret = nc_open( fname.c_str(), NC_NOWRITE, &ncid );
      if( iret != NC_NOERR )
      {
        fprintf(stderr, "NC error code = %d\n", iret);
        return NULL;
      }
#else
      return NULL;
#endif
    }
    else
    {
      fp = fopen( fname.c_str(),"rb" );
      if( !fp )
      {
        return NULL;
      }
    }

    cdm_FILE *pFile = new cdm_FILE();
    pFile->m_fp   = fp;
    pFile->m_ncid = ncid;
    pFile->m_format = format;
    pFile->m_fname = fname;

    return pFile;
  }

  /**
   * @brief 出力ファイルをオープン(binaryモード)
   * @param [in]  fname   ファイル名
   * @param [in]  format  ファイルタイプ
   * @param [in]  addmode 追記モードフラグ(true:追記)
   * @return インスタンスされたクラスのポインタ
   */
  static cdm_FILE*
  OpenWriteBinary(const std::string fname,
                  CDM::E_CDM_FORMAT format,
                  bool addMode=false)
  {
    int ncid = -1;
    FILE *fp = NULL;

    if( format == CDM::E_CDM_FMT_NETCDF4 )
    {
#ifdef _WITH_NETCDF4_
      int iret;
      if( addMode )
      {
//stmpd_printf("open write file [%s] with addMode\n", fname.c_str());
        iret = nc_open( fname.c_str(), NC_WRITE|NC_SHARE, &ncid );
      }
      else
      {
//stmpd_printf("open write file [%s] with createMode\n", fname.c_str());
//      iret = nc_create( fname.c_str(), NC_CLOBBER, &ncid );
        iret = nc_create( fname.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid );
      }
      if( iret != NC_NOERR )
      {
        fprintf(stderr, "NC error code = %d\n", iret);
        return NULL;
      }
#else
      return NULL;
#endif
    }
    else
    {
      fp = fopen( fname.c_str(),"wb" );
      if( !fp )
      {
        return NULL;
      }
    }

    cdm_FILE *pFile = new cdm_FILE();
    pFile->m_fp   = fp;
    pFile->m_ncid = ncid;
    pFile->m_format = format;
    pFile->m_addMode = addMode;
    pFile->m_fname = fname;

    return pFile;
  }

  /**
   * @brief 出力ファイルをオープン(asciiモード)
   * @param [in]  fname   ファイル名
   * @param [in]  format  ファイルタイプ
   * @param [in]  addmode 追記モードフラグ(true:追記)
   * @return インスタンスされたクラスのポインタ
   */
  static cdm_FILE*
  OpenWriteAscii(const std::string fname,
                 CDM::E_CDM_FORMAT format,
                 bool addMode=false)
  {
    int ncid = -1;
    FILE *fp = NULL;

    if( format == CDM::E_CDM_FMT_NETCDF4 )
    {
      return NULL;
    }
    else
    {
      fp = fopen( fname.c_str(),"wa" );
      if( !fp )
      {
        return NULL;
      }
    }

    cdm_FILE *pFile = new cdm_FILE();
    pFile->m_fp   = fp;
    pFile->m_ncid = ncid;
    pFile->m_format = format;
    pFile->m_addMode = addMode;
    pFile->m_fname = fname;

    return pFile;
  }

  /**
   * @brief ファイルをクローズ
   * @param [in]  fp  cdm_DFIクラスのインスタンス
   */
  static void
  CloseFile( cdm_FILE *pFile )
  {
    delete pFile;
  }

};

#endif /*_CDM_FILE_H_*/
