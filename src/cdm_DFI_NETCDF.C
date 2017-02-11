/*
###################################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */
#ifdef _WITH_NETCDF4_

/**
 * @file   cdm_DFI_NETCDF.C
 * @brief  cdm_DFI_NETCDF Class
 * @author aics
 */

#include "cdm_DFI.h"
#include "cdm_DFI_NETCDF.h"

// #################################################################
// コンストラクタ
cdm_DFI_NETCDF::cdm_DFI_NETCDF()
{
  m_writeFlag = false;
}


// #################################################################
// デストラクタ
cdm_DFI_NETCDF::~cdm_DFI_NETCDF()
{

}


// #################################################################
// NetCDFの次元変数名のセット
void
cdm_DFI_NETCDF::SetDimName( string nameX, string nameY, string nameZ, string nameT )
{
  m_nameX = nameX;
  m_nameY = nameY;
  m_nameZ = nameZ;
  m_nameT = nameT;
}

// #################################################################
// NetCDFの次元変数名の取得
void
cdm_DFI_NETCDF::GetDimName( string &nameX, string &nameY, string &nameZ, string &nameT )
{
  nameX = m_nameX;
  nameY = m_nameY;
  nameZ = m_nameZ;
  nameT = m_nameT;
}

// #################################################################
// NetCDF特有のdfiパラメータ読み込み
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::ReadAdditionalTP( cdm_TextParser tpCntl )
{
  string label, tmp;

  // TP
  TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    return CDM::E_CDM_ERROR_TEXTPARSER;
  }

  // NetCDF4要素の存在チェック
  std::string label_base = "/NetCDF4";
  if( !tpCntl.chkNode(label_base) )
  {
    // 存在しない場合、デフォルト値を使用する
    // エラーにはしない
    return CDM::E_CDM_SUCCESS;
  }

  // VariableName要素の存在チェック
  // 存在しない場合、デフォルト値を使用する
  label_base += "/VariableName";
  if( tpCntl.chkNode(label_base) )
  {
    // /NetCDF4/VariableNameに移動
    tp->changeNode(label_base);

    // xを取得
    label = "x";
    if ( tpCntl.GetValue(label, &tmp, false) )
    {
      m_nameX = tmp;
    }

    // yを取得
    label = "y";
    if ( tpCntl.GetValue(label, &tmp, false) )
    {
      m_nameY = tmp;
    }

    // zを取得
    label = "z";
    if ( tpCntl.GetValue(label, &tmp, false) )
    {
      m_nameZ = tmp;
    }

    // timeを取得
    label = "time";
    if ( tpCntl.GetValue(label, &tmp, false) )
    {
      m_nameT = tmp;
    }
  }

//stmpd_printf( "vname=%s %s %s %s\n", m_nameX.c_str(), m_nameY.c_str(), m_nameZ.c_str(), m_nameT.c_str() );

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// NetCDF特有のdfiパラメータ出力
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::WriteAdditionalTP( FILE *fp, int tab )
{
  fprintf(fp, "NetCDF4 {\n");
  fprintf(fp, "\n");

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "VariableName {\n");

  _CDM_WRITE_TAB(fp, tab+2);
  fprintf(fp, "x = \"%s\"\n", m_nameX.c_str());
  _CDM_WRITE_TAB(fp, tab+2);
  fprintf(fp, "y = \"%s\"\n", m_nameY.c_str());
  _CDM_WRITE_TAB(fp, tab+2);
  fprintf(fp, "z = \"%s\"\n", m_nameZ.c_str());
  _CDM_WRITE_TAB(fp, tab+2);
  fprintf(fp, "time = \"%s\"\n", m_nameT.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "}\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// tp文字列のチェック
void CheckTpChar( char &c )
{
#if 0
  for( char i='a';i<='z';i++ )
  {
    if( i==c ) return;
  }
  for( char i='A';i<='Z';i++ )
  {
    if( i==c ) return;
  }
  for( char i='0';i<='9';i++ )
  {
    if( i==c ) return;
  }
  if( c=='_' ) return;
  if( c=='-' ) return;
  if( c=='/' ) return;
#else
  if( c==' ' )
  {
    c = '_';
  }
  else if( c==':' )
  {
    c = '-';
  }
#endif
  return;
}

// #################################################################
// 変数の情報を取得する
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::GetVarInfo( int ncid, string var_name, stVarInfo &varInfo, int nDimCheck )
{
  // clear
  varInfo.clear();

  // name
  varInfo.name = var_name;

  // var id
  if( nc_inq_varid( ncid, var_name.c_str(), &varInfo.id ) != NC_NOERR )
  {
stmpd_printf( "error nc_inq_varid\n" );
    return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
  }

  //varのdimensionの数を取得
  // ndims
  int ndims;
  if( nc_inq_varndims(ncid, varInfo.id, &ndims) != NC_NOERR )
  {
stmpd_printf( "error nc_inq_varndims\n" );
    return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
  }

  // dimensionの数をチェック
  if( nDimCheck > 0 )
  {
    if( ndims != nDimCheck )
    {
stmpd_printf( "error ndim != %d [%s]\n", nDimCheck, var_name.c_str() );
      return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
    }
  }

  //varの情報を取得
  char    name[NC_MAX_NAME];
  nc_type type;
  int*    dimids = new int[ndims+1];
  int     natts;
  if( nc_inq_var( ncid, varInfo.id, name, &varInfo.type, 0, dimids, &natts ) != NC_NOERR )
  {
stmpd_printf( "error nc_inq_var\n" );
    delete [] dimids;
    return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
  }

  //varのendian情報を取得
  if( nc_inq_var_endian( ncid, varInfo.id, &varInfo.endian) != NC_NOERR )
  {
stmpd_printf( "error nc_inq_var_endian\n" );
    delete [] dimids;
    return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
  }

  // dimidsからdimension情報を読み込み
  for( int i=0;i<ndims;i++ )
  {
    stDimInfo dimInfo;
    dimInfo.id = dimids[i];

    // name , size
    char name[NC_MAX_NAME];
    if( nc_inq_dim(ncid, dimInfo.id, name, &dimInfo.len) != NC_NOERR )
    {
stmpd_printf( "error nc_inq_dim\n" );
      delete [] dimids;
      return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
    }

    dimInfo.name = string(name);
    varInfo.dims.push_back(dimInfo);
  }
  delete [] dimids;

  //varのunitを取得
  for( int i=0;i<natts;i++ )
  {
    char attrname[NC_MAX_NAME+1] = "\0";
    if( nc_inq_attname( ncid, varInfo.id, i, attrname ) != NC_NOERR )
    {
stmpd_printf( "error nc_inq_attrname\n" );
      return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
    }

    if( !strcasecmp(attrname,"units") )
    {
      char unitname[NC_MAX_NAME+1] = "\0";
      if( nc_get_att_text( ncid, varInfo.id, attrname, unitname ) != NC_NOERR )
      {
stmpd_printf( "error nc_inq_attrname\n" );
        return CDM::E_CDM_ERROR_READ_NETCDF_FUNC;
      }
      for(int m=0;m<strlen(unitname);m++ )
      {
//        if( unitname[m] == ' ' ) unitname[m] = '_';
//        if( unitname[m] == ':' ) unitname[m] = '_';
        CheckTpChar(unitname[m]);
      }
      varInfo.unit = string(unitname);
//stmpd_printf("var[%s] unit[%s]\n", varInfo.name.c_str(), varInfo.unit.c_str());
    }
  }

  // 1次元の場合、配列の値も取得する
  if( nDimCheck == 1 )
  {
    // 領域確保
    if( varInfo.data )
    {
      delete [] varInfo.data;
    }
    varInfo.data = new double[varInfo.dims[0].len];

    // 読み込み
    if( varInfo.type == NC_FLOAT )
    {
      float *fdata = new float[varInfo.dims[0].len];
      nc_get_var_float( ncid, varInfo.id, fdata );
      for( int i=0;i<varInfo.dims[0].len;i++ )
      {
        varInfo.data[i] = fdata[i];
      }
      delete [] fdata;
    }
    else if( varInfo.type == NC_DOUBLE )
    {
      nc_get_var_double( ncid, varInfo.id, varInfo.data );
    }
    else
    {
      delete [] varInfo.data;
      varInfo.data = NULL;
    }

#if 0
if( varInfo.data )
{
  stmpd_printf("%s:\n", varInfo.name.c_str());
  for( int i=0;i<varInfo.dims[0].len;i++ )
  {
    stmpd_printf( "%4d : %e\n", i, varInfo.data[i] );
  }
}
#endif
  }

#if 0
  stmpd_printf( "var name[%s] id[%d] type[%d] ndim[%d]\n", varInfo.name.c_str(), varInfo.id, varInfo.type, (int)varInfo.dims.size() );
  for( int i=0;i<varInfo.dims.size();i++ )
  {
    stmpd_printf( "  dim name[%s] id[%d] len[%d]\n", varInfo.dims[i].name.c_str(), varInfo.dims[i].id, (int)varInfo.dims[i].len );
  }
#endif

  return CDM::E_CDM_SUCCESS;
}


// #################################################################
// ncデータ型に対応するCDMデータ型を取得
CDM::E_CDM_DTYPE
cdm_DFI_NETCDF::GetCdmType( nc_type type )
{
  CDM::E_CDM_DTYPE rtype = CDM::E_CDM_DTYPE_UNKNOWN;

  if( type == NC_BYTE )
  {
    rtype = CDM::E_CDM_INT8;
  }
  else if( type == NC_SHORT )
  {
    rtype = CDM::E_CDM_INT16;
  }
  else if( type == NC_INT )
  {
    rtype = CDM::E_CDM_INT32;
  }
//  else if( type == NC_LONG )
//  {
//    rtype = CDM::E_CDM_INT32;
//  }
  else if( type == NC_FLOAT )
  {
    rtype = CDM::E_CDM_FLOAT32;
  }
  else if( type == NC_DOUBLE )
  {
    rtype = CDM::E_CDM_FLOAT64;
  }
  else if( type == NC_UBYTE )
  {
    rtype = CDM::E_CDM_UINT8;
  }
  else if( type == NC_USHORT )
  {
    rtype = CDM::E_CDM_UINT16;
  }
  else if( type == NC_UINT )
  {
    rtype = CDM::E_CDM_UINT32;
  }
  else if( type == NC_INT64 )
  {
    rtype = CDM::E_CDM_INT64;
  }
  else if( type == NC_UINT64 )
  {
    rtype = CDM::E_CDM_UINT64;
  }

  return rtype;
}

// #################################################################
//  CDMデータ型に対応するncデータ型を取得
nc_type
cdm_DFI_NETCDF::GetNcType( CDM::E_CDM_DTYPE type )
{
  nc_type rtype = NC_NAT;

  if( type == CDM::E_CDM_INT8 )
  {
    rtype = NC_BYTE;
  }
  else if( type == CDM::E_CDM_INT16 )
  {
    rtype = NC_SHORT;
  }
  else if( type == CDM::E_CDM_INT32 )
  {
    rtype = NC_INT;
  }
  else if( type == CDM::E_CDM_FLOAT32 )
  {
    rtype = NC_FLOAT;
  }
  else if( type == CDM::E_CDM_FLOAT64 )
  {
    rtype = NC_DOUBLE;
  }
  else if( type == CDM::E_CDM_UINT8 )
  {
    rtype = NC_UBYTE;
  }
  else if( type == CDM::E_CDM_UINT16 )
  {
    rtype = NC_USHORT;
  }
  else if( type == CDM::E_CDM_UINT32 )
  {
    rtype = NC_UINT;
  }
  else if( type == CDM::E_CDM_INT64 )
  {
    rtype = NC_INT64;
  }
  else if( type == CDM::E_CDM_UINT64 )
  {
    rtype = NC_UINT64;
  }

  return rtype;
}


// #################################################################
// nc変数配列のdimensionをチェック
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::CheckArrayVarInfo( stVarInfo &varInfo )
{
  // dimensionのチェック
  if( varInfo.dims.size() == 3 )
  {
    // timeが無いとき(z,y,x)
    if( varInfo.dims[0].id != m_varInfoZ.dims[0].id ||
        varInfo.dims[1].id != m_varInfoY.dims[0].id ||
        varInfo.dims[2].id != m_varInfoX.dims[0].id )
    {
stmpd_printf( "error dimension is invalid %s(%s,%s,%s), necessary(%s,%s,%s)\n"
  , varInfo.name.c_str()
  , varInfo.dims[0].name.c_str(), varInfo.dims[1].name.c_str(), varInfo.dims[2].name.c_str()
  , m_varInfoZ.dims[0].name.c_str(), m_varInfoY.dims[0].name.c_str(), m_varInfoX.dims[0].name.c_str()
);
      return CDM::E_CDM_ERROR;
    }
  }
  else if( varInfo.dims.size() == 4 )
  {
    // timeがあるとき(time,z,y,x)
    if( varInfo.dims[0].id != m_varInfoT.dims[0].id ||
        varInfo.dims[1].id != m_varInfoZ.dims[0].id ||
        varInfo.dims[2].id != m_varInfoY.dims[0].id ||
        varInfo.dims[3].id != m_varInfoX.dims[0].id )
    {
stmpd_printf( "error dimension is invalid %s(%s,%s,%s,%s), necessary(%s,%s,%s,%s)\n"
  , varInfo.name.c_str()
  , varInfo.dims[0].name.c_str(), varInfo.dims[1].name.c_str(), varInfo.dims[2].name.c_str(), varInfo.dims[3].name.c_str()
  , m_varInfoT.dims[0].name.c_str(), m_varInfoZ.dims[0].name.c_str(), m_varInfoY.dims[0].name.c_str(), m_varInfoX.dims[0].name.c_str()
);
      return CDM::E_CDM_ERROR;
    }
  }
  else
  {
    // 3次元でも4次元でもない場合エラー
stmpd_printf( "error dimension is invalid, necessary (%s,%s,%s) or (%s,%s,%s,%s)\n"
  , m_varInfoZ.dims[0].name.c_str(), m_varInfoY.dims[0].name.c_str(), m_varInfoX.dims[0].name.c_str()
  , m_varInfoT.dims[0].name.c_str(), m_varInfoZ.dims[0].name.c_str(), m_varInfoY.dims[0].name.c_str(), m_varInfoX.dims[0].name.c_str()
);
    return CDM::E_CDM_ERROR;
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// UnitListに変数名のunitが存在する場合、var属性として出力する
void
cdm_DFI_NETCDF::write_AttUnits( int ncid, int varid, string var_name, cdm_Unit &unit )
{
  cdm_UnitElem uElem;
  if( unit.GetUnitElem(var_name, uElem) == CDM::E_CDM_SUCCESS )
  {
    const char *ustr = uElem.Unit.c_str();
    size_t len = strlen(ustr);
    nc_put_att_text(ncid, varid, "units", len, ustr);
  }
}

// #################################################################
// ファイルのヘッダーレコード読込み
CDM::E_CDM_ERRORCODE
//cdm_DFI_NETCDF::read_HeaderRecord(int ncid,
cdm_DFI_NETCDF::read_HeaderRecord(cdm_FILE *pFile,
                                  bool matchEndian,
                                  unsigned step,
                                  const int head[3],
                                  const int tail[3],
                                  int gc,
                                  int voxsize[3],
                                  double &time)
{
//stmpd_printf( "***** function start *****\n" );
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;

  int ncid = pFile->m_ncid;

  // read header
  if( (ret = read_HeaderRecord(ncid, m_nameX, m_nameY, m_nameZ, m_nameT, m_varInfoX, m_varInfoY, m_varInfoZ, m_varInfoT)) != CDM::E_CDM_SUCCESS )
  {
    return ret;
  }

  // check voxel size
  voxsize[0] = m_varInfoX.dims[0].len;
  voxsize[1] = m_varInfoY.dims[0].len;
  voxsize[2] = m_varInfoZ.dims[0].len;
  for(int i=0; i<3; i++)
  {
    if( voxsize[i] != (tail[i]-head[i]+1+2*gc) )
    {
stmpd_printf("ERROR!!!!!!!!!!!! g(%d) h(%d,%d,%d) t(%d,%d,%d) v(%d,%d,%d)\n", gc, head[0], head[1], head[2], tail[0], tail[1], tail[2], voxsize[0], voxsize[1], voxsize[2]);
      return CDM::E_CDM_ERROR_UNMATCH_VOXELSIZE;
    }
  }

  // time
  time=0.0;
  for( int i=0; i<DFI_TimeSlice.SliceList.size(); i++ )
  {
     if( DFI_TimeSlice.SliceList[i].step == step )
     {
       time=(double)DFI_TimeSlice.SliceList[i].time;
     }
  }

  // UnitListに登録
  this->AddNcUnit( m_varInfoX.name, m_varInfoX.unit, 0.0);
  this->AddNcUnit( m_varInfoY.name, m_varInfoY.unit, 0.0);
  this->AddNcUnit( m_varInfoZ.name, m_varInfoZ.unit, 0.0);
  this->AddNcUnit( m_varInfoT.name, m_varInfoT.unit, 0.0);

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ファイルのヘッダーレコード読込み(static関数)
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::read_HeaderRecord(int ncid,
                                  std::string nameX,
                                  std::string nameY,
                                  std::string nameZ,
                                  std::string nameT,
                                  stVarInfo &varInfoX,
                                  stVarInfo &varInfoY,
                                  stVarInfo &varInfoZ,
                                  stVarInfo &varInfoT)
{
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;

  // x,y,z,timeのdimensionを取得
  if( (ret = GetVarInfo( ncid, nameX, varInfoX, 1 )) != CDM::E_CDM_SUCCESS )
  {
    return ret;
  }
  if( (ret = GetVarInfo( ncid, nameY, varInfoY, 1 )) != CDM::E_CDM_SUCCESS )
  {
    return ret;
  }
  if( (ret = GetVarInfo( ncid, nameZ, varInfoZ, 1 )) != CDM::E_CDM_SUCCESS )
  {
    return ret;
  }
  if( (ret = GetVarInfo( ncid, nameT, varInfoT, 1 )) != CDM::E_CDM_SUCCESS )
  {
    return ret;
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ファイルのデーターレコード読込み
CDM::E_CDM_ERRORCODE
//cdm_DFI_NETCDF::read_Datarecord(int ncid,
cdm_DFI_NETCDF::read_Datarecord(cdm_FILE *pFile,
                                bool matchEndian,
                                unsigned step,
                                cdm_Array* buf,
                                int head[3],
                                int nz,
                                cdm_Array* &src)
{
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;

  int ncid = pFile->m_ncid;

  // 読み込み時刻の取得
  // (FieldFilenameFormatがRANKの場合のみ、何番目に格納されている時刻かを決定)
  int step_index = 0;
  if( DFI_Finfo.FieldFilenameFormat == CDM::E_CDM_FNAME_RANK )
  {
    for( int i=0; i<DFI_TimeSlice.SliceList.size(); i++ )
    {
       if( DFI_TimeSlice.SliceList[i].step == step )
       {
         step_index = i;
       }
    }
  }

  //１層ずつ読み込み
  int hzB = head[2];

  //成分数を取得
  int nVar = GetNumVariables();
//stmpd_printf("dfi nver= %d\n", nVar );

  // 成分数回ループ
  for( int n=0;n<nVar;n++ )
  {
    // 成分名
    string var_name = getVariableName(n);

    // var情報を取得
    stVarInfo varInfo;
    if( (ret = GetVarInfo( ncid, var_name, varInfo )) != CDM::E_CDM_SUCCESS )
    {
      return ret;
    }

    // varのデータ型をチェック
    CDM::E_CDM_DTYPE cdm_type = GetCdmType( varInfo.type );
    if( cdm_type != DFI_Finfo.DataType )
    {
      return CDM::E_CDM_ERROR_READ_NETCDF_MISMATCH_TYPE;
    }

    // varのdimensionをチェック
    if( (ret = CheckArrayVarInfo( varInfo )) != CDM::E_CDM_SUCCESS )
    {
      return ret;
    }

    // start & count for read nc
    size_t start[4];
    size_t count[4];

    // data array
    if( varInfo.dims.size() == 3 )
    {
      start[0] = 0; //z
      start[1] = 0; //y
      start[2] = 0; //x
      count[0] = 1;                   //z
      count[1] = varInfo.dims[1].len; //y
      count[2] = varInfo.dims[2].len; //x
    }
    else if( varInfo.dims.size() == 4 )
    {
      start[0] = step_index;          //step
      start[1] = 0;                   //z
      start[2] = 0;                   //y
      start[3] = 0;                   //x
      count[0] = 1;                   //time
      count[1] = 1;                   //z
      count[2] = varInfo.dims[2].len; //y
      count[3] = varInfo.dims[3].len; //x
    }

    // １層ずつ読み込み
    for( int k=0;k<nz;k++ )
    {
      // headインデクスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      // Z方向の開始位置を変更
      start[1] = k;

      // １層読み込み
      if( varInfo.type == NC_BYTE )
      {
        signed char *ptr = (signed char*)buf->getData();
        nc_get_vara_schar( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_SHORT )
      {
        short *ptr = (short*)buf->getData();
        nc_get_vara_short( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_INT )
      {
        int *ptr = (int*)buf->getData();
        nc_get_vara_int( ncid, varInfo.id, start, count, ptr );
      }
//      else if( varInfo.type == NC_LONG )
//      { // same as NC_INT?
//        int *ptr = (int*)buf->getData();
//        nc_get_vara_long( ncid, varInfo.id, start, count, ptr );
//      }
      else if( varInfo.type == NC_FLOAT )
      {
        float *ptr = (float*)buf->getData();
        nc_get_vara_float( ncid, varInfo.id, start, count, ptr );
#if 0
if( k==0 && n==0 )
{
  stmpd_printf( "step=%d fname=%s\n", step, pFile->m_fname.c_str() );
  stmpd_printf( "  start = %d %d %d %d\n", start[0], start[1], start[2], start[3] );
  stmpd_printf( "  count = %d %d %d %d\n", count[0], count[1], count[2], count[3]  );
  stmpd_printf( "  val0 = %g\n", ptr[0] );
  stmpd_printf( "  val1 = %g\n", ptr[1] );
  stmpd_printf( "  val2 = %g\n", ptr[2] );
  stmpd_printf( "  val3 = %g\n", ptr[3] );
  stmpd_printf( "  val4 = %g\n", ptr[4] );
}
//MPI_Abort(MPI_COMM_WORLD,9999);
#endif
      }
      else if( varInfo.type == NC_DOUBLE )
      {
        double *ptr = (double*)buf->getData();
        nc_get_vara_double( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_UBYTE )
      {
        unsigned char *ptr = (unsigned char*)buf->getData();
        nc_get_vara_uchar( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_USHORT )
      {
        unsigned short *ptr = (unsigned short*)buf->getData();
        nc_get_vara_ushort( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_UINT )
      {
        unsigned int *ptr = (unsigned int*)buf->getData();
        nc_get_vara_uint( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_INT64 )
      {
        long long *ptr = (long long*)buf->getData();
        nc_get_vara_longlong( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_UINT64 )
      {
        unsigned long long *ptr = (unsigned long long*)buf->getData();
        nc_get_vara_ulonglong( ncid, varInfo.id, start, count, ptr );
      }
      else
      {
stmpd_printf( "ERROR!!!! not implemented data type [%d]\n", varInfo.type );
      }

      // copy
      buf->copyArrayNvari(src,n);
    }

    // 読み込みが正常終了したので、UnitListをセットする
    if( varInfo.unit.length() > 0 )
    {
      this->AddNcUnit( varInfo.name, varInfo.unit, 0.0);
//stmpd_printf("  add UnitList %s %s\n", varInfo.name.c_str(), varInfo.unit.c_str());
    }
  } // for(nVar)

  return CDM::E_CDM_SUCCESS;
}


// #################################################################
// ファイルのデーターレコード読込み(static関数)
cdm_Array*
cdm_DFI_NETCDF::read_Datarecord(int ncid,
                                unsigned step_index,
                                std::vector<std::string> &VariableName,
                                CDM::E_CDM_DTYPE DataType,
                                stVarInfo &varInfoX,
                                stVarInfo &varInfoY,
                                stVarInfo &varInfoZ,
                                stVarInfo &varInfoT,
                                std::vector<stVarInfo> &vecVarInfo,
                                CDM::E_CDM_ERRORCODE &ret)
{
  ret = CDM::E_CDM_SUCCESS;

  cdm_Array *pArray = NULL;
  cdm_Array *buf = NULL;
  int head[3] = {1,1,1};

  // 成分数回ループ
  for( int n=0;n<VariableName.size();n++ )
  {
    // 成分名
    string var_name = VariableName[n];

    // var情報を取得
    stVarInfo varInfo;
    if( (ret = GetVarInfo( ncid, var_name, varInfo )) != CDM::E_CDM_SUCCESS )
    {
      delete buf;
      delete pArray;
      return NULL;
    }
    vecVarInfo.push_back(varInfo);

    // varのデータ型をチェック
    CDM::E_CDM_DTYPE cdm_type = GetCdmType( varInfo.type );
    if( DataType != CDM::E_CDM_DTYPE_UNKNOWN )
    {
      if( cdm_type != DataType )
      {
        ret = CDM::E_CDM_ERROR_READ_NETCDF_MISMATCH_TYPE;
        delete buf;
        delete pArray;
        return NULL;
      }
    }
    else
    {
      DataType = cdm_type;
    }
    if( pArray )
    {
      if( DataType != pArray->getDataType() )
      {
        stmpd_printf("data type is not match : %s\n", varInfo.name.c_str());
        delete buf;
        delete pArray;
        return NULL;
      }
    }

    // start & count for read nc
    size_t start[4];
    size_t count[4];

    // data array size
    int nx, ny, nz;
    if( varInfo.dims.size() == 3 )
    {
      start[0] = 0; //z
      start[1] = 0; //y
      start[2] = 0; //x
      count[0] = 1;                   //z
      count[1] = varInfo.dims[1].len; //y
      count[2] = varInfo.dims[2].len; //x
      nx = varInfo.dims[2].len;
      ny = varInfo.dims[1].len;
      nz = varInfo.dims[0].len;
    }
    else if( varInfo.dims.size() == 4 )
    {
      start[0] = step_index;          //step
      start[1] = 0;                   //z
      start[2] = 0;                   //y
      start[3] = 0;                   //x
      count[0] = 1;                   //time
      count[1] = 1;                   //z
      count[2] = varInfo.dims[2].len; //y
      count[3] = varInfo.dims[3].len; //x
      nx = varInfo.dims[3].len;
      ny = varInfo.dims[2].len;
      nz = varInfo.dims[1].len;
    }
    else
    {
      stmpd_printf("not supported dimension : %d\n", (int)varInfo.dims.size());
      delete buf;
      delete pArray;
      ret = CDM::E_CDM_ERROR_READ_DFI_NETCDF;
      return NULL;
    }

    // data array
    if( !pArray )
    {
      pArray = cdm_Array::instanceArray
                   ( DataType
                   , CDM::E_CDM_IJKN
                   , nx
                   , ny
                   , nz
                   , 0
                   , (int)VariableName.size() );
      pArray->setHeadIndex(head);

    }
    else
    {
      const int *sz = pArray->getArraySizeInt();
      if( sz[0] != nx || sz[1] != ny || sz[2] != nz )
      {
        stmpd_printf("buffer array dimension size is not match : %s\n", varInfo.name.c_str());
        delete buf;
        delete pArray;
        return NULL;
      }
    }
    if( !buf )
    {
      buf = cdm_Array::instanceArray
                   ( DataType
                   , CDM::E_CDM_IJKN
                   , nx
                   , ny
                   , 1
                   , 0
                   , 1 );
      buf->setHeadIndex(head);
    }

    // １層ずつ読み込み
    for( int k=0;k<nz;k++ )
    {
      // headインデクスをずらす
      head[2] = k + 1;
      buf->setHeadIndex(head);

      // Z方向の開始位置を変更
      start[1] = k;

      // １層読み込み
      if( varInfo.type == NC_BYTE )
      {
        signed char *ptr = (signed char*)buf->getData();
        nc_get_vara_schar( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_SHORT )
      {
        short *ptr = (short*)buf->getData();
        nc_get_vara_short( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_INT )
      {
        int *ptr = (int*)buf->getData();
        nc_get_vara_int( ncid, varInfo.id, start, count, ptr );
      }
//      else if( varInfo.type == NC_LONG )
//      { // same as NC_INT?
//        int *ptr = (int*)buf->getData();
//        nc_get_vara_long( ncid, varInfo.id, start, count, ptr );
//      }
      else if( varInfo.type == NC_FLOAT )
      {
        float *ptr = (float*)buf->getData();
        nc_get_vara_float( ncid, varInfo.id, start, count, ptr );
//stmpd_printf( "start = %d %d %d %d\n", start[0], start[1], start[2], start[3] );
//stmpd_printf( "count = %d %d %d %d\n", count[0], count[1], count[2], count[3]  );
//stmpd_printf( "val0 = %g\n", ptr[0] );
//stmpd_printf( "val1 = %g\n", ptr[1] );
//stmpd_printf( "val2 = %g\n", ptr[2] );
//stmpd_printf( "val3 = %g\n", ptr[3] );
//stmpd_printf( "val4 = %g\n", ptr[4] );
//MPI_Abort(MPI_COMM_WORLD,9999);
      }
      else if( varInfo.type == NC_DOUBLE )
      {
        double *ptr = (double*)buf->getData();
        nc_get_vara_double( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_UBYTE )
      {
        unsigned char *ptr = (unsigned char*)buf->getData();
        nc_get_vara_uchar( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_USHORT )
      {
        unsigned short *ptr = (unsigned short*)buf->getData();
        nc_get_vara_ushort( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_UINT )
      {
        unsigned int *ptr = (unsigned int*)buf->getData();
        nc_get_vara_uint( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_INT64 )
      {
        long long *ptr = (long long*)buf->getData();
        nc_get_vara_longlong( ncid, varInfo.id, start, count, ptr );
      }
      else if( varInfo.type == NC_UINT64 )
      {
        unsigned long long *ptr = (unsigned long long*)buf->getData();
        nc_get_vara_ulonglong( ncid, varInfo.id, start, count, ptr );
      }
      else
      {
stmpd_printf( "ERROR!!!! not implemented data type [%d]\n", varInfo.type );
        ret = CDM::E_CDM_ERROR_READ_DFI_NETCDF;
        delete buf;
        delete pArray;
        return NULL;
      }

      // copy
      buf->copyArrayNvari(pArray,n);
    }
  } // for(nVar)

  delete buf;
//FILE *fp = fopen("test.data", "wt");
//pArray->writeAscii(fp);
//fclose(fp);
//stmpd_printf("write ascii\n");
//MPI_Abort(MPI_COMM_WORLD,9999);
  return pArray;
}


// #################################################################
// Averaged レコードの読込み
CDM::E_CDM_ERRORCODE
//cdm_DFI_NETCDF::read_averaged(int ncid,
cdm_DFI_NETCDF::read_averaged(cdm_FILE *pFile,
                              bool matchEndian,
                              unsigned dummy,
                              unsigned &step_avr,
                              double &time_avr)
{
  int ncid = pFile->m_ncid;

  return CDM::E_CDM_SUCCESS;
}


// #################################################################
// NETCDFヘッダーレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::write_HeaderRecord(cdm_FILE *pFile,
                                   const unsigned step,
                                   const double time,
                                   const int n)
{
  // 変数名
  m_varInfoX.name = m_nameX;
  m_varInfoY.name = m_nameY;
  m_varInfoZ.name = m_nameZ;
  m_varInfoT.name = m_nameT;

  // セル座標値
  int gc = DFI_Finfo.GuideCell;
  int head[3];
  for(int i=0; i<3; i++) {
    head[i] = DFI_Process.RankList[n].HeadIndex[i];
  }
  if( !m_varInfoX.data )
  {
    int nn = DFI_Process.RankList[n].VoxelSize[0];
    double *cod = new double[nn+2*gc];
    for( int i=0;i<nn+2*gc;i++ )
    {
      cod[i] = DFI_Domain->CellX(i+head[0]-1-gc);
    }
    m_varInfoX.data = cod;
  }
  if( !m_varInfoY.data )
  {
    int nn = DFI_Process.RankList[n].VoxelSize[1];
    double *cod = new double[nn+2*gc];
    for( int i=0;i<nn+2*gc;i++ )
    {
      cod[i] = DFI_Domain->CellY(i+head[1]-1-gc);
    }
    m_varInfoY.data = cod;
  }
  if( !m_varInfoZ.data )
  {
    int nn = DFI_Process.RankList[n].VoxelSize[2];
    double *cod = new double[nn+2*gc];
    for( int i=0;i<nn+2*gc;i++ )
    {
      cod[i] = DFI_Domain->CellZ(i+head[2]-1-gc);
    }
    m_varInfoZ.data = cod;
  }

  // 出力
  CDM::E_CDM_ERRORCODE ret = write_HeaderRecord(pFile, step, time, DFI_Process.RankList[n].VoxelSize, DFI_Finfo.GuideCell,
                                                DFI_Finfo.DataType, DFI_Finfo.VariableName, DFI_Unit,
                                                m_varInfoX, m_varInfoY, m_varInfoZ, m_varInfoT, m_vecVarInfo);

  return ret;
}


// #################################################################
// NETCDFヘッダーレコードの出力(static関数)
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::write_HeaderRecord(cdm_FILE *pFile,
                                   const unsigned step,
                                   const double time,
                                   int VoxelSize[3],
                                   int GuideCell,
                                   CDM::E_CDM_DTYPE DataType,
                                   vector<string> &vecVariable,
                                   cdm_Unit &Unit,
                                   stVarInfo &varInfoX,
                                   stVarInfo &varInfoY,
                                   stVarInfo &varInfoZ,
                                   stVarInfo &varInfoT,
                                   vector<stVarInfo> &vecVarInfo)
{
  // ncid
  int ncid = pFile->m_ncid;

  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;

  // name
  string nameX = varInfoX.name;
  string nameY = varInfoY.name;
  string nameZ = varInfoZ.name;
  string nameT = varInfoT.name;

  // cell coord
  double *cellX = varInfoX.data;
  double *cellY = varInfoY.data;
  double *cellZ = varInfoZ.data;

  // 変数名
  int nVar = vecVariable.size();

  // 変数情報を削除
  vecVarInfo.clear();

  //*** 新規作成モード ***//
  if( !pFile->m_addMode )
  {
    // 定義開始
    nc_redef(ncid);

    // dimensionの定義
    int dimids[4];
    int size[3];
    for( int i=0;i<3;i++ )
    {
      size[i] = VoxelSize[i]+(2*GuideCell);
    }
    nc_def_dim(ncid, nameX.c_str(), size[0], &dimids[3]);
    nc_def_dim(ncid, nameY.c_str(), size[1], &dimids[2]);
    nc_def_dim(ncid, nameZ.c_str(), size[2], &dimids[1]);
    nc_def_dim(ncid, nameT.c_str(), NC_UNLIMITED, &dimids[0]);

    // x,y,z,time配列の定義
    int varids[4];
    nc_def_var(ncid, nameX.c_str(), NC_DOUBLE, 1, &dimids[3], &varids[3]);
    nc_def_var(ncid, nameY.c_str(), NC_DOUBLE, 1, &dimids[2], &varids[2]);
    nc_def_var(ncid, nameZ.c_str(), NC_DOUBLE, 1, &dimids[1], &varids[1]);
    nc_def_var(ncid, nameT.c_str(), NC_DOUBLE, 1, &dimids[0], &varids[0]);

    // x,y,z,time配列の属性(units)
    write_AttUnits( ncid, varids[3], nameX, Unit );
    write_AttUnits( ncid, varids[2], nameY, Unit );
    write_AttUnits( ncid, varids[1], nameZ, Unit );
    write_AttUnits( ncid, varids[0], nameT, Unit );

    // 物理量配列の定義
    for( int i=0;i<nVar;i++ )
    {
      // 配列
      int varid;
      string var_name = vecVariable[i];
      nc_type ntype = GetNcType( DataType );
      nc_def_var(ncid, var_name.c_str(), ntype, 4, dimids, &varid);

      // 属性(units)
      write_AttUnits( ncid, varid, var_name, Unit );
    }

    // 定義終了
    nc_enddef(ncid);

    // x,y,z,time配列の出力
    nc_put_var_double(ncid, varids[3], cellX);
    nc_put_var_double(ncid, varids[2], cellY);
    nc_put_var_double(ncid, varids[1], cellZ);

    //// timeの出力
    {
      size_t start = 0;
      size_t count = 1;
      int rr = nc_put_vara_double( ncid, varids[0], &start, &count, &time );
    }
  }

  // varinfoの再取得
  GetVarInfo( ncid, nameX, varInfoX, 1 );
  GetVarInfo( ncid, nameY, varInfoY, 1 );
  GetVarInfo( ncid, nameZ, varInfoZ, 1 );
  GetVarInfo( ncid, nameT, varInfoT, 1 );
  for( int i=0;i<nVar;i++ )
  {
    stVarInfo varInfo;
    string var_name = vecVariable[i];
    GetVarInfo( ncid, var_name, varInfo );
    vecVarInfo.push_back(varInfo);
  }

  //*** 追記モード ***//
  if( pFile->m_addMode )
  {
    // time配列の追記
    size_t start = varInfoT.dims[0].len;
    size_t count = 1;
    nc_put_vara_double( ncid, varInfoT.id, &start, &count, &time );

    // timeのvarInfoを取得しなおし
    GetVarInfo( ncid, nameT, varInfoT, 1 );
  }

  return CDM::E_CDM_SUCCESS;
}


// #################################################################
// NETCDFデータレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::write_DataRecord(cdm_FILE *pFile,
                                 cdm_Array* val,
                                 const int gc,
                                 const int n)
{
  // 出力
  CDM::E_CDM_ERRORCODE ret = write_DataRecord(pFile, val, gc, DFI_Process.RankList[n].VoxelSize, m_varInfoT, m_vecVarInfo);

  // 出力フラグ
  m_writeFlag = true;

  return ret;
}

// #################################################################
// NETCDFデータレコードの出力(static関数)
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::write_DataRecord(cdm_FILE *pFile,
                                 cdm_Array* val,
                                 const int gc,
                                 int VoxelSize[3],
                                 stVarInfo &varInfoT,
                                 vector<stVarInfo> &vecVarInfo)
{
  // ncid
  int ncid = pFile->m_ncid;

  // 成分スカラーデータ格納領域
  CDM::E_CDM_DTYPE      dtype = val->getDataType();
  CDM::E_CDM_ARRAYSHAPE shape = CDM::E_CDM_IJKN;
  size_t sz[3];
  for( int i=0;i<3;i++ )
  {
    sz[i] = val->getArraySize()[i];
  }
  cdm_Array *wrk = cdm_Array::instanceArray( dtype, shape, sz, gc, 1 );
  int head[3];
  head[0] = val->getHeadIndex()[0];
  head[1] = val->getHeadIndex()[1];
  head[2] = val->getHeadIndex()[2];
  wrk->setHeadIndex( head );

  // ncデータ型
  nc_type type = GetNcType( dtype );

  // 仮想セルを含めた格子数
  int size[3];
  for( int i=0;i<3;i++ )
  {
    size[i] = VoxelSize[i]+(2*gc);
  }

  // timeのシフト量を取得
  int step = varInfoT.dims[0].len - 1;

  // start, count
  size_t start[4] = {step,0,0,0};
  size_t count[4] = {1, size[2], size[1], size[0]};

  // 成分ごとにスカラーデータを抽出して出力
  int nVar = vecVarInfo.size();
  for( int i=0;i<nVar;i++ )
  {
    // 抽出
    int iret = val->copyArrayNvari_to_ijk( wrk, i );

    // 配列の出力
    int varid = vecVarInfo[i].id;

    if( type == NC_BYTE )
    {
      signed char *ptr = (signed char*)wrk->getData();
      nc_put_vara_schar( ncid, varid, start, count, ptr );
    }
    else if( type == NC_SHORT )
    {
      short *ptr = (short*)wrk->getData();
      nc_put_vara_short( ncid, varid, start, count, ptr );
    }
    else if( type == NC_INT )
    {
      int *ptr = (int*)wrk->getData();
       nc_put_vara_int( ncid, varid, start, count, ptr );
    }
    else if( type == NC_FLOAT )
    {
      float *ptr = (float*)wrk->getData();
      nc_put_vara_float( ncid, varid, start, count, ptr );
    }
    else if( type == NC_DOUBLE )
    {
      double *ptr = (double*)wrk->getData();
      nc_put_vara_double( ncid, varid, start, count, ptr );
    }
    else if( type == NC_UBYTE )
    {
      unsigned char *ptr = (unsigned char*)wrk->getData();
      nc_put_vara_uchar( ncid, varid, start, count, ptr );
    }
    else if( type == NC_USHORT )
    {
      unsigned short *ptr = (unsigned short*)wrk->getData();
      nc_put_vara_ushort( ncid, varid, start, count, ptr );
    }
    else if( type == NC_UINT )
    {
      unsigned int *ptr = (unsigned int*)wrk->getData();
      nc_put_vara_uint( ncid, varid, start, count, ptr );
    }
    else if( type == NC_INT64 )
    {
      long long *ptr = (long long*)wrk->getData();
      nc_put_vara_longlong( ncid, varid, start, count, ptr );
    }
    else if( type == NC_UINT64 )
    {
      unsigned long long *ptr = (unsigned long long*)wrk->getData();
      nc_put_vara_ulonglong( ncid, varid, start, count, ptr );
    }
  }

  // ワーク配列を破棄
  delete wrk;

  // 出力フラグ
//  m_writeFlag = true;

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// 平均の出力
CDM::E_CDM_ERRORCODE
cdm_DFI_NETCDF::write_averaged(cdm_FILE *pFile,
                               const unsigned step_avr,
                               const double time_avr)
{
  int ncid = pFile->m_ncid;
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// Uuitをセットする
void
cdm_DFI_NETCDF::AddNcUnit(const std::string Name,
                          const std::string Unit,
                          const double reference,
                          const double difference,
                          const bool BsetDiff)
{
  if( Unit.size() == 0 ) return;

  /** UnitElemの生成 */
  cdm_UnitElem unit =  cdm_UnitElem(Name,Unit,reference,difference,BsetDiff);

  /** UnilListへのセット */
  m_NC_Unit.UnitList.insert(map<std::string,cdm_UnitElem>::value_type(Name,unit));
}

// #################################################################
// cdm_Unitクラスのポインタを取得
const cdm_Unit*
cdm_DFI_NETCDF::GetNcUnit()
{
  return &m_NC_Unit;
}

// #################################################################
// 出力処理を追記モードにするかどうかをチェック(NetCDF用)
bool
cdm_DFI_NETCDF::CheckAddWriteMode()
{
  if( m_writeFlag && m_output_fname == CDM::E_CDM_FNAME_RANK )
  {
    return true;
  }
  return false;
}

#endif /* _WITH_NETCDF4_ */
