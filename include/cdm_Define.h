#ifndef _CDM_DEFINE_H_
#define _CDM_DEFINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   cdm_Define.h
 * @brief  CDMの定義マクロ記述ヘッダーファイル
 * @author aics
 */

#ifdef _CDM_WITHOUT_MPI_
 #include "mpi_stubs.h"
#else
 #include "mpi.h"
#endif

#define D_CDM_DFITYPE_CARTESIAN "Cartesian"
#define D_CDM_DFITYPE_NON_UNIFORM_CARTESIAN "Non_Uniform_Cartesian"

#define D_CDM_EXT_SPH "sph"
#define D_CDM_EXT_BOV "bov"
#define D_CDM_EXT_FUNC "fun"
#define D_CDM_EXT_VTK  "vtk"
//20150918.NetCDF.s
#define D_CDM_EXT_NC  "nc"
//20150918.NetCDF.e
//20160329.fub.s
#define D_CDM_EXT_FUB "fub"
#define D_CDM_EXT_XYZ "xyz"
//20160329.fub.e


#define D_CDM_EXT_BOV_DATAFILE "dat"

#define D_CDM_ON  "on"
#define D_CDM_OFF "off"

#define D_CDM_INT8    "Int8"
#define D_CDM_INT16   "Int16"
#define D_CDM_INT32   "Int32"
#define D_CDM_INT64   "Int64"
#define D_CDM_UINT8   "UInt8"
#define D_CDM_UINT16  "UInt16"
#define D_CDM_UINT32  "UInt32"
#define D_CDM_UINT64  "UInt64"
#define D_CDM_FLOAT32 "Float32"
#define D_CDM_FLOAT64 "Float64"

#define D_CDM_BYTE   "BYTE"
#define D_CDM_INT    "INT"
#define D_CDM_FLOAT  "FLOAT"
#define D_CDM_DOUBLE "DOUBLE"

#define D_CDM_IJNK "ijkn"
#define D_CDM_NIJK "nijk"

#define D_CDM_LITTLE "little"
#define D_CDM_BIG    "big"

#define _CDM_TAB_STR "  "


/** namespace の設定 */
namespace CDM
{
  /// ファイル名のランク番号前の文字列
  static const char *C_CDM_RANKNOPREFIX = "_id";

  enum E_CDM_DFITYPE
  {
    E_CDM_DFITYPE_UNKNOWN = -1,          ///< 未定
    E_CDM_DFITYPE_CARTESIAN,             ///< Cartesian
    E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN, ///< Non Uniform Cartesian
  };

/** File 形式 */
  enum E_CDM_FORMAT
  {
    E_CDM_FMT_UNKNOWN = -1,  ///< 未定
    E_CDM_FMT_SPH,           ///< sph format
    //E_CDM_FMT_BOV            ///< bov format
    E_CDM_FMT_BOV,           ///< bov format
    E_CDM_FMT_AVS,           ///< avs format
    E_CDM_FMT_PLOT3D,        ///< plot3d format
    E_CDM_FMT_VTK,           ///< vtk format
//20150918.NetCDF.s
    E_CDM_FMT_NETCDF4,       ///< NetCDF4(/w HDF5) format
//20150918.NetCDF.e
//20160328.fub.s
    E_CDM_FMT_FUB,           ///< upacs fub field data format
    E_CDM_FMT_FUB_COD,       ///< upacs fub coordinate data format
//20160328.fub.e
  };

/** スイッチ on or off */
  enum E_CDM_ONOFF
  {
    E_CDM_OFF = 0,           ///< off
    E_CDM_ON                 ///< on
  };

/** データ形式 */
  enum E_CDM_DTYPE
  {
    E_CDM_DTYPE_UNKNOWN = 0, ///< 未定
    E_CDM_INT8,              ///< char
    E_CDM_INT16,             ///< short
    E_CDM_INT32,             ///< int
    E_CDM_INT64,             ///< long long
    E_CDM_UINT8,             ///< unsigned char
    E_CDM_UINT16,            ///< unsigned short
    E_CDM_UINT32,            ///< unsigned int
    E_CDM_UINT64,            ///< unsigned long long
    E_CDM_FLOAT32,           ///< float
    E_CDM_FLOAT64            ///< double
  };

/** 配列形式 */
  enum E_CDM_ARRAYSHAPE
  {
    E_CDM_ARRAYSHAPE_UNKNOWN=-1, ///< 未定
    E_CDM_IJKN=0,                ///< ijkn
    E_CDM_NIJK                   ///< nijk
  };

/** Endian形式 */
  enum E_CDM_ENDIANTYPE
  {
    E_CDM_ENDIANTYPE_UNKNOWN=-1,
    E_CDM_LITTLE=0,
    E_CDM_BIG
  };

/** 読込みタイプコード */
  enum E_CDM_READTYPE
  {
    E_CDM_SAMEDIV_SAMERES=1,     ///<同一分割＆同一密度
    E_CDM_SAMEDIV_REFINEMENT,    ///<同一分割＆粗密
    E_CDM_DIFFDIV_SAMERES,       ///<MxN＆同一密度
    E_CDM_DIFFDIV_REFINEMENT,    ///<MxN＆粗密
    E_CDM_READTYPE_UNKNOWN,      ///<error
  };

/** ファイルタイプ */
  enum E_CDM_FILE_TYPE
  {
    E_CDM_FILE_TYPE_DEFAULT=-1, ///<デフォルト（binary)
    E_CDM_FILE_TYPE_ASCII=0,    ///<ascii
    E_CDM_FILE_TYPE_BINARY,     ///<binary
    E_CDM_FILE_TYPE_FBINARY     ///<Fortran Binary
  };

/** フィールドファイルの命名規約タイプ */
  enum E_CDM_OUTPUT_FNAME
  {
    E_CDM_FNAME_DEFAULT=-1,  ///<出力ファイル命名規約デフォルト(step_rank)
    E_CDM_FNAME_STEP_RANK=0, ///<step_rank
    E_CDM_FNAME_RANK_STEP,   ///<rank_step
//20160418.fub.s
  //E_CDM_FNAME_RANK         ///<rank(NetCDFのみ)
    E_CDM_FNAME_RANK,        ///<rank(NetCDF,fubのみ)
    E_CDM_FNAME_CUSTOM       ///<FieldFileNameFormatに指示あり
//20160418.fub.e
  };

/** CDMのエラーコード */
  enum E_CDM_ERRORCODE
  {
    E_CDM_SUCCESS                                = 1    ///< 正常終了
,   E_CDM_ERROR                                  = -1   ///< エラー終了
,   E_CDM_ERROR_READ_DFI_GLOBALORIGIN            = 1000 ///< DFI GlobalOrigin 読込みエラー
,   E_CDM_ERROR_READ_DFI_GLOBALREGION            = 1001 ///< DFI GlobalRegion 読込みエラー
,   E_CDM_ERROR_READ_DFI_GLOBALVOXEL             = 1002 ///< DFI GlobalVoxel 読込みエラー
,   E_CDM_ERROR_READ_DFI_GLOBALDIVISION          = 1003 ///< DFI GlobalDivision 読込みエラー
,   E_CDM_ERROR_READ_DFI_DIRECTORYPATH           = 1004 ///< DFI DirectoryPath 読込みエラー
,   E_CDM_ERROR_READ_DFI_TIMESLICEDIRECTORY      = 1005 ///< DFI TimeSliceDirectoryPath 読込みエラー
,   E_CDM_ERROR_READ_DFI_PREFIX                  = 1006 ///< DFI Prefix 読込みエラー
,   E_CDM_ERROR_READ_DFI_FILEFORMAT              = 1007 ///< DFI FileFormat 読込みエラー
,   E_CDM_ERROR_READ_DFI_GUIDECELL               = 1008 ///< DFI GuideCell 読込みエラー
,   E_CDM_ERROR_READ_DFI_DATATYPE                = 1009 ///< DFI DataType 読込みエラー
,   E_CDM_ERROR_READ_DFI_ENDIAN                  = 1010 ///< DFI Endian 読込みエラー
,   E_CDM_ERROR_READ_DFI_NUMVARIABLES            = 1012 ///< DFI NumVariables 読込みエラー
,   E_CDM_ERROR_READ_DFI_FILEPATH_PROCESS        = 1013 ///< DFI FilePath/Process 読込みエラー
,   E_CDM_ERROR_READ_DFI_NO_RANK                 = 1014 ///< DFI Rank要素なし
,   E_CDM_ERROR_READ_DFI_ID                      = 1015 ///< DFI ID 読込みエラー
,   E_CDM_ERROR_READ_DFI_HOSTNAME                = 1016 ///< DFI HoatName 読込みエラー
,   E_CDM_ERROR_READ_DFI_VOXELSIZE               = 1017 ///< DFI VoxelSize 読込みエラー
,   E_CDM_ERROR_READ_DFI_HEADINDEX               = 1018 ///< DFI HeadIndex 読込みエラー
,   E_CDM_ERROR_READ_DFI_TAILINDEX               = 1019 ///< DFI TailIndex 読込みエラー
,   E_CDM_ERROR_READ_DFI_CELLID                  = 1020 ///< DFI CellID 読込みエラー
,   E_CDM_ERROR_READ_DFI_BCFLAGID                = 1021 ///< DFI BCflagID 読込みエラー
,   E_CDM_ERROR_READ_DFI_NO_SLICE                = 1022 ///< DFI TimeSlice要素なし
,   E_CDM_ERROR_READ_DFI_STEP                    = 1023 ///< DFI Step 読込みエラー
,   E_CDM_ERROR_READ_DFI_TIME                    = 1024 ///< DFI Time 読込みエラー
,   E_CDM_ERROR_READ_DFI_NO_MINMAX               = 1025 ///< DFI MinMax要素なし
,   E_CDM_ERROR_READ_DFI_MIN                     = 1026 ///< DFI Min 読込みエラー
,   E_CDM_ERROR_READ_DFI_MAX                     = 1027 ///< DFI Max 読込みエラー
,   E_CDM_ERROR_READ_DFI_DFITYPE                 = 1028 ///< DFI DFIType 読込みエラー
,   E_CDM_ERROR_READ_DFI_FIELDFILENAMEFORMAT     = 1029 ///< DFI FieldfilenameFormat 読込みエラー
,   E_CDM_ERROR_READ_DFI_COORDINATEFILE          = 1030 ///< DFI Coordinate File 読込みエラー
,   E_CDM_ERROR_READ_DFI_COORDINATEFILETYPE      = 1031 ///< DFI Coordinate File Type 読込みエラー
,   E_CDM_ERROR_READ_DFI_COORDINATEFILEPRECISION = 1032 ///< DFI Coordinate File Precision 読込みエラー
,   E_CDM_ERROR_READ_DFI_COORDINATEFILEENDIAN    = 1033 ///< DFI Coordinate File Endian 読込みエラー
,   E_CDM_ERROR_OPEN_COORDINATEFILE              = 1034 ///< Coordinate File オープンに失敗
,   E_CDM_ERROR_READ_COORDINATEFILE              = 1035 ///< Coordinate File 読込みエラー
,   E_CDM_ERROR_READ_INDEXFILE_OPENERROR         = 1050 ///< Indexファイルオープンエラー
,   E_CDM_ERROR_TEXTPARSER                       = 1051 ///< TextParserエラー
,   E_CDM_ERROR_READ_FILEINFO                    = 1052 ///< FileInfo読込みエラー
,   E_CDM_ERROR_READ_FILEPATH                    = 1053 ///< FilePath読込みエラー
,   E_CDM_ERROR_READ_UNIT                        = 1054 ///< UNIT読込みエラー
,   E_CDM_ERROR_READ_TIMESLICE                   = 1055 ///< TimeSlice読込みエラー
,   E_CDM_ERROR_READ_PROCFILE_OPENERROR          = 1056 ///< Procファイルオープンエラー
,   E_CDM_ERROR_READ_DOMAIN                      = 1057 ///< Domain読込みエラー
,   E_CDM_ERROR_READ_MPI                         = 1058 ///< MPI読込みエラー
,   E_CDM_ERROR_READ_PROCESS                     = 1059 ///< Process読込みエラー
//20150918.NetCDF.s
,   E_CDM_ERROR_READ_DFI_NETCDF                  = 1060 ///< NetCDF読込みエラー
,   E_CDM_ERROR_READ_NETCDF_MISMATCH_TYPE        = 1061 ///< DFIとNetCDFのデータ型の不一致エラー
//20150918.NetCDF.e
//20160329.fub.s
,   E_CDM_ERROR_READ_FILELIST_ID                 = 1062 ///< FileListのID読み込みエラー
,   E_CDM_ERROR_READ_FILELIST_CCORDINATEFILENAME = 1063 ///< FileListのCoordinateFileName読み込みエラー
,   E_CDM_ERROR_READ_FILELIST_FIELDDATAFILENAME  = 1064 ///< FileListのFieldDataFileName読み込みエラー
//20160329.fub.e
//20160419.fub.s
,   E_CDM_ERROR_UNDEFINED_FIELDFILENAMEFORMAT    = 1065 ///< FieldFileNameFormatがない
,   E_CDM_ERROR_READ_FIELDDATA_FILE              = 1900 ///< フィールドデータファイル読込みエラー
,   E_CDM_ERROR_READ_SPH_FILE                    = 2000 ///< SPHファイル読込みエラー
,   E_CDM_ERROR_READ_SPH_REC1                    = 2001 ///< SPHファイルレコード1読込みエラー
,   E_CDM_ERROR_READ_SPH_REC2                    = 2002 ///< SPHファイルレコード2読込みエラー
,   E_CDM_ERROR_READ_SPH_REC3                    = 2003 ///< SPHファイルレコード3読込みエラー
,   E_CDM_ERROR_READ_SPH_REC4                    = 2004 ///< SPHファイルレコード4読込みエラー
,   E_CDM_ERROR_READ_SPH_REC5                    = 2005 ///< SPHファイルレコード5読込みエラー
,   E_CDM_ERROR_READ_SPH_REC6                    = 2006 ///< SPHファイルレコード6読込みエラー
,   E_CDM_ERROR_READ_SPH_REC7                    = 2007 ///< SPHファイルレコード7読込みエラー
,   E_CDM_ERROR_UNMATCH_VOXELSIZE                = 2050 ///< SPHのボクセルサイズとDFIのボクセルサイズが合致しない
,   E_CDM_ERROR_NOMATCH_ENDIAN                   = 2051 ///< 出力Fornatが合致しない（Endian形式がBig,Little以外）
,   E_CDM_ERROR_UNMATCH_NUM_OF_VARIABLES         = 2052 ///< フィールドデータの変数の個数と登録された変数名の個数が合致しない
,   E_CDM_ERROR_NUM_OF_GUIDECELLS                = 2053 ///< フィールドデータのガイドセル値エラー
,   E_CDM_ERROR_READ_BOV_FILE                    = 2100 ///< BOVファイル読込みエラー
,   E_CDM_ERROR_READ_FIELD_HEADER_RECORD         = 2102 ///< フィールドヘッダーレコード読込み失敗
,   E_CDM_ERROR_READ_FIELD_DATA_RECORD           = 2103 ///< フィールドデータレコード読込み失敗
,   E_CDM_ERROR_READ_FIELD_AVERAGED_RECORD       = 2104 ///< フィールドAverage読込み失敗
,   E_CDM_ERROR_READ_DFI_FILELIST                = 2105 ///< FileList読み込み失敗
//,   E_CDM_ERROR_DATATYPE                         = 2500 ///< DataType error
//20150918.NetCDF.s
,   E_CDM_ERROR_READ_NETCDF_FUNC                 = 2200 ///< NetCDFのnc関数でエラー
,   E_CDM_ERROR_READ_NETCDF_VAR_1D               = 2201 ///< NetCDFの1次元配列として読み込むvariableが1次元で無い
//20150918.NetCDF.e
//20160329.fub.s
,   E_CDM_ERROR_READ_FUB_REC1                    = 2210 ///< fubの第1レコード読み込みエラー
,   E_CDM_ERROR_READ_FUB_REC2                    = 2211 ///< fubの第2レコード読み込みエラー
//20160329.fub.e
,   E_CDM_ERROR_MISMATCH_NP_SUBDOMAIN            = 3003 ///< 並列数とサブドメイン数が一致していない
,   E_CDM_ERROR_INVALID_DIVNUM                   = 3011 ///< 領域分割数が不正
,   E_CDM_ERROR_OPEN_SBDM                        = 3012 ///< ActiveSubdomainファイルのオープンに失敗
,   E_CDM_ERROR_READ_SBDM_HEADER                 = 3013 ///< ActiveSubdomainファイルのヘッダー読み込みに失敗
,   E_CDM_ERROR_READ_SBDM_FORMAT                 = 3014 ///< ActiveSubdomainファイルのフォーマットエラー
,   E_CDM_ERROR_READ_SBDM_DIV                    = 3015 ///< ActiveSubdomainファイルの領域分割数読み込みに失敗
,   E_CDM_ERROR_READ_SBDM_CONTENTS               = 3016 ///< ActiveSubdomainファイルのContents読み込みに失敗
,   E_CDM_ERROR_SBDM_NUMDOMAIN_ZERO              = 3017 ///< ActiveSubdomainファイルの活性ドメイン数が0
,   E_CDM_ERROR_MAKEDIRECTORY                    = 3100 ///< Directory生成で失敗
,   E_CDM_ERROR_OPEN_FIELDDATA                   = 3101 ///< フィールドデータのオープンに失敗
,   E_CDM_ERROR_WRITE_FIELD_HEADER_RECORD        = 3102 ///< フィールドヘッダーレコード出力失敗
,   E_CDM_ERROR_WRITE_FIELD_DATA_RECORD          = 3103 ///< フィールドデータレコード出力失敗
,   E_CDM_ERROR_WRITE_FIELD_AVERAGED_RECORD      = 3104 ///< フィールドAverage出力失敗
,   E_CDM_ERROR_WRITE_SPH_REC1                   = 3201 ///< SPHファイルレコード1出力エラー
,   E_CDM_ERROR_WRITE_SPH_REC2                   = 3202 ///< SPHファイルレコード2出力エラー
,   E_CDM_ERROR_WRITE_SPH_REC3                   = 3203 ///< SPHファイルレコード3出力エラー
,   E_CDM_ERROR_WRITE_SPH_REC4                   = 3204 ///< SPHファイルレコード4出力エラー
,   E_CDM_ERROR_WRITE_SPH_REC5                   = 3205 ///< SPHファイルレコード5出力エラー
,   E_CDM_ERROR_WRITE_SPH_REC6                   = 3206 ///< SPHファイルレコード6出力エラー
,   E_CDM_ERROR_WRITE_SPH_REC7                   = 3207 ///< SPHファイルレコード7出力エラー
,   E_CDM_ERROR_WRITE_PROCFILENAME_EMPTY         = 3500 ///< proc dfi ファイル名が未定義
,   E_CDM_ERROR_WRITE_PROCFILE_OPENERROR         = 3501 ///< proc dfi ファイルオープン失敗
,   E_CDM_ERROR_WRITE_DOMAIN                     = 3502 ///< Domain出力失敗
,   E_CDM_ERROR_WRITE_MPI                        = 3503 ///< MPI出力失敗
,   E_CDM_ERROR_WRITE_PROCESS                    = 3504 ///< Process出力失敗
,   E_CDM_ERROR_WRITE_RANKID                     = 3505 ///< 出力ランク以外
,   E_CDM_ERROR_WRITE_INDEXFILENAME_EMPTY        = 3510 ///< index dfi ファイル名が未定義
,   E_CDM_ERROR_WRITE_PREFIX_EMPTY               = 3511 ///< Prefixが未定義
,   E_CDM_ERROR_WRITE_INDEXFILE_OPENERROR        = 3512 ///< proc dfi ファイルオープン失敗
,   E_CDM_ERROR_WRITE_FILEINFO                   = 3513 ///< FileInfo出力失敗
,   E_CDM_ERROR_WRITE_UNIT                       = 3514 ///< Unit出力失敗
,   E_CDM_ERROR_WRITE_TIMESLICE                  = 3515 ///< TimeSlice出力失敗
,   E_CDM_ERROR_WRITE_FILEPATH                   = 3516 ///< FilePath出力失敗
,   E_CDM_ERROR_WRITE_VISIT                      = 3517 ///< Visit出力失敗
,   E_CDM_ERROR_WRITE_GRIDFILE                   = 3518 ///< gridファイル出力失敗
//20150918.NetCDF.s
,   E_CDM_ERROR_WRITE_DFI_NETCDF                 = 3600 ///< NetCDFのDFI出力エラー
//20150918.NetCDF.e
,   E_CDM_WARN_GETUNIT                           = 4000 ///< Unitの単位がない
  };

};

/** 3次元（スカラー）インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _CDM_IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (long long)((_K)+(_VC)) * (long long)((_NI)+2*(_VC)) * (long long)((_NJ)+2*(_VC)) \
+ (long long)((_J)+(_VC)) * (long long)((_NI)+2*(_VC)) \
+ (long long)((_I)+(_VC)) \
)

/** 2次元（スカラー）インデクス(i,j) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _CDM_IDX_IJ(_I,_J,_NI,_NJ,_VC) \
( (long long)((_J)+(_VC)) * (long long)((_NI)+2*(_VC)) \
+ (long long)((_I)+(_VC)) \
)

/** 2次元（スカラー）インデクス(n,i,j) -> 1次元インデクス変換マクロ
 *  @param[in] _N  変数インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NN 変数の個数
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _CDM_IDX_NIJ(_N,_I,_J,_NI,_NJ,_NN,_VC) \
( (long long)(_NN)*_CDM_IDX_IJ(_I,_J,_NI,_NJ,_VC) \
 + (long long)(_N) \
)


/** 3次元（ベクトル）インデクス(i,j,k,n) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  変数インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _CDM_IDX_IJKN(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
( (long long)(_N) * (long long)((_NI)+2*(_VC)) * (long long)((_NJ)+2*(_VC)) \
* (long long)((_NK)+2*(_VC)) \
+ _CDM_IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
)

/** 3次元（ベクトル）インデクス(n,i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _N  変数インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NN 変数の個数
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _CDM_IDX_NIJK(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
( (long long)(_NN) * _CDM_IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
+ (long long)(_N) )

/** DFIファイルのTab出力
 *  @param[in] _FP  ファイルポインタ
 *  @param[in] _NTAB インデント数
 */
#define _CDM_WRITE_TAB(_FP,_NTAB) {\
 for(int _NTCNT=0; _NTCNT<_NTAB; _NTCNT++) fprintf(_FP,_CDM_TAB_STR); \
}

#ifndef stmpd_printf
#define stmpd_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#endif

#endif /* _CDM_DEFINE_H_ */
