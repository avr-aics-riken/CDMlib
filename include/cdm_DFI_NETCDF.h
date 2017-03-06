#ifndef _CDM_DFI_NETCDF_H_
#define _CDM_DFI_NETCDF_H_
#ifdef _WITH_NETCDF4_

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

/**
 * @file   cdm_DFI_NETCDF.h
 * @brief  cdm_DFI_NETCDF Class Header
 * @author aics
 */

#include "cdm_DFI.h"
#include "netcdf.h"

class cdm_DFI_NETCDF : public cdm_DFI
{
public:
  struct stDimInfo
  {
    int    id;   ///< dimension ID
    string name; ///< dimension name
    size_t len;  ///< dimension
  };

  struct stVarInfo
  {
    int               id;     ///< variable ID
    string            name;   ///< variable name
    nc_type           type;   ///< data type
    vector<stDimInfo> dims;   ///< array dimension
    string            unit;   ///< variable unit name
    double*           data;   ///< data array for 1D array
    int               endian; ///< endian type

    stVarInfo()
    {
      data = NULL;
      clear();
    }

    ~stVarInfo()
    {
      clear();
    }

    void clear()
    {
      id = 0;
      name = "";
      type = 0;
      dims.clear();
      unit = "";
      if( data ) delete [] data;
      data = NULL;
      endian = NC_ENDIAN_NATIVE;
    }
  };

protected:

  string m_nameX; ///< ncファイルのX座標値配列の配列名
  string m_nameY; ///< ncファイルのY座標値配列の配列名
  string m_nameZ; ///< ncファイルのZ座標値配列の配列名
  string m_nameT; ///< ncファイルの時刻値配列の配列名

  stVarInfo m_varInfoX; ///< xの変数情報
  stVarInfo m_varInfoY; ///< yの変数情報
  stVarInfo m_varInfoZ; ///< zの変数情報
  stVarInfo m_varInfoT; ///< tの変数情報

  cdm_Unit m_NC_Unit; ///< ncファイル内のunit情報

  /** データレコード出力用のvar id(変数定義順に格納)
   *  write_HeaderRecordで格納し、write_DataRecordで配列出力時に使用
   */
  vector<stVarInfo> m_vecVarInfo;

  bool m_writeFlag; ///< 書き込み済みフラグ


public:

  /** コンストラクタ */
  cdm_DFI_NETCDF();

  /**
   * @brief コンストラクタ
   * @param [in] F_Info  FileInfo
   * @param [in] F_Path  FilePath
   * @param [in] visit   VisIt option
   * @param [in] unit    Unit
   * @param [in] domain  Domain
   * @param [in] mpi     MPI
   * @param [in] TSlice  TimeSlice
   * @param [in] process Process
   */
  cdm_DFI_NETCDF(const cdm_FileInfo F_Info,
                 const cdm_FilePath F_Path,
                 const cdm_VisIt visit,
                 const cdm_Unit unit,
                 const cdm_Domain* domain,
                 const cdm_MPI mpi,
                 const cdm_TimeSlice TSlice,
                 const cdm_Process process)
  {
    DFI_Finfo      = F_Info;
    DFI_Fpath      = F_Path;
    DFI_VisIt      = visit;
    DFI_Unit       = unit;
    DFI_Domain     = domain;
    DFI_MPI        = mpi;
    DFI_TimeSlice  = TSlice;
    DFI_Process    = process;
    m_bgrid_interp_flag = false;

    m_nameX = "x";
    m_nameY = "y";
    m_nameZ = "z";
    m_nameT = "time";

    m_writeFlag = false;
  };

  /** デストラクタ */
  ~cdm_DFI_NETCDF();

public:


  /**
   * @brief NetCDFの次元変数名のセット
   * @param[in] nameX 座標値Xの次元名
   * @param[in] nameY 座標値Yの次元名
   * @param[in] nameZ 座標値Zの次元名
   * @param[in] nameT 時刻の次元名
   */
  void
  SetDimName( string nameX, string nameY, string nameZ, string nameT );

  /**
   * @brief NetCDFの次元変数名の取得
   * @param[out] nameX 座標値Xの次元名
   * @param[out] nameY 座標値Yの次元名
   * @param[out] nameZ 座標値Zの次元名
   * @param[out] nameT 時刻の次元名
   */
  void
  GetDimName( string &nameX, string &nameY, string &nameZ, string &nameT );

  /**
   * @brief NetCDF特有のdfiパラメータ読み込み
   * @param [in] tpCntl cdm_TextParserクラス
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  ReadAdditionalTP( cdm_TextParser tpCntl );

  /**
   * @brief NetCDF特有のdfiパラメータ出力
   * @param [in] fp ファイルポインタ
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  WriteAdditionalTP( FILE *fp, int tab );

  /**
   * @brief Uuitをセットする
   * @param [in] Name       追加する単位系("Length","Velocity",,,,)
   * @param [in] Unit       単位ラベル("M","CM","MM","M/S",,,)
   * @param [in] reference  規格化したスケール値
   * @param [in] difference 差の値
   * @param [in] BsetDiff   differenceの有無
   */
  void
  AddNcUnit(const std::string Name,
            const std::string Unit,
            const double reference,
            const double difference= 0.0,
            const bool BsetDiff=false);

  /**
   * @brief cdm_Unitクラスのポインタを取得
   * @return cdm_Unitクラスポインタ
   */
  const cdm_Unit* GetNcUnit();

  /**
   * @brief 出力処理を追記モードにするかどうかをチェック(NetCDF用)
   * @return モード(true:追記モード、false:新規作成モード)
   */
  virtual
  bool CheckAddWriteMode();

  /**
   * @brief 書き込み済みフラグのセット(FCONV用)
   * @param[in] writeFlag 書き込み済みフラグ
   */
  void
  SetWriteFlag( bool writeFlag )
  {
    m_writeFlag = writeFlag;
  }

  /**
   * @brief ファイルのヘッダーレコード読込み(static関数)
   * @param[in]  ncid  ncid
   * @param[in]  nameX xの次元,変数名
   * @param[in]  nameY yの次元,変数名
   * @param[in]  nameZ zの次元,変数名
   * @param[in]  nameT timeの次元,変数名
   * @param[out] varInfoX x座標情報
   * @param[out] varInfoY y座標情報
   * @param[out] varInfoZ z座標情報
   * @param[out] varInfoT time情報
   * @return error code
   */
  static CDM::E_CDM_ERRORCODE
  read_HeaderRecord(int ncid,
                    std::string nameX,
                    std::string nameY,
                    std::string nameZ,
                    std::string nameT,
                    stVarInfo &varInfoX,
                    stVarInfo &varInfoY,
                    stVarInfo &varInfoZ,
                    stVarInfo &varInfoT);

  /**
   * @brief ncデータファイルのデータレコード読込み(static関数)
   * @param[in]    ncid         ncid
   * @param[in]    step_index   読み込むステップのインデクス
   * @param[in]    VariableName 読み込む変数名リスト
   * @param[in]    DataType     読み込む変数のデータ型
   * @param[in]    varInfoX     x座標情報
   * @param[in]    varInfoY     y座標情報
   * @param[in]    varInfoZ     z座標情報
   * @param[in]    varInfoT     time情報
   * @param[out]   vecVarInfo   読み込んだ変数のvarInfo
   * @param[out]   ret          error code
   * @return cdm_Arrayのポインタ
   */
  static cdm_Array*
  read_Datarecord(int ncid,
                  unsigned step_index,
                  std::vector<std::string> &VariableName,
                  CDM::E_CDM_DTYPE DataType,
                  stVarInfo &varInfoX,
                  stVarInfo &varInfoY,
                  stVarInfo &varInfoZ,
                  stVarInfo &varInfoT,
                  std::vector<stVarInfo> &vecVarInfo,
                  CDM::E_CDM_ERRORCODE &ret);

  /**
   * @brief ncファイルのヘッダレコードの出力(static関数)
   * @param[in] pFile  ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return error code
   */
  static CDM::E_CDM_ERRORCODE
  write_HeaderRecord(cdm_FILE* pFile,
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
                     vector<stVarInfo> &vecVarInfo);

  /**
   * @brief ncファイルのデータレコードの出力(static関数)
   * @param[in]  pFile ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  gc ガイドセル
   * @param[in]  RankID ランク番号
   * @return error code
   */
  static CDM::E_CDM_ERRORCODE
  write_DataRecord(cdm_FILE* pFile,
                   cdm_Array* val,
                   const int gc,
                   int VoxelSize[3],
                   stVarInfo &varInfoT,
                   vector<stVarInfo> &vecVarInfo);

  /**
   * @brief CDMデータ型に対応するncデータ型を取得
   * @param[in] type ncデータ型
   * @return    CDMデータ型
   */
  static nc_type
  GetNcType( CDM::E_CDM_DTYPE type );



protected:

  /**
   * @brief nc変数の情報を取得する
   * @param[in]  ncid      ファイルID(netcdf library)
   * @param[in]  var_name  情報取得対象の変数名
   * @param[out] varInfo   変数の情報
   * @param[in]  nDimCheck 変数の配列次元数のチェック(0以下:しない、1以上:その数で無いときエラー)
   * @return error code
   */
  static CDM::E_CDM_ERRORCODE
  GetVarInfo( int ncid, string var_name, stVarInfo &varInfo, int nDimCheck=0 );

  /**
   * @brief ncデータ型に対応するCDMデータ型を取得
   * @param[in] type ncデータ型
   * @return    CDMデータ型
   */
  static CDM::E_CDM_DTYPE
  GetCdmType( nc_type type );

  /**
   * @brief nc変数配列のdimensionをチェック
   *        - x,y,z,timeと同じdimensionかどうか
   *        - (z,y,x)もしくは(time,z,y,x)の形状かどうか
   * @param[in] verInfo
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  CheckArrayVarInfo( stVarInfo &varInfo );

  /**
   * @brief UnitListに変数名のunitが存在する場合、var属性として出力する
   * @param[in] ncid     ncid
   * @param[in] varid    var ID
   * @param[in] var_name var name
   * @param[in] unit     unitlist
   */
  static void
  write_AttUnits( int ncid, int varid, string var_name, cdm_Unit &unit );

  /**
   * @brief ncファイルのヘッダーレコード読込み
   * @param[in]  pFile       ファイルポインタ
   * @param[in]  matchEndian エンディアンチェックフラグ true:合致
   * @param[in]  step        ステップ番号
   * @param[in]  head        dfiのHeadIndex
   * @param[in]  tail        dfiのTailIndex
   * @param[in]  gc          dfiのガイドセル数
   * @param[out] voxsize     voxsize
   * @param[out] time        時刻
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
//read_HeaderRecord(int ncid,
  read_HeaderRecord(cdm_FILE *pFile,
                    bool matchEndian,
                    unsigned step,
                    const int head[3],
                    const int tail[3],
                    int gc,
                    int voxsize[3],
                    double &time);

  /**
   * @brief ncデータファイルのデータレコード読込み
   * @param[in]  pFile       ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        ステップ番号
   * @param[in]  buf         読込み用バッファ
   * @param[in]  head        読込みバッファHeadIndex
   * @param[in]  nz          z方向のボクセルサイズ（実セル＋ガイドセル＊２）
   * @param[out] src         読み込んだデータを格納した配列のポインタ
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
//read_Datarecord(int ncid,
  read_Datarecord(cdm_FILE *pFile,
                  bool matchEndian,
                  unsigned step,
                  cdm_Array* buf,
                  int head[3],
                  int nz,
                  cdm_Array* &src);

  /**
   * @brief ncファイルのAverageデータレコードの読込み
   * @param[in]  pFile       ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        読込みstep番号
   * @param[out] avr_step    平均ステップ
   * @param[out] avr_time    平均タイム
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
//read_averaged(int ncid,
  read_averaged(cdm_FILE *pFile,
                bool matchEndian,
                unsigned step,
                unsigned &avr_step,
                double &avr_time);

  /**
   * @brief ncファイルのヘッダレコードの出力
   * @param[in] pFile  ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
//write_HeaderRecord(FILE* fp,
  write_HeaderRecord(cdm_FILE* pFile,
                     const unsigned step,
                     const double time,
                     const int RankID);

  /**
   * @brief ncファイルのデータレコードの出力
   * @param[in]  pFile ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  gc ガイドセル
   * @param[in]  RankID ランク番号
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
//write_DataRecord(FILE* fp,
  write_DataRecord(cdm_FILE* pFile,
                   cdm_Array* val,
                   const int gc,
                   const int RankID);

  /**
   * @brief Averageレコードの出力
   * @param[in] pFile    ファイルポインタ
   * @param[in] step_avr 平均ステップ番号
   * @param[in] time_avr 平均時刻
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
//write_averaged(FILE* fp,
  write_averaged(cdm_FILE* pFile,
                 const unsigned step_avr,
                 const double time_avr);

};

#endif /* _WITH_NETCDF4_ */
#endif /* _CDM_DFI_NETCDF_H_ */
