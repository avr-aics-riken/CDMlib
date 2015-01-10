#ifndef _CDM_DFI_H_
#define _CDM_DFI_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI.h
 * @brief  cdm_DFI Class Header
 * @author aics
 */

#include "cdm_Define.h"
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <set>
#include <map>
#include <string>

#include "cdm_Version.h"


#include "cdm_PathUtil.h"
#include "cdm_TextParser.h"
#include "cdm_ActiveSubDomain.h"
#include "cdm_endianUtil.h"
#include "cdm_TypeArray.h"

#include "cdm_FileInfo.h"
#include "cdm_FilePath.h"
#include "cdm_VisIt.h"
#include "cdm_Unit.h"
#include "cdm_TimeSlice.h"
#include "cdm_Domain.h"
#include "cdm_MPI.h"
#include "cdm_Process.h"

/** CDM main class */
class cdm_DFI {

public:

protected :
  MPI_Comm          m_comm;          ///< MPI コミュニケータ
  std::string       m_directoryPath; ///< index dfi ファイルのディレクトリパス
  std::string       m_indexDfiName;  ///< index dfi ファイル名
  CDM::E_CDM_READTYPE m_read_type;   ///< 読込みタイプ

  int               m_RankID;        ///< ランク番号

  cdm_FileInfo      DFI_Finfo;       ///< FileInfo class
  cdm_FilePath      DFI_Fpath;       ///< FilePath class
  cdm_VisIt         DFI_VisIt;       ///< VisIt class
  cdm_Unit          DFI_Unit;        ///< Unit class
  const cdm_Domain* DFI_Domain;      ///< Domain class
  cdm_MPI           DFI_MPI;         ///< MPI class
  cdm_TimeSlice     DFI_TimeSlice;   ///< TimeSlice class
  cdm_Process       DFI_Process;     ///< Process class

  vector<int>m_readRankList;         ///< 読込みランクリスト

  bool m_bgrid_interp_flag;                  ///< 節点への補間フラグ
  CDM::E_CDM_FILE_TYPE  m_input_type;        ///< 入力形式(ascii,binary,FortarnBinary)
  CDM::E_CDM_FILE_TYPE  m_output_type;       ///< 出力形式(ascii,binary,FortarnBinary)
  CDM::E_CDM_FILE_TYPE  m_output_type_coord; ///< 座標データの出力形式(ascii,binary)
  CDM::E_CDM_OUTPUT_FNAME m_output_fname;    ///< 出力ファイル命名規約(step_rank,rank_step)   

public:
  /** コンストラクタ */
  cdm_DFI();
  
  /**　デストラクタ */
  virtual ~cdm_DFI();

  /**
   * @brief read インスタンス(BOVもしくはPLOT3Dの場合にインスタンス生成)
   * @param [in]  comm    MPIコミュニケータ
   * @param [in]  dfifile DFIファイル名
   * @param [in]  G_Voxel 計算空間全体のボクセルサイズ
   * @param [in]  G_Div   計算空間の領域分割数
   * @param [out] ret     終了コード
   * @return インスタンスされたクラスのポインタ
   */
  static cdm_DFI*
  ReadInit(const MPI_Comm comm, 
           const std::string dfifile,
           const int G_Voxel[3],
           const int G_Div[3],
           CDM::E_CDM_ERRORCODE &ret); 

  /**
   * @brief cdmFileInfoクラスのポインタを取得
   * @return cdm_FileInfoクラスポインタ
   */
  const cdm_FileInfo* GetcdmFileInfo();

  /**
   * @brief cdm_FilePathクラスのポインタを取得
   * @return cdm_FilePathクラスポインタ
   */
  const cdm_FilePath* GetcdmFilePath();  

  /**
   * @brief cdm_FilePathクラスのセット
   */
  void SetcdmFilePath(cdm_FilePath FPath); 

  /**
   * @brief cdm_VisItクラスのポインタを取得
   * @return cdm_VisItクラスポインタ
   */
  const cdm_VisIt* GetcdmVisIt();
  
  /**
   * @brief cdm_VisItクラスのセット
   */
  void SetcdmVisIt(cdm_VisIt Visit);
  
  /**
   * @brief cdm_Unitクラスのポインタを取得
   * @return cdm_Unitクラスポインタ
   */
  const cdm_Unit* GetcdmUnit(); 

  /**
   * @brief cdm_Unitクラスのセット
   */
  void SetcdmUnit(cdm_Unit unit); 


  /**
   * @brief cdm_Domainクラスのポインタ取得
   * @return cdm_Domainクラスポインタ
   */
  const cdm_Domain* GetcdmDomain(); 

  /**
   * @brief cdm_Domainクラスのセット
   */
  void SetcdmDomain(cdm_Domain* domain); 

  /**
   * @brief cdm_MPIクラスのポインタ取得
   * @return cdm_MPIクラスポインタ 
   */
  const cdm_MPI* GetcdmMPI();

  /**
   * @brief cdm_MPIクラスセット
   */
  void SetcdmMPI(cdm_MPI mpi);

  /**
   * @brief cdm_TimeSliceクラスのポインタ取得
   * @return cdm_TimeSliceクラスポインタ
   */
  const cdm_TimeSlice* GetcdmTimeSlice(); 

  /**
   * @brief cdm_TimeSlice クラスセット
   */
  void SetcdmTimeSlice(cdm_TimeSlice TSlice);
 
  /**
   * @brief cdm_Processクラスのポインタ取得
   * @return cdm_Processクラスポインタ
   */
  const cdm_Process* GetcdmProcess(); 

  /**
   * @brief cdm_Processクラスセット
   */
  void SetcdmProcess(cdm_Process Process);

  /**
   * @brief 出力DFIファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @return DFIファイル名
   */ 
  static std::string
  Generate_DFI_Name(const std::string prefix);

  /**
   * @brief フィールドデータ（SPH,BOV)ファイル名の作成(ディレクトリパスが付加されている）
   * @param [in] RankID ランク番号
   * @param [in] step   読込みステップ番号
   * @param [in] mio    並列判定フラグ（逐次or並列の判定用）
   * @return 生成されたファイル名　　　　　　　
   */
  std::string Generate_FieldFileName(int RankID,
                                int step, 
                                const bool mio);

  /**
   * @brief ファイル名生成
   * @param [in] prefix ベースファイル名
   * @param [in] RankID ランク番号  
   * @param [in] step   出力ステップ番号（負のとき、ステップ番号が付加されない）
   * @param [in] ext    拡張子
   * @param [in] output_fname step_rank,rank_step指示
   * @param [in] mio    並列判定フラグ
   * @param [in] TimeSliceDirFlag Time Slice 毎の出力指示
   * @return 生成されたファイル名
   */
  static
  std::string Generate_FileName(std::string prefix,
                                int RankID,
                                int step,
                                std::string ext,
                                CDM::E_CDM_OUTPUT_FNAME output_fname,
                                bool mio,
                                CDM::E_CDM_ONOFF TimeSliceDirFlag);

  /**
   * @brief write インスタンス template function (等間隔格子用)
   * @param [in] comm        MPIコミュニケータ
   * @param [in] DfiName     DFIファイル名
   * @param [in] Path        フィールドデータのディレクトリ
   * @param [in] prefix      ベースファイル名
   * @param [in] format      ファイルフォーマット
   * @param [in] GCell       出力仮想セル数　　　
   * @param [in] DataType    データタイプ　　　　
   * @param [in] nVari       変数の個数　　　　　　　
   * @param [in] proc_fname  proc.dfiファイル名
   * @param [in] G_size      グローバルボクセルサイズ　
   * @param [in] pitch       ピッチ　　　　　　　　　　
   * @param [in] G_origin    原点座標値　　　　　　　　
   * @param [in] division    領域分割数　　　　　　　　
   * @param [in] head        計算領域の開始位置　　　　
   * @param [in] tail        計算領域の終了位置　　　　
   * @param [in] hostname    ホスト名
   * @param [in] TSliceOnOff TimeSliceフラグ
   * @return インスタンスされたクラスのポインタ
   */
  template<typename T>
  static cdm_DFI*
  WriteInit(const MPI_Comm comm,
            const std::string DfiName,
            const std::string Path,
            const std::string prefix,
            const CDM::E_CDM_FORMAT format,
            const int GCell,
            const CDM::E_CDM_DTYPE DataType,
            const int nVari,
            const std::string proc_fname,
            const int G_size[3],
            const T pitch[3],
            const T G_origin[3],
            const int division[3],
            const int head[3],
            const int tail[3],
            const std::string hostname,
            const CDM::E_CDM_ONOFF TSliceOnOff);

  /**
   * @brief write インスタンス template function (不等間隔格子用)
   * @details templateの型より、座標ファイルのデータ精度を指定
   * @param [in] comm        MPIコミュニケータ
   * @param [in] DfiName     DFIファイル名
   * @param [in] Path        フィールドデータのディレクトリ
   * @param [in] prefix      ベースファイル名
   * @param [in] format      ファイルフォーマット
   * @param [in] GCell       出力仮想セル数　　　
   * @param [in] DataType    データタイプ　　　　
   * @param [in] nVari       変数の個数　　　　　　　
   * @param [in] proc_fname  proc.dfiファイル名
   * @param [in] G_size      グローバルボクセルサイズ　
   * @param [in] coord_X     X座標データポインタ
   * @param [in] coord_Y     Y座標データポインタ
   * @param [in] coord_Z     Z座標データポインタ
   * @param [in] coord_file          座標ファイル名
   * @param [in] coord_filetype      座標ファイルのファイルタイプ
   * @param [in] coord_fileEndian    座標ファイルのエンディアンタイプ
   * @param [in] division    領域分割数　　　　　　　　
   * @param [in] head        計算領域の開始位置　　　　
   * @param [in] tail        計算領域の終了位置　　　　
   * @param [in] hostname    ホスト名　　　　　　　　　
   * @param [in] TSliceOnOff TimeSliceフラグ
   * @return インスタンスされたクラスのポインタ
   */
  template<typename T>
  static cdm_DFI*
  WriteInit(const MPI_Comm comm,
            const std::string DfiName,
            const std::string Path,
            const std::string prefix,
            const CDM::E_CDM_FORMAT format,
            const int GCell,
            const CDM::E_CDM_DTYPE DataType,
            const int nVari,
            const std::string proc_fname,
            const int G_size[3],
            const T* coord_X,
            const T* coord_Y,
            const T* coord_Z,
            const std::string coord_file,
            const CDM::E_CDM_FILE_TYPE coord_filetype,
            const CDM::E_CDM_ENDIANTYPE coord_fileEndian,
            const int division[3],
            const int head[3],
            const int tail[3],
            const std::string hostname,
            const CDM::E_CDM_ONOFF TSliceOnOff);

  /**
   * @brief write インスタンス template function (等間隔格子・不等間隔格子の共通処理部分)
   * @param [in] comm        MPIコミュニケータ
   * @param [in] DfiName     DFIファイル名
   * @param [in] Path        フィールドデータのディレクトリ
   * @param [in] prefix      ベースファイル名
   * @param [in] format      ファイルフォーマット
   * @param [in] GCell       出力仮想セル数　　　
   * @param [in] DataType    データタイプ　　　　
   * @param [in] nVari       変数の個数　　　　　　　
   * @param [in] proc_fname  proc.dfiファイル名
   * @param [in] out_domain  domainインスタンス　　
   * @param [in] head        計算領域の開始位置　　　　
   * @param [in] tail        計算領域の終了位置　　　　
   * @param [in] hostname    ホスト名
   * @param [in] TSliceOnOff TimeSliceフラグ
   * @return インスタンスされたクラスのポインタ
   */
  template<typename T>
  static cdm_DFI*
  WriteInit(const MPI_Comm comm,
            const std::string DfiName,
            const std::string Path,
            const std::string prefix,
            const CDM::E_CDM_FORMAT format,
            const int GCell,
            const CDM::E_CDM_DTYPE DataType,
            const int nVari,
            const std::string proc_fname,
            const cdm_Domain* out_domain,
            const int head[3],
            const int tail[3],
            const std::string hostname,
            const CDM::E_CDM_ONOFF TSliceOnOff);

  /**
   * @brief RankIDをセットする
   * @param[in] rankID RankID
   */
  void set_RankID(const int rankID)
  { m_RankID = rankID; };

  /**
   * @brief 節点への補間フラグをセット
   * @param [in] bgrid_interp_flag 節点への補間フラグ
   */
  void set_interp_flag(bool interp_flag)
  {  m_bgrid_interp_flag = interp_flag; };

  /**
   * @brief 入力形式(ascii,binary,FortranBinary)をセット
   * @param [in] input_type 出力形式
   */
  void set_input_type(CDM::E_CDM_FILE_TYPE input_type)
  {  m_input_type = input_type; };

  /**
   * @brief 出力形式(ascii,binary,FortranBinary)をセット
   * @param [in] output_type 出力形式
   */
  void set_output_type(CDM::E_CDM_FILE_TYPE output_type)
  {  m_output_type = output_type; };

  /**
   * @brief 座標データの出力形式(ascii,binary)をセット
   * @details AVSおよびVTK形式で利用
   * @param [in] output_type_coord 座標データの出力形式
   */
  void set_output_type_coord(CDM::E_CDM_FILE_TYPE output_type_coord)
  {  m_output_type_coord = output_type_coord; };

  /**
   * @brief 出力ファイル命名規約(step_rank,rank_step)をセット
   * @param [in] output_fname 出力ファイル命名規約
   */
  void set_output_fname(CDM::E_CDM_OUTPUT_FNAME output_fname)
  { m_output_fname = output_fname; };

  /**
   * @brief DFIファイル名の取り出し
   * @return dfiファイル名
   */ 
  std::string get_dfi_fname()
  { return m_indexDfiName; };


  /**
   * @brief read field data record (template function)
   * @details 読み込んだデータのポインタを戻り値として返す
   * @param [out] ret       終了コード 1:正常、1以外：エラー  
   * @param [in] step       入力ステップ番号
   * @param [in] gc         仮想セル数　　　
   * @param [in] Gvoxel     グローバルボクセルサイズ　
   * @param [in] Gdivision  領域分割数　　　　　　　　
   * @param [in] head       計算領域の開始位置　　　　
   * @param [in] tail       計算領域の終了位置　　　　
   * @param [out] time      読み込んだ時間
   * @param [in]  mode      平均ステップ＆時間読込みフラグ　false : 読込み
   *                                                           true  : 読み込まない
   * @param [out] step_avr  平均ステップ
   * @param [out] time_avr  平均時間
   * @return 読みんだフィールドデータのポンタ
   */
//  template<class T, class TimeT, class TimeAvrT> T*
  template<class TimeT, class TimeAvrT> void*
  ReadData(CDM::E_CDM_ERRORCODE &ret,
           const unsigned step, 
           const int gc, 
           const int Gvoxel[3], 
           const int Gdivision[3], 
           const int head[3], 
           const int tail[3],
           TimeT &time,
           const bool mode, 
           unsigned &step_avr, 
           TimeAvrT &time_avr);
  /**
   * @brief read field data record (template function)
   * @details 引数で渡された配列ポインタにデータを読込む
   * @param [out] val        読み込んだデータポインタ　
   * @param [in]  step       入力ステップ番号
   * @param [in]  gc         仮想セル数　　　
   * @param [in]  Gvoxel     グローバルボクセルサイズ　
   * @param [in]  Gdivision  領域分割数　　　　　　　　
   * @param [in]  head       計算領域の開始位置　　　　
   * @param [in]  tail       計算領域の終了位置　　　　
   * @param [out] time       読み込んだ時間
   * @param [in]  mode       平均ステップ＆時間読込みフラグ　false : 読込み
   *                                                           true  : 読み込まない
   * @param [out] step_avr      平均ステップ
   * @param [out] time_avr      平均時間
   * @return 終了コード 1:正常 1以外:エラー
   */
  template<class T, class TimeT, class TimeAvrT>
  CDM::E_CDM_ERRORCODE 
  ReadData(T *val,
           const unsigned step,
           const int gc,
           const int Gvoxel[3],
           const int Gdivision[3],
           const int head[3],
           const int tail[3],
           TimeT &time,
           const bool mode,
           unsigned &step_avr,
           TimeAvrT &time_avr);

  /**
   * @brief read field data record 
   * @details template ReadData関数で型に応じた配列を確保した後、呼び出される
   * @param [out] val        読み込み先の配列をポインタで渡す　
   * @param [in]  step       読み込むステップ番号
   * @param [in]  gc         仮想セル数　　　
   * @param [in]  Gvoxel     グローバルボクセルサイズ　
   * @param [in]  Gdivision  領域分割数　　　　　　　　
   * @param [in]  head       計算領域の開始位置　　　　
   * @param [in]  tail       計算領域の終了位置　　　　
   * @param [out] time       読み込んだ時間
   * @param [in]  mode       平均ステップ＆時間読込みフラグ　false : 読込み
   *                                                           true  : 読み込まない
   * @param [out] step_avr      平均ステップ
   * @param [out] time_avr      平均時間
   * @return 終了コード 1:正常 1以外:エラー
   */
  CDM::E_CDM_ERRORCODE 
  ReadData(cdm_Array *val,
           const unsigned step,
           const int gc,
           const int Gvoxel[3],
           const int Gdivision[3],
           const int head[3],
           const int tail[3],
           double &time,
           const bool mode,
           unsigned &step_avr,
           double &time_avr);

  /**
   * @brief write field data record (template function)
   * @details                minmax[0]   =変数1のminX
   *                          minmax[1]   =変数1のmaxX
   *                               ...
   *                          minmax[2n-2]=変数nのminX
   *                          minmax[2n-1]=変数nのmaxX
   *       SPHでnVari=3のとき、minmax[2n  ]=合成値のmin
   *                          minmax[2n+1]=合成値のmax
   * @param [in] step     出力ステップ番号
   * @param [in] time     出力時刻　　　　
   * @param [in] sz       valの実ボクセルサイズ
   * @param [in] nVari    valの変数の個数
   * @param [in] gc       valの仮想セル数　　　
   * @param [in] val      出力データポインタ
   * @param [in] minmax   フィールデータのMinMax
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力 true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */ 
  template<class T, class TimeT, class TimeAvrT>
  CDM::E_CDM_ERRORCODE
  WriteData(const unsigned step, 
            TimeT time,
            const int sz[3], 
            const int nVari,
            const int gc, 
            T* val, 
            T* minmax=NULL, 
            bool avr_mode=true, 
            unsigned step_avr=0, 
            TimeAvrT time_avr=0.0);

  /**
   * @brief write field data record (template function)
   * @details フィールドデータのみ出力。dfiファイルの出力はなし。
   * @param [in] step     出力ステップ番号
   * @param [in] time     出力時刻　　　　
   * @param [in] sz       valの実ボクセルサイズ
   * @param [in] nVari    valの変数の個数
   * @param [in] gc       valの仮想セル数　　　
   * @param [in] val      出力データポインタ
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力 true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */ 
  template<class T, class TimeT, class TimeAvrT>
  CDM::E_CDM_ERRORCODE
  WriteFieldDataFile(const unsigned step,
                     TimeT time,
                     const int sz[3], 
                     const int nVari,
                     const int gc, 
                     T* val, 
                     bool avr_mode=true, 
                     unsigned step_avr=0, 
                     TimeAvrT time_avr=0.0);

  /**
   * @brief write field data record
   * @details template WriteData関数で方に応じた配列を確保した後、呼び出される 
   * @param [in] step     出力ステップ番号
   * @param [in] gc       仮想セル数　　　
   * @param [in] time     出力時刻　　　　
   * @param [in] val      出力データポインタ
   * @param [in] minmax   フィールデータのMinMax
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力
   *                                              true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */ 
  CDM::E_CDM_ERRORCODE
  WriteData(const unsigned step, 
            const int gc, 
            double time, 
            cdm_Array* val, 
            double* minmax, 
            const bool avr_mode, 
            const unsigned step_avr, 
            double time_avr);

  /**
   * @brief write field data record (not output dfi file)
   * @details template WriteFieldDataFile関数で方に応じた配列を確保した後、呼び出される 
   * @param [in] step     出力ステップ番号
   * @param [in] gc       仮想セル数　　　
   * @param [in] time     出力時刻　　　　
   * @param [in] val      出力データポインタ
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力
   *                                              true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */ 
  CDM::E_CDM_ERRORCODE
  WriteFieldDataFile(const unsigned step, 
                     const int gc, 
                     double time, 
                     cdm_Array* val, 
                     const bool avr_mode, 
                     const unsigned step_avr, 
                     double time_avr);

  /**
   * @brief proc DFIファイル出力コントロール (float)
   * @param [in] comm      MPIコミュニケータ
   * @param [in] out_host  ホスト名出力フラグ　　　　
   * @return true:出力成功 false:出力失敗
   */
/*
  CDM::E_CDM_ERRORCODE
  WriteProcDfiFile(const MPI_Comm comm, 
                   bool out_host=false);
                   float* org=NULL);
*/
  /**
   * @brief proc DFIファイル出力コントロール
   * @param [in] comm          MPIコミュニケータ
   * @param [in] out_host      ホスト名出力フラグ　
   * @param [in] cell_id       cell id
   * @param [in] bcf_id        境界ID
   * @return 終了コード 1:正常 1以外:エラー
   */
  CDM::E_CDM_ERRORCODE
  WriteProcDfiFile(const MPI_Comm comm, 
                   const bool out_host,
                   const int cell_id,
                   const int bcf_id);
                   //double* org=NULL);

  /**
   * @brief index DFIファイル出力 (API関数)
   */
  CDM::E_CDM_ERRORCODE
  WriteIndexDfiFile();

  /**
   * @brief grid ファイル出力コントロール
   * @param [in] iblank      iblankデータポインタ(PLOT3Dのxyzファイル用)　　　　
   */
  CDM::E_CDM_ERRORCODE
  WriteGridFile(const int* iblank=NULL);

  /**
   * @brief 配列形状を文字列で返す
   * @return 配列形状（文字列)
   */
  std::string 
  GetArrayShapeString();

  /**
   * @brief 配列形状を返す
   * @return 配列形状（e_num番号)
   */
  CDM::E_CDM_ARRAYSHAPE 
  GetArrayShape();

  /**
   * @brief get DataType （データタイプの取り出し関数）
   * @return データタイプ（文字列)
   */
  std::string 
  GetDataTypeString();

  /**
   * @brief get DataType （データタイプの取り出し関数）
   * @return データタイプ(e_num番号)
   */
  CDM::E_CDM_DTYPE 
  GetDataType();

  /**
   * @brief get DFIType （dfi種別の取り出し関数）
   * @return dfi種別（文字列)
   */
  std::string 
  GetDFITypeString();

  /**
   * @brief get DFIType （dfi種別の取り出し関数）
   * @return dfi種別(e_num番号)
   */
  CDM::E_CDM_DFITYPE 
  GetDFIType();

  /** 
   * @brief get FileFormat （FileFormatの取り出し関数） 
   * @return FileFormat(文字列)
   */
  std::string
  GetFileFormatString();

  /** @brief get FileFormat (FileFormatの取り出し関数）
   *  @return FileFormat(e_num番号)
   */
  CDM::E_CDM_FORMAT
  GetFileFormat();

  /**
   * @brief get Number of Variables （変数の個数の取り出し関数）
   * @return 変数の個数
   */
  int 
  GetNumVariables();

  /*
   * @brief get Number of GuideCell (仮想セル数の取り出し関数)
   * @return 仮想セル数
   */
  int
  GetNumGuideCell();

  /**
   * @brief データタイプを文字列からe_num番号に変換 
   * @param [in] datatype dfiから取得したデータタイプ
   * @return データタイプ(E_CDM_DTYPE)
   */
  static CDM::E_CDM_DTYPE 
  ConvDatatypeS2E(const std::string datatype); 

  /**
   * @brief データタイプをe_num番号から文字列に変換 
   * @param [in] Dtype データタイプ
   * @return データタイプ(string)
   */
  static std::string 
  ConvDatatypeE2S(const CDM::E_CDM_DTYPE Dtype); 

  /**
   * @brief DFI DomainのGlobalVoxelの取り出し
   * @return GlobalVoxelのポインタ
   */
  const int* 
  GetDFIGlobalVoxel(); 

  /**
   * @brief DFI DomainのGlobalDivisionの取り出し
   * @return GlobalDivisionのポインタ
   */
  const int* 
  GetDFIGlobalDivision(); 

  /**
   * @brief Uuitをセットする
   * @param [in] Name       追加する単位系("Length","Velocity",,,,)
   * @param [in] Unit       単位ラベル("M","CM","MM","M/S",,,)
   * @param [in] reference  規格化したスケール値
   * @param [in] difference 差の値
   * @param [in] BsetDiff   differenceの有無
   */
  void 
  AddUnit(const std::string Name,
          const std::string Unit,
          const double reference,
          const double difference= 0.0,
          const bool BsetDiff=false);

  /**
   * @brief TimeSliceをセットする
   * @param [in] step     出力ステップ番号
   * @param [in] time     出力時刻　　　
   * @param [in] minmax   フィールデータのMinMax
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力 true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */
  template<class T, class TimeT, class TimeAvrT>
  void 
  AddTimeSlice(const unsigned step,
               TimeT time,
               T* minmax=NULL, 
               bool avr_mode=true, 
               unsigned step_avr=0, 
               TimeAvrT time_avr=0.0);

  /**
   * @brief UuitElemを取得する
   * @param[in]  Name 取得する単位系
   * @param[out] unit 取得したcdm_UnitElem
   * @return error code
   */
  CDM::E_CDM_ERRORCODE GetUnitElem(const std::string Name,
                                   cdm_UnitElem &unit); 

  /**
   * @brief UnitElemのメンバ変数毎に取得する
   * @param[in]  Name 取得する単位系
   * @param[out] unit 単位文字列
   * @param[out] ref  reference
   * @param[out] diff difference
   * @param[out] bSetDiff differenceの有無（true:あり false:なし）
   * @return error code
   */
  CDM::E_CDM_ERRORCODE GetUnit(const std::string Name,
                               std::string &unit,
                               double &ref,
                               double &diff,
                               bool &bSetDiff); 
 
  /**
   * @brief TimeSlice OnOff フラグをセットする
   * @param [in] ONOFF
   */
  void 
  SetTimeSliceFlag(const CDM::E_CDM_ONOFF ONOFF); 

  /**
   * @brief FileInfoの変数名を登録する
   * @param [in] pvari    変数位置 0:u, 1:v, 2:w
   * @param [in] variName 変数名 "u","v","w",,,
   */
  void setVariableName(int pvari, std::string variName); 

  /**
   * @brief FileInfoの変数名を取得する
   * @param [in] pvari 変数位置 0:u, 1:v, 2:w
   * @return 変数名
   */
  std::string getVariableName(int pvari);

  /**
   * @brief DFIに出力されているminmaxの合成値を取得
   * @param [in]  step 取得するステップ
   * @param [out] vec_min 取得したminmaxの合成値
   * @param [out] vec_max 取得したminmaxの合成値
   * @return error code 取得出来たときは E_CDM_SUCCESS
   */
  CDM::E_CDM_ERRORCODE getVectorMinMax(const unsigned step,
                                       double &vec_min,
                                       double &vec_max);

  /**
   *brief DFIに出力されているminmaxを取得
   * @param [in]  step 取得するステップ
   * @param [in]  variNo 変数No(0～n)
   * @param [out] min_value 取得したmin
   * @param [out] max_value 取得したmax
   * @return error code 取得出来たときは E_CDM_SUCCESS
   */
  CDM::E_CDM_ERRORCODE getMinMax(const unsigned step,
                                 const int variNo,
                                 double &min_value,
                                 double &max_value);

  /**
   * @brief 読込みランクリストの作成
   * @details RankListがあるかないか判定しないときは新規にRankListを生成し
   *          それをもとにランクマップの生成、読込みランクリストreadRankList
   *          を生成する
   * @param [in]  dfi_domain DFIのdomain情報
   * @param [in]  head       ソルバーのHeadIndex
   * @param [in]  tail       ソルバーのTailIndex
   * @param [in]  readflag   読込み方法
   * @param [out] readRankList 読込みランクリスト
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  CheckReadRank(const cdm_Domain* dfi_domain,
                const int head[3],
                const int tail[3],
                CDM::E_CDM_READTYPE readflag,
                vector<int> &readRankList);

  /**
   * @brief 出力インターバルステップの登録
   * @details 登録しない（本メソッドがコールされない）場合はCDMでのインターバル
   *          制御は行わない
   * @param [in] interval_step インターバルステップ
   * @param [in] base_step     基準となるステップ（デフォルト0ステップ）
   * @param [in] start_step    セッション開始ステップ（デフォルト0ステップ）
   * @param [in] last_step     セッション最終ステップ（デフォルト、-1：最終ステップで出力しない）
   */
  void setIntervalStep(int interval_step,
                       int base_step =0, 
                       int start_step=0,
                       int last_step =-1); 

  /**
   * @brief インターバルタイムの登録
   * @param [in] interval_time 出力インターバルタイム
   * @param [in] dt            計算の時間間隔
   * @param [in] base_time     基準となるタイム（デフォルト0.0タイム）
   * @param [in] start_time    セッション開始タイム（デフォルト0.0タイム）
   * @param [in] last_time     せっしょん最終タイム（デフォルト、-1.0：最終タイムで出力しない）
   */
  void setIntervalTime(double interval_time,
                       double dt,
                       double base_time =0.0, 
                       double start_time=0.0,
                       double last_time =-1.0); 

  /**
   * @brief インターバルの計算に使われる全ての時間をスケールで無次元化する
   * @details (base_time, interval_time, start_time, last_time)
   * @param [in] scale スケール
   * return modeがStepのときはfalseを返す、無次元化しない
   */
  bool normalizeTime(const double scale);

  /**
   * @brief インターバルのbase_timeをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeBaseTime(const double scale);
 
  /**
   * @brief インターバルのintervalをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeIntervalTime(const double scale);
 
  /**
   * @brief インターバルのstart_timeをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeStartTime(const double scale);

  /**
   * @brief インターバルのlast_timeをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeLastTime(const double scale);

  /**
   * @brief インターバルのDetlaTをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeDelteT(const double scale);



  /**
   * @brief read field data record
   * @param [in]  fname    FieldData ファイル名
   * @param [in]  step     読込みステップ番号
   * @param [out] time     読み込んだ時間
   * @param [in]  sta      読込みスタート位置
   * @param [in]  end      読込みエンド位置
   * @param [in]  DFI_head dfiのHeadIndex
   * @param [in]  DFI_tail dfiのTailIndex
   * @param [in]  avr_mode 平均ステップ＆時間読込みフラグ　false : 読込み
   * @details　                                             true  : 読み込まない
   * @param [out] avr_step 平均ステップ
   * @param [out] avr_time 平均時間
   * @param [out] ret      終了コード
   * @return 読み込んだ配列のポインタ
   */
  virtual
  cdm_Array* 
  ReadFieldData(std::string fname,
                const unsigned step,
                double &time,
                const int sta[3],
                const int end[3],
                const int DFI_head[3],
                const int DFI_tail[3],
                bool avr_mode,
                unsigned &avr_step,
                double &avr_time,
                CDM::E_CDM_ERRORCODE &ret );

  /**
   * @brief フィールドデータファイルのヘッダーレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        ステップ番号
   * @param[in]  head        dfiのHeadIndex
   * @param[in]  tail        dfiのTailIndex
   * @param[in]  gc          dfiのガイドセル数
   * @param[out] voxsize     voxsize
   * @param[out] time        時刻
   * @return true:出力成功 false:出力失敗
   */
  virtual CDM::E_CDM_ERRORCODE 
  read_HeaderRecord(FILE* fp,
                    bool matchEndian,
                    unsigned step,
                    const int head[3],
                    const int tail[3],
                    int gc,
                    int voxsize[3],
                    double &time)=0;

  /**
   * @brief フィールドデータファイルのデータレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  buf         読込み用バッファ
   * @param[in]  head        読込みバッファHeadIndex
   * @param[in]  nz          z方向のボクセルサイズ（実セル＋ガイドセル＊２）
   * @param[out] src         読み込んだデータを格納した配列のポインタ
   */
  virtual CDM::E_CDM_ERRORCODE 
  read_Datarecord(FILE* fp,
                  bool matchEndian,
                  cdm_Array* buf,
                  int head[3],
                  int nz,
                  cdm_Array* &src)=0;

  /**
   * @brief sphファイルのAverageデータレコードの読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        読込みstep番号
   * @param[out] avr_step    平均ステップ
   * @param[out] avr_time    平均タイム
   */
  virtual CDM::E_CDM_ERRORCODE
  read_averaged(FILE* fp,
                bool matchEndian,
                unsigned step,
                unsigned &avr_step,
                double &avr_time)=0;   

protected :
  /**
   * @brief write field data record (double)
   * @param [in] fname    出力フィールドファイル名
   * @param [in] step     出力ステップ番号
   * @param [in] time     出力時刻　　　　
   * @param [in] val      出力データポインタ
   * @param [in] mode     平均ステップ＆時間出力　false : 出力 true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   * @return error code
   */
  virtual
  CDM::E_CDM_ERRORCODE 
  WriteFieldData(std::string fname,
                 const unsigned step, 
                 double time, 
                 cdm_Array* val, 
                 const bool mode, 
                 const unsigned step_avr, 
                 const double time_avr);

  /**
   * @brief SPHヘッダファイルの出力
   * @param[in] fp     ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  virtual CDM::E_CDM_ERRORCODE
  write_HeaderRecord(FILE* fp,
                     const unsigned step,
                     const double time,
                     const int RankID)=0;

  /**
   * @brief SPHデータレコードの出力
   * @param[in]  fp ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  gc ガイドセル
   * @param[in]  RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  virtual CDM::E_CDM_ERRORCODE
  write_DataRecord(FILE* fp,
                   cdm_Array* val,
                   const int gc,
                   const int RankID)=0;

  /**
   * @brief Averageレコードの出力
   * @param[in] fp       ファイルポインタ
   * @param[in] step_avr 平均ステップ番号
   * @param[in] time_avr 平均時刻
   * @return true:出力成功 false:出力失敗
   */
  virtual CDM::E_CDM_ERRORCODE
  write_averaged(FILE* fp,
                 const unsigned step_avr,
                 const double time_avr)=0;

  /**
   * @brief Grid data file 出力 コントロール
   * @param [in] iblank  iblankデータポインタ(PLOT3Dのxyzファイル用)
   */
  virtual bool
  write_GridData(const int* iblank) {}; 

//FEAST 20131125.s
  /**
   * @brief ascii ヘッダーレコード出力(bov,avs)
   * @param [in] step step番号
   * @param [in] time time
   */
  virtual
  bool 
  write_ascii_header(const unsigned step,
                     const double time)
  { return true; }; 

  /**
   * @brief データタイプ毎のサイズを取得
   * @param [in] Dtype データタイプ(Int8,Int16,,,,etc)
   * @return データサイズ
   * @return 0 エラー
   */
  static int 
  get_cdm_Datasize(CDM::E_CDM_DTYPE Dtype); 

  /**
   * @brief Create Process 
   * @param [in] comm           MPIコミュニケータ
   * @param [out] G_Process     Process class　　　
   */
  void 
  cdm_Create_dfiProcessInfo(const MPI_Comm comm,
                            cdm_Process &G_Process);


  /**
   * @brief 読込み判定判定
   * @param [in] G_voxel            計算空間全体のボクセルサイズ（自）
   * @param [in] DFI_GlobalVoxel    計算空間全体のボクセルサイズ（DFI）
   * @param [in] G_Div              分割数（自）
   * @param [in] DFI_GlobalDivision 分割数（DFI）
   * @return 読込みタイプコード
   */
  //cdm_EGlobalVoxel CheckGlobalVoxel(const int Gvoxel[3], 
  CDM::E_CDM_READTYPE 
  CheckReadType(const int G_voxel[3], 
                const int DFI_GlobalVoxel[3],
                const int G_Div[3],
                const int DFI_GlobalDivision[3]); 

  /**
   * @brief フィールドデータの読込み範囲を求める
   * @param [in] isSame   粗密フラグ true:密、false:粗
   * @param [in] head     計算領域の開始位置(自)　
   * @param [in] tail     計算領域の終了位置(自)　
   * @param [in] gc       仮想セル数(自)　
   * @param [in] DFI_head 計算領域の開始位置(DFI)　　
   * @param [in] DFI_tail 計算領域の終了位置(DFI)　　
   * @param [in] DFI_gc   仮想セル数(DFI)　
   * @param [in] readflag 読込み方法
   * @param [out] copy_sta コピー開始位置
   * @param [out] copy_end コピー終了位置　　
   * @param [out] read_sta 読込み開始位置
   * @param [out] read_end 読込み終了位置　　
   */
  void 
  CreateReadStartEnd(bool isSame,
                     const int head[3], 
                     const int tail[3], 
                     const int gc, 
                     const int DFI_head[3], 
                     const int DFI_tail[3],
                     const int DFI_gc, 
                     const CDM::E_CDM_READTYPE readflag, 
                     int copy_sta[3], 
                     int copy_end[3],
                     int read_sta[3],
                     int read_end[3]);

  /**
   * @brief index DFIファイル出力
   * @param [in] dfi_name  DFIファイル名
   * @return true:出力成功 false:出力失敗
   */
  CDM::E_CDM_ERRORCODE
  WriteIndexDfiFile(const std::string dfi_name);


public:

  /**
   * @brief セル中心データを格子点に値をセット
   * @param [out]  P 格子点データ
   * @param [in]   S セル中心data
   */
  template<class T1, class T2>
  bool setGridData(
                   cdm_TypeArray<T1>* P,
                   cdm_TypeArray<T2>* S);

  /**
   * @brief 内部の格子点のデータを重み付けでで割る
   * @param[out] P  格子点data
   */
  template<class T>
  void VolumeDataDivide(cdm_TypeArray<T> *P);


public:

  /**
   * @brief ディレクトリパスの作成(MakeDirectorySubを呼出して作成)
   * @param [in] path パス
   * @return error code　　　　　　　
   */ 
  int MakeDirectory(const std::string path);

  /**
   * @brief ディレクトリパスの作成(MakeDirectory関数を呼出して作成)
   * @return error code　　　　　　　
   */ 
  int MakeDirectoryPath();

  /**
   * @brief ディレクトリパスの作成(system関数mkdirで作成)
   * @param [in] path パス
   * @return error code　　　　　　　
   */ 
  static int MakeDirectorySub( std::string path );

  /**
   * @brief dfiのパスとDirectoryPathを連結する関数
   * @return パス名
   */
  std::string Generate_Directory_Path(); 

  /** バージョンを出力する
   */
  static std::string getVersionInfo()
  {
    std::string str(CDM_VERSION_NO);
    return str;
  }

  /**
   * @brief  コンパイルオプションで与えたバッファサイズを取得
   * @return バッファサイズ
   */
  int getBufSize(); 

};

//inline 関数
#include "inline/cdm_DFI_inline.h"


#endif // _cdm_DFI_H_
