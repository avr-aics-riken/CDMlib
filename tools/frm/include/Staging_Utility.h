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
 * @file Staging_Utility.h
 * @brief Staging Class Header
 * @author aics
 * @data 2013/4/24
 */

#ifndef _STAGING_UTILITY_H
#define _STAGING_UTILITY_H

#include "cdm_DFI.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <set>
#include <map>

#include "TextParser.h"

#include "ActiveSubDomain.h"
#include "Staging_Define.h"
#include "frm_EndianUtil.h"

/** Staging main class */
class Staging {

public:

  /** FCONV用処理リスト */
  struct step_rank_info{
   cdm_DFI* dfi;
   int stepStart;
   int stepEnd;
   int rankStart;
   int rankEnd;
  };

  /** staging 用 info.txt ファイルの Process */
  struct Rank
  {
    int RankID;
    int VoxelSize[3];
    int HeadIndex[3];
    int TailIndex[3];
    int RankPosition[3];
  };

  /** コンストラクタ **/
  Staging(){
//20160407.fub.s
    m_CoordinateStep=-1;
//20160407.fub.e
    m_infofile="";
    m_fconvfile="";
    m_fconv_numproc=1;
    m_step=-1;
    m_NumberOfRank=0;
    m_inPath="";
    m_outPath="";
    m_cropIndexStart_on=false;
    m_cropIndexEnd_on=false;
    for(int i=0; i<3; i++) {
      m_GVoxel[i]=0;
      m_Gdiv[i]=0;
    }
    m_dfi_fname.clear();
  };

  /** デストラクタ **/
  ~Staging(){};

public:

  //cdm_DFI *DFI;
  vector<cdm_DFI *> DFI;

//20160408.fub.s
//const cdm_FileInfo  *dfi_Finfo;   ///< DFI cdm_FileInfoのポインタ
        cdm_FileInfo  *dfi_Finfo;   ///< DFI cdm_FileInfoのポインタ
//20160408.fub.e
  const cdm_FilePath  *dfi_Fpath;   ///< DFI cdm_FilePathのポインタ
  const cdm_VisIt     *dfi_Visit;   ///< DFI cdm_VisItのポインタ
  const cdm_Unit      *dfi_Unit;    ///< DFI cdm_Unitのポインタ
  const cdm_Domain    *dfi_Domain;  ///< DFI cdm_Domainのポインタ
  const cdm_MPI       *dfi_MPI;     ///< DFI cdm_MPIのポインタ
  const cdm_TimeSlice *dfi_TSlice;  ///< DFI cdm_TimeSliceのポインタ
  //const cdm_Process   *dfi_Process; ///< DFI cdm_Processのポインタ
        cdm_Process   *dfi_Process; ///< DFI cdm_Processのポインタ

  int m_step;                   ///< Restart step番号
//20160407.fub.s
  int m_CoordinateStep;         ///< コピーする座標値ファイルのステップ番号
//20160407.fub.e
  string m_inPath;              ///< DFIディレクトリ
  string m_outPath;             ///< 出力ディレクトリ
  string m_infofile;            ///< Staging用procファイル名
  int m_GVoxel[3];              ///< Staging用ボクセルサイズ
  int m_Gdiv[3];                ///< Staging用分割数
  string m_ActiveSubdomain;     ///< ActiveSubdomainファイル名
  int m_NumberOfRank;           ///< Staging用並列数
  vector<Rank> m_GRankInfo;     ///< Staging用並列情報テーブル
  vector<string> m_dfi_fname;   ///<STaging用DFIファイル名
//FCONV 20140116.s
  string m_fconvfile;           ///< FCONV 入力ファイル名
  bool   m_fconv_inputfile;     ///< FCONV 入力ファイルフラグ
  int    m_fconv_numproc;       ///< FCONV 実行並列数
  stg_E_OUTPUT_CONV m_ConvType; ///< FCONV コンバートタイプ
  int    m_outGc;               ///< FCONV 出力ガイドセル数
  bool   m_cropIndexStart_on;   ///< FCONV 入力領域のスタート指示フラグ
  bool   m_cropIndexEnd_on;     ///< FCONV 入力領域のエンド位置
  int    m_CropStart[3];        ///< FCONV 入力領域スタート位置
  int    m_CropEnd[3];          ///< FCONV 入力領域エンド位置
  int    m_outList;             ///< FCONV ファイル割振り方法
  std::vector<step_rank_info> m_StepRankList; ///< ファイル割振りリスト
//FCONV 20140116.e

  std::vector<ActiveSubDomain> m_subDomainInfo; ///< 活性サブドメイン情報

  int *m_rankMap;
  int (*m_HeadTail)[6];
  int (*m_StartEnd)[6];

  std::set<int>m_headX;
  std::set<int>m_headY;
  std::set<int>m_headZ;

  headT m_mapX,m_mapY,m_mapZ;

public:

//20160407.fub.s
  /** 引数にエラーチェック
   * @ retval 終了コード(STG_SUCCESS=正常終了)
   */
  stg_ErrorCode ArgErrorCheck();

  /** fubフィルコピー（フィールドデータファイルと座標値データファイルコピー
   *
   */
  void FubFileCopy(int readRankID, int step, const bool mio,
                   std::string path2, cdm_DFI *dfi, char *cmd);
//20160407.fub.e
  /** 初期化、ファイルの読込み
   * @param[in] infofname ステージング用procファイル名
   * @retval    true  ファイル読込み成功
   * @retval    false ファイル読込み失敗
   */
  //bool Initial(string infofname, string dfifname);
  bool Initial(string infofname);

  /** ステージング用procファイルの読込み
   * @retval     true     ファイル読込み成功
   * @retval     false    ファイル読込み失敗
   */
  //bool ReadInfo(string infofile, string &fconvfile);
  bool ReadInfo();

//FCONV 20140116.s
  /** FCONV 入力ファイルの読込み
   * @retval    true      FCONV入力成功
   * @retval    false     FCONV入力失敗
   */
  //bool ReadFconvInputFile(string fconvfile);
  bool ReadFconvInputFile();

  /** FCONV 入力ファイルのチェック
   * @retval true 正常
   * @retval false エラー
   */
  bool FconvInputCheck();

  /** step基準のリスト作成
   * @param[in] myID 処理ランク
   */
  void makeStepList(int myID);

  /** rank基準のリスト生成
   * @param[in] myID 処理ランク
   */
  void makeRankList(int myID);

//FCONV 20140116.e

  /** 領域分割数の取得
   * @return 領域分割数整数配列のポインタ
   */
  const int* GetDivNum() const;

  /** 活性サブドメイン情報の存在チェック
   * @param[in] subDomain チェックする活性サブドメイン情報
   * @retval    true      存在する
   * @retval    false     存在しない
   */
  bool IsExistSubdomain( ActiveSubDomain subDomain );

  /** 活性サブドメイン情報の追加
   * @param[in] subDomain 追加する活性サブドメイン情報
   * @retval    true      追加した
   * @retval    false     追加に失敗(同じ領域分割位置で追加済み)
   */
  bool AddSubdomain( ActiveSubDomain subDomain );

  /** 活性サブドメインの数を取得
   * 活性サブドメインの数＝活性サブドメイン情報配列のサイズだが、
   * この配列が空のとき、領域分割数でサブドメイン数を決定して返す
   * @return 活性サブドメインの数
   */
  int GetSubdomainNum() const;

  /**  活性サブドメイン情報を取得
   * @param[in]  idx 登録順番号
   * @return 活性サブドメイン情報ポインタ
   */
  const ActiveSubDomain* GetSubdomainInfo( size_t idx ) const;

  /** ランクマップを生成
   * @retval ランクマップ
   */
  //bool CreateRankMap();
  int* CreateRankMap();

  /** 有効なランクマップを生成
   * @retval 有効なランクマップ
   *
   */
  int* CreateActiveRankMap();

  /** ActiveSubdomainファイルのエンディアンチェック
   *  ActiveSubdomainファイルのエンディアンをチェック
   * @param[in]  ident               ActiveSubdomainファイルのIdentifier
   * @retval     STG_Match   一致
   * @retval     STG_UnMatch 不一致
   * @retval     STG_UnKnown フォーマットが異なる
   *
   */
  static stg_EMatchType isMatchEndianSbdmMagick( int ident );

  /** ActiveSubdomainファイルの読み込み
   * ActiveSubdomainファイルを読み込み、活性ドメイン情報を生成する
   * @param[in]  subDomainFile ActiveSubdomainファイル名
   * @return   終了コード(STG_SUCCESS=正常終了)
   */
   stg_ErrorCode ReadActiveSubdomainFile( std::string subDomainFile );

  /** ActiveSubdomainファイルの読み込み(static関数)
   * ActiveSubdomainファイルを読み込み、活性ドメイン情報を生成する
   * @param[in]  subDomainFile ActiveSubdomainファイル名
   * @param[out] subDomainInfo 活性ドメイン情報
   * @param[out] div           ActiveSubdiomainファイル中の領域分割数
   * @return   終了コード(STG_SUCCESS=正常終了)
   */
   static stg_ErrorCode ReadActiveSubdomainFile( std::string subDomainFile,
                                                 std::vector<ActiveSubDomain>& subDomainInfo,
                                                 int div[3] );

  /** VOXEL数の取得
   * @return VOXEL数情報実数配列のポインタ
   */
   const int* GetVoxNum() const;

  /** 領域情報のチェック
   * 活性サブドメイン配列が空のとき、全領域が活性サブドメインになるため
   * このチェック関数内で活性サブドメイン情報を生成する.
   * @param[in] nRank 並列プロセス数
   * @return   終了コード(STG_SUCCESS=正常終了)
   */
   stg_ErrorCode CheckData( int nRank );

  /** head&tail情報の作成
   * @param[in] rankMap ランクマップ
   * @retval true  正常終了
   * @retval false エラー
   */
   bool CreateHeadTail(int* rankMap);

  /** head&tail情報の作成
   * @param[in] rankMap ランクマップ
   * @param[in] Rankinfo ランク情報
   * @retval true  正常終了
   * @retval false エラー
   */
   bool CreateHeadTail(int* rankMap, vector<Rank>& RankInfo);

  /** head&taileのセット
   */
   void SetHeadTail();

  /** head&tailをガイドセルで更新
   * @param[in]  mst_head 更新元のヘッドインデックス
   * @param[in]  mst_tail 更新元のテイルインデックス
   * @param[out] head 更新したヘッドインデックス
   * @param[out] tail 更新したテイルインデックス
   */
  void UpdateHeadTail(const int* mst_head, const int* mst_tail, int* head, int* tail);

  /** Headマップの生成
   * @param[in] head head情報
   * @param[in] map headマップ
   */
   void CreateHeadMap(std::set<int>head, headT &map);

  /** 粗密データ判定
   * @param [in] Gvoxel     計算空間全体のボクセルサイズ（自）
   * @param [in] DFI_Gvoxel 計算空間全体のボクセルサイズ（DFI）
   * @return STG_E_GV_SAME:密 STG_E_GVX2_SAME:粗 STG_E_OTHER:その他
   */
  stg_EGlobalVoxel CheckGlobalVoxel(int Gvoxel[3], int DFI_Gvoxel[3]);


  /** 粗密用にhead,tailを作成
   * @param [out] HeadTail 粗密用heat,tail情報テーブル
   * @param [in]  numrank  dfiのランク数
   */
  void MakeGVX2HeadTeail(int (*HeadTail)[6], int numrank);

  /** 検索範囲のしぼりこみ
   * @param [in] HeadTail heat,tail情報テーブル
   * @retval true  正常終了
   * @retval false エラー
   */
  bool CreateStartEnd(int (*HeadTail)[6]);

  /** ファイルをステージング用のディレクトリにコピー
   * @param [in] readRankList 読込みランクリスト
   * @param [in] myRank 処理するランク番号
   * @param [in] 処理するindex.dfiのポインタ(fub用追加)
   * @retval true  正常終了
   * @retval false エラー
   */
// 20160407.fub.s
//bool FileCopy(vector<int> readRankList, int myRank);
  bool FileCopy(vector<int> readRankList, int myRank, cdm_DFI *dfi);
// 20160407.fub.e

  /** ファイルをステージング用のディレクトリにコピー(FCONV用)
   * @param [in] info 処理リスト
   * @param [in] myRank 処理するランク番号
   * @param [in] ndfi 処理するDFI通番
   * @retval true  正常終了
   * @retval false エラー
   */
//20160509.fub.s
//bool FileCopy(step_rank_info info, int myRank);
  bool FileCopy(step_rank_info info, int myRanki, int ndfi);
//20160509.fub.e

  /** dfi ファイル出力
   * @param [in] fname DFIファイル名
   * @param [in] rankMap ランクマップ
   * @param [in] dfi 処理しているdfiポインター
   * @retval true  正常終了
   * @retval false エラー
   */
//20160408.fub.s
//bool OutputDFI(string fname, int* rankMap);
  bool OutputDFI(string fname, int* rankMap, cdm_DFI* dfi);
//20160408.fub.e

  /** index.dfiファイルの生成、出力
   * @param [in] dfi_name index dfi ファイル名
   * @param [in] dfi 処理しているdfiポインター
   */
//20160408.fub.s
//CDM::E_CDM_ERRORCODE
//WriteIndexDfiFile(const std::string dfi_name);
  CDM::E_CDM_ERRORCODE
  WriteIndexDfiFile(const std::string dfi_name, cdm_DFI *dfi);
//20160408.fub.e

  /** index.dfiファイルの生成、出力（FCONV用）
   * @param [in] dfi_name index dfi ファイル名
   * @param [in] info 処理リスト
   */
  CDM::E_CDM_ERRORCODE
  WriteIndexDfiFile(const std::string dfi_name, const step_rank_info info);

  /** ファイル名の作成
   * @param [in] RankID ランク番号
   * @param [in] step   読込みステップ番号
   * @param [in] mio    並列判定フラグ（逐次or並列の判定用）
   * @return 生成されたファイル名
   */
   std::string Generate_FileName(int RankID, int step, const bool mio);

  /** ディレクトリパスの作成
   * @param [in] path パス
   * @return error code
   */
   int MakeDirectory(string path);

  /**
   * @brief ディレクトリパスの作成(MakeDirectory関数を呼出して作成)
   * @return error code　　　　　　　
   */
   int MakeDirectoryPath();

  /**
   * @brief ディレクトリパスの作成(system関数mkdirで作成)
   * @param
   * @return error code　　　　　　　
   */
  static int MakeDirectorySub( std::string path );

  /**
   * @brief dfiのパスとDirectoryPathを連結する関数
   * @return パス名
   */
  std::string Generate_Directory_Path();


};

#endif //_STAGING_UTILITY_H
