#ifndef _CDM_PROCESS_H_
#define _CDM_PROCESS_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_Process.h
 * @brief  cdm_RANK & cdm_Process Class Header
 * @author aics    
 */

/** proc.dfi ファイルの Rank */
class cdm_Rank {

public:


  int RankID;                           ///<ランク番号
  std::string HostName;                 ///<ホスト名 
  int VoxelSize[3];                     ///<ボクセルサイズ
  int HeadIndex[3];                     ///<始点インデックス
  int TailIndex[3];                     ///<終点インデックス
  int c_id;                             ///<cell id
  int bc_id;                            ///<境界ID

  /** コンストラクタ **/
  cdm_Rank();

  /** デストラクタ **/
  ~cdm_Rank();

  /**
   * @brief read Rank(proc.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス
   * @param [in]   label_leaf ベースとなる名前（"/Process/Rank") 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl,
       std::string label_leaf);

  /**
   * @brief DFIファイル:Rank出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @return true:出力成功 false:出力失敗
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

};

/** proc.dfi ファイルのProcess */
class cdm_Process {

public:

  typedef std::map<int,int> headT;

  vector<cdm_Rank> RankList;

  int *m_rankMap;

  /** コンストラクタ */
  cdm_Process();

  /** デストラクタ */
  ~cdm_Process();

  /**
   * @brief read Rank(proc.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl); 

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
   * @brief DFIのProcessにHeadIndex,TailIndex指定が無い場合
   * @details ActiveSubDomainがあれば、読込み、なければ全て有効で
   * subDomainを生成し、CreateRankListに渡す
   * CPMと同じ分割でhead&tail情報を作成してRankListを作成する
   * @param [in]  dfi_domain DFIのdomain情報
   * @param [out] mapHeadX   headXをキーにした位置情報マップ 
   * @param [out] mapHeadY   headXをキーにした位置情報マップ 
   * @param [out] mapHeadZ   headXをキーにした位置情報マップ 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE 
  CreateRankList(const cdm_Domain* dfi_domain,
                 map<int,int> &mapHeadX,
                 map<int,int> &mapHeadY,
                 map<int,int> &mapHeadZ);

  /**
   * @brief DFIのProcessにHeadIndex,TailIndex指定が無い場合
   * 渡された、subDomainをもとにCPM同様の分割方法でRankList
   * を生成する
   * @param [in]  div      分割数
   * @param [in]  gvox     ボクセルサイズ
   * @param [out] mapHeadX   headXをキーにした位置情報マップ 
   * @param [out] mapHeadY   headXをキーにした位置情報マップ 
   * @param [out] mapHeadZ   headXをキーにした位置情報マップ 
   * @return error code 
   */
  CDM::E_CDM_ERRORCODE 
  CreateRankList(const int div[3],
                 const int gvox[3],
                 map<int,int> &mapHeadX,
                 map<int,int> &mapHeadY,
                 map<int,int> &mapHeadZ);

  /**
   * @brief ActiveSubDomain情報を作成
   * @param [in]  dfi_domain DFIのdomain情報
   * @param [out] subDomainInfo 活性ドメイン情報
   */
  CDM::E_CDM_ERRORCODE 
  CreateSubDomainInfo(const cdm_Domain* dfi_domain,
                      vector<cdm_ActiveSubDomain> &subDomainInfo); 

  /**
   * @brief ActiveSubdomainファイルのエンディアンをチェック
   * @param [in] ident ActiveSubdomainファイルのIdentifier
   * @return     1   一致
   * @return     0   不一致
   * @return     -1  フォーマットが異なる
   */
  static int isMatchEndianSbdmMagick(int ident);

  /**
   * @brief ActiveSubdomainファイルの読み込み(static関数)
   * @param [in]  subDomainFile ActiveSubdomainファイル名
   * @param [out] subDomainInfo 活性ドメイン情報
   * @param [out] div           ActiveSubdiomainファイル中の領域分割数
   * @return   終了コード(CDM_SUCCESS=正常終了)
   */
  static 
  CDM::E_CDM_ERRORCODE 
  ReadActiveSubdomainFile(
                          std::string subDomainFile,
                          std::vector<cdm_ActiveSubDomain>& subDomainInfo,
                          int div[3]);

  /**
   * @brief subdomain情報からランクマップを生成（非活性を含む）
   * @param [in] div 領域分割数
   * @param [in] subDomainInfo 活性ドメイン情報
   * @return ランクマップ
   * @return NULL
   */
  int* CreateRankMap(const int div[3],
                     std::vector<cdm_ActiveSubDomain> &subDomainInfo);

  /**
   * @brief 生成済のRankListからランクマップを生成
   * @param [in] ndiv  領域分割数
   * @param [in] mapHeadX headXをキーにした位置情報マップ
   * @param [in] mapHeadY headYをキーにした位置情報マップ
   * @param [in] mapHeadZ headZをキーにした位置情報マップ
   * @return ランクマップ
   * @return NULL
   */
  int* CreateRankMap( const int ndiv[3],
                      headT &mapHeadX,
                      headT &mapHeadY,
                      headT &mapHeadZ);
  /**
   * @brief head mapの生成
   * @param [in]  head head インデックス
   * @param [out] map  head map
   */
  void CreateHeadMap(std::set<int>head,
                     headT &map);

  /**
   * @brief head mapの生成
   * @param [in] head head インデックス
   * @param [in] ndiv 分割数
   * @param [out] map head map
   */
  void CreateHeadMap(int* head,
                     int ndiv,
                     headT &map); 


  /**
   * @brief 読込みランクファイルリストの作成
   * @param [in]  dfi_domain   DFIのDomain情報
   * @param [in]  head         計算領域の開始位置　　　　
   * @param [in]  tail         計算領域の終了位置　　　　
   * @param [in]  readflag     粗密データ判定フラグ
   * @param [in]  mapHeadX     headXをキーにした位置情報マップ
   * @param [in]  mapHeadY     headYをキーにした位置情報マップ
   * @param [in]  mapHeadZ     headZをキーにした位置情報マップ
   * @param [out] readRankList 読込みに必要なランク番号リスト
   */
  CDM::E_CDM_ERRORCODE 
  CheckStartEnd(const cdm_Domain* dfi_domain,
                const int head[3],
                const int tail[3],
                CDM::E_CDM_READTYPE readflag,
                headT mapHeadX,
                headT mapHeadY,
                headT mapHeadZ,
                vector<int> &readRankList);


  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp          ファイルポインタ
   * @param [in] tab         インデント
   * @return true:出力成功 false:出力失敗
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);


};

#endif // _CDM_PROCESS_H_
