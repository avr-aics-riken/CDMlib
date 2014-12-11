#ifndef _CONV_H_
#define _CONV_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   conv.h
 * @brief  CONV Class Header
 * @author aics
 * @date   2013/11/7
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <errno.h>
#include <float.h>

#ifndef _WIN32
#include <dirent.h>
#else
#include "sph_win32_util.h"   // for win32
#endif

#include "cpm_ParaManager.h"
#include "cdm_DFI.h"
#include "limits.h" // for UBUNTU

#include "conv_Define.h"

#include "InputParam.h"

using namespace std;

class CONV {

public:

  //処理リスト
  struct step_rank_info{
   cdm_DFI* dfi;
   int stepStart;
   int stepEnd;
   int rankStart;
   int rankEnd;
  };

  //minmax
  struct dfi_MinMax{
   cdm_DFI* dfi;
   double* Min;
   double* Max;
   dfi_MinMax(int nstep, int nvari){
     Min = new double[nvari*nstep];
     Max = new double[nvari*nstep];
     for(int i=0; i<nvari*nstep; i++ ) {
       Min[i]=DBL_MAX;
       Max[i]=-DBL_MAX;
     }
   };
   ~dfi_MinMax(){
     if( Min ) delete [] Min;
     if( Max ) delete [] Max;
   };
 };   

public:
  cpm_ParaManager* m_paraMngr; ///< Cartesian Partition Manager

  InputParam* m_param;     ///< InputParam Class
  
public:
  int m_procGrp;           ///< プロセスグループ番号
  int m_myRank;            ///< 自ノードのランク番号
  int m_numProc;           ///< 全ランク数
  std::string m_HostName;  ///< ホスト名

private:
  
  //並列実行時のSTAGINGのON/OFF
  unsigned m_staging;
  
public:
 
  int m_pflag;
  int m_pflagv;
  int m_lflag;
  int m_lflagv;

  vector<cdm_DFI *>m_in_dfi;  //入力DFIのポインタ

public:
  /** コンストラクタ */
  CONV();
  
  /**　デストラクタ */
  ~CONV();
  
protected:
  
  FILE* m_fplog;
  
public:

  /**
   * @brief conv インスタンス
   * @param [in] InputCntl InputParamクラスポインタ
   */ 
  static CONV*
  ConvInit(InputParam* param);

  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* paraMngr)
  {
    if ( !paraMngr ) return false;
    m_paraMngr = paraMngr;
    setRankInfo();
    return true;
  }
  
  /**
   * @brief ランク情報をセットする
   */
  void setRankInfo()
  {
    m_procGrp = 0;
    m_myRank  = m_paraMngr->GetMyRankID();
    m_numProc = m_paraMngr->GetNumRank();
    m_HostName= m_paraMngr->GetHostName();
  }
 
  /**
   * @brief InputParamのポインタをコピー
   * @param [in] InputCntl InputParamクラスポインタ
   * @return  エラーコード
   */
/*
  bool importInputParam(InputParam* InputCntl)
  {
    if( !InputCntl ) return false;
    m_InputCntl = InputCntl;
    return true;
  }
*/

  /**
   * @brief dfiファイルの読み込みとDfiInfoクラスデータの作成
   */
  CDM::E_CDM_ERRORCODE ReadDfiFiles();
 
  /**
   * @brief dfi毎の変数の個数、出力ガイドセルのチェックと更新
   * @return エラーコード
   */
  bool CheckDFIdata();
 
  /**
   * @brief 出力指定ディレクトリのチェック
   * @param [in] dirstr 出力ディレクトリ
   */
  void CheckDir(string dirstr);
  
  /**
   * @brief ログファイルのオープン
   */
  void OpenLogFile();
  
  /**
   * @brief ログファイルのクローズ
   */
  //void CloseLogFile();
 
  /**
   * @brief dfiのログ出力
   */
  void PrintDFI(FILE* fp); 

 
  /**
   * @brief 所要時間の記述
   * @param [in] tt 所要時間
   */
  void WriteTime(double* tt);
  
  /**
   * @brief メモリ使用量を表示する
   * @param [in] Memory メモリ量
   * @param [in] fp     ファイルポインタ
   */
  void MemoryRequirement(const double Memory, FILE* fp);
  
  /**
   * @brief メモリ使用量を表示する
   * @param [in] TotalMemory トータルメモリ使用量最大値
   * @param [in] sphMemory sphファイル読み込みのためのwkメモリ使用量最大値
   * @param [in] plot3dMemory plot3dファイル書き込みのためのメモリ使用量最大値
   * @param [in] thinMemory 間引きオプションのためのメモリ使用量最大値
   * @param [in] fp     ファイルポインタ
   */
  void MemoryRequirement(const double TotalMemory, 
                         const double sphMemory, 
                         const double plot3dMemory, 
                         const double thinMemory, 
                         FILE* fp);

  /**
   * @brief 領域分割と出力DFIのインスタンス
   */
  virtual
  void VoxelInit() { return; }

  /**
   * @brief コーンバート処理
   */
  virtual
  bool exec()=0;  

  /**
   * @brief step番号からtimeを取得
   * @param[in] dfi dfiのポインター
   * @param[in] step step番号
   * @return time
   */
  double GetSliceTime(cdm_DFI* dfi, int step); 

  /**
   * @brief 配列のコンバート
   * @param[in]  buf         読込み用バッファ
   * @param[out] src         読み込んだデータを格納した配列のポインタ
   * @param[in]  headS       出力領域のheadインデックス
   * @param[in]  tailS       出力領域のtailインデックス
   * @param[in]  n           変数位置
   */
  bool convertXY(
                 cdm_Array* buf,
                 cdm_Array* &src,
                 int headS[3],
                 int tailS[3],
                 int n);
 
  /**
   * @brief 配列のコピー
   * @param [in]  B          コピー元の配列
   * @param [out] src        コピー先の配列
   * @param [in]  sta        コピーのスタート位置
   * @param [in]  end        コピーのエンド位置
   * @param [in]  n          変数位置
   */ 
  template<class T>
  bool copyArray(cdm_TypeArray<T> *B,
                 cdm_Array *&src,
                 int sta[3],
                 int end[3],
                 int n);
 
  /**
   * @brief 配列のコピー (template 関数)
   * @param [in]  buf        コピー元の配列
   * @param [out] src        コピー先の配列
   * @param [in]  sta        コピーのスタート位置
   * @param [in]  end        コピーのエンド位置
   * @param [in]  n          変数位置
   *
   */
  template<class T1, class T2>
  bool copyArray(cdm_TypeArray<T1> *buf,
                 cdm_TypeArray<T2> *&src,
                 int sta[3],
                 int end[3],
                 int n);

  /**
   * @brief データタイプ毎にminmaxを求める
   * @param[in]  src minmaxを求める配列データのポインタ
   * @param[out] min 求められた最小値
   * @param[out] max 求められた最大値
   */
  bool DtypeMinMax(cdm_Array* src,
                  double *min,
                  double *max);

  /**
   * @brief minmaxを計算
   * @param[in]  src minmaxを求める配列データのポインタ
   * @param[out] min 求められた最小値
   * @param[out] max 求められた最大値
   */
  template<class T>
  bool calcMinMax(cdm_TypeArray<T> *src,
                  double *min,
                  double *max);

  /**
   * @brief 出力ファイル形式から拡張子を求める
   * @param [in] file_format_type
   * @return 拡張子
   */
  static
  std::string GetFilenameExt(int file_format_type); 

  /**
   * @brief step基準のリスト生成
   * @param[out] StepRankList  step基準のリスト
   */
  void makeStepList(vector<step_rank_info> &StepRankList);

  /**
   * @brief rank基準のリスト生成
   * @param[out] StepRankList  rank基準のリスト
   */
  void makeRankList(vector<step_rank_info> &StepRankList);

  /**
   * @brief index.dfiの出力
   * @param minmaxList    minmaxのリスト
   */
  bool WriteIndexDfiFile(vector<dfi_MinMax*> minmaxList);

  /**
   * @brief proc.dfiの出力
   */
  bool WriteProcDfiFile(std::string proc_name,
                        cdm_Domain* out_domain,
                        cdm_MPI* out_mpi,
                        cdm_Process* out_process); 

  /**
   * @brief Proc情報の生成 
   */
  bool makeProcInfo(cdm_DFI* dfi,
                    cdm_Domain* &out_domain,
                    cdm_MPI* &out_mpi,
                    cdm_Process* &out_process,
                    int numProc);

};

//inline 関数
#include "inline/conv_inline.h"

#endif // _CONV_H_
