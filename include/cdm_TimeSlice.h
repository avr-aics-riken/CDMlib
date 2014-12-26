#ifndef _CDM_TIMESLICE_H_
#define _CDM_TIMESLICE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_TimeSlice.h
 * @brief  cdm_Slice & cdm_TimeSliceClass Header
 * @author aics    
 */

/** index.dfi ファイルの Slice */
class cdm_Slice {

public :
  
  int    step;                          ///<ステップ番号
  double time;                          ///<時刻
  bool   avr_mode;                      ///<Average出力フラグ true:出力なし、false:出力
  int    AveragedStep;                  ///<平均ステップ
  double AveragedTime;                  ///<平均タイム
  double VectorMin;                     ///<Vectorのとき、最小値の合成値
  double VectorMax;                     ///<Vectorのとき、最大値の合成値
  vector<double> Min;                   ///<最小値
  vector<double> Max;                   ///<最大値

  /** コンストラクタ */
  cdm_Slice();

  /** デストラクタ */
  ~cdm_Slice();

  /**
  * @brief TimeSlice要素を読込む(inde.dfi)
  * @param [in]  tpCntl     cdm_TextParserクラス
  * @param [in]  label_leaf ベースとなる名前（"/TimeSlice/Slice")   
  * @param [in]  format ファイルフォーマット
  * @return error code
  */  
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl, 
       std::string label_leaf,
       CDM::E_CDM_FORMAT format);

  /**
   * @brief DFIファイル:TimeSlice要素を出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @param [in]  format ファイルフォーマット
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab,
        CDM::E_CDM_FORMAT format);

};

/** index.dfi ファイルの TimeSlice */
class cdm_TimeSlice {

public:

  vector<cdm_Slice> SliceList;

  /** コンストラクタ */
  cdm_TimeSlice();

  /** デストラクタ */
  ~cdm_TimeSlice();

  /**
  * @brief TimeSlice要素を読込む(inde.dfi)
  * @param [in]      tpCntl cdm_TextParserクラス
  * @param [in]      format ファイルフォーマット
  * @return error code
  */  
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl,
       CDM::E_CDM_FORMAT format); 

  /**
   * @brief DFIファイル:TimeSlice要素を出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @param [in] format   ファイルフォーマット
   * @return true:出力成功 false:出力失敗
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab,
        CDM::E_CDM_FORMAT format);

  /**
   * @brief DFIに出力されているminmaxの合成値を取得 
   * @param [in]  step 取得するステップ
   * @param [out] vec_min 取得したminの合成値
   * @param [out] vec_max 取得したminの合成値
   * @return error code 取得出来たときは E_CDM_SUCCESS 
   */
  CDM::E_CDM_ERRORCODE getVectorMinMax(const unsigned step,
                                       double &vec_min,
                                       double &vec_max);

  /**
   *brief DFIに出力されているminmaxとminmaxの合成値を取得
   * @param [in]  step 取得するステップ
   * @param [in]  variNo 取得する変数番号(0～n)
   * @param [out] min_value 取得したmin
   * @param [out] max_value 取得したmax
   * @return error code 取得出来たときは E_CDM_SUCCESS 
   */
  CDM::E_CDM_ERRORCODE getMinMax(const unsigned step,
                                 const int variNo,
                                 double &min_value,
                                 double &max_value);
 

  /**
   * @brief SliceListへの追加
   * @param [in]  step      ステップ番号
   * @param [in]  time      時刻
   * @param [in]  minmax    minmax
   * @param [in]  Nvari     変数の個数
   * @param [in]  format    ファイルフォーマット
   * @param [in]  avr_mode  Averageがあるかないかのフラグ
   * @param [in]  step_avr  Average step
   * @param [in]  time_avr  Average time
   */
  void AddSlice(int step,
                double time,
                double *minmax,
                int Nvari,
                CDM::E_CDM_FORMAT format,
                bool avr_mode,
                int step_avr,
                double time_avr);
                //vector<cdm_Slice> &SliceList);


};

#endif // _CDM_TIMESLICE_H_
