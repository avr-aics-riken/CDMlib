#ifndef _CIO_TIMESLICE_H_
#define _CIO_TIMESLICE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_TimeSlice.h
 * @brief  cio_Slice & cio_TimeSliceClass Header
 * @author aics    
 */

/** index.dfi ファイルの Slice */
class cio_Slice {

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
  cio_Slice();

  /** デストラクタ */
  ~cio_Slice();

  /**
  * @brief TimeSlice要素を読込む(inde.dfi)
  * @param [in]  tpCntl     cio_TextParserクラス
  * @param [in]  label_leaf ベースとなる名前（"/TimeSlice/Slice")   
  * @return error code
  */  
  CIO::E_CIO_ERRORCODE
  Read(cio_TextParser tpCntl, 
       std::string label_leaf);

  /**
   * @brief DFIファイル:TimeSlice要素を出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

};

/** index.dfi ファイルの TimeSlice */
class cio_TimeSlice {

public:

  vector<cio_Slice> SliceList;

  /** コンストラクタ */
  cio_TimeSlice();

  /** デストラクタ */
  ~cio_TimeSlice();

  /**
  * @brief TimeSlice要素を読込む(inde.dfi)
  * @param [in]      tpCntl cio_TextParserクラス
  * @return error code
  */  
  CIO::E_CIO_ERRORCODE
  Read(cio_TextParser tpCntl); 

  /**
   * @brief DFIファイル:TimeSlice要素を出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @return true:出力成功 false:出力失敗
   */
  CIO::E_CIO_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

  /**
   * @brief DFIに出力されているminmaxの合成値を取得 
   * @param [in]  step 取得するステップ
   * @param [out] vec_min 取得したminの合成値
   * @param [out] vec_max 取得したminの合成値
   * @return error code 取得出来たときは E_CIO_SUCCESS 
   */
  CIO::E_CIO_ERRORCODE getVectorMinMax(const unsigned step,
                                       double &vec_min,
                                       double &vec_max);

  /**
   *brief DFIに出力されているminmaxとminmaxの合成値を取得
   * @param [in]  step 取得するステップ
   * @param [in]  compNo 取得する成分番号(0～n)
   * @param [out] min_value 取得したmin
   * @param [out] max_value 取得したmax
   * @return error code 取得出来たときは E_CIO_SUCCESS 
   */
  CIO::E_CIO_ERRORCODE getMinMax(const unsigned step,
                                 const int compNo,
                                 double &min_value,
                                 double &max_value);
 

  /**
   * @brief SliceListへの追加
   * @param [in]  step      ステップ番号
   * @param [in]  time      時刻
   * @param [in]  minmax    minmax
   * @param [in]  Ncomp     コンポーネント数
   * @param [in]  avr_mode  Averageがあるかないかのフラグ
   * @param [in]  step_avr  Average step
   * @param [in]  time_avr  Average time
   */
  void AddSlice(int step,
                double time,
                double *minmax,
                int Ncomp,
                bool avr_mode,
                int step_avr,
                double time_avr);
                //vector<cio_Slice> &SliceList);


};

#endif // _CIO_TIMESLICE_H_
