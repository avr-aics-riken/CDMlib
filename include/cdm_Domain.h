#ifndef _CDM_DOMAIN_H_
#define _CDM_DOMAIN_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_Domain.h
 * @brief  cdm_Domain Class Header
 * @author aics    
 */
  
/** proc.dfi ファイルの Domain */
class cdm_Domain {

public:

  double GlobalOrigin[3];             ///<計算空間の起点座標
  double GlobalRegion[3];             ///<計算空間の各軸方向の長さ
  int GlobalVoxel[3];                 ///<計算領域全体のボクセル数
  int GlobalDivision[3];              ///<計算領域の分割数
  std::string ActiveSubdomainFile;    ///<ActiveSubdomainファイル名

  /** コンストラクタ **/
  cdm_Domain();

  /** 
  * @brief コンストラクタ 
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalRegion   各軸方向の長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  */ 
  cdm_Domain(const double* _GlobalOrigin, 
             const double* _GlobalRegion, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision);

  /** デストラクタ **/
  ~cdm_Domain();

  /**
   * @brief read Domain(proc.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp         ファイルポインタ
   * @param [in] tab        インデント
   * @return true:出力成功 false:出力失敗
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

};

#endif // _CDM_DOMAIN_H_
