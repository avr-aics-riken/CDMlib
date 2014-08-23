#ifndef _CIO_DOMAIN_H_
#define _CIO_DOMAIN_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_Domain.h
 * @brief  cio_Domain Class Header
 * @author aics    
 */
  
/** proc.dfi ファイルの Domain */
class cio_Domain {

public:

  double GlobalOrigin[3];             ///<計算空間の起点座標
  double GlobalRegion[3];             ///<計算空間の各軸方向の長さ
  int GlobalVoxel[3];                 ///<計算領域全体のボクセル数
  int GlobalDivision[3];              ///<計算領域の分割数
  std::string ActiveSubdomainFile;    ///<ActiveSubdomainファイル名

  /** コンストラクタ **/
  cio_Domain();

  /** 
  * @brief コンストラクタ 
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalRegion   各軸方向の長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  */ 
  cio_Domain(const double* _GlobalOrigin, 
             const double* _GlobalRegion, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision);

  /** デストラクタ **/
  ~cio_Domain();

  /**
   * @brief read Domain(proc.dfi)
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Read(cio_TextParser tpCntl);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp         ファイルポインタ
   * @param [in] tab        インデント
   * @return true:出力成功 false:出力失敗
   */
  CIO::E_CIO_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab);

};

#endif // _CIO_DOMAIN_H_
