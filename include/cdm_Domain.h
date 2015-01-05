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

private:
  double Pitch[3];                    ///<計算空間のピッチ

protected:
  virtual void Clear();
public:
  /** コンストラクタ **/
  cdm_Domain();

  /** 
  * @brief コンストラクタ (_GlobalOrigin,_GlobalPitchは、double型とfloat型の両方あり)
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalPitch    ボクセルの長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  */ 
  cdm_Domain(const double* _GlobalOrigin,
             const double* _GlobalPitch, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision);
  cdm_Domain(const float* _GlobalOrigin, 
             const float* _GlobalPitch, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision);

  /** デストラクタ **/
  virtual ~cdm_Domain();

  /**
   * @brief セル中心のX座標を取得
   * @param [in] i X方向のセル番号
   * @return セル中心のX座標
   */
  virtual double CellX(int i) const{
    return GlobalOrigin[0] + Pitch[0]*(i+0.5);
  }

  /**
   * @brief セル中心のY座標を取得
   * @param [in] j Y方向のセル番号
   * @return セル中心のY座標
   */
  virtual double CellY(int j) const{
    return GlobalOrigin[1] + Pitch[1]*(j+0.5);
  }

  /**
   * @brief セル中心のZ座標を取得
   * @param [in] k Z方向のセル番号
   * @return セル中心のZ座標
   */
  virtual double CellZ(int k) const{
    return GlobalOrigin[2] + Pitch[2]*(k+0.5);
  }

  /**
   * @brief 格子点のX座標を取得
   * @param [in] i X方向の格子番号
   * @return 格子点のX座標
   */
  virtual double NodeX(int i) const{
    return GlobalOrigin[0] + Pitch[0]*i;
  }

  /**
   * @brief 格子点のY座標を取得
   * @param [in] j Y方向の格子番号
   * @return 格子点のY座標
   */
  virtual double NodeY(int j) const{
    return GlobalOrigin[1] + Pitch[1]*j;
  }

  /**
   * @brief 格子点のZ座標を取得
   * @param [in] k Z方向の格子番号
   * @return 格子点のZ座標
   */
  virtual double NodeZ(int k) const{
    return GlobalOrigin[2] + Pitch[2]*k;
  }

  /** 座標ファイル名を取得 */
  virtual std::string GetCoordinateFile() const{};

  /** 座標ファイルのファイルタイプを取得 */
  virtual CDM::E_CDM_FILE_TYPE GetCoordinateFileType() const{};

  /** 座標ファイルのデータ精度を取得 */
  virtual CDM::E_CDM_DTYPE GetCoordinateFilePrecision() const{};

  /** 座標ファイルのエンディアンタイプを取得 */
  virtual CDM::E_CDM_ENDIANTYPE GetCoordinateFileEndian() const{};

  /**
   * @brief read Domain(proc.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @param [in]   dirName DFIのディレクトリパス
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  virtual Read(cdm_TextParser tpCntl,
               std::string dirName);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp         ファイルポインタ
   * @param [in] tab        インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  virtual Write(FILE* fp, 
                const unsigned tab) const;

};

#endif // _CDM_DOMAIN_H_
