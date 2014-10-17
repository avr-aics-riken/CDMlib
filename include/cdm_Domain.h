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
  const int* iblank;                  ///<iblankデータポインタ(PLOT3Dのxyzファイル用)

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
  * @param [in] _iblank         iblankデータポインタ(PLOT3Dのxyzファイル用)
  */ 
  cdm_Domain(const double* _GlobalOrigin,
             const double* _GlobalPitch, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision,
             const int* _iblank);
  cdm_Domain(const float* _GlobalOrigin, 
             const float* _GlobalPitch, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision,
             const int* _iblank);

  /** デストラクタ **/
  virtual ~cdm_Domain();

  /** セル中心の座標を出力 */
  virtual double CellX(int i) const{
    return GlobalOrigin[0] + Pitch[0]*(i+0.5);
  }
  virtual double CellY(int j) const{
    return GlobalOrigin[1] + Pitch[1]*(j+0.5);
  }
  virtual double CellZ(int k) const{
    return GlobalOrigin[2] + Pitch[2]*(k+0.5);
  }
  /** 格子の座標を出力 */
  virtual double NodeX(int i) const{
    return GlobalOrigin[0] + Pitch[0]*i;
  }
  virtual double NodeY(int j) const{
    return GlobalOrigin[1] + Pitch[1]*j;
  }
  virtual double NodeZ(int k) const{
    return GlobalOrigin[2] + Pitch[2]*k;
  }

  virtual cdm_Domain& operator=(const cdm_Domain& other){
    for(int i=0;i<3;++i){
      this->GlobalOrigin[i] = other.GlobalOrigin[i];
      this->GlobalRegion[i] = other.GlobalRegion[i];
      this->GlobalVoxel[i] = other.GlobalVoxel[i];
      this->GlobalDivision[i] = other.GlobalDivision[i];
      this->Pitch[i] = other.Pitch[i];
    }
    this->iblank = other.iblank;
  }

  /**
   * @brief read Domain(proc.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  virtual Read(cdm_TextParser tpCntl);

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
