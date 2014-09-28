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
  const int* iblank;

private:
  double Pitch[3];                    ///<計算空間のピッチ

protected:
  virtual void Clear();
public:
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
  virtual ~cdm_Domain();

  double CellX(int i){
    return GlobalOrigin[0] + Pitch[0]*(i+0.5);
  }
  double CellY(int j){
    return GlobalOrigin[1] + Pitch[1]*(j+0.5);
  }
  double CellZ(int k){
    return GlobalOrigin[2] + Pitch[2]*(k+0.5);
  }
  double NodeX(int i){
    return GlobalOrigin[0] + Pitch[0]*i;
  }
  double NodeY(int j){
    return GlobalOrigin[1] + Pitch[1]*j;
  }
  double NodeZ(int k){
    return GlobalOrigin[2] + Pitch[2]*k;
  }
  template<class T>
  CDM::E_CDM_ERRORCODE CellXYZ(int i,int j,int k,T xyz[3]){
    xyz[0] = CellX(i);
    xyz[1] = CellY(j);
    xyz[2] = CellZ(k);
    return CDM::E_CDM_SUCCESS;
  }
  template<class T>
  CDM::E_CDM_ERRORCODE NodeXYZ(int i,int j,int k,T xyz[3]){
    xyz[0] = NodeX(i);
    xyz[1] = NodeY(j);
    xyz[2] = NodeZ(k);
    return CDM::E_CDM_SUCCESS;
  }

  cdm_Domain& operator=(const cdm_Domain& other){
    for(int i=0;i<3;++i){
      this->GlobalOrigin[i] = other.GlobalOrigin[i];
      this->GlobalRegion[i] = other.GlobalRegion[i];
      this->GlobalVoxel[i] = other.GlobalVoxel[i];
      this->GlobalDivision[i] = other.GlobalDivision[i];
      this->Pitch[i] = other.Pitch[i];
    }
    this->iblank = other.iblank;
  }

  static CDM::E_CDM_ERRORCODE Read(cdm_TextParser tpCntl,cdm_Domain* &domain);

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
