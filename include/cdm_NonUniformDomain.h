#ifndef _CDM_NONUNIFORMDOMAIN_H_
#define _CDM_NONUNIFORMDOMAIN_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_NonUniformDomain.h
 * @brief  cdm_NonUniformDomain Class Header
 * @author advancesoft
 */
  
/** proc.dfi ファイルの Domain */
class cdm_NonUniformDomain : public cdm_Domain {
private:
  double *XCoordinates;
  double *YCoordinates;
  double *ZCoordinates;
protected:
  virtual void Clear()
  {
    cdm_Domain::Clear();
    if( XCoordinates != NULL ){ delete[] XCoordinates; }
    if( YCoordinates != NULL ){ delete[] YCoordinates; }
    if( ZCoordinates != NULL ){ delete[] ZCoordinates; }
  }
public:

  /** コンストラクタ **/
  cdm_NonUniformDomain();

  /** 
  * @brief コンストラクタ 
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalRegion   各軸方向の長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  */ 
  cdm_NonUniformDomain(
             const double* _GlobalOrigin, 
             const double* _GlobalRegion, 
             const int* _GlobalVoxel, 
             const int* _GlobalDivision,
             const double* _XCoordinates,
             const double* _YCoordinates,
             const double* _ZCoordinates)
  : cdm_Domain(_GlobalOrigin,_GlobalRegion,_GlobalVoxel,_GlobalDivision)
  {
    XCoordinates = new double[GlobalVoxel[0]+1];
    YCoordinates = new double[GlobalVoxel[1]+1];
    ZCoordinates = new double[GlobalVoxel[2]+1];
    for(int i=0;i<GlobalVoxel[0]+1;++i){
      XCoordinates[i] = _XCoordinates[i];
    }
    for(int j=0;j<GlobalVoxel[1]+1;++j){
      YCoordinates[j] = _YCoordinates[j];
    }
    for(int k=0;k<GlobalVoxel[2]+1;++k){
      ZCoordinates[k] = _ZCoordinates[k];
    }
  }

  /** デストラクタ **/
  ~cdm_NonUniformDomain()
  {
    if( XCoordinates != NULL ){ delete[] XCoordinates; }
    if( YCoordinates != NULL ){ delete[] YCoordinates; }
    if( ZCoordinates != NULL ){ delete[] ZCoordinates; }
  }

  double CellX(int i){
    return 0.5*(XCoordinates[i]+XCoordinates[i+1]);
  }
  double CellY(int j){
    return 0.5*(YCoordinates[j]+YCoordinates[j+1]);
  }
  double CellZ(int k){
    return 0.5*(ZCoordinates[k]+ZCoordinates[k+1]);
  }
  double NodeX(int i){
    return XCoordinates[i];
  }
  double NodeY(int j){
    return YCoordinates[j];
  }
  double NodeZ(int k){
    return ZCoordinates[k];
  }

  cdm_NonUniformDomain& operator=(const cdm_NonUniformDomain& other){
    Clear();
    for(int i=0;i<3;++i){
      this->GlobalOrigin[i] = other.GlobalOrigin[i];
      this->GlobalRegion[i] = other.GlobalRegion[i];
      this->GlobalVoxel[i] = other.GlobalVoxel[i];
      this->GlobalDivision[i] = other.GlobalDivision[i];
    }
    XCoordinates = new double[GlobalVoxel[0]+1];
    for(int i=0;i<GlobalVoxel[0]+1;++i){
      this->XCoordinates[i] = other.XCoordinates[i];
    }
    YCoordinates = new double[GlobalVoxel[1]+1];
    for(int j=0;j<GlobalVoxel[1]+1;++j){
      this->YCoordinates[j] = other.YCoordinates[j];
    }
    ZCoordinates = new double[GlobalVoxel[2]+1];
    for(int k=0;k<GlobalVoxel[2]+1;++k){
      this->ZCoordinates[k] = other.ZCoordinates[k];
    }
  }

  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, const unsigned tab);

};

#endif // _CDM_DOMAIN_H_
