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
    if( XCoordinate != NULL ){ delete[] XCoordinate; }
    if( YCoordinate != NULL ){ delete[] YCoordinate; }
    if( ZCoordinate != NULL ){ delete[] ZCoordinate; }
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
             const double* _XCoordinate,
             const double* _YCoordinate,
             const double* _ZCoordinate)
  : cdm_Domain(_GlobalOrigin,_GlobalRegion,_GlobalVoxel,_GlobalDivision)
  {
    XCoordinate = new double[GlobalVoxel[0]+1];
    YCoordinate = new double[GlobalVoxel[1]+1];
    ZCoordinate = new double[GlobalVoxel[2]+1];
    for(int i=0;i<GlobalVoxel[0]+1;++i){
      XCoordinate[i] = _XCoordinate[i];
    }
    for(int j=0;j<GlobalVoxel[1]+1;++j){
      YCoordinate[j] = _YCoordinate[j];
    }
    for(int k=0;k<GlobalVoxel[2]+1;++k){
      ZCoordinate[k] = _ZCoordinate[k];
    }
  }

  /** デストラクタ **/
  ~cdm_NonUniformDomain()
  {
    if( XCoordinate != NULL ){ delete[] XCoordinate; }
    if( YCoordinate != NULL ){ delete[] YCoordinate; }
    if( ZCoordinate != NULL ){ delete[] ZCoordinate; }
  }

  double CellX(int i){
    return 0.5*(XCoordinate[i]+XCoordinate[i+1]);
  }
  double CellY(int j){
    return 0.5*(YCoordinate[j]+YCoordinate[j+1]);
  }
  double CellZ(int k){
    return 0.5*(ZCoordinate[k]+ZCoordinate[k+1]);
  }
  double NodeX(int i){
    return XCoordinate[i];
  }
  double NodeY(int j){
    return YCoordinate[j];
  }
  double NodeZ(int k){
    return ZCoordinate[k];
  }

  cdm_NonUniformDomain& operator=(const cdm_NonUniformDomain& other){
    Clear();
    for(int i=0;i<3;++i){
      this->GlobalOrigin[i] = other.GlobalOrigin[i];
      this->GlobalRegion[i] = other.GlobalRegion[i];
      this->GlobalVoxel[i] = other.GlobalVoxel[i];
      this->GlobalDivision[i] = other.GlobalDivision[i];
      this->Pitch[i] = other.Pitch[i];
    }
    XCoordinate = new double[GlobalVoxel[0]+1];
    for(int i=0;i<GlobalVoxel[0]+1;++i){
      this->XCoordinate[i] = other.XCoordinate[i];
    }
    YCoordinate = new double[GlobalVoxel[1]+1];
    for(int j=0;j<GlobalVoxel[1]+1;++j){
      this->YCoordinate[j] = other.YCoordinate[j];
    }
    ZCoordinate = new double[GlobalVoxel[2]+1];
    for(int k=0;k<GlobalVoxel[2]+1;++k){
      this->ZCoordinate[k] = other.ZCoordinate[k];
    }
  }

  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, const unsigned tab);

};

#endif // _CDM_DOMAIN_H_
