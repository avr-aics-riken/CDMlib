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
  double *XCoordinates;                        ///<X座標データポインタ(Domainの格子点)
  double *YCoordinates;                        ///<Y座標データポインタ(Domainの格子点)
  double *ZCoordinates;                        ///<Z座標データポインタ(Domainの格子点)
  std::string CoordinateFile;                  ///<CoordinateFileファイル名
  CDM::E_CDM_OUTPUT_TYPE CoordinateFileFormat; ///<座標ファイルのデータフォーマット
  CDM::E_CDM_DTYPE CoordinateFilePrecision;    ///<座標ファイルのデータタイプ

protected:
  virtual void Clear();

public:
  /** コンストラクタ **/
  cdm_NonUniformDomain();

  /** 
  * @brief コンストラクタ 
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalRegion   各軸方向の長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  * @param [in] _iblank         iblankデータポインタ(PLOT3Dのxyzファイル用)
  * @param [in] _XCoordinates   X座標データポインタ(Domainの格子点)
  * @param [in] _YCoordinates   Y座標データポインタ(Domainの格子点)
  * @param [in] _ZCoordinates   Z座標データポインタ(Domainの格子点)
  */ 
  cdm_NonUniformDomain(const double* _GlobalOrigin, 
                          const double* _GlobalRegion, 
                          const int* _GlobalVoxel, 
                          const int* _GlobalDivision,
                          const int* _iblank,
                          const double* _XCoordinates,
                          const double* _YCoordinates,
                          const double* _ZCoordinates);

  /** デストラクタ **/
  ~cdm_NonUniformDomain();

  double CellX(int i) const{
    //if( XCoordinates != NULL ){}
    return 0.5*(XCoordinates[i]+XCoordinates[i+1]);
  }
  double CellY(int j) const{
    return 0.5*(YCoordinates[j]+YCoordinates[j+1]);
  }
  double CellZ(int k) const{
    return 0.5*(ZCoordinates[k]+ZCoordinates[k+1]);
  }
  double NodeX(int i) const{
    return XCoordinates[i];
  }
  double NodeY(int j) const{
    return YCoordinates[j];
  }
  double NodeZ(int k) const{
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
  Read_CoordinateFile(FILE* fp);

};

#endif // _CDM_DOMAIN_H_
