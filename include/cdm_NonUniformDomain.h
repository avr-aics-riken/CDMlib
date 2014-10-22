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
 * @author aics
 */

#include "cdm_Domain.h"

/** proc.dfi ファイルの Domain */
template<class T>
class cdm_NonUniformDomain : public cdm_Domain {

private:
  T *XCoordinates;                             ///<X座標データポインタ(Domainの格子点)
  T *YCoordinates;                             ///<Y座標データポインタ(Domainの格子点)
  T *ZCoordinates;                             ///<Z座標データポインタ(Domainの格子点)
  std::string CoordinateFile;                  ///<CoordinateFileファイル名
  CDM::E_CDM_OUTPUT_TYPE CoordinateFileFormat; ///<座標ファイルのデータフォーマット
  CDM::E_CDM_DTYPE CoordinateFilePrecision;    ///<座標ファイルのデータタイプ

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
  cdm_NonUniformDomain() : cdm_Domain(), XCoordinates(NULL), YCoordinates(NULL), ZCoordinates(NULL)
  {
    CoordinateFile = "";
    CoordinateFileFormat = CDM::E_CDM_OUTPUT_TYPE_DEFAULT;
    CoordinateFilePrecision = CDM::E_CDM_DTYPE_UNKNOWN;
  }

  /** 
  * @brief コンストラクタ 
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalPitch    ボクセルの長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  * @param [in] _iblank         iblankデータポインタ(PLOT3Dのxyzファイル用)
  * @param [in] _XCoordinates   X座標データポインタ(Domainの格子点)
  * @param [in] _YCoordinates   Y座標データポインタ(Domainの格子点)
  * @param [in] _ZCoordinates   Z座標データポインタ(Domainの格子点)
  * @param [in] _CoordinateFile          座標データ名
  * @param [in] _CoordinateFileFormat    座標データのファイルタイプ
  * @param [in] _CoordinateFilePrecision 座標データの精度
  */ 
  cdm_NonUniformDomain(const T* _GlobalOrigin, 
                       const T* _GlobalPitch, 
                       const int* _GlobalVoxel, 
                       const int* _GlobalDivision,
                       const int* _iblank,
                       const T* _XCoordinates,
                       const T* _YCoordinates,
                       const T* _ZCoordinates,
                       const std::string _CoordinateFile,
                       const CDM::E_CDM_OUTPUT_TYPE _CoordinateFileFormat,
                       const CDM::E_CDM_DTYPE _CoordinateFilePrecision)
  : cdm_Domain(_GlobalOrigin,_GlobalPitch,_GlobalVoxel,_GlobalDivision,_iblank)
  {
    XCoordinates = new T[GlobalVoxel[0]+1];
    YCoordinates = new T[GlobalVoxel[1]+1];
    ZCoordinates = new T[GlobalVoxel[2]+1];
    for(int i=0;i<GlobalVoxel[0]+1;++i){
      XCoordinates[i] = _XCoordinates[i];
    }
    for(int j=0;j<GlobalVoxel[1]+1;++j){
      YCoordinates[j] = _YCoordinates[j];
    }
    for(int k=0;k<GlobalVoxel[2]+1;++k){
      ZCoordinates[k] = _ZCoordinates[k];
    }
    CoordinateFile = _CoordinateFile;
    CoordinateFileFormat = _CoordinateFileFormat;
    CoordinateFilePrecision = _CoordinateFilePrecision;

    //Overwrite GlobalOrigin and GlobalRegion using XCoordinates, YCoordinates, ZCoordinates
    GlobalOrigin[0] = XCoordinates[0];
    GlobalOrigin[1] = YCoordinates[0];
    GlobalOrigin[2] = ZCoordinates[0];
    GlobalRegion[0] = XCoordinates[GlobalVoxel[0]] - XCoordinates[0];
    GlobalRegion[1] = YCoordinates[GlobalVoxel[1]] - YCoordinates[0];
    GlobalRegion[2] = ZCoordinates[GlobalVoxel[2]] - ZCoordinates[0];
  }

  /** デストラクタ **/
  ~cdm_NonUniformDomain()
  {
    if( XCoordinates != NULL ){ delete[] XCoordinates; }
    if( YCoordinates != NULL ){ delete[] YCoordinates; }
    if( ZCoordinates != NULL ){ delete[] ZCoordinates; }
  }

  /** セル中心の座標を出力 */
  double CellX(int i) const{
    return (double)0.5*(XCoordinates[i]+XCoordinates[i+1]);
  }
  double CellY(int j) const{
    return (double)0.5*(YCoordinates[j]+YCoordinates[j+1]);
  }
  double CellZ(int k) const{
    return (double)0.5*(ZCoordinates[k]+ZCoordinates[k+1]);
  }
  /** 格子の座標を出力 */
  double NodeX(int i) const{
    return (double)XCoordinates[i];
  }
  double NodeY(int j) const{
    return (double)YCoordinates[j];
  }
  double NodeZ(int k) const{
    return (double)ZCoordinates[k];
  }

  cdm_NonUniformDomain& operator=(const cdm_NonUniformDomain& other){
    Clear();
    for(int i=0;i<3;++i){
      this->GlobalOrigin[i] = other.GlobalOrigin[i];
      this->GlobalRegion[i] = other.GlobalRegion[i];
      this->GlobalVoxel[i] = other.GlobalVoxel[i];
      this->GlobalDivision[i] = other.GlobalDivision[i];
    }
    XCoordinates = new T[GlobalVoxel[0]+1];
    for(int i=0;i<GlobalVoxel[0]+1;++i){
      this->XCoordinates[i] = other.XCoordinates[i];
    }
    YCoordinates = new T[GlobalVoxel[1]+1];
    for(int j=0;j<GlobalVoxel[1]+1;++j){
      this->YCoordinates[j] = other.YCoordinates[j];
    }
    ZCoordinates = new T[GlobalVoxel[2]+1];
    for(int k=0;k<GlobalVoxel[2]+1;++k){
      this->ZCoordinates[k] = other.ZCoordinates[k];
    }
  }

  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  CDM::E_CDM_ERRORCODE
  Read_CoordinateFile(FILE* fp);

  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, const unsigned tab) const;

};

//inline 関数
#include "inline/cdm_NonUniformDomain_inline.h"

#endif // _CDM_NONUNIFORMDOMAIN_H_
