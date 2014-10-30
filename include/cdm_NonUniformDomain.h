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
#include <typeinfo>

/** proc.dfi ファイルの Domain */
template<class T>
class cdm_NonUniformDomain : public cdm_Domain {

private:
  T *XCoordinates;                           ///<X座標データポインタ(Domainの格子点)
  T *YCoordinates;                           ///<Y座標データポインタ(Domainの格子点)
  T *ZCoordinates;                           ///<Z座標データポインタ(Domainの格子点)
  std::string CoordinateFile;                ///<CoordinateFileファイル名
  CDM::E_CDM_FILE_TYPE CoordinateFileType;   ///<座標ファイルのファイルタイプ
  CDM::E_CDM_DTYPE CoordinateFilePrecision;  ///<座標ファイルのデータ精度
  T pit_gcXsta;                              ///<X方向のガイドセルの格子幅(始点側)
  T pit_gcXend;                              ///<X方向のガイドセルの格子幅(終点側)
  T pit_gcYsta;                              ///<Y方向のガイドセルの格子幅(始点側)
  T pit_gcYend;                              ///<Y方向のガイドセルの格子幅(終点側)
  T pit_gcZsta;                              ///<Z方向のガイドセルの格子幅(始点側)
  T pit_gcZend;                              ///<Z方向のガイドセルの格子幅(終点側)

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
    CoordinateFileType = CDM::E_CDM_FILE_TYPE_DEFAULT;
    CoordinateFilePrecision = CDM::E_CDM_DTYPE_UNKNOWN;
  }

  /** 
  * @brief コンストラクタ 
  * @param [in] _XCoordinates   X座標データポインタ(Domainの格子点)
  * @param [in] _YCoordinates   Y座標データポインタ(Domainの格子点)
  * @param [in] _ZCoordinates   Z座標データポインタ(Domainの格子点)
  * @param [in] _CoordinateFile          座標データ名
  * @param [in] _CoordinateFileType      座標データのファイルタイプ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  * @param [in] _iblank         iblankデータポインタ(PLOT3Dのxyzファイル用)
  * @param [in] _gc             ガイドセル数
  */ 
  cdm_NonUniformDomain(const T* _XCoordinates,
                       const T* _YCoordinates,
                       const T* _ZCoordinates,
                       const std::string _CoordinateFile,
                       const CDM::E_CDM_FILE_TYPE _CoordinateFileType,
                       const int* _GlobalVoxel, 
                       const int* _GlobalDivision,
                       const int* _iblank,
                       const int _gc=0)
  {
    //GlobalVoxel,GlobalDivision,iblankの設定
    GlobalVoxel[0]=_GlobalVoxel[0];
    GlobalVoxel[1]=_GlobalVoxel[1];
    GlobalVoxel[2]=_GlobalVoxel[2];

    GlobalDivision[0]=_GlobalDivision[0];
    GlobalDivision[1]=_GlobalDivision[1];
    GlobalDivision[2]=_GlobalDivision[2];

    iblank = _iblank;

    // 座標に関する設定
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
    CoordinateFileType = _CoordinateFileType;
    if( typeid(T) == typeid(float) ){
      CoordinateFilePrecision = CDM::E_CDM_FLOAT32;
    } else if( typeid(T) == typeid(double) ) {
      CoordinateFilePrecision = CDM::E_CDM_FLOAT64;
    }

    //GlobalOrigin,GlobalRegionの設定
    GlobalOrigin[0] = XCoordinates[0];
    GlobalOrigin[1] = YCoordinates[0];
    GlobalOrigin[2] = ZCoordinates[0];

    GlobalRegion[0] = XCoordinates[GlobalVoxel[0]] - XCoordinates[0];
    GlobalRegion[1] = YCoordinates[GlobalVoxel[1]] - YCoordinates[0];
    GlobalRegion[2] = ZCoordinates[GlobalVoxel[2]] - ZCoordinates[0];

    //ガイドセルがある場合、ガイドセルのピッチ幅を算出
    if( _gc>0 ) {
      pit_gcXsta = XCoordinates[1] - XCoordinates[0];
      pit_gcXend = XCoordinates[GlobalVoxel[0]] - XCoordinates[GlobalVoxel[0]-1];
      pit_gcYsta = YCoordinates[1] - YCoordinates[0];
      pit_gcYend = YCoordinates[GlobalVoxel[1]] - YCoordinates[GlobalVoxel[1]-1];
      pit_gcZsta = ZCoordinates[1] - ZCoordinates[0];
      pit_gcZend = ZCoordinates[GlobalVoxel[2]] - ZCoordinates[GlobalVoxel[2]-1];
    }
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
    if( i < 0 ){
      return (double)0.5*(XCoordinates[0]+XCoordinates[1]) - (double)pit_gcXsta*(-i);
    }else if( i >= GlobalVoxel[0] ){
      return (double)0.5*(XCoordinates[GlobalVoxel[0]-1]+XCoordinates[GlobalVoxel[0]]) + (double)pit_gcXend*(i-(GlobalVoxel[0]-1));
    }else{
      return (double)0.5*(XCoordinates[i]+XCoordinates[i+1]);
    }
  }
  double CellY(int j) const{
    if( j < 0 ) {
      return (double)0.5*(YCoordinates[0]+YCoordinates[1]) - (double)pit_gcYsta*(-j);
    } else if( j >= GlobalVoxel[1] ) {
      return (double)0.5*(YCoordinates[GlobalVoxel[1]-1]+YCoordinates[GlobalVoxel[1]]) + (double)pit_gcYend*(j-(GlobalVoxel[1]-1));
    } else {
      return (double)0.5*(YCoordinates[j]+YCoordinates[j+1]);
    }
  }
  double CellZ(int k) const{
    if( k < 0 ) {
      return (double)0.5*(ZCoordinates[0]+ZCoordinates[1]) - (double)pit_gcZsta*(-k);
    } else if( k >= GlobalVoxel[2] ) {
      return (double)0.5*(ZCoordinates[GlobalVoxel[2]-1]+ZCoordinates[GlobalVoxel[2]]) + (double)pit_gcZend*(k-(GlobalVoxel[2]-1));
    } else {
      return (double)0.5*(ZCoordinates[k]+ZCoordinates[k+1]);
    }
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
