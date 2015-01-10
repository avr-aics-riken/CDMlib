#ifndef _CDM_NONUNIFORMDOMAIN_H_
#define _CDM_NONUNIFORMDOMAIN_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
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
  T *XCoordinates;                            ///<X座標データポインタ(Domainの格子点)
  T *YCoordinates;                            ///<Y座標データポインタ(Domainの格子点)
  T *ZCoordinates;                            ///<Z座標データポインタ(Domainの格子点)
  std::string CoordinateFile;                 ///<CoordinateFileファイル名
  CDM::E_CDM_FILE_TYPE CoordinateFileType;    ///<座標ファイルのファイルタイプ
  CDM::E_CDM_DTYPE CoordinateFilePrecision;   ///<座標ファイルのデータ精度
  CDM::E_CDM_ENDIANTYPE CoordinateFileEndian; ///<座標ファイルのエンディアンタイプ
  T pit_gcXsta;                               ///<X方向のガイドセルの格子幅(始点側)
  T pit_gcXend;                               ///<X方向のガイドセルの格子幅(終点側)
  T pit_gcYsta;                               ///<Y方向のガイドセルの格子幅(始点側)
  T pit_gcYend;                               ///<Y方向のガイドセルの格子幅(終点側)
  T pit_gcZsta;                               ///<Z方向のガイドセルの格子幅(始点側)
  T pit_gcZend;                               ///<Z方向のガイドセルの格子幅(終点側)

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
    CoordinateFileEndian = CDM::E_CDM_ENDIANTYPE_UNKNOWN;
  }

  /** 
  * @brief コンストラクタ 
  * @param [in] _XCoordinates   X座標データポインタ(Domainの格子点)
  * @param [in] _YCoordinates   Y座標データポインタ(Domainの格子点)
  * @param [in] _ZCoordinates   Z座標データポインタ(Domainの格子点)
  * @param [in] _CoordinateFile          座標データ名
  * @param [in] _CoordinateFileType      座標データのファイルタイプ
  * @param [in] _CoordinateFileEndian    座標データのエンディアンタイプ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  * @param [in] _gc             ガイドセル数
  */ 
  cdm_NonUniformDomain(const T* _XCoordinates,
                       const T* _YCoordinates,
                       const T* _ZCoordinates,
                       const std::string _CoordinateFile,
                       const CDM::E_CDM_FILE_TYPE _CoordinateFileType,
                       const CDM::E_CDM_ENDIANTYPE _CoordinateFileEndian,
                       const int* _GlobalVoxel, 
                       const int* _GlobalDivision,
                       const int _gc=0)
  {
    //GlobalVoxel,GlobalDivisionの設定
    GlobalVoxel[0]=_GlobalVoxel[0];
    GlobalVoxel[1]=_GlobalVoxel[1];
    GlobalVoxel[2]=_GlobalVoxel[2];

    GlobalDivision[0]=_GlobalDivision[0];
    GlobalDivision[1]=_GlobalDivision[1];
    GlobalDivision[2]=_GlobalDivision[2];

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
    CoordinateFileEndian = _CoordinateFileEndian;

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

  /** 
  * @brief コンストラクタ(proc.dfi出力用)
  * @param [in] _GlobalOrigin   起点座標
  * @param [in] _GlobalRegion   各軸方向の長さ
  * @param [in] _GlobalVoxel    ボクセル数
  * @param [in] _GlobalDivision 分割数
  * @param [in] _CoordinateFile          座標データ名
  * @param [in] _CoordinateFileType      座標データのファイルタイプ
  * @param [in] _CoordinateFilePrecision 座標データのデータ精度
  * @param [in] _CoordinateFileEndian    座標データのエンディアンタイプ
  */ 
  cdm_NonUniformDomain(const T* _GlobalOrigin,
                       const T* _GlobalRegion,
                       const int* _GlobalVoxel, 
                       const int* _GlobalDivision,
                       const std::string _CoordinateFile,
                       const CDM::E_CDM_FILE_TYPE _CoordinateFileType,
                       const CDM::E_CDM_DTYPE _CoordinateFilePrecision,
                       const CDM::E_CDM_ENDIANTYPE _CoordinateFileEndian)
  {
    GlobalOrigin[0]=_GlobalOrigin[0];
    GlobalOrigin[1]=_GlobalOrigin[1];
    GlobalOrigin[2]=_GlobalOrigin[2];

    GlobalRegion[0]=_GlobalRegion[0];
    GlobalRegion[1]=_GlobalRegion[1];
    GlobalRegion[2]=_GlobalRegion[2];

    GlobalVoxel[0]=_GlobalVoxel[0];
    GlobalVoxel[1]=_GlobalVoxel[1];
    GlobalVoxel[2]=_GlobalVoxel[2];

    GlobalDivision[0]=_GlobalDivision[0];
    GlobalDivision[1]=_GlobalDivision[1];
    GlobalDivision[2]=_GlobalDivision[2];

    CoordinateFile = _CoordinateFile;
    CoordinateFileType = _CoordinateFileType;
    CoordinateFilePrecision = _CoordinateFilePrecision;
    CoordinateFileEndian = _CoordinateFileEndian;

  }

  /** デストラクタ **/
  ~cdm_NonUniformDomain()
  {
    if( XCoordinates != NULL ){ delete[] XCoordinates; }
    if( YCoordinates != NULL ){ delete[] YCoordinates; }
    if( ZCoordinates != NULL ){ delete[] ZCoordinates; }
  }

  /**
   * @brief セル中心のX座標を取得
   * @param [in] i X方向のセル番号
   * @return セル中心のX座標
   */
  double CellX(int i) const{
    if( i < 0 ){
      return (double)0.5*(XCoordinates[0]+XCoordinates[1]) - (double)pit_gcXsta*(-i);
    }else if( i >= GlobalVoxel[0] ){
      return (double)0.5*(XCoordinates[GlobalVoxel[0]-1]+XCoordinates[GlobalVoxel[0]]) + (double)pit_gcXend*(i-(GlobalVoxel[0]-1));
    }else{
      return (double)0.5*(XCoordinates[i]+XCoordinates[i+1]);
    }
  }

  /**
   * @brief セル中心のY座標を取得
   * @param [in] j Y方向のセル番号
   * @return セル中心のX座標
   */
  double CellY(int j) const{
    if( j < 0 ) {
      return (double)0.5*(YCoordinates[0]+YCoordinates[1]) - (double)pit_gcYsta*(-j);
    } else if( j >= GlobalVoxel[1] ) {
      return (double)0.5*(YCoordinates[GlobalVoxel[1]-1]+YCoordinates[GlobalVoxel[1]]) + (double)pit_gcYend*(j-(GlobalVoxel[1]-1));
    } else {
      return (double)0.5*(YCoordinates[j]+YCoordinates[j+1]);
    }
  }

  /**
   * @brief セル中心のZ座標を取得
   * @param [in] k Z方向のセル番号
   * @return セル中心のX座標
   */
  double CellZ(int k) const{
    if( k < 0 ) {
      return (double)0.5*(ZCoordinates[0]+ZCoordinates[1]) - (double)pit_gcZsta*(-k);
    } else if( k >= GlobalVoxel[2] ) {
      return (double)0.5*(ZCoordinates[GlobalVoxel[2]-1]+ZCoordinates[GlobalVoxel[2]]) + (double)pit_gcZend*(k-(GlobalVoxel[2]-1));
    } else {
      return (double)0.5*(ZCoordinates[k]+ZCoordinates[k+1]);
    }
  }

  /**
   * @brief 格子点のX座標を取得
   * @param [in] i X方向の格子番号
   * @return 格子点のX座標
   */
  double NodeX(int i) const{
    return (double)XCoordinates[i];
  }

  /**
   * @brief 格子点のY座標を取得
   * @param [in] j Y方向の格子番号
   * @return 格子点のY座標
   */
  double NodeY(int j) const{
    return (double)YCoordinates[j];
  }

  /**
   * @brief 格子点のZ座標を取得
   * @param [in] k Z方向の格子番号
   * @return 格子点のZ座標
   */
  double NodeZ(int k) const{
    return (double)ZCoordinates[k];
  }

  /**
   * @brief 座標ファイル名を取得
   * @return 座標ファイル名
   */
  std::string GetCoordinateFile() const{
    return CoordinateFile;
  }

  /**
   * @brief 座標ファイルのファイルタイプを取得
   * @return 座標ファイルのファイルタイプ
   */
  CDM::E_CDM_FILE_TYPE GetCoordinateFileType() const{
    return CoordinateFileType;
  }

  /**
   * @brief 座標ファイルのデータ精度を取得
   * @return 座標ファイルのデータ精度
   */
  CDM::E_CDM_DTYPE GetCoordinateFilePrecision() const{
    return CoordinateFilePrecision;
  }

  /**
   * @brief 座標ファイルのエンディアンタイプを取得
   * @return 座標ファイルのエンディアンタイプ
   */
  CDM::E_CDM_ENDIANTYPE GetCoordinateFileEndian() const{
    return CoordinateFileEndian;
  }

  /**
   * @brief Domain(proc.dfi)を読込む
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @param [in]   dirName DFIのディレクトリパス
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl,
       std::string dirName);

  /**
   * @brief CoordinateFileを読込む
   * @param [in]   fp  ファイルポインタ
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read_CoordinateFile(FILE* fp);

  /**
   * @brief CoordinateFileの各方向の座標データ数を読込む
   * @param [in]   fp          ファイルポインタ
   * @param [in]   matchEndian true:Endian一致
   * @param [in]   globalVoxel 計算領域のボクセル数
   * @param [out]  dataCount   座標データ数
   * @param [out]  nread       読込みデータ数
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  readCoordDataCount(FILE* fp,
                     bool matchEndian,
                     int globalVoxel,
                     int* dataCount,
                     size_t* nread);

  /**
   * @brief CoordinateFileの各方向の座標データを読込む
   * @param [in]   fp           ファイルポインタ
   * @param [in]   matchEndian  true:Endian一致
   * @param [in]   dataCount    座標データ数
   * @param [in]   globalOrigin 計算空間の起点座標
   * @param [in]   globalRegion 計算空間の各軸方向の長さ
   * @param [out]  coordinates  読み込んだ座標データのポインタ
   * @param [out]  nread        読込みデータ数
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  readCoordData(FILE* fp,
                bool matchEndian,
                int* dataCount,
                double globalOrigin,
                double globalRegion,
                T* coordinates,
                size_t* nread);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in]   fp  ファイルポインタ
   * @param [in]   tab インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp,
        const unsigned tab) const;

};

//inline 関数
#include "inline/cdm_NonUniformDomain_inline.h"

#endif // _CDM_NONUNIFORMDOMAIN_H_
