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

#include "cdm_Domain.h"

/** proc.dfi ファイルの Domain */
template<class T>
class cdm_NonUniformDomain : public cdm_Domain {

private:
  T *XCoordinates;                        ///<X座標データポインタ(Domainの格子点)
  T *YCoordinates;                        ///<Y座標データポインタ(Domainの格子点)
  T *ZCoordinates;                        ///<Z座標データポインタ(Domainの格子点)
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
  */ 
  cdm_NonUniformDomain(const T* _GlobalOrigin, 
                       const T* _GlobalPitch, 
                       const int* _GlobalVoxel, 
                       const int* _GlobalDivision,
                       const int* _iblank,
                       const T* _XCoordinates,
                       const T* _YCoordinates,
                       const T* _ZCoordinates)
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
    CoordinateFile = "";
    CoordinateFileFormat = CDM::E_CDM_OUTPUT_TYPE_DEFAULT;
    CoordinateFilePrecision = CDM::E_CDM_DTYPE_UNKNOWN;
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

};

template<typename T>
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain<T>::Read(cdm_TextParser tpCntl) 
{

  std::string str;
  std::string label;
  double v[3];
  int iv[3];

  //GlobalOrign
  label = "/Domain/GlobalOrigin";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_GLOBALORIGIN;
  }
  GlobalOrigin[0]=v[0];
  GlobalOrigin[1]=v[1];
  GlobalOrigin[2]=v[2];

  //GlobalRegion
  label = "/Domain/GlobalRegion";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_GLOBALREGION;
  }
  GlobalRegion[0]=v[0];
  GlobalRegion[1]=v[1];
  GlobalRegion[2]=v[2];

  //Global_Voxel
  label = "/Domain/GlobalVoxel";
  for (int n=0; n<3; n++) iv[n]=0;
  if ( !(tpCntl.GetVector(label, iv, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_GLOBALVOXEL;
  }
  GlobalVoxel[0]=iv[0];
  GlobalVoxel[1]=iv[1];
  GlobalVoxel[2]=iv[2];

  //Global_Division
  label = "/Domain/GlobalDivision";
  for (int n=0; n<3; n++) iv[n]=0;
  if ( !(tpCntl.GetVector(label, iv, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_GLOBALDIVISION;
  }
  GlobalDivision[0]=iv[0];
  GlobalDivision[1]=iv[1];
  GlobalDivision[2]=iv[2];

  //ActiveSubdomain
  label = "/Domain/ActiveSubdomainFile";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    str="";
  }
  ActiveSubdomainFile=str;

  //CoordinateFile
  label = "/Domain/CoordinateFile";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    str="";
  }
  CoordinateFile=str;

  //CoordinateFileFormat
  label = "/Domain/CoordinateFileFormat";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEFORMAT;
  }
  if( !strcasecmp(str.c_str(),"ascii" ) ) {
    CoordinateFileFormat=CDM::E_CDM_OUTPUT_TYPE_ASCII;
  } else if( !strcasecmp(str.c_str(),"binary" ) ) {
    CoordinateFileFormat=CDM::E_CDM_OUTPUT_TYPE_BINARY;
  } else {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEFORMAT;
  }

  //CoordinateFilePrecision
  if ( CoordinateFileFormat == CDM::E_CDM_OUTPUT_TYPE_BINARY ) {
    label = "/Domain/CoordinateFilePrecision";
    if ( !(tpCntl.GetValue(label, &str )) )
    {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEPRECISION;
    }
    if( !strcasecmp(str.c_str(),"float32" ) ) {
      CoordinateFilePrecision=CDM::E_CDM_FLOAT32;
    } else if( !strcasecmp(str.c_str(),"float64" ) ) {
      CoordinateFilePrecision=CDM::E_CDM_FLOAT64;
    } else {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEPRECISION;
    }
  }

  //Read CoordinateFile
  FILE* fp;
  CDM::E_CDM_ERRORCODE ret;
  if( !(fp=fopen(CoordinateFile.c_str(),"rb")) ) {
    printf("\tCan't open file. (%s)\n",CoordinateFile.c_str());
    return CDM::E_CDM_ERROR_OPEN_COORDINATEFILE;
  } else {
    ret = Read_CoordinateFile(fp);
    if( ret != CDM::E_CDM_SUCCESS ) {
      return ret;
    }
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// Domain の読込み関数
template<typename T>
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain<T>::Read_CoordinateFile(FILE* fp)
{

  int szGrid[3];
  double eps = 1e-10;
  double reg_Coord;

  //ascii
  if( CoordinateFileFormat == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {

    //x
    fscanf(fp,"%d\n",&szGrid[0]);
    //GlobalVoxelの確認
    if ( szGrid[0] != GlobalVoxel[0]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    T *XCoordinates = new T[szGrid[0]];
    for(int i=0; i<szGrid[0]; i++) {
      fscanf(fp,"%lf\n",&(XCoordinates[i]));
    }
    //GlobalOriginとGlobalRegionの確認
    if ( fabs(XCoordinates[0]-GlobalOrigin[0]) > eps) {
      printf("\tError in Read CoordinateFile: Origin in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    reg_Coord = XCoordinates[szGrid[0]-1]-XCoordinates[0];
    if ( fabs(reg_Coord-GlobalRegion[0]) > eps) {
      printf("\tError in Read CoordinateFile: Region in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }

    //y
    fscanf(fp,"%d\n",&szGrid[1]);
    if ( szGrid[1] != GlobalVoxel[1]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    T *YCoordinates = new T[szGrid[1]];
    for(int j=0; j<szGrid[1]; j++) {
      fscanf(fp,"%lf\n",&(YCoordinates[j]));
    }
    if ( fabs(YCoordinates[0]-GlobalOrigin[1]) > eps) {
      printf("\tError in Read CoordinateFile: Origin in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    reg_Coord = YCoordinates[szGrid[1]-1]-YCoordinates[0];
    if ( fabs(reg_Coord-GlobalRegion[1]) > eps) {
      printf("\tError in Read CoordinateFile: Region in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }

    //z
    fscanf(fp,"%d\n",&szGrid[2]);
    if ( szGrid[2] != GlobalVoxel[2]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    T *ZCoordinates = new T[szGrid[2]];
    for(int k=0; k<szGrid[2]; k++) {
      fscanf(fp,"%lf\n",&(ZCoordinates[k]));
    }
    if ( fabs(ZCoordinates[0]-GlobalOrigin[2]) > eps) {
      printf("\tError in Read CoordinateFile: Origin in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    reg_Coord = ZCoordinates[szGrid[2]-1]-ZCoordinates[0];
    if ( fabs(reg_Coord-GlobalRegion[2]) > eps) {
      printf("\tError in Read CoordinateFile: Region in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }

  //binary
  } else {

    //x
    fread(&szGrid[0], sizeof(int), 1, fp);
    //GlobalVoxelの確認
    if ( szGrid[0] != GlobalVoxel[0]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    T *XCoordinates = new T[szGrid[0]];
    fread(XCoordinates, sizeof(T), szGrid[0], fp);
    //GlobalOriginとGlobalRegionの確認
    if ( fabs(XCoordinates[0]-GlobalOrigin[0]) > eps) {
      printf("\tError in Read CoordinateFile: Origin in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    reg_Coord = XCoordinates[szGrid[0]-1]-XCoordinates[0];
    if ( fabs(reg_Coord-GlobalRegion[0]) > eps) {
      printf("\tError in Read CoordinateFile: Region in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }

    //y
    fread(&szGrid[1], sizeof(int), 1, fp);
    T *YCoordinates = new T[szGrid[1]];
    if ( szGrid[1] != GlobalVoxel[1]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    fread(YCoordinates, sizeof(T), szGrid[1], fp);
    if ( fabs(YCoordinates[0]-GlobalOrigin[1]) > eps) {
      printf("\tError in Read CoordinateFile: Origin in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    reg_Coord = YCoordinates[szGrid[1]-1]-YCoordinates[0];
    if ( fabs(reg_Coord-GlobalRegion[1]) > eps) {
      printf("\tError in Read CoordinateFile: Region in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }

    //z
    fread(&szGrid[2], sizeof(int), 1, fp);
    if ( szGrid[2] != GlobalVoxel[2]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    T *ZCoordinates = new T[szGrid[2]];
    fread(ZCoordinates, sizeof(T), szGrid[2], fp);
    if ( fabs(ZCoordinates[0]-GlobalOrigin[2]) > eps) {
      printf("\tError in Read CoordinateFile: Origin in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    reg_Coord = ZCoordinates[szGrid[2]-1]-ZCoordinates[0];
    if ( fabs(reg_Coord-GlobalRegion[2]) > eps) {
      printf("\tError in Read CoordinateFile: Region in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }

  }

  return CDM::E_CDM_SUCCESS;
}



#endif // _CDM_DOMAIN_H_
