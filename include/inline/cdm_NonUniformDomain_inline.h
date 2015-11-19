#ifndef _CDM_NONUNIFORMDOMAIN_INLINE_H_
#define _CDM_NONUNIFORMDOMAIN_INLINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_NonUniformDomain_inline.h
 * @brief  cdm_NonUniformDomain template Header
 * @author aics    
 */

#ifdef CDM_INLINE
 #undef CDM_INLINE
#endif

#ifndef CDM_NO_INLINE
 #define CDM_INLINE inline
#else
 #define CDM_INLINE
#endif

// #################################################################
// Domain の読込み関数
template<typename T>
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain<T>::Read(cdm_TextParser tpCntl,
                              std::string dirName)
{

  std::string str;
  std::string label;
  double v[3];
  int iv[3];

  //GlobalOrigin
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
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILE;
  }
  CoordinateFile=str;
  //Check of extension of CoordinateFile
  string ext;
  size_t pos_dot = CoordinateFile.rfind('.');
  if(pos_dot != string::npos){
    ext = CoordinateFile.substr(pos_dot+1, CoordinateFile.size()-pos_dot);
  } else {
    printf("\tFail to get extension of '%s'\n",CoordinateFile.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILE;
  }
  if(ext != "crd"){
    printf("\tExtension of Coordinate File '%s' is not 'crd'. \n",ext.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILE;
  }

  //CoordinateFileType
  label = "/Domain/CoordinateFileType";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILETYPE;
  }
  if( !strcasecmp(str.c_str(),"ascii" ) ) {
    CoordinateFileType=CDM::E_CDM_FILE_TYPE_ASCII;
  } else if( !strcasecmp(str.c_str(),"binary" ) ) {
    CoordinateFileType=CDM::E_CDM_FILE_TYPE_BINARY;
  } else {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILETYPE;
  }

  //CoordinateFilePrecision
  if ( CoordinateFileType == CDM::E_CDM_FILE_TYPE_BINARY ) {
    label = "/Domain/CoordinateFilePrecision";
    if ( !(tpCntl.GetValue(label, &str )) )
    {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEPRECISION;
    }
    if( !strcasecmp(str.c_str(),"Float32" ) ) {
      CoordinateFilePrecision=CDM::E_CDM_FLOAT32;
    } else if( !strcasecmp(str.c_str(),"Float64" ) ) {
      CoordinateFilePrecision=CDM::E_CDM_FLOAT64;
    } else {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEPRECISION;
    }
  }

  //CoordinateFileEndian
  if ( CoordinateFileType == CDM::E_CDM_FILE_TYPE_BINARY ) {
    label = "/Domain/CoordinateFileEndian";
    if ( !(tpCntl.GetValue(label, &str )) )
    {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEENDIAN;
    }
    if( !strcasecmp(str.c_str(),"little" ) ) {
      CoordinateFileEndian=CDM::E_CDM_LITTLE;
    } else if( !strcasecmp(str.c_str(),"big" ) ) {
      CoordinateFileEndian=CDM::E_CDM_BIG;
    } else {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_COORDINATEFILEENDIAN;
    }
  }

  //Read CoordinateFile
  FILE* fp;
  std::string crdfile;
  CDM::E_CDM_ERRORCODE ret;
    if( CDM::cdmPath_isAbsolute(CoordinateFile) ){
      crdfile = CoordinateFile;
    } else {
      crdfile = CDM::cdmPath_ConnectPath(dirName,CoordinateFile);
    }

  if( !(fp=fopen(crdfile.c_str(),"rb")) ) {
    printf("\tCan't open file. (%s)\n",crdfile.c_str());
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
// CoordinateFile の読込み関数
template<typename T>
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain<T>::Read_CoordinateFile(FILE* fp)
{

  int szGrid[3];
  double eps = 1e-10;
  double reg_Coord;

  //ascii
  if( CoordinateFileType == CDM::E_CDM_FILE_TYPE_ASCII ) {

    //x
    fscanf(fp,"%d\n",&szGrid[0]);
    //GlobalVoxelの確認
    if ( szGrid[0] != GlobalVoxel[0]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    if( XCoordinates ){ delete[] XCoordinates; }
    XCoordinates = new T[szGrid[0]];
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
    //GlobalVoxelの確認
    if ( szGrid[1] != GlobalVoxel[1]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    if( YCoordinates ){ delete[] YCoordinates; }
    YCoordinates = new T[szGrid[1]];
    for(int j=0; j<szGrid[1]; j++) {
      fscanf(fp,"%lf\n",&(YCoordinates[j]));
    }

    //GlobalOriginとGlobalRegionの確認
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
    //GlobalVoxelの確認
    if ( szGrid[2] != GlobalVoxel[2]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    if( ZCoordinates ){ delete[] ZCoordinates; }
    ZCoordinates = new T[szGrid[2]];
    for(int k=0; k<szGrid[2]; k++) {
      fscanf(fp,"%lf\n",&(ZCoordinates[k]));
    }

    //GlobalOriginとGlobalRegionの確認
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

    /** Endian セット */
    int idumy = 1;
    char* cdumy = (char*)(&idumy);
    CDM::E_CDM_ENDIANTYPE Endian=CDM::E_CDM_ENDIANTYPE_UNKNOWN;
    if( cdumy[0] == 0x01 ) Endian = CDM::E_CDM_LITTLE;
    if( cdumy[0] == 0x00 ) Endian = CDM::E_CDM_BIG;

    bool matchEndian = true;
    if( Endian != CoordinateFileEndian ) matchEndian = false;

    size_t nread_sum = 0;
    size_t nread = 0;
    int xdataCount;
    int ydataCount;
    int zdataCount;
    CDM::E_CDM_ERRORCODE tmpRes;

    //x
    tmpRes = readCoordDataCount(fp, matchEndian, GlobalVoxel[0], &xdataCount, &nread);
    if( tmpRes != CDM::E_CDM_SUCCESS ) {
      return tmpRes;
    }
    nread_sum += nread;
    if( XCoordinates ){ delete[] XCoordinates; }
    XCoordinates = new T[xdataCount];
    tmpRes = readCoordData(fp, matchEndian, &xdataCount, GlobalOrigin[0], GlobalRegion[0], XCoordinates, &nread);
    if( tmpRes != CDM::E_CDM_SUCCESS) {
      return tmpRes;
    }
    nread_sum += nread;

    //y
    tmpRes = readCoordDataCount(fp, matchEndian, GlobalVoxel[1], &ydataCount, &nread);
    if( tmpRes != CDM::E_CDM_SUCCESS ) {
        return tmpRes;
    }
    nread_sum += nread;
    if( YCoordinates ){ delete[] YCoordinates; }
    YCoordinates = new T[ydataCount];
    tmpRes = readCoordData(fp, matchEndian, &ydataCount, GlobalOrigin[1], GlobalRegion[1], YCoordinates, &nread);
    if( tmpRes != CDM::E_CDM_SUCCESS ) {
        return tmpRes;
    }
    nread_sum += nread;

    //z
    tmpRes = readCoordDataCount(fp, matchEndian, GlobalVoxel[2], &zdataCount, &nread);
    if( tmpRes != CDM::E_CDM_SUCCESS ) {
        return tmpRes;
    }
    nread_sum += nread;
    if( ZCoordinates ){ delete[] ZCoordinates; }
    ZCoordinates = new T[zdataCount];
    tmpRes = readCoordData(fp, matchEndian, &zdataCount, GlobalOrigin[2], GlobalRegion[2], ZCoordinates, &nread);
    if( tmpRes != CDM::E_CDM_SUCCESS ) {
        return tmpRes;
    }
    nread_sum += nread;

    //CoordinatesFileから読み込んだデータ数のチェック
    if ( nread_sum != (xdataCount+ydataCount+zdataCount+3) ) {
      printf("\tError in Read CoordinateFile: Num of CoordinatesFile data\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// CoordinateFileの各方向の座標データ数を読込む関数
template<typename T>
CDM_INLINE
CDM::E_CDM_ERRORCODE 
cdm_NonUniformDomain<T>::readCoordDataCount(FILE* fp,
                                            bool matchEndian,
                                            int globalVoxel,
                                            int* dataCount,
                                            size_t* nread)
{

  *nread = fread(dataCount, sizeof(int), 1, fp);
  if( !matchEndian ) BSWAPVEC(dataCount,1);

  //GlobalVoxelの確認
  if ( dataCount[0] != globalVoxel+1 ) {
    printf("\tError in Read CoordinateFile: Number of grid\n");
    return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// CoordinateFileの各方向の座標データを読込む関数
template<typename T>
CDM_INLINE
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain<T>::readCoordData(FILE* fp,
                                       bool matchEndian,
                                       int* dataCount,
                                       double globalOrigin,
                                       double globalRegion,
                                       T* coordinates,
                                       size_t* nread)
{

  *nread = fread(coordinates, sizeof(T), dataCount[0], fp);
  if( !matchEndian ) {
    if( typeid(T) == typeid(float) ){
      BSWAPVEC(coordinates, (size_t)dataCount[0]);
    } else if( typeid(T) == typeid(double) ) {
      DBSWAPVEC(coordinates, (size_t)dataCount[0]);
    }
  }

  //GlobalOriginとGlobalRegionの確認
  double eps = 1e-10;
  if ( fabs(coordinates[0]-globalOrigin) > eps) {
    printf("\tError in Read CoordinateFile: Origin\n");
    cout << globalOrigin << endl;
    return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
  }
  double reg_Coord = coordinates[dataCount[0]-1]-coordinates[0];
  if ( fabs(reg_Coord-globalRegion) > eps) {
    printf("\tError in Read CoordinateFile: Region\n");
    return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// DFIファイル:Domain要素を出力する
template<typename T>
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain<T>::Write(FILE* fp, 
                               const unsigned tab) const
{

  fprintf(fp, "Domain {\n");
  fprintf(fp, "\n");  

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "GlobalOrigin        = (%e, %e, %e)\n",
          GlobalOrigin[0],
          GlobalOrigin[1],
          GlobalOrigin[2]);

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "GlobalRegion        = (%e, %e, %e)\n",
          GlobalRegion[0],
          GlobalRegion[1],
          GlobalRegion[2]);

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "GlobalVoxel         = (%d, %d, %d)\n",
          GlobalVoxel[0],
          GlobalVoxel[1],
          GlobalVoxel[2]);  

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "GlobalDivision      = (%d, %d, %d)\n",
          GlobalDivision[0],
          GlobalDivision[1],
          GlobalDivision[2]);

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp,"ActiveSubdomainFile = \"%s\"\n",ActiveSubdomainFile.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp,"CoordinateFile      = \"%s\"\n",CoordinateFile.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  if(        CoordinateFileType == CDM::E_CDM_FILE_TYPE_ASCII ) {
    fprintf(fp, "CoordinateFileType  = \"ascii\"\n");
  } else if( CoordinateFileType == CDM::E_CDM_FILE_TYPE_BINARY ) {
    fprintf(fp, "CoordinateFileType  = \"binary\"\n");
  }

  _CDM_WRITE_TAB(fp, tab+1);
  std::string Dtype = cdm_DFI::ConvDatatypeE2S((CDM::E_CDM_DTYPE)CoordinateFilePrecision);
  fprintf(fp, "CoordinateFilePrecision = \"%s\"\n",Dtype.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  if(        CoordinateFileEndian == CDM::E_CDM_LITTLE ) {
    fprintf(fp, "CoordinateFileEndian = \"little\"\n");
  } else if( CoordinateFileEndian == CDM::E_CDM_BIG ) {
    fprintf(fp, "CoordinateFileEndian = \"big\"\n");
  }

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}

#endif // _CDM_NONUNIFORMDOMAIN_INLINE_H_
