#ifndef _CDM_NONUNIFORMDOMAIN_INLINE_H_
#define _CDM_NONUNIFORMDOMAIN_INLINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
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

    size_t nread=0;

    //x
    nread += fread(&szGrid[0], sizeof(int), 1, fp);
    //GlobalVoxelの確認
    if ( szGrid[0] != GlobalVoxel[0]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in X direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    if( XCoordinates ){ delete[] XCoordinates; }
    XCoordinates = new T[szGrid[0]];
    nread += fread(XCoordinates, sizeof(T), szGrid[0], fp);

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
    nread += fread(&szGrid[1], sizeof(int), 1, fp);
    //GlobalVoxelの確認
    if ( szGrid[1] != GlobalVoxel[1]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Y direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    if( YCoordinates ){ delete[] YCoordinates; }
    YCoordinates = new T[szGrid[1]];
    nread += fread(YCoordinates, sizeof(T), szGrid[1], fp);

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
    nread += fread(&szGrid[2], sizeof(int), 1, fp);
    //GlobalVoxelの確認
    if ( szGrid[2] != GlobalVoxel[2]+1 ) {
      printf("\tError in Read CoordinateFile: Number of grid in Z direction\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
    if( ZCoordinates ){ delete[] ZCoordinates; }
    ZCoordinates = new T[szGrid[2]];
    nread += fread(ZCoordinates, sizeof(T), szGrid[2], fp);

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

    //Check of number of coordinates
    if ( nread != (szGrid[0]+szGrid[1]+szGrid[2]+3) ) {
      printf("\tError in Read CoordinateFile: Num of coordinates\n");
      return CDM::E_CDM_ERROR_READ_COORDINATEFILE;
    }
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


  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}

#endif // _CDM_NONUNIFORMDOMAIN_INLINE_H_
