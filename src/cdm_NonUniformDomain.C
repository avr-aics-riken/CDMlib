/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_NonUniformDomain.C
 * @brief  cdm_NonUniformDomain Class
 * @author advancesoft
 */

#include "cdm_DFI.h"
#include "cdm_NonUniformDomain.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cdm_NonUniformDomain::cdm_NonUniformDomain()
: cdm_Domain(), XCoordinates(NULL), YCoordinates(NULL), ZCoordinates(NULL)
{
  CoordinateFile = "";
  CoordinateFileFormat = CDM::E_CDM_OUTPUT_TYPE_DEFAULT;
  CoordinateFilePrecision = CDM::E_CDM_DTYPE_UNKNOWN;
}

// #################################################################
// コンストラクタ
cdm_NonUniformDomain::cdm_NonUniformDomain(const double* _GlobalOrigin, 
                                           const double* _GlobalRegion, 
                                           const int* _GlobalVoxel,
                                           const int* _GlobalDivision,
                                           const int* _iblank,
                                           const double* _XCoordinates,
                                           const double* _YCoordinates,
                                           const double* _ZCoordinates)
: cdm_Domain(_GlobalOrigin,_GlobalRegion,_GlobalVoxel,_GlobalDivision,_iblank)
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
  CoordinateFile = "";
  CoordinateFileFormat = CDM::E_CDM_OUTPUT_TYPE_DEFAULT;
  CoordinateFilePrecision = CDM::E_CDM_DTYPE_UNKNOWN;
}

void cdm_NonUniformDomain::Clear()
{
  cdm_Domain::Clear();
  if( XCoordinates != NULL ){ delete[] XCoordinates; }
  if( YCoordinates != NULL ){ delete[] YCoordinates; }
  if( ZCoordinates != NULL ){ delete[] ZCoordinates; }
}

// #################################################################
// デストラクタ
cdm_NonUniformDomain::~cdm_NonUniformDomain()
{
  if( XCoordinates != NULL ){ delete[] XCoordinates; }
  if( YCoordinates != NULL ){ delete[] YCoordinates; }
  if( ZCoordinates != NULL ){ delete[] ZCoordinates; }
}


// #################################################################
// Domain の読込み関数
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain::Read(cdm_TextParser tpCntl) 
{

  std::string str;
  std::string label;
  double v[3];
  int iv[3];

  printf("\tRead of cdm_NonUniformDomain\n");

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
  // if ( CoordinateFile != "")
  if( !(fp=fopen(CoordinateFile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",CoordinateFile.c_str());
    return CDM::E_CDM_ERROR_OPEN_COORDINATEFILE;
  } else {
    printf("Open file. (%s)\n",CoordinateFile.c_str());
    Read_CoordinateFile(fp);
    //ret = Read_CoordinateFile(fp) としてエラーコード返す？関数Read_CoordinateFile内でエラー設定してから
    //if( ret != CDM::E_CDM_SUCCESS ){return CDM::E_CDM_ERROR_??}
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// Domain の読込み関数
CDM::E_CDM_ERRORCODE
cdm_NonUniformDomain::Read_CoordinateFile(FILE* fp)
{
  printf("\tRead CoordinateFile\n");

  int szGrid[3];

  //ascii
  if( CoordinateFileFormat == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
    std::cout << "CoordinateFileFormat is Ascii " << std::endl;

    //x
    fscanf(fp,"%d\n",&szGrid[0]);
    std::cout << "szGrid[0] " << szGrid[0] << std::endl;
    if ( szGrid[0] != GlobalVoxel[0]+1 ) {
      std::cout << "szGrid[0] is not equal to GlobalVoxel[0]+1 " << std::endl;
      //return CDM::E_CDM_ERROR_??
    }
    double *XCoordinates = new double[szGrid[0]];
    for(int i=0; i<szGrid[0]; i++) {
      fscanf(fp,"%lf\n",&(XCoordinates[i]));  //double
      std::cout << "XCoordinates for i= " << i << " " << XCoordinates[i] << std::endl;
      //実際の座標の数とszGrid[0]と一致しているかも確認必要か。
    }
    //if ( XCoordinates[0] != GlobalOrigin[0] ) {
    //if ( (XCoordinates[szGrid[0]-1]-XCoordinates[0]) != GlobalRegion[0] ) {
    //実数比較ってこれでOKなの?数値誤差は考慮できてる？

    //y
    fscanf(fp,"%d\n",&szGrid[1]);
    cout << "szGrid[1] " << szGrid[1] << endl;
    if ( szGrid[1] != GlobalVoxel[1]+1 ) {
      std::cout << "szGrid[1] is not equal to GlobalVoxel[1]+1 " << std::endl;
      //return CDM::E_CDM_ERROR_??
    }
    double *YCoordinates = new double[szGrid[1]];
    for(int j=0; j<szGrid[1]; j++) {
      fscanf(fp,"%lf\n",&(YCoordinates[j]));  //double
      std::cout << "YCoordinates for j= " << j << " " << YCoordinates[j] << std::endl;
    }

    //z
    fscanf(fp,"%d\n",&szGrid[2]);
    cout << "szGrid[2] " << szGrid[2] << endl;
    if ( szGrid[2] != GlobalVoxel[2]+1 ) {
      std::cout << "szGrid[2] is not equal to GlobalVoxel[2]+1 " << std::endl;
      //return CDM::E_CDM_ERROR_??
    }
    double *ZCoordinates = new double[szGrid[2]];
    for(int k=0; k<szGrid[2]; k++) {
      fscanf(fp,"%lf\n",&(ZCoordinates[k]));  //double
      std::cout << "ZCoordinates for k= " << k << " " << ZCoordinates[k] << std::endl;
    }

  //binary
  } else {
    std::cout << "CoordinatesFileFormat is Binary " << std::endl;

    //x
    fread(&szGrid[0], sizeof(int), 1, fp);
    cout << "szGrid[0] " << szGrid[0] << endl;
    if ( szGrid[0] != GlobalVoxel[0]+1 ) {
      std::cout << "szGrid[0] is not equal to GlobalVoxel[0]+1 " << std::endl;
      //return CDM::E_CDM_ERROR_??
    }
    double *XCoordinates = new double[szGrid[0]];
    fread(XCoordinates, sizeof(double), szGrid[0], fp);
    for(int i=0; i<szGrid[0]; i++) {
      cout << "XCoordinates for i= " << i << " " << XCoordinates[i] << endl;
    }

    //y
    fread(&szGrid[1], sizeof(int), 1, fp);
    cout << "szGrid[1] " << szGrid[1] << endl;
    double *YCoordinates = new double[szGrid[1]];
    if ( szGrid[1] != GlobalVoxel[1]+1 ) {
      std::cout << "szGrid[1] is not equal to GlobalVoxel[1]+1 " << std::endl;
      //return CDM::E_CDM_ERROR_??
    }
    fread(YCoordinates, sizeof(double), szGrid[1], fp);
    for(int j=0; j<szGrid[1]; j++) {
      cout << "YCoordinates for j= " << j << " " << YCoordinates[j] << endl;
    }

    //z
    fread(&szGrid[2], sizeof(int), 1, fp);
    cout << "szGrid[2] " << szGrid[2] << endl;
    if ( szGrid[2] != GlobalVoxel[2]+1 ) {
      std::cout << "szGrid[2] is not equal to GlobalVoxel[2]+1 " << std::endl;
      //return CDM::E_CDM_ERROR_??
    }
    double *ZCoordinates = new double[szGrid[2]];
    fread(ZCoordinates, sizeof(double), szGrid[2], fp);
    for(int k=0; k<szGrid[2]; k++) {
      cout << "ZCoordinates for k= " << k << " " << ZCoordinates[k] << endl;
    }

  }

  //さらに、GlobalOrigin等の内容が、Coordinateファイルの内容と一致するかの確認も必要。
  //GlobalVoxel(ボクセル数)はfread,scanfの時に確認している

  return CDM::E_CDM_SUCCESS;
}
