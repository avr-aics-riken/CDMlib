/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_Domain.C
 * @brief  cdm_Domain Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cdm_Domain::cdm_Domain()
{
  for(int i=0; i<3; i++) GlobalOrigin[i]=0.0;
  for(int i=0; i<3; i++) GlobalRegion[i]=0.0;
  for(int i=0; i<3; i++) GlobalVoxel[i]=0;
  for(int i=0; i<3; i++) GlobalDivision[i]=0;
  for(int i=0; i<3; i++) Pitch[i]=0.0;
  ActiveSubdomainFile="";
}

// #################################################################
// コンストラクタ
cdm_Domain::cdm_Domain(const double* _GlobalOrigin,
                       const double* _GlobalPitch,
                       const int* _GlobalVoxel,
                       const int* _GlobalDivision)
{
  GlobalOrigin[0]=_GlobalOrigin[0];
  GlobalOrigin[1]=_GlobalOrigin[1];
  GlobalOrigin[2]=_GlobalOrigin[2];

  Pitch[0]=_GlobalPitch[0];
  Pitch[1]=_GlobalPitch[1];
  Pitch[2]=_GlobalPitch[2];

  GlobalVoxel[0]=_GlobalVoxel[0];
  GlobalVoxel[1]=_GlobalVoxel[1];
  GlobalVoxel[2]=_GlobalVoxel[2];

  GlobalDivision[0]=_GlobalDivision[0];
  GlobalDivision[1]=_GlobalDivision[1];
  GlobalDivision[2]=_GlobalDivision[2];

  GlobalRegion[0]=_GlobalPitch[0]*_GlobalVoxel[0];
  GlobalRegion[1]=_GlobalPitch[1]*_GlobalVoxel[1];
  GlobalRegion[2]=_GlobalPitch[2]*_GlobalVoxel[2];
}

cdm_Domain::cdm_Domain(const float* _GlobalOrigin,
                       const float* _GlobalPitch,
                       const int* _GlobalVoxel,
                       const int* _GlobalDivision)
{
  GlobalOrigin[0]=(double)_GlobalOrigin[0];
  GlobalOrigin[1]=(double)_GlobalOrigin[1];
  GlobalOrigin[2]=(double)_GlobalOrigin[2];

  Pitch[0]=(double)_GlobalPitch[0];
  Pitch[1]=(double)_GlobalPitch[1];
  Pitch[2]=(double)_GlobalPitch[2];

  GlobalVoxel[0]=_GlobalVoxel[0];
  GlobalVoxel[1]=_GlobalVoxel[1];
  GlobalVoxel[2]=_GlobalVoxel[2];

  GlobalDivision[0]=_GlobalDivision[0];
  GlobalDivision[1]=_GlobalDivision[1];
  GlobalDivision[2]=_GlobalDivision[2];

  GlobalRegion[0]=(double)_GlobalPitch[0]*_GlobalVoxel[0];
  GlobalRegion[1]=(double)_GlobalPitch[1]*_GlobalVoxel[1];
  GlobalRegion[2]=(double)_GlobalPitch[2]*_GlobalVoxel[2];
}

void cdm_Domain::Clear()
{
  for(int i=0; i<3; i++) GlobalOrigin[i]=0.0;
  for(int i=0; i<3; i++) GlobalRegion[i]=0.0;
  for(int i=0; i<3; i++) GlobalVoxel[i]=0;
  for(int i=0; i<3; i++) GlobalDivision[i]=0;
  for(int i=0; i<3; i++) Pitch[i]=0.0;
  ActiveSubdomainFile="";
}

// #################################################################
// デストラクタ
cdm_Domain::~cdm_Domain()
{
}

// #################################################################
// Domain の読込み関数
CDM::E_CDM_ERRORCODE
cdm_Domain::Read(cdm_TextParser tpCntl,
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

  //Pitch
  Pitch[0] = GlobalRegion[0] / GlobalVoxel[0];
  Pitch[1] = GlobalRegion[1] / GlobalVoxel[1];
  Pitch[2] = GlobalRegion[2] / GlobalVoxel[2];

  //ActiveSubdomain
  label = "/Domain/ActiveSubdomainFile";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    str="";
  }
  ActiveSubdomainFile=str;

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// DFIファイル:Domain要素を出力する
CDM::E_CDM_ERRORCODE
cdm_Domain::Write(FILE* fp, 
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

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}
