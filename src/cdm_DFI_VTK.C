/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_VTK.C
 * @brief  cdm_DFI_VTK Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include "cdm_DFI_VTK.h"

// #################################################################
// コンストラクタ
cdm_DFI_VTK::cdm_DFI_VTK()
{

}


// #################################################################
// デストラクタ
cdm_DFI_VTK::~cdm_DFI_VTK()
{

}

// #################################################################
// ヘッダーレコード出力
CDM::E_CDM_ERRORCODE
cdm_DFI_VTK::write_HeaderRecord(FILE* fp,
                                const unsigned step,
                                const double time,
                                const int n)
{

  if( !fp ) return CDM::E_CDM_ERROR;

  fprintf( fp, "# vtk DataFile Version 2.0\n" );
  fprintf( fp, "step=%d,time=%g\n", step, time );

  if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_BINARY ) {
    fprintf( fp, "BINARY\n" );
  } else if( m_output_type ==  CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
    fprintf( fp, "ASCII\n" );
  }

  fprintf( fp, "DATASET STRUCTURED_POINTS\n" );

  int imax,jmax,kmax;
  imax = (int)DFI_Process.RankList[n].VoxelSize[0]+(2*(int)DFI_Finfo.GuideCell);
  jmax = (int)DFI_Process.RankList[n].VoxelSize[1]+(2*(int)DFI_Finfo.GuideCell);
  kmax = (int)DFI_Process.RankList[n].VoxelSize[2]+(2*(int)DFI_Finfo.GuideCell);
  fprintf( fp, "DIMENSIONS %d %d %d\n", imax+1, jmax+1, kmax+1 );

  //double t_org[3];
  double t_pit[3];
  for(int i=0; i<3; i++ ) t_pit[i]=DFI_Domain->GlobalRegion[i]/
                                   (double)DFI_Domain->GlobalVoxel[i];
  //for(int i=0; i<3; i++ ) t_org[i]=DFI_Domain.GlobalOrigin[i]-(t_pit[0]*0.5);
  fprintf( fp, "ORIGIN %e %e %e\n",DFI_Domain->GlobalOrigin[0],
                                   DFI_Domain->GlobalOrigin[1],
                                   DFI_Domain->GlobalOrigin[2]);

  fprintf( fp, "ASPECT_RATIO %e %e %e\n", t_pit[0], t_pit[1], t_pit[2] );

  //int nw = imax*jmax*kmax;
  //fprintf( fp, "CELL_DATA %d\n", nw );
  int nw = (imax+1)*(jmax+1)*(kmax+1);
  fprintf( fp, "POINT_DATA %d\n", nw );

  std::string d_type;
  if(      DFI_Finfo.DataType == CDM::E_CDM_UINT8  ) d_type="unsigned_char";
  else if( DFI_Finfo.DataType == CDM::E_CDM_INT8   ) d_type="char";
  else if( DFI_Finfo.DataType == CDM::E_CDM_UINT16 ) d_type="unsigned_short";
  else if( DFI_Finfo.DataType == CDM::E_CDM_INT16  ) d_type="short";
  else if( DFI_Finfo.DataType == CDM::E_CDM_UINT32 ) d_type="unsigned_int";
  else if( DFI_Finfo.DataType == CDM::E_CDM_INT32  ) d_type="int";
  else if( DFI_Finfo.DataType == CDM::E_CDM_UINT64 ) d_type="unsigned_long";
  else if( DFI_Finfo.DataType == CDM::E_CDM_INT64  ) d_type="long";
  else if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32) d_type="float";
  else if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64) d_type="double";

  if( DFI_Finfo.Component == 1 )
  {
    fprintf( fp, "SCALARS %s %s\n", DFI_Finfo.Prefix.c_str(),d_type.c_str() );
    fprintf( fp, "LOOKUP_TABLE default\n" );
  }
  else if( DFI_Finfo.Component == 3 )
  {
    fprintf( fp, "VECTORS %s %s\n", DFI_Finfo.Prefix.c_str(),d_type.c_str() );
  }
  else
  {
    fprintf( fp, "FIELD %s 1\n", DFI_Finfo.Prefix.c_str() );
    fprintf( fp, "%s %d %d %s\n", DFI_Finfo.Prefix.c_str(), DFI_Finfo.Component, 
             nw, d_type.c_str() );
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// データレコード出力
CDM::E_CDM_ERRORCODE
cdm_DFI_VTK::write_DataRecord(FILE* fp, 
                              cdm_Array* val, 
                              const int gc, 
                              const int n)
{

  const int* sz = val->getArraySizeInt();
  size_t dLen = (size_t)sz[0]*(size_t)sz[1]*(size_t)sz[2]*val->getNcomp();

  if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_BINARY ) {

    //出力実数タイプがuint8のとき
    if( val->getDataType() == CDM::E_CDM_UINT8 ) {
      unsigned char *data = (unsigned char*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned char), dLen, fp );

    //出力実数タイプがint8のとき
    }else if( val->getDataType() == CDM::E_CDM_INT8 ) {
      char *data = (char*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(char), dLen, fp );

    //出力実数タイプがuint16のとき
    }else if( val->getDataType() == CDM::E_CDM_UINT16 ) {
      unsigned short *data = (unsigned short*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned short), dLen, fp );

    //出力実数タイプがint16のとき
    }else if( val->getDataType() == CDM::E_CDM_INT16 ) {
      short *data = (short*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(short), dLen, fp );

    //出力実数タイプがuint32のとき
    }else if( val->getDataType() == CDM::E_CDM_UINT32 ) {
      unsigned int *data = (unsigned int*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned int), dLen, fp );

    //出力実数タイプがint32のとき
    }else if( val->getDataType() == CDM::E_CDM_INT32 ) {
      int *data = (int*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(int), dLen, fp );

    //出力実数タイプがuint64のとき
    }else if( val->getDataType() == CDM::E_CDM_UINT64 ) {
      unsigned long long *data = (unsigned long long*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned long long), dLen, fp );

    //出力実数タイプがint64のとき
    }else if( val->getDataType() == CDM::E_CDM_INT64 ) {
      long long *data = (long long*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(long long), dLen, fp );

    //出力実数タイプがfloatのとき
    }else if( val->getDataType() == CDM::E_CDM_FLOAT32 ) {
      float *data = (float*)val->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(float), dLen, fp );

    //出力実数タイプがdoubleのとき
    }else if( val->getDataType() == CDM::E_CDM_FLOAT64 ) {
      double *data = (double*)val->getData();
      DBSWAPVEC(data,dLen);
      fwrite( data, sizeof(double), dLen, fp );
    }

    fprintf( fp, "\n" );
  } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {

    
    if( val->writeAscii(fp) != dLen ) {
      return CDM::E_CDM_ERROR;
    }
    fprintf( fp, "\n" );
  
  }
  return CDM::E_CDM_SUCCESS;
}

