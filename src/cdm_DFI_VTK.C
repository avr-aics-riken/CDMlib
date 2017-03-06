/*
###################################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
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
//cdm_DFI_VTK::write_HeaderRecord(FILE* fp,
cdm_DFI_VTK::write_HeaderRecord(cdm_FILE* pFile,
                                const unsigned step,
                                const double time,
                                const int n)
{
  FILE *fp = pFile->m_fp;

  if( !fp ) return CDM::E_CDM_ERROR;

  fprintf( fp, "# vtk DataFile Version 2.0\n" );
  fprintf( fp, "step=%d,time=%g\n", step, time );

  if( m_output_type ==  CDM::E_CDM_FILE_TYPE_ASCII ) {
    fprintf( fp, "ASCII\n" );
  } else {
    fprintf( fp, "BINARY\n" );
  }

  if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {
    fprintf( fp, "DATASET STRUCTURED_POINTS\n" );
  } else if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
    fprintf( fp, "DATASET RECTILINEAR_GRID\n" );
  }

  int sz[3];
  int head[3];
  for(int i=0; i<3; i++) {
    sz[i] = DFI_Process.RankList[n].VoxelSize[i];
    head[i] = DFI_Process.RankList[n].HeadIndex[i];
  }

  int gc = DFI_Finfo.GuideCell;

  int imax,jmax,kmax;
  if( m_bgrid_interp_flag ) {
    //格子点補間処理ありの場合は、配列サイズを+1。(ガイドセル出力はなし)
    imax = sz[0]+1;
    jmax = sz[1]+1;
    kmax = sz[2]+1;
  } else {
    //ガイドセルも含めた格子点の数
    imax = sz[0]+2*gc;
    jmax = sz[1]+2*gc;
    kmax = sz[2]+2*gc;
  }
  fprintf( fp, "DIMENSIONS %d %d %d\n", imax, jmax, kmax );

  if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {

    if( m_bgrid_interp_flag ) {
      fprintf( fp, "ORIGIN %e %e %e\n",DFI_Domain->NodeX(head[0]-1),
                                       DFI_Domain->NodeY(head[1]-1),
                                       DFI_Domain->NodeZ(head[2]-1));
    } else {
      fprintf( fp, "ORIGIN %e %e %e\n",DFI_Domain->CellX(head[0]-1-gc),
                                       DFI_Domain->CellY(head[1]-1-gc),
                                       DFI_Domain->CellZ(head[2]-1-gc));
    }

    double t_pit[3];
    for(int i=0; i<3; i++ ) t_pit[i]=DFI_Domain->GlobalRegion[i]/
                                     (double)DFI_Domain->GlobalVoxel[i];
    fprintf( fp, "ASPECT_RATIO %e %e %e\n", t_pit[0], t_pit[1], t_pit[2] );

  } else if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {

    std::string d_type_coord;
    if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT32 ) {
      d_type_coord = "float";
    } else if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT64 ) {
      d_type_coord = "double";
    }

    //ascii
    if( m_output_type == CDM::E_CDM_FILE_TYPE_ASCII ) {

      if( m_bgrid_interp_flag ) {
        //格子点補間する場合(ガイドセル出力はなし)
        //x
        fprintf( fp, "X_COORDINATES %d %s\n", imax, d_type_coord.c_str() );
        for (int i=0; i<sz[0]; i++ ) {
          fprintf( fp, "%e ", DFI_Domain->NodeX(i+head[0]-1) );
        }
        fprintf( fp, "%e\n", DFI_Domain->NodeX(sz[0]+head[0]-1) );

        //y
        fprintf( fp, "Y_COORDINATES %d %s\n", jmax, d_type_coord.c_str() );
        for (int j=0; j<sz[1]; j++ ) {
          fprintf( fp, "%e ", DFI_Domain->NodeY(j+head[1]-1) );
        }
        fprintf( fp, "%e\n", DFI_Domain->NodeY(sz[1]+head[1]-1) );

        //z
        fprintf( fp, "Z_COORDINATES %d %s\n", kmax, d_type_coord.c_str() );
        for (int k=0; k<sz[2]; k++ ) {
          fprintf( fp, "%e ", DFI_Domain->NodeZ(k+head[2]-1) );
        }
        fprintf( fp, "%e\n", DFI_Domain->NodeZ(sz[2]+head[2]-1) );

      } else {

        //x
        fprintf( fp, "X_COORDINATES %d %s\n", imax, d_type_coord.c_str() );
        for (int i=0; i<sz[0]+2*gc-1; i++ ) {
          fprintf( fp, "%e ", DFI_Domain->CellX(i+head[0]-1-gc) );
        }
        fprintf( fp, "%e\n", DFI_Domain->CellX(sz[0]+gc+head[0]-2) );

        //y
        fprintf( fp, "Y_COORDINATES %d %s\n", jmax, d_type_coord.c_str() );
        for (int j=0; j<sz[1]+2*gc-1; j++ ) {
          fprintf( fp, "%e ", DFI_Domain->CellY(j+head[1]-1-gc) );
        }
        fprintf( fp, "%e\n", DFI_Domain->CellY(sz[1]+gc+head[1]-2) );

        //z
        fprintf( fp, "Z_COORDINATES %d %s\n", kmax, d_type_coord.c_str() );
        for (int k=0; k<sz[2]+2*gc-1; k++ ) {
          fprintf( fp, "%e ", DFI_Domain->CellZ(k+head[2]-1-gc) );
        }
        fprintf( fp, "%e\n", DFI_Domain->CellZ(sz[2]+gc+head[2]-2) );

      }

    //binary
    } else {

      if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT32 ) {

        float *coord_X = NULL;
        float *coord_Y = NULL;
        float *coord_Z = NULL;
        coord_X = new float[imax];
        coord_Y = new float[jmax];
        coord_Z = new float[kmax];

        if( m_bgrid_interp_flag ) {
          //格子点補間する場合(ガイドセル出力はなし)
          for(int i=0; i<sz[0]+1; i++) coord_X[i] = (float)(DFI_Domain->NodeX(i+head[0]-1));
          for(int j=0; j<sz[1]+1; j++) coord_Y[j] = (float)(DFI_Domain->NodeY(j+head[1]-1));
          for(int k=0; k<sz[2]+1; k++) coord_Z[k] = (float)(DFI_Domain->NodeZ(k+head[2]-1));
        } else {
          for(int i=0; i<sz[0]+2*gc; i++) coord_X[i] = (float)(DFI_Domain->CellX(i+head[0]-1-gc));
          for(int j=0; j<sz[1]+2*gc; j++) coord_Y[j] = (float)(DFI_Domain->CellY(j+head[1]-1-gc));
          for(int k=0; k<sz[2]+2*gc; k++) coord_Z[k] = (float)(DFI_Domain->CellZ(k+head[2]-1-gc));
        }

        //x
        fprintf( fp, "X_COORDINATES %d %s\n", imax, d_type_coord.c_str() );
        BSWAPVEC(coord_X,imax);
        fwrite(coord_X, sizeof(float), imax, fp);
        fprintf( fp, "\n" );

        //y
        fprintf( fp, "Y_COORDINATES %d %s\n", jmax, d_type_coord.c_str() );
        BSWAPVEC(coord_Y,jmax);
        fwrite(coord_Y, sizeof(float), jmax, fp);
        fprintf( fp, "\n" );

        //z
        fprintf( fp, "Z_COORDINATES %d %s\n", kmax, d_type_coord.c_str() );
        BSWAPVEC(coord_Z,kmax);
        fwrite(coord_Z, sizeof(float), kmax, fp);
        fprintf( fp, "\n" );

        delete [] coord_X;
        delete [] coord_Y;
        delete [] coord_Z;

      } else if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT64 ) {

        double *coord_X = NULL;
        double *coord_Y = NULL;
        double *coord_Z = NULL;
        coord_X = new double[imax];
        coord_Y = new double[jmax];
        coord_Z = new double[kmax];

        if( m_bgrid_interp_flag ) {
          //格子点補間する場合(ガイドセル出力はなし)
          for(int i=0; i<sz[0]+1; i++) coord_X[i] = (double)(DFI_Domain->NodeX(i+head[0]-1));
          for(int j=0; j<sz[1]+1; j++) coord_Y[j] = (double)(DFI_Domain->NodeY(j+head[1]-1));
          for(int k=0; k<sz[2]+1; k++) coord_Z[k] = (double)(DFI_Domain->NodeZ(k+head[2]-1));
        } else {
          for(int i=0; i<sz[0]+2*gc; i++) coord_X[i] = (double)(DFI_Domain->CellX(i+head[0]-1-gc));
          for(int j=0; j<sz[1]+2*gc; j++) coord_Y[j] = (double)(DFI_Domain->CellY(j+head[1]-1-gc));
          for(int k=0; k<sz[2]+2*gc; k++) coord_Z[k] = (double)(DFI_Domain->CellZ(k+head[2]-1-gc));
        }

        //x
        fprintf( fp, "X_COORDINATES %d %s\n", imax, d_type_coord.c_str() );
        DBSWAPVEC(coord_X,imax);
        fwrite(coord_X, sizeof(double), imax, fp);
        fprintf( fp, "\n" );

        //y
        fprintf( fp, "Y_COORDINATES %d %s\n", jmax, d_type_coord.c_str() );
        DBSWAPVEC(coord_Y,jmax);
        fwrite(coord_Y, sizeof(double), jmax, fp);
        fprintf( fp, "\n" );

        //z
        fprintf( fp, "Z_COORDINATES %d %s\n", kmax, d_type_coord.c_str() );
        DBSWAPVEC(coord_Z,kmax);
        fwrite(coord_Z, sizeof(double), kmax, fp);
        fprintf( fp, "\n" );

        delete [] coord_X;
        delete [] coord_Y;
        delete [] coord_Z;

      }

    }

  }

  int nw = imax*jmax*kmax;
  fprintf( fp, "POINT_DATA %d\n", nw );

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// データレコード出力
CDM::E_CDM_ERRORCODE
//cdm_DFI_VTK::write_DataRecord(FILE* fp,
cdm_DFI_VTK::write_DataRecord(cdm_FILE* pFile,
                              cdm_Array* val,
                              const int gc,
                              const int n)
{
  FILE *fp = pFile->m_fp;

  const int* sz = val->getArraySizeInt();
  size_t dLen;
  if( !m_bgrid_interp_flag ) {
    dLen = (size_t)(sz[0]+2*gc)*(size_t)(sz[1]+2*gc)*(size_t)(sz[2]+2*gc);
  } else {
    dLen = (size_t)sz[0]*(size_t)sz[1]*(size_t)sz[2];
  }

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

  cdm_Array *out = cdm_Array::instanceArray
                   (val->getDataType(),
                    CDM::E_CDM_IJKN,
                    (int *)sz,
                    val->getGc(),
                    1); //変数毎に出力するので、変数の個数は1

  for(int nv=0; nv<DFI_Finfo.NumVariables; nv++) {

    fprintf( fp, "SCALARS %s %s\n", getVariableName(nv).c_str(),d_type.c_str() );
    fprintf( fp, "LOOKUP_TABLE default\n" );

    if( val->copyArrayNvari_to_ijk(out,nv) != 0 ) {
      printf("\tError : copyArrayNvari_to_ijk in class cdm_DFI_VTK\n");
      return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;
    }

    //ascii
    if( m_output_type == CDM::E_CDM_FILE_TYPE_ASCII ) {

      if( out->writeAscii(fp) != dLen ) {
        return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;
      }
      fprintf( fp, "\n" );

    //binary
    } else {

      //出力実数タイプがuint8のとき
      if( out->getDataType() == CDM::E_CDM_UINT8 ) {
        unsigned char *data = (unsigned char*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(unsigned char), dLen, fp );

      //出力実数タイプがint8のとき
      }else if( out->getDataType() == CDM::E_CDM_INT8 ) {
        char *data = (char*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(char), dLen, fp );

      //出力実数タイプがuint16のとき
      }else if( out->getDataType() == CDM::E_CDM_UINT16 ) {
        unsigned short *data = (unsigned short*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(unsigned short), dLen, fp );

      //出力実数タイプがint16のとき
      }else if( out->getDataType() == CDM::E_CDM_INT16 ) {
        short *data = (short*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(short), dLen, fp );

      //出力実数タイプがuint32のとき
      }else if( out->getDataType() == CDM::E_CDM_UINT32 ) {
        unsigned int *data = (unsigned int*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(unsigned int), dLen, fp );

      //出力実数タイプがint32のとき
      }else if( out->getDataType() == CDM::E_CDM_INT32 ) {
        int *data = (int*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(int), dLen, fp );

      //出力実数タイプがuint64のとき
      }else if( out->getDataType() == CDM::E_CDM_UINT64 ) {
        unsigned long long *data = (unsigned long long*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(unsigned long long), dLen, fp );

      //出力実数タイプがint64のとき
      }else if( out->getDataType() == CDM::E_CDM_INT64 ) {
        long long *data = (long long*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(long long), dLen, fp );

      //出力実数タイプがfloatのとき
      }else if( out->getDataType() == CDM::E_CDM_FLOAT32 ) {
        float *data = (float*)out->getData();
        BSWAPVEC(data,dLen);
        fwrite( data, sizeof(float), dLen, fp );

      //出力実数タイプがdoubleのとき
      }else if( out->getDataType() == CDM::E_CDM_FLOAT64 ) {
        double *data = (double*)out->getData();
        DBSWAPVEC(data,dLen);
        fwrite( data, sizeof(double), dLen, fp );
      }

      fprintf( fp, "\n" );

    }

  }
  return CDM::E_CDM_SUCCESS;
}
