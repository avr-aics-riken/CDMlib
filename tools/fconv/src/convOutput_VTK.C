/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   convOutput_VTK.C
 * @brief  convOutput_VTK Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_VTK.h"

// #################################################################
// コンストラクタ
convOutput_VTK::convOutput_VTK()
{


}

// #################################################################
// デストラクタ
convOutput_VTK::~convOutput_VTK()
{


}

// #################################################################
// 出力ファイルをオープンする。
FILE* convOutput_VTK::OutputFile_Open(
                                      const std::string prefix,
                                      const unsigned step,
                                      const int id,
                                      const bool mio)
{
  FILE* fp;

  //ファイル名の生成
  std::string outfile;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir()+ "/"+
            cdm_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       "vtk",
                                       fnameformat,
                                       mio,
                                       CDM::E_CDM_OFF);

  //ファイルオープン
  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_BINARY ) {
    if( (fp = fopen(outfile.c_str(), "w")) == NULL ) {
      printf("\tCan't open file.(%s)\n",outfile.c_str());
      Exit(0);
    }
  } else if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {
    if( (fp = fopen(outfile.c_str(), "wa")) == NULL ) {
      printf("\tCan't open file.(%s)\n",outfile.c_str());
      Exit(0);
    }
  }

  return fp;

}

// #################################################################
//
bool convOutput_VTK::WriteHeaderRecord(
                                        int step, 
                                        int dim, 
                                        CDM::E_CDM_DTYPE d_type, 
                                        int imax, 
                                        int jmax, 
                                        int kmax,
                                        double time, 
                                        double* org, 
                                        double* pit, 
                                        std::string prefix, 
                                        FILE *fp)
{
  if( !fp ) return false;

  fprintf( fp, "# vtk DataFile Version 2.0\n" );
  fprintf( fp, "step=%d,time=%g\n", step, time );

  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_BINARY ) {
    fprintf( fp, "BINARY\n" );
  } else if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {
    fprintf( fp, "ASCII\n" );
  }

  fprintf( fp, "DATASET STRUCTURED_POINTS\n" );

  fprintf( fp, "DIMENSIONS %d %d %d\n", imax+1, jmax+1, kmax+1 );

  double t_org[3];
  t_org[0]=org[0]-(pit[0]*0.5);
  t_org[1]=org[1]-(pit[1]*0.5);
  t_org[2]=org[2]-(pit[2]*0.5);
  fprintf( fp, "ORIGIN %e %e %e\n", t_org[0], t_org[1], t_org[2] );

  fprintf( fp, "ASPECT_RATIO %e %e %e\n", pit[0], pit[1], pit[2] );

  //int nw = imax*jmax*kmax;
  //fprintf( fp, "CELL_DATA %d\n", nw );
  int nw = (imax+1)*(jmax+1)*(kmax+1);
  fprintf( fp, "POINT_DATA %d\n", nw );

  std::string out_d_type;
  if(      d_type == CDM::E_CDM_UINT8  ) out_d_type="unsigned_char";
  else if( d_type == CDM::E_CDM_INT8   ) out_d_type="char";
  else if( d_type == CDM::E_CDM_UINT16 ) out_d_type="unsigned_short";
  else if( d_type == CDM::E_CDM_INT16  ) out_d_type="short";
  else if( d_type == CDM::E_CDM_UINT32 ) out_d_type="unsigned_int";
  else if( d_type == CDM::E_CDM_INT32  ) out_d_type="int";
  else if( d_type == CDM::E_CDM_UINT64 ) out_d_type="unsigned_long";
  else if( d_type == CDM::E_CDM_INT64  ) out_d_type="long";
  else if( d_type == CDM::E_CDM_FLOAT32) out_d_type="float";
  else if( d_type == CDM::E_CDM_FLOAT64) out_d_type="double";
  
  if( dim == 1 )
  {
    fprintf( fp, "SCALARS %s %s\n", prefix.c_str(), out_d_type.c_str() );
    fprintf( fp, "LOOKUP_TABLE default\n" );
  }
  else if( dim == 3 )
  {
    fprintf( fp, "VECTORS %s %s\n", prefix.c_str(), out_d_type.c_str() );
  }
  else
  {
    fprintf( fp, "FIELD %s 1\n", prefix.c_str() );
    fprintf( fp, "%s %d %d %s\n", prefix.c_str(), dim, nw, out_d_type.c_str() );
  }

  return true;
}

// #################################################################
//
bool convOutput_VTK::WriteHeaderRecord(int step,
                                       double time,
                                       int dim, 
                                       CDM::E_CDM_DTYPE d_type, 
                                       CDM::E_CDM_DFITYPE dfi_type,
                                       cdm_Domain* out_domain,
                                       cdm_Process* out_process,
                                       int gc,
                                       std::string prefix,
                                       FILE *fp)
{
  if( !fp ) return false;

  fprintf( fp, "# vtk DataFile Version 2.0\n" );
  fprintf( fp, "step=%d,time=%g\n", step, time );

  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_BINARY ) {
    fprintf( fp, "BINARY\n" );
  } else if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {
    fprintf( fp, "ASCII\n" );
  }

  if( dfi_type == CDM::E_CDM_DFITYPE_CARTESIAN ) {
    fprintf( fp, "DATASET STRUCTURED_POINTS\n" );
  } else if( dfi_type == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
    fprintf( fp, "DATASET RECTILINEAR_GRID\n" );
  }

  int sz[3];
  int head[3];
  for(int i=0; i<3; i++) {
    sz[i] = out_process->RankList[0].VoxelSize[i];
    head[i] = out_process->RankList[0].HeadIndex[i];
  }

  int imax,jmax,kmax; 
  if( m_InputCntl->Get_Interp_flag() ) {
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

  if( dfi_type == CDM::E_CDM_DFITYPE_CARTESIAN ) {

    if( m_InputCntl->Get_Interp_flag() ) {
      fprintf( fp, "ORIGIN %e %e %e\n",out_domain->NodeX(head[0]-1),
                                       out_domain->NodeY(head[1]-1),
                                       out_domain->NodeZ(head[2]-1));
    } else {
      fprintf( fp, "ORIGIN %e %e %e\n",out_domain->CellX(head[0]-1-gc),
                                       out_domain->CellY(head[1]-1-gc),
                                       out_domain->CellZ(head[2]-1-gc));
    }

    double t_pit[3];
    for(int i=0; i<3; i++ ) t_pit[i]=out_domain->GlobalRegion[i]/
                                     (double)out_domain->GlobalVoxel[i];
    fprintf( fp, "ASPECT_RATIO %e %e %e\n", t_pit[0], t_pit[1], t_pit[2] );

  } else if( dfi_type == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {

    std::string d_type_coord;
    if( out_domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT32 ) {
      d_type_coord = "float";
    } else if( out_domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT64 ) {
      d_type_coord = "double";
    }

    //ascii
    if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {

      if( m_InputCntl->Get_Interp_flag() ) {
        //格子点補間する場合(ガイドセル出力はなし)
        //x
        fprintf( fp, "X_COORDINATES %d %s\n", imax, d_type_coord.c_str() );
        for (int i=0; i<sz[0]; i++ ) {
          fprintf( fp, "%e ", out_domain->NodeX(i+head[0]-1) );
        }
        fprintf( fp, "%e\n", out_domain->NodeX(sz[0]+head[0]-1) );

        //y
        fprintf( fp, "Y_COORDINATES %d %s\n", jmax, d_type_coord.c_str() );
        for (int j=0; j<sz[1]; j++ ) {
          fprintf( fp, "%e ", out_domain->NodeY(j+head[1]-1) );
        }
        fprintf( fp, "%e\n", out_domain->NodeY(sz[1]+head[1]-1) );

        //z
        fprintf( fp, "Z_COORDINATES %d %s\n", kmax, d_type_coord.c_str() );
        for (int k=0; k<sz[2]; k++ ) {
          fprintf( fp, "%e ", out_domain->NodeZ(k+head[2]-1) );
        }
        fprintf( fp, "%e\n", out_domain->NodeZ(sz[2]+head[2]-1) );

      } else {

        //x
        fprintf( fp, "X_COORDINATES %d %s\n", imax, d_type_coord.c_str() );
        for (int i=0; i<sz[0]+2*gc-1; i++ ) {
          fprintf( fp, "%e ", out_domain->CellX(i+head[0]-1-gc) );
        }
        fprintf( fp, "%e\n", out_domain->CellX(sz[0]+gc+head[0]-2) );

        //y
        fprintf( fp, "Y_COORDINATES %d %s\n", jmax, d_type_coord.c_str() );
        for (int j=0; j<sz[1]+2*gc-1; j++ ) {
          fprintf( fp, "%e ", out_domain->CellY(j+head[1]-1-gc) );
        }
        fprintf( fp, "%e\n", out_domain->CellY(sz[1]+gc+head[1]-2) );

        //z
        fprintf( fp, "Z_COORDINATES %d %s\n", kmax, d_type_coord.c_str() );
        for (int k=0; k<sz[2]+2*gc-1; k++ ) {
          fprintf( fp, "%e ", out_domain->CellZ(k+head[2]-1-gc) );
        }
        fprintf( fp, "%e\n", out_domain->CellZ(sz[2]+gc+head[2]-2) );

      }

    //binary
    } else {

      if( out_domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT32 ) {

        float *coord_X = NULL;
        float *coord_Y = NULL;
        float *coord_Z = NULL;
        coord_X = new float[imax];
        coord_Y = new float[jmax];
        coord_Z = new float[kmax];

        if( m_InputCntl->Get_Interp_flag() ) {
          //格子点補間する場合(ガイドセル出力はなし)
          for(int i=0; i<sz[0]+1; i++) coord_X[i] = (float)(out_domain->NodeX(i+head[0]-1));
          for(int j=0; j<sz[1]+1; j++) coord_Y[j] = (float)(out_domain->NodeY(j+head[1]-1));
          for(int k=0; k<sz[2]+1; k++) coord_Z[k] = (float)(out_domain->NodeZ(k+head[2]-1));
        } else {
          for(int i=0; i<sz[0]+2*gc; i++) coord_X[i] = (float)(out_domain->CellX(i+head[0]-1-gc));
          for(int j=0; j<sz[1]+2*gc; j++) coord_Y[j] = (float)(out_domain->CellY(j+head[1]-1-gc));
          for(int k=0; k<sz[2]+2*gc; k++) coord_Z[k] = (float)(out_domain->CellZ(k+head[2]-1-gc));
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

      } else if( out_domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT64 ) {

        double *coord_X = NULL;
        double *coord_Y = NULL;
        double *coord_Z = NULL;
        coord_X = new double[imax];
        coord_Y = new double[jmax];
        coord_Z = new double[kmax];

        if( m_InputCntl->Get_Interp_flag() ) {
          //格子点補間する場合(ガイドセル出力はなし)
          for(int i=0; i<sz[0]+1; i++) coord_X[i] = (double)(out_domain->NodeX(i+head[0]-1));
          for(int j=0; j<sz[1]+1; j++) coord_Y[j] = (double)(out_domain->NodeY(j+head[1]-1));
          for(int k=0; k<sz[2]+1; k++) coord_Z[k] = (double)(out_domain->NodeZ(k+head[2]-1));
        } else {
          for(int i=0; i<sz[0]+2*gc; i++) coord_X[i] = (double)(out_domain->CellX(i+head[0]-1-gc));
          for(int j=0; j<sz[1]+2*gc; j++) coord_Y[j] = (double)(out_domain->CellY(j+head[1]-1-gc));
          for(int k=0; k<sz[2]+2*gc; k++) coord_Z[k] = (double)(out_domain->CellZ(k+head[2]-1-gc));
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

  return true;
}

// #################################################################
//
bool convOutput_VTK::WriteFieldData(FILE* fp, cdm_Array* src, size_t dLen)
{

  const int* sz = src->getArraySizeInt();
  cdm_Array *out = cdm_Array::instanceArray
                   (src->getDataType(),
                    src->getArrayShape(),
                    (int *)sz,
                    0,
                    src->getNvari());

  int ret = src->copyArray(out);

  //バイナリー出力のとき
  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_BINARY ) {

    //出力タイプごとにポインターを取得して出力
    if( out->getDataType() == CDM::E_CDM_UINT8 ) {
      //UINT8
      unsigned char *data = (unsigned char*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned char), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT8 ) {
      //INT8
      char *data = (char*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(char), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_UINT16 ) {
      //UINT16
      unsigned short *data = (unsigned short*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned short), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT16 ) {
      //INT16
      short *data = (short*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(short), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_UINT32 ) {
      //UINT32
      unsigned int *data = (unsigned int*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned int), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT32 ) {
      //INT32
      int *data = (int*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(int), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_UINT64 ) {
      //UINT64
      unsigned long long *data = (unsigned long long*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned long long), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT64 ) {
      //INT64
      long long *data = (long long*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(long long), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_FLOAT32 ) {
      //FLOAT32
      float *data = (float*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(float), dLen, fp );

    } else if( out->getDataType() == CDM::E_CDM_FLOAT64 ) {
      //FLOAT 64
      double *data = (double*)out->getData();
      DBSWAPVEC(data,dLen);
      fwrite( data, sizeof(double), dLen, fp );

    } else {
      printf("\tIllegal datatype\n");
      delete out;
      Exit(0);
    }

  //アスキー出力のとき
  } else if ( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {

    if( out->writeAscii(fp) != dLen ) {
      delete out;
      Exit(0);
    }

  }

  delete out;
  return true;
}

// #################################################################
//
bool convOutput_VTK::WriteFieldData(FILE* fp,
                                    cdm_Array* src,
                                    size_t dLen,
                                    CDM::E_CDM_DTYPE d_type,
                                    bool flag_variname,
                                    std::string variname)
{

  if( flag_variname ) {
    std::string out_d_type;
    if(      d_type == CDM::E_CDM_UINT8  ) out_d_type="unsigned_char";
    else if( d_type == CDM::E_CDM_INT8   ) out_d_type="char";
    else if( d_type == CDM::E_CDM_UINT16 ) out_d_type="unsigned_short";
    else if( d_type == CDM::E_CDM_INT16  ) out_d_type="short";
    else if( d_type == CDM::E_CDM_UINT32 ) out_d_type="unsigned_int";
    else if( d_type == CDM::E_CDM_INT32  ) out_d_type="int";
    else if( d_type == CDM::E_CDM_UINT64 ) out_d_type="unsigned_long";
    else if( d_type == CDM::E_CDM_INT64  ) out_d_type="long";
    else if( d_type == CDM::E_CDM_FLOAT32) out_d_type="float";
    else if( d_type == CDM::E_CDM_FLOAT64) out_d_type="double";

    fprintf( fp, "\n" );
    fprintf( fp, "SCALARS %s %s\n", variname.c_str(), out_d_type.c_str() );
    fprintf( fp, "LOOKUP_TABLE default\n" );
  }

  const int* sz = src->getArraySizeInt();
  cdm_Array *out = cdm_Array::instanceArray
                   (src->getDataType(),
                    src->getArrayShape(),
                    (int *)sz,
                    0,
                    src->getNvari());

  int ret = src->copyArray(out);

  //バイナリー出力のとき
  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_BINARY ) {

    //出力タイプごとにポインターを取得して出力
    if( out->getDataType() == CDM::E_CDM_UINT8 ) {
      //UINT8
      unsigned char *data = (unsigned char*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned char), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT8 ) {
      //INT8
      char *data = (char*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(char), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_UINT16 ) {
      //UINT16
      unsigned short *data = (unsigned short*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned short), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT16 ) {
      //INT16
      short *data = (short*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(short), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_UINT32 ) {
      //UINT32
      unsigned int *data = (unsigned int*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned int), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT32 ) {
      //INT32
      int *data = (int*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(int), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_UINT64 ) {
      //UINT64
      unsigned long long *data = (unsigned long long*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned long long), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_INT64 ) {
      //INT64
      long long *data = (long long*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(long long), dLen, fp);

    } else if( out->getDataType() == CDM::E_CDM_FLOAT32 ) {
      //FLOAT32
      float *data = (float*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(float), dLen, fp );

    } else if( out->getDataType() == CDM::E_CDM_FLOAT64 ) {
      //FLOAT 64
      double *data = (double*)out->getData();
      DBSWAPVEC(data,dLen);
      fwrite( data, sizeof(double), dLen, fp );

    } else {
      printf("\tIllegal datatype\n");
      delete out;
      Exit(0);
    }

  //アスキー出力のとき
  } else if ( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {

    if( out->writeAscii(fp) != dLen ) {
      delete out;
      Exit(0);
    }

  }

  delete out;
  return true;
}

// #################################################################
//
bool convOutput_VTK::WriteDataMarker(int dmy, FILE* fp, bool out)
{
  if( !out ) return true;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  return true;
}

// #################################################################
//
void convOutput_VTK::OutputFile_Close(FILE* fp)
{
  fprintf( fp, "\n" );
  fclose(fp);
}
