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
  CIO::E_CIO_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir()+ "/"+
            cio_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       "vtk",
                                       fnameformat,
                                       mio,
                                       CIO::E_CIO_OFF);

  //ファイルオープン
  if( m_InputCntl->Get_OutputFormatType() == CIO::E_CIO_OUTPUT_TYPE_BINARY ) {
    if( (fp = fopen(outfile.c_str(), "w")) == NULL ) {
      printf("\tCan't open file.(%s)\n",outfile.c_str());
      Exit(0);
    }
  } else if( m_InputCntl->Get_OutputFormatType() == CIO::E_CIO_OUTPUT_TYPE_ASCII ) {
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
                                        CIO::E_CIO_DTYPE d_type, 
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

  if( m_InputCntl->Get_OutputFormatType() == CIO::E_CIO_OUTPUT_TYPE_BINARY ) {
    fprintf( fp, "BINARY\n" );
  } else if( m_InputCntl->Get_OutputFormatType() == CIO::E_CIO_OUTPUT_TYPE_ASCII ) {
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
  if(      d_type == CIO::E_CIO_UINT8  ) out_d_type="unsigned_char";
  else if( d_type == CIO::E_CIO_INT8   ) out_d_type="char";
  else if( d_type == CIO::E_CIO_UINT16 ) out_d_type="unsigned_short";
  else if( d_type == CIO::E_CIO_INT16  ) out_d_type="short";
  else if( d_type == CIO::E_CIO_UINT32 ) out_d_type="unsigned_int";
  else if( d_type == CIO::E_CIO_INT32  ) out_d_type="int";
  else if( d_type == CIO::E_CIO_UINT64 ) out_d_type="unsigned_long";
  else if( d_type == CIO::E_CIO_INT64  ) out_d_type="long";
  else if( d_type == CIO::E_CIO_FLOAT32) out_d_type="float";
  else if( d_type == CIO::E_CIO_FLOAT64) out_d_type="double";
  
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
bool convOutput_VTK::WriteFieldData(FILE* fp, cio_Array* src, size_t dLen)
{

  const int* sz = src->getArraySizeInt();
  cio_Array *out = cio_Array::instanceArray
                   (src->getDataType(),
                    src->getArrayShape(),
                    (int *)sz,
                    0,
                    src->getNcomp());

  int ret = src->copyArray(out);

  //バイナリー出力のとき
  if( m_InputCntl->Get_OutputFormatType() == CIO::E_CIO_OUTPUT_TYPE_BINARY ) {

    //出力タイプごとにポインターを取得して出力
    if( out->getDataType() == CIO::E_CIO_UINT8 ) {
      //UINT8
      unsigned char *data = (unsigned char*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned char), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_INT8 ) {
      //INT8
      char *data = (char*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(char), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_UINT16 ) {
      //UINT16
      unsigned short *data = (unsigned short*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned short), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_INT16 ) {
      //INT16
      short *data = (short*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(short), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_UINT32 ) {
      //UINT32
      unsigned int *data = (unsigned int*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned int), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_INT32 ) {
      //INT32
      int *data = (int*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(int), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_UINT64 ) {
      //UINT64
      unsigned long long *data = (unsigned long long*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(unsigned long long), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_INT64 ) {
      //INT64
      long long *data = (long long*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(long long), dLen, fp);

    } else if( out->getDataType() == CIO::E_CIO_FLOAT32 ) {
      //FLOAT32
      float *data = (float*)out->getData();
      BSWAPVEC(data,dLen);
      fwrite( data, sizeof(float), dLen, fp );

    } else if( out->getDataType() == CIO::E_CIO_FLOAT64 ) {
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
  } else if ( m_InputCntl->Get_OutputFormatType() == CIO::E_CIO_OUTPUT_TYPE_ASCII ) {

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
