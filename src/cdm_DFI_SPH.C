/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_SPH.C
 * @brief  cdm_DFI_SPH Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include "cdm_DFI_SPH.h"

// #################################################################
// コンストラクタ
cdm_DFI_SPH::cdm_DFI_SPH()
{

}


// #################################################################
// デストラクタ
cdm_DFI_SPH::~cdm_DFI_SPH()
{

}

// #################################################################
// ファイルのヘッダーレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_SPH::read_HeaderRecord(FILE* fp, 
                               bool matchEndian, 
                               unsigned step,
                               const int head[3],
                               const int tail[3],
                               int gc, 
                               int voxsize[3],  
                               double &time) 
{

  unsigned int dmy,type_dmy;

//REC1  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != 8 ) {
    BSWAP32(dmy);
    if( dmy != 8 ) { 
      fclose(fp); 
      return CDM::E_CDM_ERROR_READ_SPH_REC1; 
    } else {
      fclose(fp);
      return CDM::E_CDM_ERROR_NOMATCH_ENDIAN;
    }
  }

  DataDims data_dims;
  if( fread(&data_dims, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; }
  if( !matchEndian ) BSWAP32(data_dims); 
  if( data_dims == _SCALAR && DFI_Finfo.NumVariables != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; } 
  if( data_dims == _VECTOR && DFI_Finfo.NumVariables <= 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; } 

  int real_type;

  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; }
  if( !matchEndian ) BSWAP32(real_type); 
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != 8 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC1; }

  if( real_type == _FLOAT ) {
    if( DFI_Finfo.DataType != CDM::E_CDM_FLOAT32 ) return CDM::E_CDM_ERROR_READ_SPH_REC1;
    type_dmy=12;
  } else if( real_type == _DOUBLE) {
    if( DFI_Finfo.DataType != CDM::E_CDM_FLOAT64 ) return CDM::E_CDM_ERROR_READ_SPH_REC1;
    type_dmy=24;
  }

//REC2
//ボクセルサイズ
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC2; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC2; }
  if( real_type == _FLOAT ) {
    if( fread(voxsize, sizeof(int), 3, fp) != 3 ){fclose(fp);return CDM::E_CDM_ERROR_READ_SPH_REC2;}
    if( !matchEndian ) {
      BSWAP32(voxsize[0]); 
      BSWAP32(voxsize[1]); 
      BSWAP32(voxsize[2]);
    }
  } else if( real_type == _DOUBLE ) {
    long long tmp[3];
    if( fread(tmp, sizeof(long long), 3, fp) != 3 ){fclose(fp);return CDM::E_CDM_ERROR_READ_SPH_REC2;}
    if( !matchEndian ) {
      BSWAP64(tmp[0]); 
      BSWAP64(tmp[1]); 
      BSWAP64(tmp[2]);
    }
    voxsize[0]=(int)tmp[0];
    voxsize[1]=(int)tmp[1];
    voxsize[2]=(int)tmp[2];
  } 
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC2; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_FILE; }

//REC3
//原点座標
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC3; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC3; }
  if( real_type == _FLOAT ) {
    float voxorg[3];
    if( fread(voxorg, sizeof(float), 3, fp) != 3 ){fclose(fp);return CDM::E_CDM_ERROR_READ_SPH_REC3;}
    if( !matchEndian ) {
      BSWAP32(voxorg[0]); 
      BSWAP32(voxorg[1]); 
      BSWAP32(voxorg[2]); 
    }
  } else if( real_type == _DOUBLE ) {
    double voxorg[3];
    if( fread(voxorg, sizeof(double), 3, fp) != 3 ){fclose(fp);return CDM::E_CDM_ERROR_READ_SPH_REC3;}
    if( !matchEndian ) {
      BSWAP64(voxorg[0]); 
      BSWAP64(voxorg[1]); 
      BSWAP64(voxorg[2]); 
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC3; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC3; }

//REC4  
//pit
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC4; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC4; }
  if( real_type == _FLOAT ) {
    float voxpit[3];
    if( fread(voxpit, sizeof(float), 3, fp) != 3 ){fclose(fp);return CDM::E_CDM_ERROR_READ_SPH_REC4;}
    if( !matchEndian ) {
      BSWAP32(voxpit[0]); 
      BSWAP32(voxpit[1]); 
      BSWAP32(voxpit[2]); 
    }
  } else if( real_type == _DOUBLE ) {
    double voxpit[3];
    if( fread(voxpit, sizeof(double), 3, fp) != 3 ){fclose(fp);return CDM::E_CDM_ERROR_READ_SPH_REC4;}
    if( !matchEndian ) {
      BSWAP64(voxpit[0]); 
      BSWAP64(voxpit[1]); 
      BSWAP64(voxpit[2]); 
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC4; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_FILE; }

//REC5
//step,time
  if( real_type == _FLOAT ) type_dmy = 8;
  if( real_type == _DOUBLE) type_dmy = 16;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
  if( real_type == _FLOAT ) { 
    int r_step;
    if( fread(&r_step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
    if( !matchEndian ) BSWAP32(r_step);
    if( r_step != step ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
  } else if( real_type == _DOUBLE ) {
    long long r_step;
    if( fread(&r_step, sizeof(long long), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
    if( !matchEndian ) BSWAP64(r_step);
    if( r_step != step ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
  }
  if( real_type == _FLOAT ) {
    float r_time;
    if( fread(&r_time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
    if( !matchEndian ) BSWAP32(r_time); 
    time = r_time;
  } else if( real_type == _DOUBLE ) {
    double r_time;
    if( fread(&r_time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
    if( !matchEndian ) BSWAP64(r_time); 
    time = r_time;
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }
  if( !matchEndian ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC5; }

  for(int i=0; i<3; i++) {
    if( voxsize[i] != (tail[i]-head[i]+1+2*gc) ) return CDM::E_CDM_ERROR_UNMATCH_VOXELSIZE;
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// ファイルのデーターレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_SPH::read_Datarecord(FILE* fp,
                             bool matchEndian, 
                             cdm_Array* buf, 
                             int head[3],
                             int nz,
                             cdm_Array* &src)
{

  //１層ずつ読み込み
  int hzB = head[2];

  // fortran record の読込み
  int idmy;
  if( fread(&idmy,sizeof(int),1,fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC6; }
  if( !matchEndian ) BSWAP32(idmy); 

  for( int k=0; k<nz; k++ ) {
    //headインデクスをずらす
    head[2]=hzB+k;
    buf->setHeadIndex(head);

    //１層読み込み
    size_t ndata = buf->getArrayLength();
    if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_SPH_REC6;

    // コピー
    buf->copyArray(src);
  }

  // fortran record の読込み
  if( fread(&idmy,sizeof(int),1,fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC6; }
  if( !matchEndian ) BSWAP32(idmy); 

  return CDM::E_CDM_SUCCESS;

}


// #################################################################
// Averaged レコードの読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_SPH::read_averaged(FILE* fp,
                           bool matchEndian,
                           unsigned dummy,
                           unsigned &step_avr,
                           double &time_avr)
{

  unsigned int dmy,type_dmy;

  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) type_dmy = 8;
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) type_dmy = 16;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
  if( !matchEndian ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) {
    int r_step;
    if( fread(&r_step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
    if( !matchEndian ) BSWAP32(r_step);
    step_avr=(unsigned)r_step;
  } else if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) {
    long long r_step;
    if( fread(&r_step, sizeof(long long), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
    if( !matchEndian ) BSWAP64(r_step);
    step_avr=(unsigned)r_step;
  }
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) { 
    float r_time;
    if( fread(&r_time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
    if( !matchEndian ) BSWAP32(r_time);
    time_avr = (double)r_time;
  } else if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) {
    double r_time;
    if( fread(&r_time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
    if( !matchEndian ) BSWAP64(r_time);
    time_avr = r_time;
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }
  if( !matchEndian ) BSWAP32(dmy);
  if( dmy != type_dmy ) { fclose(fp); return CDM::E_CDM_ERROR_READ_SPH_REC7; }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// SPHヘッダーレコードの出力
CDM::E_CDM_ERRORCODE 
cdm_DFI_SPH::write_HeaderRecord(FILE* fp, 
                                     const unsigned step, 
                                     const double time, 
                                     const int n)
{

  //REC1
  int svType = 0;
  if( DFI_Finfo.NumVariables == 1 ) svType = 1;
  if( DFI_Finfo.NumVariables > 1  ) svType = 2;
  if( svType == 0 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC1;

  int dType = 0;
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) dType = 1; 
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) dType = 2; 
  if( dType == 0 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC1;

  unsigned int dmy;
  dmy = 8;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC1;
  if( fwrite(&svType, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC1;
  if( fwrite(&dType, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC1;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC1;


  if( dType == 1 ) dmy = 12; //float
  else             dmy = 24; //double

  //REC2
  //voxel size
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC2;
  if( dType == 1 ) {
    int size[3];
    for(int i=0; i<3; i++ ) size[i] = (int)DFI_Process.RankList[n].VoxelSize[i]+(int)(2*DFI_Finfo.GuideCell);
    if( fwrite(size, sizeof(int), 3, fp) !=3 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC2;
  } else {
    long long  size[3];
    for(int i=0; i<3; i++ ) size[i] = (long long)DFI_Process.RankList[n].VoxelSize[i]+(long long)(2*DFI_Finfo.GuideCell);
    if( fwrite(size, sizeof(long long), 3, fp) !=3 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC2;
  }

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC2;

  //REC3
  //origin
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC3;
  if( dType == 1 ) {
    float pch[3];
    for(int i=0; i<3; i++ ) pch[i]=(float)DFI_Domain->GlobalRegion[i]/DFI_Domain->GlobalVoxel[i];
    float org[3];
    //for(int i=0; i<3; i++ ) org[i]=(float)DFI_Domain.GlobalOrigin[i];
    for(int i=0; i<3; i++ ) org[i]=(float)DFI_Domain->GlobalOrigin[i]+0.5*pch[i];
    if( DFI_Finfo.GuideCell>1 ) for(int i=0; i<3; i++) org[i]=org[i]-pch[i]*(float)DFI_Finfo.GuideCell;
    if( fwrite(org, sizeof(float), 3, fp) !=3 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC3;
  } else {
    double pch[3];
    for(int i=0; i<3; i++ ) pch[i]=(double)DFI_Domain->GlobalRegion[i]/DFI_Domain->GlobalVoxel[i];
    double org[3];
    for(int i=0; i<3; i++ ) org[i]=(double)DFI_Domain->GlobalOrigin[i]+0.5*pch[i];
    if( DFI_Finfo.GuideCell>1 ) for(int i=0; i<3; i++) org[i]=org[i]-pch[i]*(double)DFI_Finfo.GuideCell;
    if( fwrite(org, sizeof(double), 3, fp) !=3 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC3;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC3;

  //REC4
  //pitch
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC4;
  if( dType == 1 ) {
    float pch[3];
    for(int i=0; i<3; i++ ) pch[i]=(float)DFI_Domain->GlobalRegion[i]/DFI_Domain->GlobalVoxel[i];
    if( fwrite(pch, sizeof(float), 3, fp) !=3 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC4;
  } else {
    double pch[3];
    for(int i=0; i<3; i++ ) pch[i]=(double)DFI_Domain->GlobalRegion[i]/DFI_Domain->GlobalVoxel[i];
    if( fwrite(pch, sizeof(double), 3, fp) !=3 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC4;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC4;

  //REC5
  //step&time
  int Int_size,Real_size;
  if ( dType == 1 ) {
    dmy = 8;
    Int_size  = sizeof(int);
    Real_size = sizeof(float);
  }else{
    dmy = 16;
    Int_size  = sizeof(long long);
    Real_size = sizeof(double);
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC5;
  if( dType == 1 ){
    float ttime = (float)time;
    if( fwrite(&step, Int_size, 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC5;
    if( fwrite(&ttime, Real_size, 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC5;
  } else {
    long long dstep = (long long)step;
    double ttime = (double)time;
    if( fwrite(&dstep, Int_size, 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC5;
    if( fwrite(&ttime, Real_size, 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC5;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC5;
  
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// SPHデータレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_SPH::write_DataRecord(FILE* fp, 
                              cdm_Array* val, 
                              const int gc, 
                              const int n)
{

  CDM::E_CDM_DTYPE Dtype = (CDM::E_CDM_DTYPE)DFI_Finfo.DataType;
  int Real_size = get_cdm_Datasize(Dtype);

  int size[3];
  for(int i=0; i<3; i++ ) size[i] = (int)DFI_Process.RankList[n].VoxelSize[i]+(int)(2*gc);

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.NumVariables > 1 ) dLen *= 3;  

  unsigned int dmy = dLen * Real_size;

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC6;
  if( val->writeBinary(fp) != dLen ) return CDM::E_CDM_ERROR_WRITE_SPH_REC6;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_SPH_REC6;
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// 平均の出力
CDM::E_CDM_ERRORCODE
cdm_DFI_SPH::write_averaged(FILE* fp,
                            const unsigned step_avr,
                            const double time_avr)
{
  int dType = 0;
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT32 ) dType = 1;
  if( DFI_Finfo.DataType == CDM::E_CDM_FLOAT64 ) dType = 2;

  unsigned int dmy;
  int Int_size,Real_size;
  if ( dType == 1 ) {
    dmy = 8;
    Int_size  = sizeof(int);
    Real_size = sizeof(float);
  }else{
    dmy = 16;
    Int_size  = sizeof(long long);
    Real_size = sizeof(double);
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_SPH_REC7;
  }

  //averaged step time の出力
  if( dType == 1 ){
    //float 型
    int istep = (int)step_avr;
    float ttime = (float)time_avr;
    if( fwrite(&istep, Int_size, 1, fp) != 1 ) {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_SPH_REC7;
    }
    if( fwrite(&ttime, Real_size, 1, fp) != 1 ) {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_SPH_REC7;
    }
  } else {
    //doublet 型
    long long dstep = (long long)step_avr;
    double ttime = (double)time_avr;
    if( fwrite(&dstep, Int_size, 1, fp) != 1 ) {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_SPH_REC7;
    }
    if( fwrite(&ttime, Real_size, 1, fp) != 1 ) {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_SPH_REC7;
    }
  }

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_SPH_REC7;
  }

  return CDM::E_CDM_SUCCESS;

}


