/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_AVS.C
 * @brief  cdm_DFI_AVS Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include "cdm_DFI_AVS.h"

// #################################################################
// コンストラクタ
cdm_DFI_AVS::cdm_DFI_AVS()
{

}


// #################################################################
// デストラクタ
cdm_DFI_AVS::~cdm_DFI_AVS()
{

}

// #################################################################
// ヘッダーレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_AVS::write_HeaderRecord(FILE* fp,
                                const unsigned step,
                                const double time,
                                const int n)
{

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// AVSデータレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_AVS::write_DataRecord(FILE* fp,
                              cdm_Array* val,
                              const int gc,
                              const int n)
{

  CDM::E_CDM_DTYPE Dtype = (CDM::E_CDM_DTYPE)DFI_Finfo.DataType;
  int Real_size = get_cdm_Datasize(Dtype);

  const int* sz_without_gc = val->getArraySizeInt();
  int sz[3];
  for(int i=0; i<3; i++) sz[i] = sz_without_gc[i];
  if( !m_bgrid_interp_flag ) for(int i=0; i<3; i++) sz[i] += (int)(2*gc);
  size_t dLen = (size_t)sz[0]*(size_t)sz[1]*(size_t)sz[2]*val->getNvari();

  if( val->writeBinary(fp) != dLen ) return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ヘッダーレコードの出力
bool cdm_DFI_AVS::write_ascii_header(const unsigned step,
                                     const double time)
{

  //格子点の数をセット
  int dims[3];
  for(int i=0; i<3; i++) dims[i] = DFI_Process.RankList[m_RankID].VoxelSize[i];
  int gc = DFI_Finfo.GuideCell;
  if( m_bgrid_interp_flag ) {
    for(int i=0; i<3; i++) dims[i] += 1;
  } else {
    for(int i=0; i<3; i++) dims[i] += 2*gc;
  }

  //座標値データファイルの出力
  if( !write_avs_cord(dims,gc) ) return false;

  //ヘッダーデータファイルの出力
  if( !write_avs_header(dims) ) return false; 

  return true;

}

// #################################################################
// 座標値データファイルの出力
bool cdm_DFI_AVS::write_avs_cord(int dims[3],
                                 int gc)
{

  FILE* fp=NULL;

  //ファイル名の作成
  bool mio = false;
  if( DFI_MPI.NumberOfRank > 1 ) mio = true;
  
  std::string fname,tmp;
  tmp = Generate_FileName("cord",m_RankID,-1,"cod",m_output_fname,mio,
                          DFI_Finfo.TimeSliceDirFlag);
  if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
    fname = DFI_Finfo.DirectoryPath + "/" + tmp;
  } else {
    fname = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
  }

  //printf("cord file name : %s\n",fname.c_str());

  //座標値データファイルオープン
  if( (fp = fopen(fname.c_str(),"w"))  == NULL ) {
    printf("\tCan't open file.(%s)\n",fname.c_str());
    return false;
  }

  int head[3],tail[3];
  for(int i=0; i<3; i++) {
    head[i] = DFI_Process.RankList[m_RankID].HeadIndex[i];
    tail[i] = DFI_Process.RankList[m_RankID].TailIndex[i];
  }

  if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {

    double min_ext[3],max_ext[3];
    if( m_bgrid_interp_flag ) {
      min_ext[0] = DFI_Domain->NodeX(head[0]-1);
      min_ext[1] = DFI_Domain->NodeY(head[1]-1);
      min_ext[2] = DFI_Domain->NodeZ(head[2]-1);
      max_ext[0] = DFI_Domain->NodeX(tail[0]);
      max_ext[1] = DFI_Domain->NodeY(tail[1]);
      max_ext[2] = DFI_Domain->NodeZ(tail[2]);
    } else {
      min_ext[0] = DFI_Domain->CellX(head[0]-1-gc);
      min_ext[1] = DFI_Domain->CellY(head[1]-1-gc);
      min_ext[2] = DFI_Domain->CellZ(head[2]-1-gc);
      max_ext[0] = DFI_Domain->CellX(tail[0]-1+gc);
      max_ext[1] = DFI_Domain->CellY(tail[1]-1+gc);
      max_ext[2] = DFI_Domain->CellZ(tail[2]-1+gc);
    }
    //座標値データ（min,max)の出力
    fprintf(fp,"#### X #####\n");
    fprintf(fp,"%.6f\n",min_ext[0]);
    fprintf(fp,"%.6f\n",max_ext[0]);
    fprintf(fp,"#### Y #####\n");
    fprintf(fp,"%.6f\n",min_ext[1]);
    fprintf(fp,"%.6f\n",max_ext[1]);
    fprintf(fp,"#### Z #####\n");
    fprintf(fp,"%.6f\n",min_ext[2]);
    fprintf(fp,"%.6f\n",max_ext[2]);

  } else if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {

    //ascii
    if( m_output_type_coord == CDM::E_CDM_FILE_TYPE_ASCII ) {

      if( m_bgrid_interp_flag ) {
        //格子点補間する場合(ガイドセル出力はなし)
        fprintf(fp,"#### X #####\n");
        for (int i=0; i<dims[0]; i++ ) {
          fprintf( fp, "%.6f\n", DFI_Domain->NodeX(i+head[0]-1) );
        }
        fprintf(fp,"#### Y #####\n");
        for (int j=0; j<dims[1]; j++ ) {
          fprintf( fp, "%.6f\n", DFI_Domain->NodeY(j+head[1]-1) );
        }
        fprintf(fp,"#### Z #####\n");
        for (int k=0; k<dims[2]; k++ ) {
          fprintf( fp, "%.6f\n", DFI_Domain->NodeZ(k+head[2]-1) );
        }
      } else {
        fprintf(fp,"#### X #####\n");
        for (int i=0; i<dims[0]; i++ ) {
          fprintf( fp, "%.6f\n", DFI_Domain->CellX(i+head[0]-1-gc) );
        }
        fprintf(fp,"#### Y #####\n");
        for (int j=0; j<dims[1]; j++ ) {
          fprintf( fp, "%.6f\n", DFI_Domain->CellY(j+head[1]-1-gc) );
        }
        fprintf(fp,"#### Z #####\n");
        for (int k=0; k<dims[2]; k++ ) {
          fprintf( fp, "%.6f\n", DFI_Domain->CellZ(k+head[2]-1-gc) );
        }
      }

    //binary
    } else {

      //AVS形式の座標データはfloat型のみサポート
      float *coord_X = NULL;
      float *coord_Y = NULL;
      float *coord_Z = NULL;
      coord_X = new float[dims[0]];
      coord_Y = new float[dims[1]];
      coord_Z = new float[dims[2]];

      if( m_bgrid_interp_flag ) {
        //格子点補間する場合(ガイドセル出力はなし)
        for(int i=0; i<dims[0]; i++) coord_X[i] = (float)(DFI_Domain->NodeX(i+head[0]-1));
        for(int j=0; j<dims[1]; j++) coord_Y[j] = (float)(DFI_Domain->NodeY(j+head[1]-1));
        for(int k=0; k<dims[2]; k++) coord_Z[k] = (float)(DFI_Domain->NodeZ(k+head[2]-1));
      } else {
        for(int i=0; i<dims[0]; i++) coord_X[i] = (float)(DFI_Domain->CellX(i+head[0]-1-gc));
        for(int j=0; j<dims[1]; j++) coord_Y[j] = (float)(DFI_Domain->CellY(j+head[1]-1-gc));
        for(int k=0; k<dims[2]; k++) coord_Z[k] = (float)(DFI_Domain->CellZ(k+head[2]-1-gc));
      }

      fwrite(coord_X, sizeof(float), dims[0], fp);
      fwrite(coord_Y, sizeof(float), dims[1], fp);
      fwrite(coord_Z, sizeof(float), dims[2], fp);

      delete [] coord_X;
      delete [] coord_Y;
      delete [] coord_Z;

    }

  }

  //座標値データファイルクローズ
  fclose(fp);

  return true;

}

// #################################################################
// ヘッダーデータファイルの出力
bool cdm_DFI_AVS::write_avs_header(int dims[3])
{
  FILE* fp=NULL;
  std::string dType;
  std::string out_fname;

  bool mio=false;

  //データタイプのセット
  if( GetDataType() == CDM::E_CDM_INT8 ) {
    dType = "byte";
  } else if( GetDataType() == CDM::E_CDM_INT16 ) {
    dType = "short";
  } else if( GetDataType() == CDM::E_CDM_INT32 ) {
    dType = "integer";
  } else if( GetDataType() == CDM::E_CDM_FLOAT32 ) {
    dType = "float";
  } else if( GetDataType() == CDM::E_CDM_FLOAT64 ) {
    dType = "double";
  } else {
    dType = GetDataTypeString();
    printf("\tillergal data type.(%s)\n",dType.c_str());
    return false;
  }

  //ファイル名生成

  if( DFI_MPI.NumberOfRank > 1 ) mio = true;
  std::string fname,tmp;
  tmp = Generate_FileName(DFI_Finfo.Prefix,m_RankID,-1,"fld",m_output_fname,mio,
                          DFI_Finfo.TimeSliceDirFlag);
  if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
    fname = DFI_Finfo.DirectoryPath +"/"+ tmp;
  } else {
    fname = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
  }

  //printf("fld file name : %s\n",fname.c_str());

  //出力ヘッダーファイルオープン
  if( (fp = fopen(fname.c_str(),"w"))  == NULL ) {
    printf("\tCan't open file.(%s)\n",fname.c_str());
    return false;
  }

  int ndim = 3;
  int nspace = 3;
  //dims = DFI_Process.RankList[m_RankID].VoxelSize[0]

  //先頭レコードの出力
  fprintf(fp,"# AVS field file\n");

  //計算空間の次元数を出力
  fprintf(fp,"ndim=%d\n",ndim);

  //計算空間サイズを出力
  fprintf(fp,"dim1=%d\n",dims[0]);
  fprintf(fp,"dim2=%d\n",dims[1]);
  fprintf(fp,"dim3=%d\n",dims[2]);

  //物理空間の次元数を出力
  fprintf(fp,"nspace=%d\n",nspace);

  //変数の個数の出力
  fprintf(fp,"veclen=%d\n",DFI_Finfo.NumVariables);

  //データのタイプ出力
  fprintf(fp,"data=%s\n",dType.c_str());

  //座標定義情報の出力
  if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {
    fprintf(fp,"field=uniform\n");
  } else if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
    fprintf( fp, "field=rectilinear\n" );
  }

  //labelの出力
  for(int i=0; i<DFI_Finfo.NumVariables; i++) {
    std::string label=getVariableName(i);
    if( label == "" ) continue;
    fprintf(fp,"label=%s\n",label.c_str());
  }

  //step毎の出力
  if( DFI_TimeSlice.SliceList.size()>1 ) {
    fprintf(fp,"nstep=%d\n",(int)DFI_TimeSlice.SliceList.size());
  }
  for(int i=0; i<DFI_TimeSlice.SliceList.size(); i++) {
    fprintf(fp,"time value=%.6f\n",DFI_TimeSlice.SliceList[i].time);

    //field data file name 出力
    for(int j=1; j<=DFI_Finfo.NumVariables; j++) {
      int skip;
      if( dType == "float" ) {
        skip=(j-1)*4;
        //skip=96+(j-1)*4;
      } else {
        skip=(j-1)*8;
        //skip=140+(j-1)*8;
      }
      out_fname=Generate_FileName(DFI_Finfo.Prefix,
                                  m_RankID,
                                  DFI_TimeSlice.SliceList[i].step,
                                  "dat",
                                  m_output_fname,
                                  mio,
                                  DFI_Finfo.TimeSliceDirFlag);
      //std::string xxx = CDM::cdmPath_FileName(out_fname,"sph");
      fprintf(fp,"variable %d file=%s filetype=binary skip=%d stride=%d\n",
              j,out_fname.c_str(),skip,DFI_Finfo.NumVariables);
    }

    //coord data file name 出力
    tmp = Generate_FileName("cord",m_RankID,-1,"cod",m_output_fname,mio,
                            DFI_Finfo.TimeSliceDirFlag);
    if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {
      fprintf(fp,"coord 1 file=%s filetype=ascii skip=1\n",tmp.c_str());
      fprintf(fp,"coord 2 file=%s filetype=ascii skip=4\n",tmp.c_str());
      fprintf(fp,"coord 3 file=%s filetype=ascii skip=7\n",tmp.c_str());
    } else if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
      if( m_output_type_coord == CDM::E_CDM_FILE_TYPE_ASCII ) {
        fprintf(fp,"coord 1 file=%s filetype=ascii skip=%d\n",tmp.c_str(),1);
        fprintf(fp,"coord 2 file=%s filetype=ascii skip=%d\n",tmp.c_str(),dims[0]+2);
        fprintf(fp,"coord 3 file=%s filetype=ascii skip=%d\n",tmp.c_str(),dims[0]+dims[1]+3);
      } else {
        //AVS形式の座標データはfloat型(4バイト)のみサポート
        fprintf(fp,"coord 1 file=%s filetype=binary skip=%d\n",tmp.c_str(),0);
        fprintf(fp,"coord 2 file=%s filetype=binary skip=%d\n",tmp.c_str(),4*dims[0]);
        fprintf(fp,"coord 3 file=%s filetype=binary skip=%d\n",tmp.c_str(),4*(dims[0]+dims[1]));
      }
    }
    fprintf(fp,"EOT\n");
    
  }

  //出力ヘッダーファイルクローズ
  fclose(fp);

  //if( tmp ) delete tmp;

  return true;
}
