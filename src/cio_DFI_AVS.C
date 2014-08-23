/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_AVS.C
 * @brief  cio_DFI_AVS Class
 * @author aics    
 */

#include "cio_DFI.h"
#include "cio_DFI_AVS.h"

// #################################################################
// コンストラクタ
cio_DFI_AVS::cio_DFI_AVS()
{

}


// #################################################################
// デストラクタ
cio_DFI_AVS::~cio_DFI_AVS()
{

}

// #################################################################
// ヘッダーレコードの出力
CIO::E_CIO_ERRORCODE
cio_DFI_AVS::write_HeaderRecord(FILE* fp,
                                const unsigned step,
                                const double time,
                                const int n)
{

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// SPHデータレコードの出力
CIO::E_CIO_ERRORCODE
cio_DFI_AVS::write_DataRecord(FILE* fp,
                              cio_Array* val,
                              const int gc,
                              const int n)
{

  CIO::E_CIO_DTYPE Dtype = (CIO::E_CIO_DTYPE)DFI_Finfo.DataType;
  int Real_size = get_cio_Datasize(Dtype);

  const int *size = val->getArraySizeInt();
  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.Component > 1 ) dLen *= 3;

  if( val->writeBinary(fp) != dLen ) return CIO::E_CIO_ERROR_WRITE_FIELD_DATA_RECORD;

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// ヘッダーレコードの出力
bool cio_DFI_AVS::write_ascii_header(const unsigned step,
                                     const double time)
{

  int ndim,nspace;
  int dims[3];
  double min_ext[3],max_ext[3];
  double pit[3];

  //ピッチを計算
  for(int i=0; i<3; i++) {
    pit[i]=(DFI_Domain.GlobalRegion[i]/DFI_Domain.GlobalVoxel[i]);
  }

  //座標値の最小値、最大値をセット
  for(int i=0; i<3; i++) {
    min_ext[i]=DFI_Domain.GlobalOrigin[i]-pit[i]*0.5;
    max_ext[i]=min_ext[i]+((double)DFI_Process.RankList[m_RankID].VoxelSize[i])*pit[i];
  }

  //座標値データファイルの出力
  if( !write_avs_cord(min_ext,max_ext) ) return false;

  //ヘッダーデータファイルの出力
  if( !write_avs_header() ) return false; 

  return true;

}

// #################################################################
// 座標値データファイルの出力
bool cio_DFI_AVS::write_avs_cord(double min_ext[3],
                                 double max_ext[3])
{

  FILE* fp=NULL;

  //ファイル名の作成
  bool mio = false;
  if( DFI_MPI.NumberOfRank > 1 ) mio = true;
  
  std::string fname,tmp;
  tmp = Generate_FileName("cord",m_RankID,-1,"cod",m_output_fname,mio,
                          DFI_Finfo.TimeSliceDirFlag);
  if( CIO::cioPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
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

  //座標値データファイルクローズ
  fclose(fp);

  return true;

}

// #################################################################
// ヘッダーデータファイルの出力
bool cio_DFI_AVS::write_avs_header()
{
  FILE* fp=NULL;
  std::string dType;
  std::string out_fname;

  bool mio=false;

  //データタイプのセット
  if( GetDataType() == CIO::E_CIO_INT8 ) {
    dType = "byte";
  } else if( GetDataType() == CIO::E_CIO_INT16 ) {
    dType = "short";
  } else if( GetDataType() == CIO::E_CIO_INT32 ) {
    dType = "integer";
  } else if( GetDataType() == CIO::E_CIO_FLOAT32 ) {
    dType = "float";
  } else if( GetDataType() == CIO::E_CIO_FLOAT64 ) {
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
  if( CIO::cioPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
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
  fprintf(fp,"dim1=%d\n",DFI_Process.RankList[m_RankID].VoxelSize[0]+1);
  fprintf(fp,"dim2=%d\n",DFI_Process.RankList[m_RankID].VoxelSize[1]+1);
  fprintf(fp,"dim3=%d\n",DFI_Process.RankList[m_RankID].VoxelSize[2]+1);

  //物理空間の次元数を出力
  fprintf(fp,"nspace=%d\n",nspace);

  //成分数の出力
  fprintf(fp,"veclen=%d\n",DFI_Finfo.Component);

  //データのタイプ出力
  fprintf(fp,"data=%s\n",dType.c_str());

  //座標定義情報の出力
  fprintf(fp,"field=uniform\n");

  //labelの出力
  for(int i=0; i<DFI_Finfo.Component; i++) {
    std::string label=getComponentVariable(i);
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
    for(int j=1; j<=DFI_Finfo.Component; j++) {
      int skip;
      if( dType == "float" ) {
        skip=96+(j-1)*4;
      } else {
        skip=140+(j-1)*8;
      }
      out_fname=Generate_FileName(DFI_Finfo.Prefix,
                                  m_RankID,
                                  DFI_TimeSlice.SliceList[i].step,
                                  "sph",
                                  m_output_fname,
                                  mio,
                                  DFI_Finfo.TimeSliceDirFlag);
      //std::string xxx = CIO::cioPath_FileName(out_fname,"sph");
      fprintf(fp,"variable %d file=%s filetype=binary skip=%d stride=%d\n",
              j,out_fname.c_str(),skip,DFI_Finfo.Component);
    }

    //coord data file name 出力
    tmp = Generate_FileName("cord",m_RankID,-1,"cod",m_output_fname,mio,
                            DFI_Finfo.TimeSliceDirFlag);
    fprintf(fp,"coord 1 file=%s filetype=ascii skip=1\n",tmp.c_str());
    fprintf(fp,"coord 2 file=%s filetype=ascii skip=4\n",tmp.c_str());
    fprintf(fp,"coord 3 file=%s filetype=ascii skip=7\n",tmp.c_str());
    fprintf(fp,"EOT\n");
    
  }

  //出力ヘッダーファイルクローズ
  fclose(fp);

  //if( tmp ) delete tmp;

  return true;
}
