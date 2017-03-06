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
 * @file   convOutput_AVS.C
 * @brief  convOutput_AVS Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_AVS.h"

// #################################################################
// コンストラクタ
convOutput_AVS::convOutput_AVS()
{


}

// #################################################################
// デストラクタ
convOutput_AVS::~convOutput_AVS()
{


}

// #################################################################
// 出力ファイルをオープンする。
cdm_FILE* convOutput_AVS::OutputFile_Open(const std::string prefix,
                                const unsigned step,
                                const int id,
                                const bool mio)
{
  cdm_FILE* pFile;

  //ファイル名の生成
  std::string outfile;
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  outfile = m_InputCntl->Get_OutputDir()+"/"+
            cdm_DFI::Generate_FileName(prefix,
                                       id,
                                       step,
                                       "dat",
                                       fnameformat,
                                       mio,
                                       CDM::E_CDM_OFF);

  //ファイルオープン
  if( (pFile = cdm_FILE::OpenWriteBinary(outfile, CDM::E_CDM_FMT_AVS)) == NULL ) {
    printf("\tCan't open file.(%s)\n",outfile.c_str());
    Exit(0);
  }

  return pFile;
}

// #################################################################
bool
convOutput_AVS::WriteFieldData(cdm_FILE* pFile, cdm_Array* src, size_t dLen)
{
  FILE *fp = pFile->m_fp;

  const int* sz = src->getArraySizeInt();

  cdm_Array *out = cdm_Array::instanceArray
                   (src->getDataType(),
                    src->getArrayShape(),
                    (int *)sz,
                    0,
                    src->getNvari());
  int ret = src->copyArray(out);
  if( out->writeBinary(fp) != dLen ) {
    delete out;
    Exit(0);
  }

  delete out;
  return true;
}

// #################################################################
// output_avs
/*
void
convOutput_AVS::output_avs(
                           int myRank,
                           vector<cdm_DFI *>in_dfi,
                           cpm_ParaManager* paraMngr,
                           int *head)
{

  if      ( m_InputCntl->Get_ConvType() == E_OUTPUT_Mx1 ) {
    output_avs_Mx1(myRank, in_dfi);
  }else if( m_InputCntl->Get_ConvType() == E_OUTPUT_MxM ) {
    output_avs_MxM(myRank, in_dfi);
  }else if( m_InputCntl->Get_ConvType() == E_OUTPUT_MxN ) {
    output_avs_MxN(myRank, in_dfi, paraMngr, head);
  }

}
*/
// #################################################################
// output avs file
void convOutput_AVS::output_avs(int myRank,
                                vector<cdm_DFI *>in_dfi)
{

  if( myRank != 0 ) return; //myRank==0のときのみヘッダーレコードを出力

  //間引き数のセット
  int thin_count = m_InputCntl->Get_ThinOut();

  int ndim,nspace;

  int dims[3];
  double min_ext[3],max_ext[3];

  for(int i=0; i<in_dfi.size(); i++) {

    //cdm_Domainクラスポインタの取得
    const cdm_Domain* DFI_Domain  = in_dfi[i]->GetcdmDomain();

    //間引きを考慮しての計算空間サイズをセット
    dims[0]=DFI_Domain->GlobalVoxel[0]/thin_count;
    dims[1]=DFI_Domain->GlobalVoxel[1]/thin_count;
    dims[2]=DFI_Domain->GlobalVoxel[2]/thin_count;
    if(DFI_Domain->GlobalVoxel[0]%thin_count != 0 ) dims[0]++;
    if(DFI_Domain->GlobalVoxel[1]%thin_count != 0 ) dims[1]++;
    if(DFI_Domain->GlobalVoxel[2]%thin_count != 0 ) dims[2]++;

    if( i==0 ) {
      //座標値の最小値、最大値をセット
      for(int j=0; j<3; j++) {
        double pit=(DFI_Domain->GlobalRegion[j])/(double)(DFI_Domain->GlobalVoxel[j])*thin_count;
        //min_ext[j]=DFI_Domain->GlobalOrigin[j]+0.5*(pit*(double)thin_count);
        //max_ext[j]=(DFI_Domain->GlobalOrigin[j]+DFI_Domain->GlobalRegion[j])-0.5*(pit*(double)thin_count);
        min_ext[j]=DFI_Domain->GlobalOrigin[j];
        max_ext[j]=DFI_Domain->GlobalOrigin[j]+(pit*(double)dims[j]);
      }

      //coord データファイル出力
      output_avs_coord(myRank, false, min_ext, max_ext);
    }

    //avsのヘッダーファイル出力
    ndim=3;
    nspace=3;
    output_avs_header(in_dfi[i], myRank, false, ndim, nspace, dims);

  }

  return;
}
// #################################################################
// output avs file (不等間隔格子対応版)
void convOutput_AVS::output_avs(int myRank,
                                vector<cdm_DFI *>in_dfi,
                                cdm_Domain* out_domain,
                                cdm_Process* out_process,
                                int gc)
{

  if( myRank != 0 ) return; //myRank==0のときのみヘッダーレコードを出力

  //格子点の数,headをセット
  int dims[3];
  for(int i=0; i<3; i++) dims[i] = out_process->RankList[myRank].VoxelSize[i];
  if( m_InputCntl->Get_Interp_flag() ) {
    for(int i=0; i<3; i++) dims[i] += 1;
  } else {
    for(int i=0; i<3; i++) dims[i] += 2*gc;
  }

  for(int i=0; i<in_dfi.size(); i++) {

    CDM::E_CDM_DFITYPE dfi_type = in_dfi[i]->GetDFIType();
    //coord データファイル出力
    if( i==0 ) output_avs_coord(myRank, false, out_domain, out_process, dfi_type, dims, gc);

    //avsのヘッダーファイル出力
    output_avs_header(in_dfi[i], myRank, false, dfi_type, dims);

  }

  return;
}
// #################################################################
// output avs filei MxM
void convOutput_AVS::output_avs_MxM(int myRank,
                                    vector<cdm_DFI *>in_dfi)
{
  if( myRank != 0 ) return; //myRank==0のときのみヘッダーレコードを出力

  //間引き数のセット
  int thin_count = m_InputCntl->Get_ThinOut();

  int ndim,nspace;

  int dims[3];
  double min_ext[3],max_ext[3];

  double pit[3];

  for(int i=0; i<in_dfi.size(); i++) {
    //cdm_Domainクラスポインタの取得
    const cdm_Domain* DFI_Domain  = in_dfi[i]->GetcdmDomain();
    //ピッチのセット
    pit[0]=(DFI_Domain->GlobalRegion[0]/(double)DFI_Domain->GlobalVoxel[0])*(double)thin_count;
    pit[1]=(DFI_Domain->GlobalRegion[1]/(double)DFI_Domain->GlobalVoxel[1])*(double)thin_count;
    pit[2]=(DFI_Domain->GlobalRegion[2]/(double)DFI_Domain->GlobalVoxel[2])*(double)thin_count;

    //cdm_Processクラスポインタの取得
    const cdm_Process* DFI_Process = in_dfi[i]->GetcdmProcess();

    for(int j=0; j<DFI_Process->RankList.size(); j++) {
      //間引きを考慮しての計算空間サイズをセット
      dims[0]=DFI_Process->RankList[j].VoxelSize[0]/thin_count;
      dims[1]=DFI_Process->RankList[j].VoxelSize[1]/thin_count;
      dims[2]=DFI_Process->RankList[j].VoxelSize[2]/thin_count;
      if(DFI_Process->RankList[j].VoxelSize[0]%thin_count != 0) dims[0]++;
      if(DFI_Process->RankList[j].VoxelSize[1]%thin_count != 0) dims[1]++;
      if(DFI_Process->RankList[j].VoxelSize[2]%thin_count != 0) dims[2]++;

      if( i==0 ) {
        //座標値の最小値、最大値をセット
        for(int k=0; k<3; k++) {
          int head = (DFI_Process->RankList[j].HeadIndex[k]-1)/thin_count;
          if( (DFI_Process->RankList[j].HeadIndex[k]-1)%thin_count != 0 ) head++;
          min_ext[k]=DFI_Domain->GlobalOrigin[k]+0.5*pit[k]+double(head)*pit[k];
          max_ext[k]=min_ext[k]+(double(dims[k]-1))*pit[k];
        }

        //coord データファイル出力
        output_avs_coord(j, true, min_ext, max_ext);
      }

      //avsのヘッダーファイル出力
      ndim=3;
      nspace=3;
      output_avs_header(in_dfi[i], j, true, ndim, nspace, dims);

    }
  }

  return;

}

// #################################################################
// output avs filei MxN
void convOutput_AVS::output_avs_MxN(int myRank,
                                    vector<cdm_DFI *>in_dfi,
                                    cpm_ParaManager* paraMngr,
                                    int *mHead)
{
  int ndim,nspace;
  int dims[3];
  double min_ext[3],max_ext[3];

  //間引き数のセット
  int thin_count = m_InputCntl->Get_ThinOut();

  //自ノードでのボクセルサイズ取得
  const int* tmp = paraMngr->GetLocalVoxelSize();

  //間引きを考慮しての計算空間サイズをセット
  for(int i=0; i<3; i++) {
    dims[i]=tmp[i]/thin_count;
    if( tmp[i]%thin_count != 0 ) dims[i]++;
  }

  //pitを取得
  const double* dtmp;
  double pit[3];
  dtmp = paraMngr->GetPitch();
  for(int i=0; i<3; i++) pit[i]=dtmp[i]*double(thin_count);

  //座標値の最小値、最大値をセット
  dtmp = paraMngr->GetGlobalOrigin();
  for(int i=0; i<3; i++) {
    int head=(mHead[i]-1)/thin_count;
    if( (mHead[i]-1)%thin_count ) head++;
    min_ext[i]=dtmp[i]+0.5*pit[i]+double(head)*pit[i];
    max_ext[i]=min_ext[i]+(double(dims[i]-1))*pit[i];
  }

  //coord データファイル出力
  output_avs_coord(myRank, true, min_ext, max_ext);

  for(int i=0; i<in_dfi.size(); i++) {
    //avsのヘッダーファイル出力
    ndim=3;
    nspace=3;
    output_avs_header(in_dfi[i], myRank, true, ndim, nspace, dims);
  }

  return;

}

// #################################################################
// output coord data (cod)
void convOutput_AVS::output_avs_coord(int RankID,
                                      bool mio,
                                      double min_ext[3],
                                      double max_ext[3])
{

  FILE* fp;
  std::string cod_fname;

  //座標値データファイルオープン
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  cod_fname = m_InputCntl->Get_OutputDir() +"/"+
              cdm_DFI::Generate_FileName("cord",
                                         RankID,
                                         -1,
                                         "cod",
                                         fnameformat,
                                         mio,
                                         CDM::E_CDM_OFF);

  if( (fp = fopen(cod_fname.c_str(),"w")) == NULL ) {
    printf("\tCan't open file.(%s)\n",cod_fname.c_str());
    Exit(0);
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

}
// #################################################################
// output coord data (cod) (不等間隔格子対応版)
void convOutput_AVS::output_avs_coord(int RankID,
                                      bool mio,
                                      cdm_Domain* out_domain,
                                      cdm_Process* out_process,
                                      CDM::E_CDM_DFITYPE dfi_type,
                                      int dims[3],
                                      int gc)
{

  FILE* fp;
  std::string cod_fname;

  //座標値データファイルオープン
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  cod_fname = m_InputCntl->Get_OutputDir() +"/"+
              cdm_DFI::Generate_FileName("cord",
                                         RankID,
                                         -1,
                                         "cod",
                                         fnameformat,
                                         mio,
                                         CDM::E_CDM_OFF);

  if( (fp = fopen(cod_fname.c_str(),"w")) == NULL ) {
    printf("\tCan't open file.(%s)\n",cod_fname.c_str());
    Exit(0);
  }

  int head[3],tail[3];
  for(int i=0; i<3; i++) {
    head[i] = out_process->RankList[RankID].HeadIndex[i];
    tail[i] = out_process->RankList[RankID].TailIndex[i];
  }

  if( dfi_type == CDM::E_CDM_DFITYPE_CARTESIAN ) {

    double min_ext[3],max_ext[3];
    if( m_InputCntl->Get_Interp_flag() ) {
      min_ext[0] = out_domain->NodeX(head[0]-1);
      min_ext[1] = out_domain->NodeY(head[1]-1);
      min_ext[2] = out_domain->NodeZ(head[2]-1);
      max_ext[0] = out_domain->NodeX(tail[0]);
      max_ext[1] = out_domain->NodeY(tail[1]);
      max_ext[2] = out_domain->NodeZ(tail[2]);
    } else {
      min_ext[0] = out_domain->CellX(head[0]-1-gc);
      min_ext[1] = out_domain->CellY(head[1]-1-gc);
      min_ext[2] = out_domain->CellZ(head[2]-1-gc);
      max_ext[0] = out_domain->CellX(tail[0]-1+gc);
      max_ext[1] = out_domain->CellY(tail[1]-1+gc);
      max_ext[2] = out_domain->CellZ(tail[2]-1+gc);
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

  } else if( dfi_type == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {

    //ascii
    if( m_InputCntl->Get_OutputFileTypeCoord() == CDM::E_CDM_FILE_TYPE_ASCII ) {

      if( m_InputCntl->Get_Interp_flag() ) {
        //格子点補間する場合(ガイドセル出力はなし)
        fprintf(fp,"#### X #####\n");
        for (int i=0; i<dims[0]; i++ ) {
          fprintf( fp, "%.6f\n", out_domain->NodeX(i+head[0]-1) );
        }
        fprintf(fp,"#### Y #####\n");
        for (int j=0; j<dims[1]; j++ ) {
          fprintf( fp, "%.6f\n", out_domain->NodeY(j+head[1]-1) );
        }
        fprintf(fp,"#### Z #####\n");
        for (int k=0; k<dims[2]; k++ ) {
          fprintf( fp, "%.6f\n", out_domain->NodeZ(k+head[2]-1) );
        }
      } else {
        fprintf(fp,"#### X #####\n");
        for (int i=0; i<dims[0]; i++ ) {
          fprintf( fp, "%.6f\n", out_domain->CellX(i+head[0]-1-gc) );
        }
        fprintf(fp,"#### Y #####\n");
        for (int j=0; j<dims[1]; j++ ) {
          fprintf( fp, "%.6f\n", out_domain->CellY(j+head[1]-1-gc) );
        }
        fprintf(fp,"#### Z #####\n");
        for (int k=0; k<dims[2]; k++ ) {
          fprintf( fp, "%.6f\n", out_domain->CellZ(k+head[2]-1-gc) );
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

      if( m_InputCntl->Get_Interp_flag() ) {
        //格子点補間する場合(ガイドセル出力はなし)
        for(int i=0; i<dims[0]; i++) coord_X[i] = (float)(out_domain->NodeX(i+head[0]-1));
        for(int j=0; j<dims[1]; j++) coord_Y[j] = (float)(out_domain->NodeY(j+head[1]-1));
        for(int k=0; k<dims[2]; k++) coord_Z[k] = (float)(out_domain->NodeZ(k+head[2]-1));
      } else {
        for(int i=0; i<dims[0]; i++) coord_X[i] = (float)(out_domain->CellX(i+head[0]-1-gc));
        for(int j=0; j<dims[1]; j++) coord_Y[j] = (float)(out_domain->CellY(j+head[1]-1-gc));
        for(int k=0; k<dims[2]; k++) coord_Z[k] = (float)(out_domain->CellZ(k+head[2]-1-gc));
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

}
// #################################################################
// output avc Header (fld)
void convOutput_AVS::output_avs_header(cdm_DFI* dfi,
                                       int RankID,
                                       bool mio,
                                       int ndim,
                                       int nspace,
                                       int dims[3])
{

  FILE* fp;
  std::string dType;
  std::string fld_fname, out_fname;
  std::string cod_fname;

  //cdm_FileInfoクラスポインタの取得
  const cdm_FileInfo* DFI_FInfo = dfi->GetcdmFileInfo();
  //cdm_TimeSliceクラスポインタの取得
  const cdm_TimeSlice* TSlice   = dfi->GetcdmTimeSlice();

  //データタイプのセット
  int out_dtype = m_InputCntl->Get_OutputDataType();
  if( out_dtype == CDM::E_CDM_DTYPE_UNKNOWN ) out_dtype = dfi->GetDataType();
  if(      out_dtype == CDM::E_CDM_INT8    ) {
    dType="byte";
  } else if( out_dtype == CDM::E_CDM_INT16   ) {
    dType="short";
  } else if( out_dtype == CDM::E_CDM_INT32   ) {
    dType="integer";
  } else if( out_dtype == CDM::E_CDM_FLOAT32 ) {
    dType="float";
  } else if( out_dtype == CDM::E_CDM_FLOAT64 ) {
    dType="double";
  } else {
    dType = m_InputCntl->Get_OutputDataType_string();
    printf("\tillergal data type.(%s)\n",dType.c_str());
    Exit(0);
  }

  //出力ヘッダーファイルオープン
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  fld_fname = m_InputCntl->Get_OutputDir() +"/"+
              cdm_DFI::Generate_FileName(DFI_FInfo->Prefix,
                                         RankID,
                                         -1,
                                         "fld",
                                         fnameformat,
                                         mio,
                                         CDM::E_CDM_OFF);

  if( (fp = fopen(fld_fname.c_str(),"w")) == NULL ) {
    printf("\tCan't open file.(%s)\n",fld_fname.c_str());
    Exit(0);
  }

  //先頭レコードの出力
  fprintf(fp,"# AVS field file\n");

  //計算空間の次元数を出力
  fprintf(fp,"ndim=%d\n",ndim);

  //計算空間サイズを出力
  fprintf(fp,"dim1=%d\n",dims[0]+1);
  fprintf(fp,"dim2=%d\n",dims[1]+1);
  fprintf(fp,"dim3=%d\n",dims[2]+1);

  //物理空間の次元数を出力
  fprintf(fp,"nspace=%d\n",nspace);

  //変数の個数の出力
  fprintf(fp,"veclen=%d\n",DFI_FInfo->NumVariables);

  //データのタイプ出力
  fprintf(fp,"data=%s\n",dType.c_str());

  //座標定義情報の出力
  fprintf(fp,"field=uniform\n");

  //labelの出力
  for(int j=0; j<DFI_FInfo->NumVariables; j++) {
    std::string label=dfi->getVariableName(j);
    if( label == "" ) continue;
    fprintf(fp,"label=%s\n",label.c_str());
  }

  //step毎の出力
  if( TSlice->SliceList.size()>1 ) {
    fprintf(fp,"nstep=%d\n",TSlice->SliceList.size());
  }
  for(int j=0; j<TSlice->SliceList.size(); j++ ) {
    fprintf(fp,"time value=%.6f\n",TSlice->SliceList[j].time);
    for(int n=1; n<=DFI_FInfo->NumVariables; n++) {
      int skip;
      /*
      if( dType == "float" ) {
        skip=96+(n-1)*4;
      } else {
        skip=140+(n-1)*8;
      }
      */
      skip=0;
      out_fname = cdm_DFI::Generate_FileName(DFI_FInfo->Prefix,
                                             RankID,
                                             TSlice->SliceList[j].step,
                                             "sph",
                                             fnameformat,
                                             mio,
                                             CDM::E_CDM_OFF);

      fprintf(fp,"variable %d file=%s filetype=binary skip=%d stride=%d\n",
              n,out_fname.c_str(),skip,DFI_FInfo->NumVariables);
    }
    //coord data
    cod_fname = cdm_DFI::Generate_FileName("cord",
                                           RankID,
                                           -1,
                                           "cod",
                                           fnameformat,
                                           mio,
                                           CDM::E_CDM_OFF);

    fprintf(fp,"coord 1 file=%s filetype=ascii skip=1\n",cod_fname.c_str());
    fprintf(fp,"coord 2 file=%s filetype=ascii skip=4\n",cod_fname.c_str());
    fprintf(fp,"coord 3 file=%s filetype=ascii skip=7\n",cod_fname.c_str());
    fprintf(fp,"EOT\n");
  }


  //出力ヘッダーファイルクローズ
  fclose(fp);

}
// #################################################################
// output avc Header (fld) (不等間隔格子対応版)
void convOutput_AVS::output_avs_header(cdm_DFI* dfi,
                                       int RankID,
                                       bool mio,
                                       CDM::E_CDM_DFITYPE dfi_type,
                                       int dims[3])
{

  FILE* fp;
  std::string dType;
  std::string fld_fname, out_fname;
  std::string cod_fname;

  //cdm_FileInfoクラスポインタの取得
  const cdm_FileInfo* DFI_FInfo = dfi->GetcdmFileInfo();
  //cdm_TimeSliceクラスポインタの取得
  const cdm_TimeSlice* TSlice   = dfi->GetcdmTimeSlice();

  //データタイプのセット
  int out_dtype = m_InputCntl->Get_OutputDataType();
  if( out_dtype == CDM::E_CDM_DTYPE_UNKNOWN ) out_dtype = dfi->GetDataType();
  if(      out_dtype == CDM::E_CDM_INT8    ) {
    dType="byte";
  } else if( out_dtype == CDM::E_CDM_INT16   ) {
    dType="short";
  } else if( out_dtype == CDM::E_CDM_INT32   ) {
    dType="integer";
  } else if( out_dtype == CDM::E_CDM_FLOAT32 ) {
    dType="float";
  } else if( out_dtype == CDM::E_CDM_FLOAT64 ) {
    dType="double";
  } else {
    dType = m_InputCntl->Get_OutputDataType_string();
    printf("\tillergal data type.(%s)\n",dType.c_str());
    Exit(0);
  }

  //出力ヘッダーファイルオープン
  CDM::E_CDM_OUTPUT_FNAME fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  fld_fname = m_InputCntl->Get_OutputDir() +"/"+
              cdm_DFI::Generate_FileName(DFI_FInfo->Prefix,
                                         RankID,
                                         -1,
                                         "fld",
                                         fnameformat,
                                         mio,
                                         CDM::E_CDM_OFF);

  if( (fp = fopen(fld_fname.c_str(),"w")) == NULL ) {
    printf("\tCan't open file.(%s)\n",fld_fname.c_str());
    Exit(0);
  }

  int ndim = 3;
  int nspace = 3;

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
  fprintf(fp,"veclen=%d\n",DFI_FInfo->NumVariables);

  //データのタイプ出力
  fprintf(fp,"data=%s\n",dType.c_str());

  //座標定義情報の出力
  if( dfi_type == CDM::E_CDM_DFITYPE_CARTESIAN ) {
    fprintf(fp,"field=uniform\n");
  } else if( dfi_type == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
    fprintf( fp, "field=rectilinear\n" );
  }

  //labelの出力
  for(int j=0; j<DFI_FInfo->NumVariables; j++) {
    std::string label=dfi->getVariableName(j);
    if( label == "" ) continue;
    fprintf(fp,"label=%s\n",label.c_str());
  }

  //step毎の出力
  if( TSlice->SliceList.size()>1 ) {
    fprintf(fp,"nstep=%d\n",TSlice->SliceList.size());
  }
  for(int j=0; j<TSlice->SliceList.size(); j++ ) {
    fprintf(fp,"time value=%.6f\n",TSlice->SliceList[j].time);

    //field data file name 出力
    for(int n=1; n<=DFI_FInfo->NumVariables; n++) {
      int skip;
      if( dType == "float" ) {
        skip=(n-1)*4;
      } else {
        skip=(n-1)*8;
      }
      //skip=0;
      out_fname = cdm_DFI::Generate_FileName(DFI_FInfo->Prefix,
                                             RankID,
                                             TSlice->SliceList[j].step,
                                             "dat",
                                             fnameformat,
                                             mio,
                                             CDM::E_CDM_OFF);

      fprintf(fp,"variable %d file=%s filetype=binary skip=%d stride=%d\n",
              n,out_fname.c_str(),skip,DFI_FInfo->NumVariables);
    }
    //coord data
    cod_fname = cdm_DFI::Generate_FileName("cord",
                                           RankID,
                                           -1,
                                           "cod",
                                           fnameformat,
                                           mio,
                                           CDM::E_CDM_OFF);

    if( dfi_type == CDM::E_CDM_DFITYPE_CARTESIAN ) {
      fprintf(fp,"coord 1 file=%s filetype=ascii skip=1\n",cod_fname.c_str());
      fprintf(fp,"coord 2 file=%s filetype=ascii skip=4\n",cod_fname.c_str());
      fprintf(fp,"coord 3 file=%s filetype=ascii skip=7\n",cod_fname.c_str());
    } else if( dfi_type == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
      if( m_InputCntl->Get_OutputFileTypeCoord() == CDM::E_CDM_FILE_TYPE_ASCII ) {
        fprintf(fp,"coord 1 file=%s filetype=ascii skip=%d\n",cod_fname.c_str(),1);
        fprintf(fp,"coord 2 file=%s filetype=ascii skip=%d\n",cod_fname.c_str(),dims[0]+2);
        fprintf(fp,"coord 3 file=%s filetype=ascii skip=%d\n",cod_fname.c_str(),dims[0]+dims[1]+3);
      } else {
        //AVS形式の座標データはfloat型(4バイト)のみサポート
        fprintf(fp,"coord 1 file=%s filetype=binary skip=%d\n",cod_fname.c_str(),0);
        fprintf(fp,"coord 2 file=%s filetype=binary skip=%d\n",cod_fname.c_str(),4*dims[0]);
        fprintf(fp,"coord 3 file=%s filetype=binary skip=%d\n",cod_fname.c_str(),4*(dims[0]+dims[1]));
      }
    }
    fprintf(fp,"EOT\n");
  }

  //出力ヘッダーファイルクローズ
  fclose(fp);

}
