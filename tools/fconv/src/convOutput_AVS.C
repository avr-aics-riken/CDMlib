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
FILE* convOutput_AVS::OutputFile_Open(const std::string prefix,
                                const unsigned step,
                                const int id,
                                const bool mio)
{
  FILE* fp;

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
  if( (fp = fopen(outfile.c_str(), "wb")) == NULL ) {
    printf("\tCan't open file.(%s)\n",outfile.c_str());
    Exit(0);
  }

  return fp;
}

// #################################################################
bool 
convOutput_AVS::WriteFieldData(FILE* fp, cdm_Array* src, size_t dLen)
{

  const int* sz = src->getArraySizeInt();
  
  cdm_Array *out = cdm_Array::instanceArray
                   (src->getDataType(),
                    src->getArrayShape(),
                    (int *)sz,
                    0,
                    src->getNcomp());
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

  //成分数の出力
  fprintf(fp,"veclen=%d\n",DFI_FInfo->Component);

  //データのタイプ出力
  fprintf(fp,"data=%s\n",dType.c_str());

  //座標定義情報の出力
  fprintf(fp,"field=uniform\n");

  //labelの出力
  for(int j=0; j<DFI_FInfo->Component; j++) {
    std::string label=dfi->getComponentVariable(j);
    if( label == "" ) continue;
    fprintf(fp,"label=%s\n",label.c_str());
  }

  //step毎の出力
  if( TSlice->SliceList.size()>1 ) {
    fprintf(fp,"nstep=%d\n",TSlice->SliceList.size());
  }
  for(int j=0; j<TSlice->SliceList.size(); j++ ) {
    fprintf(fp,"time value=%.6f\n",TSlice->SliceList[j].time);
    for(int n=1; n<=DFI_FInfo->Component; n++) {
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
              n,out_fname.c_str(),skip,DFI_FInfo->Component);
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
