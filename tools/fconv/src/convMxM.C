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
 * @file   convMxM.C 
 * @brief  convMxM Class
 * @author aics
 */

#include "convMxM.h"

// #################################################################
// コンストラクタ
convMxM::convMxM()
{

  m_StepRankList.clear();

}

// #################################################################
// デストラクタ
convMxM::~convMxM()
{

}

// #################################################################
//
bool convMxM::exec()
{

  if( m_myRank == 0 ) {
    printf("Convert M x M\n");
  }

  //並列実行時のファイル割振り方法の取得
  E_CONV_OUTPUT_MULTI_FILE_CAST outlist = m_param->Get_MultiFileCasting();

  //並列実行時のファイル割振り方法でのリスト生成
  if( outlist == E_CONV_OUTPUT_STEP || outlist == E_CONV_OUTPUT_CAST_UNKNOWN ) {
    makeStepList(m_StepRankList);
  } else if( outlist == E_CONV_OUTPUT_RANK ) {
    makeRankList(m_StepRankList);
  }

  //出力dfiファイル名の取得
  //vector<std::string> out_dfi_name = m_InputCntl->Get_OutdfiNameList();
  //vector<std::string> out_proc_name = m_InputCntl->Get_OutprocNameList();

  //minmaxの格納構造体のインスタンス
  vector<dfi_MinMax*> minmaxList;

  for(int i=0; i<m_in_dfi.size(); i++){
    const cdm_TimeSlice* TSlice = m_in_dfi[i]->GetcdmTimeSlice();
    int nVari = m_in_dfi[i]->GetNumVariables();

    dfi_MinMax *MinMax;
    if( nVari == 1 ) MinMax = new dfi_MinMax(TSlice->SliceList.size(),nVari);
    else             MinMax = new dfi_MinMax(TSlice->SliceList.size(),nVari+1);

    MinMax->dfi = m_in_dfi[i];
    minmaxList.push_back(MinMax);
  }

  //Listのループ
  for (int i=0; i<m_StepRankList.size(); i++) {

    //dfiのstepリストの取得
    const cdm_TimeSlice* TSlice = m_StepRankList[i].dfi->GetcdmTimeSlice();

    //変数の個数の取得
    int nVari = m_StepRankList[i].dfi->GetNumVariables();

    //stepのループ
    for(int j=m_StepRankList[i].stepStart; j<=m_StepRankList[i].stepEnd; j++) {

      //minmaxの初期化
      int nsize = nVari;
      if( nVari > 1 ) nsize++;
      double *min = new double[nsize]; 
      double *max = new double[nsize]; 
      for(int n=0; n<nsize; n++) {
         min[n]=DBL_MAX;
         max[n]=-DBL_MAX;
      }

      //rankのループ
      for(int k=m_StepRankList[i].rankStart; k<=m_StepRankList[i].rankEnd; k++) {
        //MxMの読込みコンバート出力
        if( !mxmsolv(m_StepRankList[i].dfi->get_dfi_fname(),
                     m_StepRankList[i].dfi, 
                     TSlice->SliceList[j].step, 
                     (float)TSlice->SliceList[j].time, 
                     k,
                     min,
                     max) ) return false;
      } 

      //dfiごとに登録
      for(int ndfi = 0; ndfi<minmaxList.size(); ndfi++) {
        if( minmaxList[ndfi]->dfi != m_StepRankList[i].dfi ) continue;
        for(int n=0; n<nsize; n++) {
          if( minmaxList[ndfi]->Min[j*nsize+n] > min[n] ) minmaxList[ndfi]->Min[j*nsize+n] = min[n];
          if( minmaxList[ndfi]->Max[j*nsize+n] < max[n] ) minmaxList[ndfi]->Max[j*nsize+n] = max[n]; 
        }
      }
    }
  }


  //出力dfiファイルがないときはreturn
  //if( out_dfi_name.size() < 1 || out_proc_name.size() < 1 ) return true;
  if( !m_param->Get_Outputdfi_on() ) return true;

  //ランク間で通信してMINMAXを求めてランク０に送信
  for(int i=0; i<minmaxList.size(); i++) {
   int nVari = minmaxList[i]->dfi->GetNumVariables();
   const cdm_TimeSlice* TSlice = minmaxList[i]->dfi->GetcdmTimeSlice();
   int nStep = TSlice->SliceList.size();

   int n = nVari*nStep;
   if( nVari > 1 ) n = (nVari+1)*nStep;

   //minの通信
   double *send1 =  minmaxList[i]->Min;
   double *recv1 = new double[n];
   MPI_Reduce(send1, recv1, n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
   minmaxList[i]->Min = recv1;

   //maxの通信
   double *send2 = minmaxList[i]->Max;
   double *recv2 = new double[n];
   MPI_Reduce(send2, recv2, n, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   minmaxList[i]->Max = recv2;

  }

  //出力dfiファイルの出力
  if( m_myRank == 0 ) {
    WriteIndexDfiFile(minmaxList);

    for(int i=0; i<m_in_dfi.size(); i++) {
      cdm_Domain* out_domain = NULL;
      cdm_MPI* out_mpi = NULL;
      cdm_Process* out_process = NULL; 
      const cdm_MPI* dfi_mpi = m_in_dfi[i]->GetcdmMPI();
      int numProc = dfi_mpi->NumberOfRank; 

      //Proc情報の生成 
      makeProcInfo(m_in_dfi[i],out_domain,out_mpi,out_process,numProc);

      //Procファイル出力
      //WriteProcDfiFile(out_proc_name[i],out_domain,out_mpi,out_process);
      WriteProcDfiFile(m_param->m_dfiList[i].out_proc_name,
                       out_domain,out_mpi,out_process);
    }

  }

  return true;

}

// #################################################################
//
bool convMxM::mxmsolv(std::string dfiname,
                      cdm_DFI* dfi,
                      int l_step,
                      double l_time,
                      int RankID,
                      double* min,
                      double* max)
{

  const cdm_Process* DFI_Process = dfi->GetcdmProcess();
  cdm_Domain* DFI_Domain = (cdm_Domain *)dfi->GetcdmDomain();
  const cdm_MPI* DFI_MPI = dfi->GetcdmMPI();
  const cdm_FileInfo* DFI_FInfo = dfi->GetcdmFileInfo();
  const cdm_TimeSlice* TSlice = dfi->GetcdmTimeSlice();

  bool mio = false;
  if( DFI_MPI->NumberOfRank > 1) mio=true;

  //間引き数のセット
  int thin_count = m_param->Get_ThinOut();

  //出力ガイドセルの設定
  int outGc=0;
  if( m_param->Get_OutputGuideCell() > 0 ) outGc = m_param->Get_OutputGuideCell();
  if( outGc > 0 ) {
    const cdm_FileInfo* DFI_FInfo = dfi->GetcdmFileInfo();
    if( outGc > DFI_FInfo->GuideCell ) outGc=DFI_FInfo->GuideCell;
  }
  //間引きありのとき、出力ガイドセルを0に設定
  if( thin_count > 1 ) outGc=0;
  //格子点出力のときガイドセルを0に設定
  if( m_param->Get_Interp_flag() ) outGc=0;

  //ピッチのセット
  double l_dpit[3];
  l_dpit[0]=DFI_Domain->GlobalRegion[0]/(double)DFI_Domain->GlobalVoxel[0];
  l_dpit[1]=DFI_Domain->GlobalRegion[1]/(double)DFI_Domain->GlobalVoxel[1];
  l_dpit[2]=DFI_Domain->GlobalRegion[2]/(double)DFI_Domain->GlobalVoxel[2];
  double out_dpit[3];
  for (int i=0;i<3;i++) out_dpit[i]=l_dpit[i]*double(thin_count);

  //全体のボクセルサイズを求める
  int voxel[3];
  for(int i=0; i<3; i++) voxel[i]=DFI_Domain->GlobalVoxel[i];

  //間引きを考慮したサイズのセット
  int l_imax_th = DFI_Process->RankList[RankID].VoxelSize[0]/thin_count;
  int l_jmax_th = DFI_Process->RankList[RankID].VoxelSize[1]/thin_count;
  int l_kmax_th = DFI_Process->RankList[RankID].VoxelSize[2]/thin_count;

  //間引き後のサイズが1つも無い領域のときエラー
  if( l_imax_th < 1 || l_jmax_th < 1 || l_kmax_th < 1 ) {
    printf("\toutput domain size error\n");
    return false;
  }

  if(DFI_Process->RankList[RankID].VoxelSize[0]%thin_count != 0) l_imax_th++;
  if(DFI_Process->RankList[RankID].VoxelSize[1]%thin_count != 0) l_jmax_th++;
  if(DFI_Process->RankList[RankID].VoxelSize[2]%thin_count != 0) l_kmax_th++;
  
  //間引き後のheadインデックスを求める
  int head[3];
  head[0] = (DFI_Process->RankList[RankID].HeadIndex[0]-1)/thin_count;
  head[1] = (DFI_Process->RankList[RankID].HeadIndex[1]-1)/thin_count;
  head[2] = (DFI_Process->RankList[RankID].HeadIndex[2]-1)/thin_count;
  if( (DFI_Process->RankList[RankID].HeadIndex[0]-1)%thin_count != 0 ) head[0]++;
  if( (DFI_Process->RankList[RankID].HeadIndex[1]-1)%thin_count != 0 ) head[1]++;
  if( (DFI_Process->RankList[RankID].HeadIndex[2]-1)%thin_count != 0 ) head[2]++;
  //間引き後のオリジンを求める
  double l_dorg[3];
  l_dorg[0]= DFI_Domain->GlobalOrigin[0];
  l_dorg[1]= DFI_Domain->GlobalOrigin[1];
  l_dorg[2]= DFI_Domain->GlobalOrigin[2];
  //自ノードの計算領域に合わせて原点シフト
  //(SPH,BOV形式のみ。他の形式については、各形式のクラス内で実施。)
  if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_SPH || 
      m_param->Get_OutputFormat() == CDM::E_CDM_FMT_BOV ) {
    l_dorg[0] += head[0]*out_dpit[0];
    l_dorg[1] += head[1]*out_dpit[1];
    l_dorg[2] += head[2]*out_dpit[2];
  }

  //出力タイプのセット
  CDM::E_CDM_DTYPE d_type;
  if( m_param->Get_OutputDataType() == CDM::E_CDM_DTYPE_UNKNOWN )
  {
    d_type = dfi->GetDataType();
  } else {
    d_type = m_param->Get_OutputDataType();
  }
  
  //出力バッファのインスタンス
  int szS[3];
  szS[0]=l_imax_th;
  szS[1]=l_jmax_th;
  szS[2]=l_kmax_th;
  cdm_Array* src = cdm_Array::instanceArray
                 ( d_type
                 , m_param->Get_OutputArrayShape()
                 , szS
                 , outGc
                 , dfi->GetNumVariables() );

  //読込みファイル名の生成
  std::string inPath = CDM::cdmPath_DirName(dfiname);
  std::string infile =  CDM::cdmPath_ConnectPath(inPath,dfi->Generate_FieldFileName(
                        RankID,l_step,mio));

  //ファイルの読込み
  unsigned int avr_step;
  double l_dtime, avr_time;
  CDM::E_CDM_ERRORCODE ret;
  int read_sta[3],read_end[3];
  for(int i=0; i<3; i++) {
    read_sta[i]=DFI_Process->RankList[RankID].HeadIndex[i]-outGc;
    read_end[i]=DFI_Process->RankList[RankID].TailIndex[i]+outGc;
  }

  cdm_Array* buf = dfi->ReadFieldData(infile, l_step, l_dtime,
                                      read_sta,
                                      read_end,
                                      DFI_Process->RankList[RankID].HeadIndex,
                                      DFI_Process->RankList[RankID].TailIndex,
                                      true, avr_step, avr_time, ret);
  if( ret != CDM::E_CDM_SUCCESS ) return false;

  //間引き及び型変換がない場合
  if( thin_count == 1 && buf->getDataType() == src->getDataType() && 
      buf->getArrayShape() == src->getArrayShape() ) {
    src = buf;
  } else {
  //間引きまたは型変換がある場合
    int headS[3],tailS[3];
    for(int i=0; i<3; i++) {
      headS[i]=DFI_Process->RankList[RankID].HeadIndex[i]-1-outGc;
      tailS[i]=DFI_Process->RankList[RankID].TailIndex[i]-1+outGc;
    }
    buf->setHeadIndex( headS );
    src->setHeadIndex( head );

    for(int n=0; n<dfi->GetNumVariables(); n++) convertXY(buf,src,headS,tailS,n);
    //delete buf;
  }

  //出力DFIのインスタンス
  int tail[3];
  head[0]=head[0]+1;
  head[1]=head[1]+1;
  head[2]=head[2]+1;
  tail[0]=head[0]+l_imax_th-1;
  tail[1]=head[1]+l_jmax_th-1;
  tail[2]=head[2]+l_kmax_th-1;

  cdm_DFI *out_dfi = NULL;
  if( DFI_FInfo->DFIType == CDM::E_CDM_DFITYPE_CARTESIAN )
  {
    //等間隔格子の場合
    out_dfi = cdm_DFI::WriteInit<double>(MPI_COMM_WORLD,
                                         "",
                                         m_param->Get_OutputDir(),
                                         DFI_FInfo->Prefix,
                                         m_param->Get_OutputFormat(),
                                         outGc,
                                         d_type,
                                         DFI_FInfo->NumVariables,
                                         "",
                                         voxel,
                                         out_dpit,
                                         l_dorg,
                                         DFI_Domain->GlobalDivision,
                                         head,
                                         tail,
                                         m_HostName,
                                         CDM::E_CDM_OFF);
  }
  else if( DFI_FInfo->DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN )
  {
    //不等間隔格子の場合
    int head_org[3],tail_org[3]; //間引く前のhead,tail
    for(int i=0; i<3; i++) {
      head_org[i] = DFI_Process->RankList[RankID].HeadIndex[i];
      tail_org[i] = DFI_Process->RankList[RankID].TailIndex[i];
    }

    if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT32 )
    {
      float *coord_X = NULL;
      float *coord_Y = NULL;
      float *coord_Z = NULL;

      //全計算領域のサイズの配列を用意し、head,tailで自ランクの座標データをセット
      coord_X = new float[voxel[0]+1]; //+1はセル数ではなく格子数のため。
      coord_Y = new float[voxel[1]+1];
      coord_Z = new float[voxel[2]+1];

      //配列(coord_X,coord_Y,coord_Z)に値をセット
      /* ＜注意点＞
       * 間引き時(thin_count>1)は、ガイドセルなし(outGc=0)
       * 間引かない時(thin_count=1)は、ガイドセルあり(outGc>0)で、head,tailは、head_org,tail_orgと一致
       */
      //x
      for(int ni=-outGc; ni<(l_imax_th+outGc); ni++) {
        if ( (ni+head[0]-1 >= 0) && (ni+head[0]-1 < voxel[0]+1) ) {
          coord_X[ni+head[0]-1] = (float)(DFI_Domain->NodeX(head_org[0]-1+ni*thin_count));
        }
      }
      if ( tail[0]+outGc < voxel[0]+1 ) {
        coord_X[tail[0]+outGc] = (float)(DFI_Domain->NodeX(tail_org[0]+outGc));
      }
      //y
      for(int nj=-outGc; nj<(l_jmax_th+outGc); nj++) {
        if ( (nj+head[1]-1 >= 0) && (nj+head[1]-1 < voxel[1]+1) ) {
          coord_Y[nj+head[1]-1] = (float)(DFI_Domain->NodeY(head_org[1]-1+nj*thin_count));
        }
      }
      if ( tail[1]+outGc < voxel[1]+1 ) {
        coord_Y[tail[1]+outGc] = (float)(DFI_Domain->NodeY(tail_org[1]+outGc));
      }
      //z
      for(int nk=-outGc; nk<(l_kmax_th+outGc); nk++) {
        if ( (nk+head[2]-1 >= 0) && (nk+head[2]-1 < voxel[2]+1) ) {
          coord_Z[nk+head[2]-1] = (float)(DFI_Domain->NodeZ(head_org[2]-1+nk*thin_count));
        }
      }
      if ( tail[2]+outGc < voxel[2]+1 ) {
        coord_Z[tail[2]+outGc] = (float)(DFI_Domain->NodeZ(tail_org[2]+outGc));
      }

      out_dfi = cdm_DFI::WriteInit<float>(MPI_COMM_WORLD,
                                          "",
                                          m_param->Get_OutputDir(),
                                          DFI_FInfo->Prefix,
                                          m_param->Get_OutputFormat(),
                                          outGc,
                                          d_type,
                                          DFI_FInfo->NumVariables,
                                          "",
                                          voxel,
                                          coord_X,
                                          coord_Y,
                                          coord_Z,
                                          DFI_Domain->GetCoordinateFile(),
                                          DFI_Domain->GetCoordinateFileType(),
                                          DFI_Domain->GetCoordinateFileEndian(),
                                          DFI_Domain->GlobalDivision,
                                          head,
                                          tail,
                                          m_HostName,
                                          CDM::E_CDM_OFF);
    }
    else if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT64 )
    {
      double *coord_X = NULL;
      double *coord_Y = NULL;
      double *coord_Z = NULL;

      //全計算領域のサイズの配列を用意し、head,tailで自ランクの座標データをセット
      coord_X = new double[voxel[0]+1]; //+1はセル数ではなく格子数のため。
      coord_Y = new double[voxel[1]+1];
      coord_Z = new double[voxel[2]+1];

      //配列(coord_X,coord_Y,coord_Z)に値をセット
      /* ＜注意点＞
       * 間引き時(thin_count>1)は、ガイドセルなし(outGc=0)
       * 間引かない時(thin_count=1)は、ガイドセルあり(outGc>0)で、head,tailは、head_org,tail_orgと一致
       */
      //x
      for(int ni=-outGc; ni<(l_imax_th+outGc); ni++) {
        if ( (ni+head[0]-1 >= 0) && (ni+head[0]-1 < voxel[0]+1) ) {
          coord_X[ni+head[0]-1] = (double)(DFI_Domain->NodeX(head_org[0]-1+ni*thin_count));
        }
      }
      if ( tail[0]+outGc < voxel[0]+1 ) {
        coord_X[tail[0]+outGc] = (double)(DFI_Domain->NodeX(tail_org[0]+outGc));
      }
      //y
      for(int nj=-outGc; nj<(l_jmax_th+outGc); nj++) {
        if ( (nj+head[1]-1 >= 0) && (nj+head[1]-1 < voxel[1]+1) ) {
          coord_Y[nj+head[1]-1] = (double)(DFI_Domain->NodeY(head_org[1]-1+nj*thin_count));
        }
      }
      if ( tail[1]+outGc < voxel[1]+1 ) {
        coord_Y[tail[1]+outGc] = (double)(DFI_Domain->NodeY(tail_org[1]+outGc));
      }
      //z
      for(int nk=-outGc; nk<(l_kmax_th+outGc); nk++) {
        if ( (nk+head[2]-1 >= 0) && (nk+head[2]-1 < voxel[2]+1) ) {
          coord_Z[nk+head[2]-1] = (double)(DFI_Domain->NodeZ(head_org[2]-1+nk*thin_count));
        }
      }
      if ( tail[2]+outGc < voxel[2]+1 ) {
        coord_Z[tail[2]+outGc] = (double)(DFI_Domain->NodeZ(tail_org[2]+outGc));
      }

      out_dfi = cdm_DFI::WriteInit<double>(MPI_COMM_WORLD,
                                           "",
                                           m_param->Get_OutputDir(),
                                           DFI_FInfo->Prefix,
                                           m_param->Get_OutputFormat(),
                                           outGc,
                                           d_type,
                                           DFI_FInfo->NumVariables,
                                           "",
                                           voxel,
                                           coord_X,
                                           coord_Y,
                                           coord_Z,
                                           DFI_Domain->GetCoordinateFile(),
                                           DFI_Domain->GetCoordinateFileType(),
                                           DFI_Domain->GetCoordinateFileEndian(),
                                           DFI_Domain->GlobalDivision,
                                           head,
                                           tail,
                                           m_HostName,
                                           CDM::E_CDM_OFF);
    }
  }
  if( out_dfi == NULL ) {
    printf("\tFails to instance dfi\n");
    return false;
  }

  out_dfi->set_RankID(RankID);
  out_dfi->SetcdmMPI(*DFI_MPI);

  //cdm_Processの作成＆更新
  cdm_Process out_Process;
  cdm_Rank out_Rank;
  for(int i=0; i<DFI_Process->RankList.size(); i++) {
    out_Rank.RankID = DFI_Process->RankList[i].RankID;
    out_Rank.HostName = "";
    for(int j=0; j<3; j++) {
      out_Rank.HeadIndex[j]=head[j];
      out_Rank.TailIndex[j]=tail[j];
      out_Rank.VoxelSize[j]=tail[j]-head[j]+1;
    }
    out_Process.RankList.push_back(out_Rank);
  }

  //out_dfi->SetcdmProcess(*DFI_Process);
  out_dfi->SetcdmProcess(out_Process);
  out_dfi->SetcdmTimeSlice(*TSlice);

  //出力形式（ascii,binary,Fbinary)のセット
  out_dfi->set_output_type(m_param->Get_OutputFileType());

  //節点への補間フラグのセット(AVSおよびVTK形式)
  if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_AVS || 
      m_param->Get_OutputFormat() == CDM::E_CDM_FMT_VTK ) {
    out_dfi->set_interp_flag(m_param->Get_Interp_flag());
  }

  //座標データの出力形式のセット(AVS形式)
  if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_AVS ) {
    out_dfi->set_output_type_coord(m_param->Get_OutputFileTypeCoord());
  }

  //gridファイルを出力(PLOT3D形式，iblankはすべて1にセット)
  if (m_param->Get_OutputFormat() == CDM::E_CDM_FMT_PLOT3D) {
    size_t size_ib=(l_imax_th+2*outGc)*(l_jmax_th+2*outGc)*(l_kmax_th+2*outGc);
    int *iblank;
    iblank = new int[size_ib];
    for(int i=0; i<size_ib; i++) {
      iblank[i] = 1;
    }
    out_dfi->WriteGridFile(iblank);
    delete [] iblank;
  }

  //フィールドデータ出力
  CDM::E_CDM_OUTPUT_FNAME output_fname = m_param->Get_OutputFilenameFormat();
  out_dfi->set_output_fname(output_fname);
  double tmp_minmax[8];
  unsigned idummy=0;
  double ddummy=0.0;
  ret = out_dfi->WriteData(
                           (unsigned)l_step,
                           //0,
                           outGc,
                           l_time,
                           src,
                           tmp_minmax,
                           true,
                           idummy,
                           ddummy);

  if( ret != CDM::E_CDM_SUCCESS ) return false;
 
  //minmaxを求める
  if( !DtypeMinMax(src,min,max) ) return false;

  //delete
  delete out_dfi;
  delete src;

  return true;
}
