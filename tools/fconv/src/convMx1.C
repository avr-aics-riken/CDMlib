/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   convMx1.C 
 * @brief  convMx1 Class
 * @author aics
 */

#include "convMx1.h"

// #################################################################
// コンストラクタ
convMx1::convMx1()
{

  m_StepRankList.clear();

}

// #################################################################
// デストラクタ
convMx1::~convMx1()
{

}

// #################################################################
//
bool convMx1::exec()
{

  if( m_myRank == 0 ) {
    printf("Convert M x 1\n");
  } 

  // 出力ファイル形式クラスのインスタンス
  ConvOut = convOutput::OutputInit(m_param->Get_OutputFormat());

  // InputParamのインスタンス
  if( !ConvOut->importInputParam(m_param) ) {
    return false;
  }

  string prefix,outfile,infile,inPath;
  cdm_FILE *pFile;
  int dummy;
  
  int l_rank;
  int l_d_type;
  CDM::E_CDM_DTYPE d_type;
  int l_step, l_imax, l_jmax, l_kmax;
  float l_time;
  double l_dorg[3], l_dpit[3];
  double l_dtime;
  int xsize,ysize,zsize,asize,vsize;
  int dim;
  
  int div[3];
  int l_imax_th, l_jmax_th, l_kmax_th;

  //間引き数のセット
  int thin_count = m_param->Get_ThinOut();
 
  // 入力モード
  bool mio = false;
  //const cdm_MPI* DFI_MPI = m_in_dfi[0]->GetcdmMPI();
  //if( DFI_MPI->NumberOfRank > 1) mio=true;

  headT mapHeadX;
  headT mapHeadY;
  headT mapHeadZ;

  //step基準のリスト生成
  makeStepList(m_StepRankList);

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

  //入力領域指示を考慮
  int IndexStart[3];
  int IndexEnd[3];
  const int *cropIndexStart = m_param->Get_CropIndexStart();
  const int *cropIndexEnd =   m_param->Get_CropIndexEnd();
  for(int i=0; i<3; i++ ) {
    IndexStart[i]=cropIndexStart[i];
    IndexEnd[i]=cropIndexEnd[i];
  }
  //dfi*stepのループ 
  for (int i=0;i<m_StepRankList.size();i++) {

    const cdm_Domain* DFI_Domain = m_StepRankList[i].dfi->GetcdmDomain();
    cdm_Process* DFI_Process = (cdm_Process *)m_StepRankList[i].dfi->GetcdmProcess();
    //全体サイズのキープ
    l_imax= DFI_Domain->GlobalVoxel[0];
    l_jmax= DFI_Domain->GlobalVoxel[1];
    l_kmax= DFI_Domain->GlobalVoxel[2];

    //オリジナルのピッチを計算
    double o_pit[3];
    for(int j=0; j<3; j++) o_pit[j]=DFI_Domain->GlobalRegion[j]/DFI_Domain->GlobalVoxel[j];

    //入力領域指示を考慮
    if( !m_param->Get_CropIndexStart_on() ) {
      IndexStart[0] = 1;
      IndexStart[1] = 1;
      IndexStart[2] = 1;
    }
    if( !m_param->Get_CropIndexEnd_on() ) {
      IndexEnd[0] = l_imax;
      IndexEnd[1] = l_jmax;
      IndexEnd[2] = l_kmax;
    }
    l_imax = IndexEnd[0]-IndexStart[0]+1;
    l_jmax = IndexEnd[1]-IndexStart[1]+1;
    l_kmax = IndexEnd[2]-IndexStart[2]+1;

    //入力領域を考慮したリージョンの作成
    double region[3];
    region[0]=l_imax*o_pit[0];
    region[1]=l_jmax*o_pit[1];
    region[2]=l_kmax*o_pit[2];
      
    //間引きを考慮
    l_imax_th=l_imax/thin_count;//間引き後のxサイズ
    l_jmax_th=l_jmax/thin_count;//間引き後のyサイズ
    l_kmax_th=l_kmax/thin_count;//間引き後のzサイズ
    if(l_imax%thin_count != 0) l_imax_th++;
    if(l_jmax%thin_count != 0) l_jmax_th++;
    if(l_kmax%thin_count != 0) l_kmax_th++;

    //オリジンとピッチを作成
    l_dpit[0]=region[0]/(double)l_imax_th;
    l_dpit[1]=region[1]/(double)l_jmax_th;
    l_dpit[2]=region[2]/(double)l_kmax_th;
    l_dorg[0]=DFI_Domain->GlobalOrigin[0];
    l_dorg[1]=DFI_Domain->GlobalOrigin[1];
    l_dorg[2]=DFI_Domain->GlobalOrigin[2];
    //格子点の原点からセル中心の原点へのシフト
    //(SPH,BOV形式のみ。他の形式については、各形式の出力時に実施。)
    if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_SPH || 
        m_param->Get_OutputFormat() == CDM::E_CDM_FMT_BOV ) {
      l_dorg[0] += 0.5*l_dpit[0];
      l_dorg[1] += 0.5*l_dpit[1];
      l_dorg[2] += 0.5*l_dpit[2];
    }

    //GRID データ 出力
    const cdm_FileInfo* DFI_FInfo = m_StepRankList[i].dfi->GetcdmFileInfo();
    int sz[3];
    sz[0]=l_imax;
    sz[1]=l_jmax;
    sz[2]=l_kmax;

    //出力ガイドセルの設定
    int outGc=0;
    if( m_param->Get_OutputGuideCell() > 0 ) outGc = m_param->Get_OutputGuideCell();
    if( outGc > 0 ) {
      if( outGc > DFI_FInfo->GuideCell ) outGc=DFI_FInfo->GuideCell;
    }
    if( thin_count > 1 || m_param->Get_Interp_flag() ) outGc=0;

    //格子データ出力用の情報を作成
    cdm_Domain* out_domain = NULL;
    cdm_Process* out_process = NULL;

    if( !makeGridInfo(m_StepRankList[i].dfi,out_domain,out_process,outGc) ) {
      printf("\tError : makeProcInfo\n");
      return false;
    }

    /*
    ConvOut->WriteGridData(DFI_FInfo->Prefix,0, m_myRank, m_in_dfi[0]->GetDataType(),
                           DFI_FInfo->GuideCell, l_dorg, l_dpit, sz);
    */
    ConvOut->WriteGridData(DFI_FInfo->Prefix,0, m_myRank, m_in_dfi[0]->GetDataType(),
                           outGc, out_domain, out_process);

    //dfiファイルのディレクトリの取得
    inPath = CDM::cdmPath_DirName(m_StepRankList[i].dfi->get_dfi_fname());

    //const cdm_FileInfo* DFI_FInfo = m_StepRankList[i].dfi->GetcdmFileInfo();
    prefix=DFI_FInfo->Prefix;
    LOG_OUTV_ fprintf(m_fplog,"  COMBINE SPH START : %s\n", prefix.c_str());
    STD_OUTV_ printf("  COMBINE SPH START : %s\n", prefix.c_str());
    
    //Number of Variables
    dim=m_StepRankList[i].dfi->GetNumVariables();

    const cdm_TimeSlice* TSlice = m_StepRankList[i].dfi->GetcdmTimeSlice();

    // convOutputに情報をセット(NetCDF4用)
    ConvOut->m_pTSlice = TSlice;
    ConvOut->m_pFinfo  = DFI_FInfo;
    ConvOut->m_pUnit   = m_in_dfi[i]->GetcdmUnit();

    div[0]=DFI_Domain->GlobalDivision[0];
    div[1]=DFI_Domain->GlobalDivision[1];
    div[2]=DFI_Domain->GlobalDivision[2];
    //mapHeadX,Y,Zの生成
    DFI_Process->CreateRankList(DFI_Domain,
                                mapHeadX,
                                mapHeadY,
                                mapHeadZ);

    //入力ファイルの並列フラグセット
    const cdm_MPI* DFI_MPI = m_StepRankList[i].dfi->GetcdmMPI();
    if( DFI_MPI->NumberOfRank > 1) mio=true;
    else mio=false;

    //step Loop
    for(int j=m_StepRankList[i].stepStart; j<=m_StepRankList[i].stepEnd; j++) {

      //minmaxの初期化
      int nsize = dim;
      if( dim > 1 ) nsize++;
      double *min = new double[nsize];
      double *max = new double[nsize];
      for(int n=0; n<nsize; n++) {
         min[n]=DBL_MAX;
         max[n]=-DBL_MAX;
      } 

      l_step=TSlice->SliceList[j].step;
      l_time=(float)TSlice->SliceList[j].time;
      
      LOG_OUTV_ fprintf(m_fplog,"\tstep = %d\n", l_step);
      STD_OUTV_ printf("\tstep = %d\n", l_step);
      
      //連結出力ファイルオープン
      pFile = ConvOut->OutputFile_Open(prefix, l_step, 0, false);

      //m_d_typeのセット (float or double)
      if( m_StepRankList[i].dfi->GetDataType() == CDM::E_CDM_FLOAT32 ) {
        l_d_type = SPH_FLOAT;
      } else if( m_StepRankList[i].dfi->GetDataType() == CDM::E_CDM_FLOAT64 ) {
        l_d_type = SPH_DOUBLE;
      }
      
      if( m_param->Get_OutputDataType() == CDM::E_CDM_DTYPE_UNKNOWN )
      {
        d_type = m_StepRankList[i].dfi->GetDataType();
      } else {
        d_type = m_param->Get_OutputDataType();
      }

      //出力ガイドセルによるオリジンの更新
      double t_org[3];
      for(int n=0; n<3; n++) t_org[n]=l_dorg[n];
      if( outGc > 0 ) {
        for(int n=0; n<3; n++) t_org[n]=t_org[n]-(double)outGc*l_dpit[n];
      }

      //ヘッダーレコードを出力
      if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_VTK ) {
        if( !(ConvOut->WriteHeaderRecord(l_step,
                                         l_time,
                                         dim,
                                         d_type,
                                         DFI_FInfo->DFIType,
                                         out_domain,
                                         out_process,
                                         outGc,
                                         prefix,
                                         pFile)) ) {
          printf("\twrite header error (VTK format)\n");
          return false;
        }
      } else {
        if( !(ConvOut->WriteHeaderRecord(l_step, dim, d_type, 
                                         l_imax_th+2*outGc, l_jmax_th+2*outGc, l_kmax_th+2*outGc,
                                         l_time, t_org, l_dpit, prefix, pFile)) ) {
          printf("\twrite header error\n");
          return false;
        }
      }

      //全体の大きさの計算とデータのヘッダ書き込み
      size_t dLen;
      bool flag_Mark = false;

      //dLen = size_t(l_imax_th) * size_t(l_jmax_th) * size_t(l_kmax_th);
      dLen = size_t(l_imax_th+2*outGc) * size_t(l_jmax_th+2*outGc) * size_t(l_kmax_th+2*outGc);
      if( dim >1 ) dLen *= (size_t)dim;
      //if( m_param->Get_OutputDataType() == CDM::E_CDM_FLOAT32 ) {
      if( d_type == CDM::E_CDM_FLOAT32 ) {
         dummy = dLen * sizeof(float);
      } else {
         dummy = dLen * sizeof(double);
      }
      if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_SPH ||
          m_param->Get_OutputFormat() == CDM::E_CDM_FMT_PLOT3D) {
        flag_Mark = true;
      }
      if( !(ConvOut->WriteDataMarker(dummy, pFile, flag_Mark)) ) {
        printf("\twrite data header error\n");
        return false;
      }
      
      //書き込みworkareaのサイズ決め
      xsize=l_imax_th;
      ysize=l_jmax_th;
      asize=xsize*ysize;
      vsize=0;
      for(int n=0; n< DFI_Process->RankList.size(); n++ ) {
        int szx,szy,szz;
        szx=DFI_Process->RankList[n].VoxelSize[0];
        szy=DFI_Process->RankList[n].VoxelSize[1];
        szz=DFI_Process->RankList[n].VoxelSize[2];
        int vdum=szx*szy*szz;
        if(vsize < vdum) vsize=vdum;
      }
      // メモリチェック
      LOG_OUTV_ fprintf(m_fplog,"\tNode %4d - Node %4d\n", 0,
                        DFI_Domain->GlobalVoxel[2] -1);
      STD_OUTV_ printf("\tNode %4d - Node %4d\n", 0,
                       DFI_Domain->GlobalVoxel[2] -1);
      double mc1 = (double)asize*(double)dim;
      double mc2 = (double)vsize*(double)dim;
      if(mc1>(double)INT_MAX){// 整数値あふれ出しチェック //参考 894*894*894*3=2143550952 INT_MAX 2147483647
        printf("\tsize error : mc1>INT_MAX\n");
        return false;
      }
      if(mc2>(double)INT_MAX){// 整数値あふれ出しチェック //参考 894*894*894*3=2143550952 INT_MAX 2147483647
        printf("\tsize error : mc2>INT_MAX\n");
        return false;
      }
      double TotalMemory=0.0; // = mc * (double)sizeof(REAL_TYPE);
      if( m_param->Get_OutputDataType() == CDM::E_CDM_FLOAT32 ) {
        TotalMemory = TotalMemory + mc1 * (double)sizeof(float);
      } else {
        TotalMemory = TotalMemory + mc1 * (double)sizeof(double);
      }

      if( l_d_type == SPH_FLOAT ) {
        TotalMemory = TotalMemory + mc2 * (double)sizeof(float);
      } else {
        TotalMemory = TotalMemory + mc2 * (double)sizeof(double);
      }
      LOG_OUT_ MemoryRequirement(TotalMemory,m_fplog);
      STD_OUT_ MemoryRequirement(TotalMemory,stdout);
      
      int szS[3];
      szS[0]=l_imax_th;
      szS[1]=l_jmax_th;
      szS[2]=1;
     
      CDM::E_CDM_ARRAYSHAPE output_AShape = m_param->Get_OutputArrayShape();

      if( (output_AShape == CDM::E_CDM_NIJK || DFI_FInfo->NumVariables == 1) &&
          (m_param->Get_OutputFormat() != CDM::E_CDM_FMT_VTK) ){  //VTKの出力はijknの方を使う

        //output nijk
        if( !convMx1_out_nijk(pFile,
                             inPath,
                             l_step,
                             l_dtime,
                             d_type,
                             mio,
                             div,
                             szS,
                             //m_stepList[i].dfi,
                             m_StepRankList[i].dfi,
                             DFI_Process,
                             mapHeadX,mapHeadY,mapHeadZ,
                             min,max
                             ) ) return false;

      } else {

        //output IJKN
        if( !convMx1_out_ijkn(pFile,
                             inPath,
                             l_step,
                             l_dtime,
                             d_type,
                             mio,
                             div,
                             szS,
                             //m_stepList[i].dfi,
                             m_StepRankList[i].dfi,
                             DFI_Process,
                             mapHeadX,mapHeadY,mapHeadZ,
                             min,max
                             ) ) return false;

      } 

      //dfiごとにminmaxを登録
      for(int ndfi = 0; ndfi<minmaxList.size(); ndfi++) {
        if( minmaxList[ndfi]->dfi != m_StepRankList[i].dfi ) continue;
        for(int n=0; n<nsize; n++) {
          if( minmaxList[ndfi]->Min[j*nsize+n] > min[n] ) minmaxList[ndfi]->Min[j*nsize+n] = min[n];
          if( minmaxList[ndfi]->Max[j*nsize+n] < max[n] ) minmaxList[ndfi]->Max[j*nsize+n] = max[n];
        }
      }

      //データのフッタ書き込み
      if( !(ConvOut->WriteDataMarker(dummy, pFile, flag_Mark)) ) {
        printf("\twrite data error\n");
        return false;
      }
      
      //出力ファイルクローズ
      ConvOut->OutputFile_Close(pFile);

    }

    //avsのヘッダーファイル出力
    if( i==0 ){
      ConvOut->output_avs(m_myRank,
                          m_in_dfi,
                          out_domain,
                          out_process,
                          outGc);
    }

    //クラスout_domain,out_processの解放
    delete out_domain;
    delete out_process;

  }

  //avsのヘッダーファイル出力
  //ConvOut->output_avs(m_myRank, m_in_dfi);

  //出力dfiファイル名の取得
  //vector<std::string> out_dfi_name = m_InputCntl->Get_OutdfiNameList();
  //vector<std::string> out_proc_name = m_InputCntl->Get_OutprocNameList();

  //出力dfiファイルの出力
  //if( out_dfi_name.size() == 0 || out_proc_name.size() == 0 ) return true;
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

  if( m_myRank == 0 ) {

    WriteIndexDfiFile(minmaxList);

    for(int i=0; i<m_in_dfi.size(); i++) {
      cdm_Domain* out_domain = NULL;
      cdm_MPI* out_mpi = NULL;
      cdm_Process* out_process = NULL;
      const cdm_MPI* dfi_mpi = m_in_dfi[i]->GetcdmMPI();
      int numProc = dfi_mpi->NumberOfRank;

      //Proc情報の生成
      makeProcInfo(m_in_dfi[i],out_domain,out_mpi,out_process,1);

      //Procファイル出力
      //WriteProcDfiFile(out_proc_name[i],out_domain,out_mpi,out_process);
      WriteProcDfiFile(m_param->m_dfiList[i].out_proc_name,
                       out_domain,out_mpi,out_process);
     }
  }
  
  return true;

}

// #################################################################
// NIJK to NIJK
bool
convMx1::convMx1_out_nijk(cdm_FILE* pFile,
                          std::string inPath,
                          int l_step,
                          double l_dtime,
                          CDM::E_CDM_DTYPE d_type,
                          bool mio,
                          int div[3],
                          int sz[3],
                          cdm_DFI* dfi, 
                          cdm_Process* DFI_Process,
                          headT mapHeadX, headT mapHeadY, headT mapHeadZ,
                          double* min, double* max)
{

  //cdm_Domain* DFI_Domain = (cdm_Domain *)m_in_dfi[0]->GetcdmDomain();
  cdm_Domain* DFI_Domain = (cdm_Domain *)dfi->GetcdmDomain();

  int thin_count = m_param->Get_ThinOut();

  int outGc=0;
  int interp_Gc=0;

  //出力ガイドセルの設定
  if( m_param->Get_OutputGuideCell() > 0 ) outGc = m_param->Get_OutputGuideCell();
  if( outGc > 0 ) {
    const cdm_FileInfo* DFI_FInfo = dfi->GetcdmFileInfo();
    if( outGc > DFI_FInfo->GuideCell ) outGc=DFI_FInfo->GuideCell;
  }

  //間引きありのとき、出力ガイドセルを0に設定
  if( thin_count > 1 ) outGc=0;
  interp_Gc = outGc;

  //格子点出力のときガイドセルが0のとき1にセット
  if( m_param->Get_Interp_flag() && outGc==0 ) interp_Gc=1;

  //cell出力のとき、出力ガイドセルを0に設定
  if( !m_param->Get_Interp_flag() ) interp_Gc=0;

  //出力のヘッダー、フッターをセット
  int headS[3],tailS[3];
  const int* CorpIndexStart = m_param->Get_CropIndexStart();
  const int* CorpIndexEnd = m_param->Get_CropIndexEnd();
  int IndexStart[3],IndexEnd[3];
  for(int i=0;i<3; i++) {
    IndexStart[i]=CorpIndexStart[i]-outGc;
    IndexEnd[i]=CorpIndexEnd[i]+outGc;
  }
  if( !m_param->Get_CropIndexStart_on() ) {
    IndexStart[0]=1-outGc;
    IndexStart[1]=1-outGc;
    IndexStart[2]=1-outGc;
  }
  if( !m_param->Get_CropIndexEnd_on() ) {
    IndexEnd[0]=DFI_Domain->GlobalVoxel[0]+outGc;
    IndexEnd[1]=DFI_Domain->GlobalVoxel[1]+outGc;
    IndexEnd[2]=DFI_Domain->GlobalVoxel[2]+outGc;
  }
  
  headS[0]=IndexStart[0]-1;
  headS[1]=IndexStart[1]-1;
  tailS[0]=IndexEnd[0]-1;
  tailS[1]=IndexEnd[1]-1;

  //変数の個数の取り出し
  int nVari = dfi->GetNumVariables();

  //配列形状の設定
  CDM::E_CDM_ARRAYSHAPE out_shape;
  if( nVari == 1 ) out_shape = CDM::E_CDM_IJKN;
  else if( nVari > 1 ) out_shape = CDM::E_CDM_NIJK;

  //セル中心出力のときガイドセル数を考慮してサイズ更新
  if( !m_param->Get_Interp_flag() ) {
    sz[0]=sz[0]+2*outGc;
    sz[1]=sz[1]+2*outGc;
  }

  //出力バッファのインスタンス
  cdm_Array* src = cdm_Array::instanceArray
                   ( d_type
                   //, dfi->GetArrayShape()
                   , out_shape
                   , sz
                   , interp_Gc
                   , nVari );

  //補間用バッファ,格子点出力バッファのインスタンス
  cdm_Array* src_old = NULL;
  cdm_Array* outArray = NULL;
  if( m_param->Get_Interp_flag() ) {
    src_old = cdm_Array::instanceArray
              ( d_type
              //, dfi->GetArrayShape()
              , out_shape
              , sz
              , interp_Gc
              , nVari );

    int szOut[3];
    for(int i=0; i<2; i++) szOut[i]=sz[i]+1;
    szOut[2]=sz[2];
    outArray = cdm_Array::instanceArray
             ( d_type
              //, dfi->GetArrayShape()
              , out_shape
              , szOut
              , interp_Gc
              , nVari );
  }

  int kdiv,jdiv,idiv;
  int l_rank;
  std::string infile;

  //z方向の分割数回のループ
  for( headT::iterator itz=mapHeadZ.begin(); itz!= mapHeadZ.end(); itz++ ) {

    //z層のスタートエンドを設定
    kdiv = itz->second;
    int kp_sta,kp_end;
    kp_sta = itz->first;
    int nrank = _CDM_IDX_IJK(0,0,kdiv,div[0],div[1],div[2],0);
    kp_end = kp_sta + DFI_Process->RankList[nrank].VoxelSize[2];

    //z層のスタートエンドをガイドセルの考慮
    if( kdiv == 0 ) kp_sta = kp_sta-outGc;
    if( kdiv == div[2]-1 ) kp_end = kp_end+outGc;

    //同一Z面のループ
    for(int kp=kp_sta; kp< kp_end; kp++) {

      //入力領域外のときスキップ
      if( kp < IndexStart[2] || kp > IndexEnd[2] ) continue; 

      int kk = kp-1;
      //間引きの層のときスキップ
      //if( kk%thin_count != 0 ) continue;
      if( (kp-IndexStart[2])%thin_count != 0 ) continue;

      //y方向の分割数のループ
      for( headT::iterator ity=mapHeadY.begin(); ity!= mapHeadY.end(); ity++ ) {

        //yのスタートエンドの設定
        jdiv = ity->second;
        int jp_sta,jp_end;
        jp_sta = ity->first;
        int nrank = _CDM_IDX_IJK(0,jdiv,kdiv,div[0],div[1],div[2],0);
        jp_end = jp_sta + DFI_Process->RankList[nrank].VoxelSize[1];

        //x方向の分割数のループ
        for( headT::iterator itx=mapHeadX.begin(); itx!= mapHeadX.end(); itx++ ) {

          //xのスタートエンドの設定
          idiv = itx->second;
          int ip_sta,ip_end;
          ip_sta = itx->first;
          int nrank = _CDM_IDX_IJK(idiv,jdiv,kdiv,div[0],div[1],div[2],0);
          ip_end = ip_sta + DFI_Process->RankList[nrank].VoxelSize[0];

          int RankID = _CDM_IDX_IJK(idiv,jdiv,kdiv,div[0],div[1],div[2],0);

          //読込み範囲の設定
          if( IndexStart[0] > ip_end || IndexEnd[0] < ip_sta ) continue;
          if( IndexStart[1] > jp_end || IndexEnd[1] < jp_sta ) continue;

          int read_sta[3],read_end[3];
          read_sta[0]=ip_sta;
          if( IndexStart[0] > ip_sta ) read_sta[0] = IndexStart[0];
          read_sta[1]=jp_sta;
          if( IndexStart[1] > jp_sta ) read_sta[1] = IndexStart[1];
          read_sta[2]=kp;
          read_end[0]=ip_end-1;
          if( IndexEnd[0] < ip_end ) read_end[0] = IndexEnd[0];
          read_end[1]=jp_end-1;
          if( IndexEnd[1] < jp_end ) read_end[1] = IndexEnd[1];
          read_end[2]=kp;

          //ガイドセルを考慮して読込み範囲を更新
          if( idiv == 0 ) read_sta[0] = read_sta[0]-outGc;
          if( idiv == div[0]-1 ) read_end[0] = read_end[0]+outGc;
          if( jdiv == 0 ) read_sta[1] = read_sta[1]-outGc;
          if( jdiv == div[1]-1 ) read_end[1] = read_end[1]+outGc;

          l_rank=DFI_Process->RankList[RankID].RankID;
          //連結対象ファイル名の生成
          infile = CDM::cdmPath_ConnectPath(inPath,dfi->Generate_FieldFileName(l_rank,l_step,mio));

          unsigned int avr_step;
          double avr_time;
          CDM::E_CDM_ERRORCODE ret;
          //連結対象ファイルの読込み
          cdm_Array* buf = dfi->ReadFieldData(infile, l_step, l_dtime,
                                              read_sta, read_end,
                                              DFI_Process->RankList[RankID].HeadIndex,
                                              DFI_Process->RankList[RankID].TailIndex,
                                              true, avr_step, avr_time, ret);
          if( ret != CDM::E_CDM_SUCCESS ) {
            printf("\tCan't Read Field Data Record %s\n",infile.c_str());
            return false;
          } 
          //headIndexを０スタートにしてセット
          int headB[3];
          headB[0]=read_sta[0]-1;
          headB[1]=read_sta[1]-1;
          headB[2]=read_sta[2]-1;
          buf->setHeadIndex( headB );

          int headS0[3];
          //headS0[0]=headS[0];
          //headS0[1]=headS[1];
          headS0[0]=headS[0]/thin_count;
          headS0[1]=headS[1]/thin_count;
          headS0[2]=kk/thin_count;
          src->setHeadIndex( headS0 );

          headS0[0]=headS[0];
          headS0[1]=headS[1];
          headS0[2]=kk;
          tailS[2]=headS0[2];

          //出力配列へのコンバイン
          for(int n=0; n<nVari; n++) convertXY(buf,src,headS0,tailS,n);
          delete buf;

        } /// Loop itx
      } /// Loop ity
      //補間処理
      if( m_param->Get_Interp_flag() ) {
        if( kp == kp_sta ) {
          for(int n=0; n<nVari; n++) {
            if( !InterPolate(src,src,outArray,n,n) ) return false;
          }
        } else {
          for(int n=0; n<nVari; n++) {
            if( !InterPolate(src_old,src,outArray,n,n) ) return false;
          }
        }
      } else {
        outArray = src;
        //outArrayのz方向のheadを0にセット
        const int* head_src = src->getHeadIndex();
        int head_out[3];
        head_out[0] = head_src[0];
        head_out[1] = head_src[1];
        head_out[2] = 0;
        outArray->setHeadIndex( head_out );
      }

      //一層分出力
      if( outArray ) {
        const int* szOutArray = outArray->getArraySizeInt();
        size_t dLen = szOutArray[0]*szOutArray[1]*szOutArray[2]*outArray->getNvari();
        if( ConvOut->WriteFieldData(pFile,
                                    outArray,
                                    dLen ) != true ) return false;

      }

      //minmaxを求める
      if( !DtypeMinMax(outArray,min,max) ) return false;  

      //補間ありのとき、読込んだ層の配列ポインタをsrc_oldにコピー
      if( m_param->Get_Interp_flag() ) {
        cdm_Array* tmp = src;
        src = src_old;
        src_old = tmp;
      }
    } /// Loop kp
  } /// Loop itz

  if( m_param->Get_Interp_flag() ) {
    for(int n=0; n<nVari; n++) {
      if( !InterPolate(src_old,src_old,outArray,n,n) ) return false;
    }
    if( outArray ) {
      const int* szOutArray = outArray->getArraySizeInt();
      size_t dLen = szOutArray[0]*szOutArray[1]*szOutArray[2]*outArray->getNvari();
      if( ConvOut->WriteFieldData(pFile,
                                  outArray,
                                  dLen ) != true ) return false;

       //minmaxを求める
       if( !DtypeMinMax(src,min,max) ) return false;

    } else return false;
  }
  delete src;

  if( m_param->Get_Interp_flag() ) {
    delete src_old;
    delete outArray;
  }

  return true;

}

// #################################################################
// NIJK to IJKN または IJKN to IJKN
bool
convMx1::convMx1_out_ijkn(cdm_FILE* pFile,
                          std::string inPath,
                          int l_step,
                          double l_dtime,
                          CDM::E_CDM_DTYPE d_type,
                          bool mio,
                          int div[3],
                          int sz[3],
                          cdm_DFI* dfi,
                          cdm_Process* DFI_Process,
                          headT mapHeadX, headT mapHeadY, headT mapHeadZ,
                          double* min, double* max)
{

  //cdm_Domain* DFI_Domain = (cdm_Domain *)m_in_dfi[0]->GetcdmDomain();
  cdm_Domain* DFI_Domain = (cdm_Domain *)dfi->GetcdmDomain();
  const cdm_FileInfo* DFI_FInfo = dfi->GetcdmFileInfo();

  int thin_count = m_param->Get_ThinOut();

  int outGc=0;
  int interp_Gc=0;

  //出力ガイドセルの設定
  if( m_param->Get_OutputGuideCell() > 0 ) outGc = m_param->Get_OutputGuideCell();
  if( outGc > 0 ) {
    if( outGc > DFI_FInfo->GuideCell ) outGc=DFI_FInfo->GuideCell;
  }

  //間引きありのとき、出力ガイドセルを0に設定
  if( thin_count > 1 ) outGc=0;
  interp_Gc = outGc;

  //格子点出力のときガイドセルが0のとき1にセット
  if( m_param->Get_Interp_flag() && outGc==0 ) interp_Gc=1;

  //cell出力のとき、出力ガイドセルを0に設定
  if( !m_param->Get_Interp_flag() ) interp_Gc=0;

  //出力のヘッダー、フッターをセット
  int headS[3],tailS[3];
  const int* CorpIndexStart = m_param->Get_CropIndexStart();
  const int* CorpIndexEnd = m_param->Get_CropIndexEnd();
  int IndexStart[3],IndexEnd[3];
  for(int i=0;i<3; i++) {
    IndexStart[i]=CorpIndexStart[i]-outGc;
    IndexEnd[i]=CorpIndexEnd[i]+outGc;
  }
  if( !m_param->Get_CropIndexStart_on() ) {
    IndexStart[0]=1-outGc;
    IndexStart[1]=1-outGc;
    IndexStart[2]=1-outGc;
  }
  if( !m_param->Get_CropIndexEnd_on() ) {
    IndexEnd[0]=DFI_Domain->GlobalVoxel[0]+outGc;
    IndexEnd[1]=DFI_Domain->GlobalVoxel[1]+outGc;
    IndexEnd[2]=DFI_Domain->GlobalVoxel[2]+outGc;
  }

  headS[0]=IndexStart[0]-1;
  headS[1]=IndexStart[1]-1;
  tailS[0]=IndexEnd[0]-1;
  tailS[1]=IndexEnd[1]-1;

  //変数の個数の取り出し
  int nVari = dfi->GetNumVariables();

  //セル中心出力のときガイドセル数を考慮してサイズ更新
  if( !m_param->Get_Interp_flag() ) {
    sz[0]=sz[0]+2*outGc;
    sz[1]=sz[1]+2*outGc;
  }

  //出力バッファのインスタンス(読込み配列形状でのDFIでインスタンス）
  cdm_Array* src = cdm_Array::instanceArray
                   ( d_type
                   , dfi->GetArrayShape()
                   , sz
                   , interp_Gc
                   , nVari );

  //補間用バッファ（読込み配列形状でのDFIでインスタンス）
  cdm_Array* src_old = NULL;
  cdm_Array* outArray = NULL;
  if( m_param->Get_Interp_flag() ) {
    src_old = cdm_Array::instanceArray
              ( d_type
              , dfi->GetArrayShape()
              , sz
              , interp_Gc
              , nVari );

    int szOut[3];
    for(int i=0; i<2; i++) szOut[i]=sz[i]+1;
    szOut[2]=sz[2];
    outArray = cdm_Array::instanceArray
             ( d_type
              , dfi->GetArrayShape()
              //, CDM::E_CDM_IJKN,
              , szOut
              , interp_Gc
              , 1 );
  }

  int kdiv,jdiv,idiv;
  int l_rank;
  std::string infile;

  //変数の個数のループ
  for(int n=0; n<nVari; n++) {

    bool flag_variname = true; //VTK形式における変数名出力フラグ
    std::string variname = DFI_FInfo->VariableName[n];

    //z方向の分割数回のループ
    for( headT::iterator itz=mapHeadZ.begin(); itz!= mapHeadZ.end(); itz++ ) {

      //z層のスタートエンドを設定
      kdiv = itz->second;
      int kp_sta,kp_end;
      kp_sta = itz->first;
      int nrank = _CDM_IDX_IJK(0,0,kdiv,div[0],div[1],div[2],0);
      kp_end = kp_sta + DFI_Process->RankList[nrank].VoxelSize[2];

      //z層のスタートエンドをガイドセルの考慮
      if( kdiv == 0 ) kp_sta = kp_sta-outGc;
      if( kdiv == div[2]-1 ) kp_end = kp_end+outGc;

      //同一Z面のループ
      for(int kp=kp_sta; kp< kp_end; kp++) {

        //入力領域外のときスキップ
        if( kp < IndexStart[2] || kp > IndexEnd[2] ) continue;

        int kk = kp-1;
        //間引きの層のときスキップ
        //if( kk%thin_count != 0 ) continue;
        if( (kp-IndexStart[2])%thin_count != 0 ) continue;

        //y方向の分割数のループ
        for( headT::iterator ity=mapHeadY.begin(); ity!= mapHeadY.end(); ity++ ) {

          //yのスタートエンドの設定
          jdiv = ity->second;
          int jp_sta,jp_end;
          jp_sta = ity->first;
          int nrank = _CDM_IDX_IJK(0,jdiv,kdiv,div[0],div[1],div[2],0);
          jp_end = jp_sta + DFI_Process->RankList[nrank].VoxelSize[1];

          //x方向の分割数のループ
          for( headT::iterator itx=mapHeadX.begin(); itx!= mapHeadX.end(); itx++ ) {

            //xのスタートエンドの設定
            idiv = itx->second;
            int ip_sta,ip_end;
            ip_sta = itx->first;
            int nrank = _CDM_IDX_IJK(idiv,jdiv,kdiv,div[0],div[1],div[2],0);
            ip_end = ip_sta + DFI_Process->RankList[nrank].VoxelSize[0];

            int RankID = _CDM_IDX_IJK(idiv,jdiv,kdiv,div[0],div[1],div[2],0);

            //読込み範囲の設定
            if( IndexStart[0] > ip_end || IndexEnd[0] < ip_sta ) continue;
            if( IndexStart[1] > jp_end || IndexEnd[1] < jp_sta ) continue;

            int read_sta[3],read_end[3];
            read_sta[0]=ip_sta;
            if( IndexStart[0] > ip_sta ) read_sta[0] = IndexStart[0];
            read_sta[1]=jp_sta;
            if( IndexStart[1] > jp_sta ) read_sta[1] = IndexStart[1];
            read_sta[2]=kp;
            read_end[0]=ip_end-1;
            if( IndexEnd[0] < ip_end ) read_end[0] = IndexEnd[0];
            read_end[1]=jp_end-1;
            if( IndexEnd[1] < jp_end ) read_end[1] = IndexEnd[1];
            read_end[2]=kp;

            //ガイドセルを考慮して読込み範囲を更新
            if( idiv == 0 ) read_sta[0] = read_sta[0]-outGc;
            if( idiv == div[0]-1 ) read_end[0] = read_end[0]+outGc;
            if( jdiv == 0 ) read_sta[1] = read_sta[1]-outGc;
            if( jdiv == div[1]-1 ) read_end[1] = read_end[1]+outGc;

            l_rank=DFI_Process->RankList[RankID].RankID;
            //連結対象ファイル名の生成
            infile = CDM::cdmPath_ConnectPath(inPath,dfi->Generate_FieldFileName(l_rank,l_step,mio));
            unsigned int avr_step;
            double avr_time;
            CDM::E_CDM_ERRORCODE ret;
            //連結対象ファイルの読込み
            cdm_Array* buf = dfi->ReadFieldData(infile, l_step, l_dtime,
                                                read_sta, read_end,
                                                DFI_Process->RankList[RankID].HeadIndex,
                                                DFI_Process->RankList[RankID].TailIndex,
                                                true, avr_step, avr_time, ret);

            if( ret != CDM::E_CDM_SUCCESS ) {
              printf("\tCan't Read Field Data Record %s\n",infile.c_str());
              return false;
            } 
            //headIndexを０スタートにしてセット
            int headB[3];
            headB[0]=read_sta[0]-1;
            headB[1]=read_sta[1]-1;
            headB[2]=read_sta[2]-1;
            buf->setHeadIndex( headB );

            int headS0[3];
            headS0[0]=headS[0];
            headS0[1]=headS[1];
            headS0[2]=kk/thin_count;
            src->setHeadIndex( headS0 );

            headS0[2]=kk;
            tailS[2]=headS0[2];

            //出力配列へのコンバイン
            convertXY(buf,src,headS0,tailS,n);

            //minmaxを求める
            if( n==0 ) if( !DtypeMinMax(buf,min,max) ) return false;

            delete buf;

          } /// Loop itx
        } /// Loop ity 
        //補間処理
        if( m_param->Get_Interp_flag() ) {
          if( kp == kp_sta ) {
            if( !InterPolate(src,src,outArray,n,0) ) return false;
          } else {
            if( !InterPolate(src_old,src,outArray,n,0) ) return false;
          }
        } else {
          //NIJK レコードをIJKにコピー
          outArray = nijk_to_ijk(src,n);
          outArray->setHeadIndex((int*)src->getHeadIndex());
        }

        //一層分出力
        if( outArray ) {
          const int* szOutArray = outArray->getArraySizeInt();
          size_t dLen = szOutArray[0]*szOutArray[1]*szOutArray[2]*outArray->getNvari();
          if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_VTK ||
              m_param->Get_OutputFormat() == CDM::E_CDM_FMT_NETCDF4 ) {
            if( ConvOut->WriteFieldData(pFile,
                                        outArray,
                                        dLen,
                                        d_type,
                                        flag_variname,
                                        variname) != true ) return false;
            flag_variname = false;
          } else {
            if( ConvOut->WriteFieldData(pFile,
                                        outArray,
                                        dLen ) != true ) return false;
          }
        }
        //補間ありのとき、読込んだ層の配列ポインタをsrc_oldにコピー
        if( m_param->Get_Interp_flag() ) {
          cdm_Array* tmp = src;
          src = src_old;
          src_old = tmp;
        }
      } ///Loop kp
    } ///Loop itz

    if( m_param->Get_Interp_flag() ) {
      for(int n=0; n<nVari; n++) {
        if( !InterPolate(src_old,src_old,outArray,n,0) ) return false;
      }
      if( outArray ) {
        const int* szOutArray = outArray->getArraySizeInt();
        size_t dLen = szOutArray[0]*szOutArray[1]*szOutArray[2]*outArray->getNvari();
        if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_VTK ) {
          if( ConvOut->WriteFieldData(pFile,
                                      outArray,
                                      dLen,
                                      d_type,
                                      flag_variname,
                                      variname) != true ) return false;
          flag_variname = false;
        } else {
          if( ConvOut->WriteFieldData(pFile,
                                      outArray,
                                      dLen ) != true ) return false;
        }
      } else return false;
    }
  } ///Loop n

  delete src;
  if( m_param->Get_Interp_flag() ) {
    delete src_old;
    delete outArray;
  }

  return true;

}

// #################################################################
// 補間処理
bool
convMx1::InterPolate(cdm_Array* src_old, cdm_Array* src, cdm_Array* outArray,
                     int ivar_src, int ivar_out )
{

  if( !src_old || !src || !outArray ) return false;
  //if( !src_old || !src ) return NULL;

  //データタイプの取得
  //int nVari = src->getNvari();
  CDM::E_CDM_DTYPE dtype = src->getDataType();

  //char
  if( dtype == CDM::E_CDM_INT8 ) {
    cdm_TypeArray<char> *O = dynamic_cast<cdm_TypeArray<char>*>(outArray);
    cdm_TypeArray<char> *S = dynamic_cast<cdm_TypeArray<char>*>(src);
    cdm_TypeArray<char> *S_old = dynamic_cast<cdm_TypeArray<char>*>(src_old);

    //足しこみ領域のゼロクリア
    zeroClearArray(O,ivar_out);
    //srcの足しこみ
    setGridData_XY(O,S,    ivar_out,ivar_src);
    //src_oldの足しこみ
    setGridData_XY(O,S_old,ivar_out,ivar_src);
    //平均化（８で割る）
    VolumeDataDivide8(O,ivar_out);
  }
  //short
  else if( dtype == CDM::E_CDM_INT16 ) {
    cdm_TypeArray<short> *O = dynamic_cast<cdm_TypeArray<short>*>(outArray);
    cdm_TypeArray<short> *S = dynamic_cast<cdm_TypeArray<short>*>(src);
    cdm_TypeArray<short> *S_old = dynamic_cast<cdm_TypeArray<short>*>(src_old);

    //足しこみ領域のゼロクリア
    zeroClearArray(O,ivar_out);
    //srcの足しこみ
    setGridData_XY(O,S,    ivar_out,ivar_src);
    //src_oldの足しこみ
    setGridData_XY(O,S_old,ivar_out,ivar_src);
    //平均化（８で割る）
    VolumeDataDivide8(O,ivar_out);
  }
  //int
  else if( dtype == CDM::E_CDM_INT32 ) {
    cdm_TypeArray<int> *O = dynamic_cast<cdm_TypeArray<int>*>(outArray);
    cdm_TypeArray<int> *S = dynamic_cast<cdm_TypeArray<int>*>(src);
    cdm_TypeArray<int> *S_old = dynamic_cast<cdm_TypeArray<int>*>(src_old);

    //足しこみ領域のゼロクリア
    zeroClearArray(O,ivar_out);
    //srcの足しこみ
    setGridData_XY(O,S,    ivar_out,ivar_src);
    //src_oldの足しこみ
    setGridData_XY(O,S_old,ivar_out,ivar_src);
    //平均化（８で割る）
    VolumeDataDivide8(O,ivar_out);
  }
  //float
  else if( dtype == CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *O = dynamic_cast<cdm_TypeArray<float>*>(outArray);
    cdm_TypeArray<float> *S = dynamic_cast<cdm_TypeArray<float>*>(src);
    cdm_TypeArray<float> *S_old = dynamic_cast<cdm_TypeArray<float>*>(src_old);

    //足しこみ領域のゼロクリア
    zeroClearArray(O,ivar_out);
    //srcの足しこみ
    setGridData_XY(O,S,    ivar_out,ivar_src);
    //src_oldの足しこみ
    setGridData_XY(O,S_old,ivar_out,ivar_src);

    //平均化（８で割る）
    VolumeDataDivide8(O,ivar_out);
  }
  //double
  else if( dtype == CDM::E_CDM_FLOAT64 ) {
    cdm_TypeArray<double> *O = dynamic_cast<cdm_TypeArray<double>*>(outArray);
    cdm_TypeArray<double> *S = dynamic_cast<cdm_TypeArray<double>*>(src);
    cdm_TypeArray<double> *S_old = dynamic_cast<cdm_TypeArray<double>*>(src_old);

    //足しこみ領域のゼロクリア
    zeroClearArray(O,ivar_out);
    //srcの足しこみ
    setGridData_XY(O,S,    ivar_out,ivar_src);
    //src_oldの足しこみ
    setGridData_XY(O,S_old,ivar_out,ivar_src);
    //平均化（８で割る）
    VolumeDataDivide8(O,ivar_out);
  }

  return outArray;

}

// #################################################################
// NIJK配列をスカラーのIJK配列にコピー
cdm_Array*
convMx1::nijk_to_ijk(cdm_Array* src, int ivar)
{

  //コピー元配列のサイズとデータタイプの取得
  const int *sz = src->getArraySizeInt();
  CDM::E_CDM_DTYPE d_type = src->getDataType();

  cdm_Array* outArray = cdm_Array::instanceArray
                       ( d_type
                       , CDM::E_CDM_IJKN
                       , (int *)sz
                       , 0 
                       , 1 );
  //unsigned char
  if( d_type == CDM::E_CDM_UINT8 ) {
    cdm_TypeArray<unsigned char> *S = dynamic_cast<cdm_TypeArray<unsigned char>*>(src);
    cdm_TypeArray<unsigned char> *O = dynamic_cast<cdm_TypeArray<unsigned char>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //char
  else if( d_type == CDM::E_CDM_INT8 ) {
    cdm_TypeArray<char> *S = dynamic_cast<cdm_TypeArray<char>*>(src);
    cdm_TypeArray<char> *O = dynamic_cast<cdm_TypeArray<char>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //unsigned short 
  else if( d_type == CDM::E_CDM_UINT16 ) {
    cdm_TypeArray<unsigned short> *S = dynamic_cast<cdm_TypeArray<unsigned short>*>(src);
    cdm_TypeArray<unsigned short> *O = dynamic_cast<cdm_TypeArray<unsigned short>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //short
  else if( d_type == CDM::E_CDM_INT16 ) {
    cdm_TypeArray<short> *S = dynamic_cast<cdm_TypeArray<short>*>(src);
    cdm_TypeArray<short> *O = dynamic_cast<cdm_TypeArray<short>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //unsigned int 
  else if( d_type == CDM::E_CDM_UINT32 ) {
    cdm_TypeArray<unsigned int> *S = dynamic_cast<cdm_TypeArray<unsigned int>*>(src);
    cdm_TypeArray<unsigned int> *O = dynamic_cast<cdm_TypeArray<unsigned int>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //int
  else if( d_type == CDM::E_CDM_INT32 ) {
    cdm_TypeArray<int> *S = dynamic_cast<cdm_TypeArray<int>*>(src);
    cdm_TypeArray<int> *O = dynamic_cast<cdm_TypeArray<int>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //unsigned long 
  else if( d_type == CDM::E_CDM_UINT64 ) {
    cdm_TypeArray<unsigned long long> *S = dynamic_cast<cdm_TypeArray<unsigned long long>*>(src);
    cdm_TypeArray<unsigned long long> *O = dynamic_cast<cdm_TypeArray<unsigned long long>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //long
  else if( d_type == CDM::E_CDM_INT64 ) {
    cdm_TypeArray<long long> *S = dynamic_cast<cdm_TypeArray<long long>*>(src);
    cdm_TypeArray<long long> *O = dynamic_cast<cdm_TypeArray<long long>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //float
  else if( d_type == CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *S = dynamic_cast<cdm_TypeArray<float>*>(src);
    cdm_TypeArray<float> *O = dynamic_cast<cdm_TypeArray<float>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }
  //double
  else if( d_type == CDM::E_CDM_FLOAT64 ) {
    cdm_TypeArray<double> *S = dynamic_cast<cdm_TypeArray<double>*>(src);
    cdm_TypeArray<double> *O = dynamic_cast<cdm_TypeArray<double>*>(outArray);
    copyArray_nijk_ijk(S,O,ivar);
  }

  return outArray;
}

