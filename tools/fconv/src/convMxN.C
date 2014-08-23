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
 * @file   convMxN.C 
 * @brief  convMxN Class
 * @author aics
 */

#include "convMxN.h"

// #################################################################
// コンストラクタ
convMxN::convMxN()
{

}

// #################################################################
// デストラクタ
convMxN::~convMxN()
{
  for(int i=0; i<m_out_dfi.size(); i++ ) if( !m_out_dfi[i] ) delete m_out_dfi[i];
}

// #################################################################
// 領域分割と出力DFIのインスタンス
void convMxN::VoxelInit()
{
  std::string outdfiname;
  int iret=0;
  const cio_Domain *DFI_Domain = m_in_dfi[0]->GetcioDomain();

  //ピッチのセット
  double dfi_pit[3];
  for(int i=0; i<3; i++) dfi_pit[i]=DFI_Domain->GlobalRegion[i]/(double)DFI_Domain->GlobalVoxel[i];

  //入力領域指示のセット
  int IndexStart[3];
  int IndexEnd[3];
  for(int i=0; i<3; i++) IndexStart[i]=1;
  for(int i=0; i<3; i++) IndexEnd[i]=DFI_Domain->GlobalVoxel[i];
  if( m_param->Get_CropIndexStart_on() ) {
    const int *cropIndexStart = m_param->Get_CropIndexStart();
    for(int i=0; i<3; i++) IndexStart[i]=cropIndexStart[i];
  }
  if( m_param->Get_CropIndexEnd_on() ) {
    const int *cropIndexEnd = m_param->Get_CropIndexEnd();
    for(int i=0; i<3; i++) IndexEnd[i]=cropIndexEnd[i];
  }

  //全体のサイズをセット
  int voxel[3];
  for(int i=0; i<3; i++) voxel[i]=IndexEnd[i]-IndexStart[i]+1;

  //リージョンのセット
  double region[3];
  for(int i=0; i<3; i++) region[i]=voxel[i]*dfi_pit[i];

  //出力領域の分割数の取得
  int* Gdiv = m_param->Get_OutputDivision();

  if( Gdiv[0]>0 && Gdiv[1]>0 && Gdiv[2]>0 ) {
    //分割数が指示されている場合
    /*
    iret = m_paraMngr->VoxelInit(Gdiv, (int *)DFI_Domain->GlobalVoxel,
                              (double *)DFI_Domain->GlobalOrigin,
                              (double *)DFI_Domain->GlobalRegion, 0, 0);
    */
    iret = m_paraMngr->VoxelInit(Gdiv, voxel,
                                 (double *)DFI_Domain->GlobalOrigin,
                                 region, 0, 0);
    if( iret != 0 ) {
       printf("\tVoxelInit Error cpm_ErrorCode : %d\n",iret);
       Exit(0);
    }
  } else {
    //分割数が指示されていない場合
    /*
    iret =  m_paraMngr->VoxelInit((int *)DFI_Domain->GlobalVoxel,
                               (double *)DFI_Domain->GlobalOrigin,
                               (double *)DFI_Domain->GlobalRegion, 0, 0);
    */
    iret =  m_paraMngr->VoxelInit(voxel,
                                  (double *)DFI_Domain->GlobalOrigin,
                                  region, 0, 0);
    if( iret != 0 ) {
       printf("\tVoxelInit Error cpm_ErrorCode : %d\n",iret);
       Exit(0);
    }
  }


  const int* tmp;
  /*
  tmp = m_paraMngr->GetGlobalVoxelSize();
  for(int i=0; i<3; i++) m_Gvoxel[i]=tmp[i];
  */
  for(int i=0; i<3; i++) m_Gvoxel[i]=DFI_Domain->GlobalVoxel[i];
  
  tmp = m_paraMngr->GetVoxelHeadIndex();
  for(int i=0; i<3; i++) m_Head[i]=tmp[i]+1;
  tmp = m_paraMngr->GetVoxelTailIndex();
  for(int i=0; i<3; i++) m_Tail[i]=tmp[i]+1;
  tmp = m_paraMngr->GetDivNum();
  for(int i=0; i<3; i++) m_Gdiv[i]=tmp[i];

  //入力指示を考慮したヘッド、テイルに更新
 
  tmp = m_paraMngr->GetLocalVoxelSize();
 
  for(int i=0; i<3; i++) {
    //if( m_Head[i] < IndexStart[i] ) m_Head[i]=IndexStart[i];
    if( m_Head[i] < IndexStart[i] ) {
      m_Head[i]=m_Head[i]+IndexStart[i]-1;
      m_Tail[i]=m_Head[i]+tmp[i]-1;
    }
    if( m_Tail[i] > IndexEnd[i]   ) m_Tail[i]=IndexEnd[i];
  }

  //間引き数のセット
  int thin_count = m_param->Get_ThinOut();

  //間引きを考慮した全体サイズのセット
  int voxel_thin[3];
  for(int i=0; i<3; i++) {
    voxel_thin[i]=voxel[i]/thin_count;
    if( voxel[i]%thin_count != 0 ) voxel_thin[i]++;
  }


  //間引きを考慮したヘッド、テイルインデックスの作成
  int head[3],tail[3];
  for(int i=0; i<3; i++) {
    head[i]=(m_Head[i]-1)/thin_count;
    if( (m_Head[i]-1)%thin_count != 0 ) head[i]++;
    tail[i]=(m_Tail[i]-1)/thin_count;
  }

  const double* dtmp;
  double pit[3],org[3];
  //dtmp = m_paraMngr->GetPitch();
  //for(int i=0; i<3; i++) pit[i]=dtmp[i];
  //for(int i=0; i<3; i++) pit[i]=dtmp[i]*double(thin_count);
  pit[0]=region[0]/voxel_thin[0];
  pit[1]=region[1]/voxel_thin[1];
  pit[2]=region[2]/voxel_thin[2];

  dtmp = m_paraMngr->GetGlobalOrigin();
  for(int i=0; i<3; i++) org[i]=dtmp[i];

  const int *tmp_head = m_paraMngr->GetVoxelHeadIndex();
  const int *tmp_tail = m_paraMngr->GetVoxelTailIndex();
  for(int i=0; i<3; i++) {
    head[i]=tmp_head[i]/thin_count;
    if( tmp_head[i]%thin_count != 0 ) head[i]++;
    tail[i]=tmp_tail[i]/thin_count;
  }

  for(int i=0; i<3; i++) org[i]+=double(head[i])*pit[i];

  for(int i=0; i<3; i++) {
    head[i]=head[i]+1;
    tail[i]=tail[i]+1;
  }

  //出力DFIの初期化
  for(int i=0; i<m_in_dfi.size(); i++) {
    const cio_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcioFileInfo();

    std::string outdfifname="";
    std::string outprocfname="";
    if( m_param->Get_Outputdfi_on() ) {
       outdfifname =m_param->m_dfiList[i].out_dfi_name;
       outprocfname=m_param->m_dfiList[i].out_proc_name;
    }

    //出力タイプのセット
    CIO::E_CIO_DTYPE d_type;
    if( m_param->Get_OutputDataType() == CIO::E_CIO_DTYPE_UNKNOWN )
    {
      d_type = m_in_dfi[i]->GetDataType();
    } else {
      d_type = m_param->Get_OutputDataType();
    }

    //出力ガイドセルの設定
    int outGc=0;
    if( m_param->Get_OutputGuideCell() > 1 ) outGc = m_param->Get_OutputGuideCell();
    if( outGc > 1 ) {
      const cio_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcioFileInfo();
      if( outGc > DFI_FInfo->GuideCell ) outGc=DFI_FInfo->GuideCell;
    }
    //間引きありのとき、出力ガイドセルを0に設定
    if( thin_count > 1 ) outGc=0;
    //格子点出力のとき、出力ガイドセルを0に設定
    if( m_bgrid_interp_flag ) outGc=0; 

    cio_DFI* dfi=cio_DFI::WriteInit(MPI_COMM_WORLD,
                          outdfifname,
                          m_param->Get_OutputDir(),
                          DFI_FInfo->Prefix,
                          m_param->Get_OutputFormat(),
                          //0,
                          outGc,
                          d_type,
                          m_param->Get_OutputArrayShape(),
                          DFI_FInfo->Component,
                          outprocfname,
                          voxel_thin,
                          pit,
                          org,
                          m_Gdiv,
                          head,
                          tail,
                          m_HostName,
                          CIO::E_CIO_OFF);
    if( dfi == NULL ) {
      printf("\tFails to instance dfi\n");
      Exit(0);
    }

    //Procファイル出力
    if( m_param->Get_Outputdfi_on() ) dfi->WriteProcDfiFile(MPI_COMM_WORLD,false);

    //出力形式（ascii,binary,Fbinary)のセット
    dfi->set_output_type(m_param->Get_OutputFormatType());

    //Unitのセット
    std::string unit;
    double ref;
    double diff;
    bool bdiff;
    m_in_dfi[i]->GetUnit("Length",unit,ref,diff,bdiff);
    dfi->AddUnit("Length",unit,ref,diff,bdiff);
    m_in_dfi[i]->GetUnit("Velocity",unit,ref,diff,bdiff);
    dfi->AddUnit("Velocity",unit,ref,diff,bdiff);
    m_in_dfi[i]->GetUnit("Pressure",unit,ref,diff,bdiff);
    dfi->AddUnit("Pressure",unit,ref,diff,bdiff);

    //成分名の取り出しとセット
    for(int n=0; n<DFI_FInfo->Component; n++) {
       std::string variable = m_in_dfi[i]->getComponentVariable(n);
       if( variable != "" ) dfi->setComponentVariable(n,variable);
    } 

    m_out_dfi.push_back(dfi);
  }

}

// #################################################################
//
bool convMxN::exec()
{

  if( m_myRank == 0 ) {
    printf("Convert M x N\n");
  }


  // 出力ファイル形式クラスのインスタンス
  convOutput *ConvOut = convOutput::OutputInit(m_param->Get_OutputFormat());

  // InputParamのインスタンス
  if( !ConvOut->importInputParam(m_param) ) {
    //Exit(0);
    return false;
  }

  //出力ファイル名の取得
  //vector<std::string> out_dfi_name = m_InputCntl->Get_OutdfiNameList();
  std::string prefix,outfile;

  FILE *fp;
  int dummy;

  CIO::E_CIO_DTYPE d_type;

  CIO::E_CIO_ERRORCODE ret; 
  double rtime;
  unsigned idummy;
  double ddummy;
  float fminmax[8];
  double dminmax[8];

  bool mio;
  mio = false;
  if( m_numProc > 1 ) mio=true;

  //間引き数のセット
  int thin_count = m_param->Get_ThinOut();

  //入力領域指示のセット
  int IndexStart[3];
  int IndexEnd[3];
  for(int i=0; i<3; i++) IndexStart[i]=m_Head[i];
  for(int i=0; i<3; i++) IndexEnd[i]=m_Tail[i];
  if( m_param->Get_CropIndexStart_on() ) {
    const int *cropIndexStart = m_param->Get_CropIndexStart();
    for(int i=0; i<3; i++) {
      if( IndexStart[i] < cropIndexStart[i] ) IndexStart[i]=cropIndexStart[i];
    }
  }
  if( m_param->Get_CropIndexEnd_on() ) {
    const int *cropIndexEnd = m_param->Get_CropIndexEnd();
    for(int i=0; i<3; i++) {
      if( IndexEnd[i] > cropIndexEnd[i] ) IndexEnd[i]=cropIndexEnd[i];
    }
  }

  //自ノードのボクセルサイズの取得
  int sz[3];
  const int* tmp = m_paraMngr->GetLocalVoxelSize();
  for(int i=0; i<3; i++) sz[i]=tmp[i];

  //自ノードのボクセルサイズを入力指示を考慮して更新
  for(int i=0; i<3; i++) {
    if( sz[i] > (IndexEnd[i]-IndexStart[i]+1) ) sz[i]=(IndexEnd[i]-IndexStart[i]+1);
  }

  //出力workareaのサイズ
  int szS[3];
  const int *cropIndexStart = m_param->Get_CropIndexStart();
  for(int i=0; i<3; i++) {
    szS[i]=sz[i]/thin_count; 
    if( szS[i] < 1 ) {
      printf("\toutput domain size error\n");
      return false;
    }
    if( m_param->Get_CropIndexStart_on() ) {
      if( IndexStart[i] == cropIndexStart[i] ) {
        if( sz[i]%thin_count != 0 ) szS[i]++;
      }
    } else {
      if( sz[i]%thin_count != 0 ) szS[i]++;
    }
  }

  int head[3],tail[3];
  for(int i=0; i<3; i++) {
    head[i]=(m_Head[i]-1)/thin_count;
    if( (m_Head[i]-1)%thin_count != 0 ) head[i]++;
    tail[i]=(m_Tail[i]-1)/thin_count;
  }
  const double* dtmp;
  double pit[3],org[3];
  dtmp = m_paraMngr->GetPitch();
  for(int i=0; i<3; i++) pit[i]=dtmp[i]*double(thin_count);
  dtmp = m_paraMngr->GetGlobalOrigin();
  for(int i=0; i<3; i++) org[i]=dtmp[i]+0.5*pit[i];
  for(int i=0; i<3; i++) org[i]+=double(head[i])*pit[i];

  const cio_FileInfo* DFI_FInfo = m_in_dfi[0]->GetcioFileInfo();

  //dfiのループ
  for (int i=0; i<m_in_dfi.size(); i++) {


    int nComp = m_in_dfi[i]->GetNumComponent();

    int outGc=0;
    if( m_param->Get_OutputGuideCell() > 1 ) outGc = m_param->Get_OutputGuideCell();
    if( outGc > 0 ) {
      const cio_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcioFileInfo();
      if( outGc > DFI_FInfo->GuideCell ) outGc = DFI_FInfo->GuideCell;
    }

    if( thin_count > 1 ) outGc=0;
    if( m_bgrid_interp_flag ) outGc=0; 

    //読込みバッファのインスタンス
    cio_Array* buf = cio_Array::instanceArray
    ( m_in_dfi[i]->GetDataType(),
      m_in_dfi[i]->GetArrayShape(),
      sz,
      //0,
      outGc,
      //m_in_dfi[i]->GetNumComponent());
      nComp);

    //出力タイプのセット
    if( m_param->Get_OutputDataType() == CIO::E_CIO_DTYPE_UNKNOWN )
    {
      d_type = m_in_dfi[i]->GetDataType();
    } else {
      d_type = m_param->Get_OutputDataType();
    }

    //出力バッファのインスタンス
    cio_Array* src = cio_Array::instanceArray
    ( d_type,
      //m_in_dfi[i]->GetArrayShape(),
      m_param->Get_OutputArrayShape(),
      szS,
      //0,
      outGc,
      //m_in_dfi[i]->GetNumComponent());
      nComp);
   
    //DFI_FInfoクラスの取得
    const cio_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcioFileInfo();
    prefix=DFI_FInfo->Prefix; 

    //TimeSliceクラスの取得
    const cio_TimeSlice* TSlice = m_in_dfi[i]->GetcioTimeSlice();


    //ステップ数のループ
    for ( int j=0; j<TSlice->SliceList.size(); j++ ) {

      //MxNの読込み
      ret = m_in_dfi[i]->ReadData(buf,
                                (unsigned)TSlice->SliceList[j].step,
                                //0,
                                outGc,
                                m_Gvoxel,
                                m_Gdiv,
                                m_Head,
                                m_Tail,
                                rtime,
                                true,
                                idummy,
                                ddummy);
      if( ret != CIO::E_CIO_SUCCESS ) {
        printf("ReadData Error\n");
        return false;
      }

      //読込みバッファのheadIndexのセット
      int headB[3];
      for(int k=0; k<3; k++) headB[k]=m_Head[k]-1;
      buf->setHeadIndex( headB );

      //間引き及び型変換がない場合
      if( thin_count == 1 && buf->getDataType() == src->getDataType() &&
          buf->getArrayShape() == src->getArrayShape() ) {
        src=buf;
      } else {
      //間引きまたは型変換がある場合
        //出力バッファの間引きなしでのHeadIndex,TailIndex
        int headS[3],tailS[3];
        for(int k=0; k<3; k++) {
          headS[k]=m_Head[k]-1;
          tailS[k]=m_Tail[k]-1;
        }
        //出力バッファのHeadIndexセット
        int headS0[3];
        if( m_param->Get_CropIndexStart_on() ) {
          for(int k=0; k<3; k++) {
            headS0[k]=headS[k]/thin_count;
          }
        } else {
          for(int k=0; k<3; k++) {
            headS0[k]=headS[k]/thin_count;
            if( headS[k]%thin_count != 0 ) headS0[k]++;
          }
        }

        src->setHeadIndex( headS0 );

        for(int n=0; n<nComp; n++) convertXY(buf,src,headS,tailS,n);
      }

      CIO::E_CIO_OUTPUT_FNAME output_fname = m_param->Get_OutputFilenameFormat();
      m_out_dfi[i]->set_output_fname(output_fname);

      //minmaxの初期化
      int nsize = nComp;
      if( nComp > 1 ) nsize++;
      double *min = new double[nsize];
      double *max = new double[nsize];
      for(int n=0; n<nsize; n++) {
        min[n]=DBL_MAX;
        max[n]=-DBL_MAX;
      }
      //minmaxを求める
      if( !DtypeMinMax(src,min,max) ) return false;

      //if( out_dfi_name.size() > 1 ) {
      if( m_param->Get_Outputdfi_on() ) {
        //ランク間で通信してMINMAXを求めてランク０に送信
        int nbuff = nsize*1;
        //minの通信
        double *send1 = min;
        double *recv1 = new double[nbuff];
        MPI_Reduce(send1, recv1, nbuff, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        min = recv1;
        //maxの通信
        double *send2 = max;
        double *recv2 = new double[nbuff];
        MPI_Reduce(send2, recv2, nbuff, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        max = recv2;
      } 

      //出力処理
      double *tmp_minmax = new double[nsize*2];
      for(int n=0; n<nsize; n++ ) {
        tmp_minmax[n*2+0] = min[n];
        tmp_minmax[n*2+1] = max[n];
      }
  
      m_out_dfi[i]->SetcioTimeSlice(*TSlice);
 
      ret = m_out_dfi[i]->WriteData(
                                    (unsigned)TSlice->SliceList[j].step,
                                    //0,
                                    outGc,
                                    rtime,
                                    src,
                                    tmp_minmax,
                                    true,
                                    idummy,
                                    ddummy);


    } 
    delete src;
  }

  return true;

}

