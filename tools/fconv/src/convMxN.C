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
 * @file   convMxN.C 
 * @brief  convMxN Class
 * @author aics
 */

#include "convMxN.h"
#include "cdm_NonUniformDomain.h"

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
  const cdm_Domain *DFI_Domain = m_in_dfi[0]->GetcdmDomain();

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
//20160512.fub.s
/*
  for(int i=0; i<3; i++) {
    region[i] = DFI_Domain->NodeX(IndexEnd[i]) - DFI_Domain->NodeX(IndexStart[i]-1);
  }
*/
  region[0] = DFI_Domain->NodeX(IndexEnd[0]) - DFI_Domain->NodeX(IndexStart[0]-1);
  region[1] = DFI_Domain->NodeY(IndexEnd[1]) - DFI_Domain->NodeY(IndexStart[1]-1);
  region[2] = DFI_Domain->NodeZ(IndexEnd[2]) - DFI_Domain->NodeZ(IndexStart[2]-1);
//20160512.fub.e

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

  //自ノードの計算領域に合わせて原点シフト
  for(int i=0; i<3; i++) org[i]+=double(head[i])*pit[i];

  for(int i=0; i<3; i++) {
    head[i]=head[i]+1;
    tail[i]=tail[i]+1;
  }

  //出力DFIの初期化
  for(int i=0; i<m_in_dfi.size(); i++) {
//20160425.fub.s
    SetPrefixFileInfo(m_in_dfi[i],i);
//20160425.fub.e
    const cdm_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcdmFileInfo();

    std::string outdfifname="";
    std::string outprocfname="";
    if( m_param->Get_Outputdfi_on() ) {
       outdfifname =m_param->m_dfiList[i].out_dfi_name;
       outprocfname=m_param->m_dfiList[i].out_proc_name;
    }

    //出力タイプのセット
    CDM::E_CDM_DTYPE d_type;
    if( m_param->Get_OutputDataType() == CDM::E_CDM_DTYPE_UNKNOWN )
    {
      d_type = m_in_dfi[i]->GetDataType();
    } else {
      d_type = m_param->Get_OutputDataType();
    }

    //出力ガイドセルの設定
    int outGc=0;
    if( m_param->Get_OutputGuideCell() > 0 ) outGc = m_param->Get_OutputGuideCell();
    if( outGc > 0 ) {
      const cdm_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcdmFileInfo();
      if( outGc > DFI_FInfo->GuideCell ) outGc=DFI_FInfo->GuideCell;
    }
    //間引きありのとき、出力ガイドセルを0に設定
    if( thin_count > 1 ) outGc=0;
    //格子点出力のとき、出力ガイドセルを0に設定
    if( m_param->Get_Interp_flag() ) outGc=0; 

    cdm_DFI *dfi = NULL;
    if( DFI_FInfo->DFIType == CDM::E_CDM_DFITYPE_CARTESIAN )
    {
      //等間隔格子の場合
      dfi=cdm_DFI::WriteInit<double>(MPI_COMM_WORLD,
                                     outdfifname,
                                     m_param->Get_OutputDir(),
                                     DFI_FInfo->Prefix,
                                     m_param->Get_OutputFormat(),
                                     outGc,
                                     d_type,
                                     DFI_FInfo->NumVariables,
                                     outprocfname,
                                     voxel_thin,
                                     pit,
                                     org,
                                     m_Gdiv,
                                     head,
                                     tail,
                                     m_HostName,
                                     CDM::E_CDM_OFF);
    }
    else if( DFI_FInfo->DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN )
    {
      //不等間隔格子の場合
      if( DFI_Domain->GetCoordinateFilePrecision() == CDM::E_CDM_FLOAT32 )
      {
        float *coord_X = NULL;
        float *coord_Y = NULL;
        float *coord_Z = NULL;

        //全計算領域の座標をWriteInitに渡す
        coord_X = new float[voxel_thin[0]+1]; //+1はセル数ではなく格子数のため。
        coord_Y = new float[voxel_thin[1]+1];
        coord_Z = new float[voxel_thin[2]+1];

        //配列(coord_X,coord_Y,coord_Z)に値をセット
        //x
        for(int ni=0; ni<voxel_thin[0]; ni++) {
          coord_X[ni] = (float)(DFI_Domain->NodeX(ni*thin_count));
        }
        coord_X[voxel_thin[0]] = (float)(DFI_Domain->NodeX(DFI_Domain->GlobalVoxel[0]));
        //y
        for(int nj=0; nj<voxel_thin[1]; nj++) {
          coord_Y[nj] = (float)(DFI_Domain->NodeY(nj*thin_count));
        }
        coord_Y[voxel_thin[1]] = (float)(DFI_Domain->NodeY(DFI_Domain->GlobalVoxel[1]));
        //z
        for(int nk=0; nk<voxel_thin[2]; nk++) {
          coord_Z[nk] = (float)(DFI_Domain->NodeZ(nk*thin_count));
        }
        coord_Z[voxel_thin[2]] = (float)(DFI_Domain->NodeZ(DFI_Domain->GlobalVoxel[2]));

        dfi=cdm_DFI::WriteInit<float>(MPI_COMM_WORLD,
                                      outdfifname,
                                      m_param->Get_OutputDir(),
                                      DFI_FInfo->Prefix,
                                      m_param->Get_OutputFormat(),
                                      outGc,
                                      d_type,
                                      DFI_FInfo->NumVariables,
                                      outprocfname,
                                      voxel_thin,
                                      coord_X,
                                      coord_Y,
                                      coord_Z,
                                      DFI_Domain->GetCoordinateFile(),
                                      DFI_Domain->GetCoordinateFileType(),
                                      DFI_Domain->GetCoordinateFileEndian(),
                                      m_Gdiv,
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

        //全計算領域の座標をWriteInitに渡す
        coord_X = new double[voxel_thin[0]+1]; //+1はセル数ではなく格子数のため。
        coord_Y = new double[voxel_thin[1]+1];
        coord_Z = new double[voxel_thin[2]+1];

        //配列(coord_X,coord_Y,coord_Z)に値をセット
        //x
        for(int ni=0; ni<voxel_thin[0]; ni++) {
          coord_X[ni] = (double)(DFI_Domain->NodeX(ni*thin_count));
        }
        coord_X[voxel_thin[0]] = (double)(DFI_Domain->NodeX(DFI_Domain->GlobalVoxel[0]));
        //y
        for(int nj=0; nj<voxel_thin[1]; nj++) {
          coord_Y[nj] = (double)(DFI_Domain->NodeY(nj*thin_count));
        }
        coord_Y[voxel_thin[1]] = (double)(DFI_Domain->NodeY(DFI_Domain->GlobalVoxel[1]));
        //z
        for(int nk=0; nk<voxel_thin[2]; nk++) {
          coord_Z[nk] = (double)(DFI_Domain->NodeZ(nk*thin_count));
        }
        coord_Z[voxel_thin[2]] = (double)(DFI_Domain->NodeZ(DFI_Domain->GlobalVoxel[2]));

        dfi=cdm_DFI::WriteInit<double>(MPI_COMM_WORLD,
                                       outdfifname,
                                       m_param->Get_OutputDir(),
                                       DFI_FInfo->Prefix,
                                       m_param->Get_OutputFormat(),
                                       outGc,
                                       d_type,
                                       DFI_FInfo->NumVariables,
                                       outprocfname,
                                       voxel_thin,
                                       coord_X,
                                       coord_Y,
                                       coord_Z,
                                       DFI_Domain->GetCoordinateFile(),
                                       DFI_Domain->GetCoordinateFileType(),
                                       DFI_Domain->GetCoordinateFileEndian(),
                                       m_Gdiv,
                                       head,
                                       tail,
                                       m_HostName,
                                       CDM::E_CDM_OFF);
      }
    }
    if( dfi == NULL ) {
      printf("\tFails to instance dfi\n");
      Exit(0);
    }

    //Procファイル出力
    const cdm_Process *DFI_Process = m_in_dfi[i]->GetcdmProcess();
    if( m_param->Get_Outputdfi_on() ) {
      dfi->WriteProcDfiFile(MPI_COMM_WORLD,
                            false,
                            DFI_Process->RankList[0].c_id,   //RankID=0のCellIDをセット
                            DFI_Process->RankList[0].bc_id); //RankID=0の境界IDをセット
      printf("CellID and BCflagID of all ranks were converted into those of rank 0.\n");
    }

    //出力形式（ascii,binary,Fbinary)のセット
    dfi->set_output_type(m_param->Get_OutputFileType());

    //節点への補間フラグのセット(AVSおよびVTK形式)
    if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_AVS || 
        m_param->Get_OutputFormat() == CDM::E_CDM_FMT_VTK ) {
      dfi->set_interp_flag(m_param->Get_Interp_flag());
    }

    //座標データの出力形式のセット(AVS形式)
    if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_AVS ) {
      dfi->set_output_type_coord(m_param->Get_OutputFileTypeCoord());
    }

    //gridファイルを出力(PLOT3D形式，iblankはすべて1にセット)
    if (m_param->Get_OutputFormat() == CDM::E_CDM_FMT_PLOT3D) {
      int voxel_ib[3];
      for(int i=0; i<3; i++) voxel_ib[i] = tail[i]-head[i]+1;
      size_t size_ib=(voxel_ib[0]+2*outGc)*(voxel_ib[1]+2*outGc)*(voxel_ib[2]+2*outGc);
      int *iblank;
      iblank = new int[size_ib];
      for(int i=0; i<size_ib; i++) {
        iblank[i] = 1;
      }
      dfi->WriteGridFile(iblank);
      delete [] iblank;
    }

    //Unitのセット
#if 0
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
#else
    const cdm_Unit *pUnitIn = m_in_dfi[i]->GetcdmUnit();
    map<string,cdm_UnitElem>::const_iterator it = pUnitIn->UnitList.begin();
    for( ;it!=pUnitIn->UnitList.end();it++ )
    {
      string Name = (it->second).Name;
      string Unit = (it->second).Unit;
      double ref  = (it->second).reference;
      double diff = (it->second).difference;
      bool bdiff  = (it->second).BsetDiff;
      dfi->AddUnit(Name, Unit, ref, diff, bdiff);
    }
#endif

    //変数名の取り出しとセット
    for(int n=0; n<DFI_FInfo->NumVariables; n++) {
       std::string variable = m_in_dfi[i]->getVariableName(n);
       if( variable != "" ) dfi->setVariableName(n,variable);
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

  CDM::E_CDM_DTYPE d_type;

  CDM::E_CDM_ERRORCODE ret; 
  double rtime;
  unsigned idummy = 0;
  double ddummy = 0.0;
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

  const cdm_FileInfo* DFI_FInfo = m_in_dfi[0]->GetcdmFileInfo();

  //dfiのループ
  for (int i=0; i<m_in_dfi.size(); i++) {


    int nVari = m_in_dfi[i]->GetNumVariables();

    int outGc=0;
    if( m_param->Get_OutputGuideCell() > 0 ) outGc = m_param->Get_OutputGuideCell();
    if( outGc > 0 ) {
      const cdm_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcdmFileInfo();
      if( outGc > DFI_FInfo->GuideCell ) outGc = DFI_FInfo->GuideCell;
    }

    if( thin_count > 1 ) outGc=0;
    if( m_param->Get_Interp_flag() ) outGc=0; 

    //読込みバッファのインスタンス
    cdm_Array* buf = cdm_Array::instanceArray
    ( m_in_dfi[i]->GetDataType(),
      m_in_dfi[i]->GetArrayShape(),
      sz,
      //0,
      outGc,
      //m_in_dfi[i]->GetNumVariables());
      nVari);

//20160425.fub.s
    //fubファイルのとき読み込み座標値バッファのインスタンス
    cdm_FieldFileNameFormat * Ffformat =
       (cdm_FieldFileNameFormat *)m_in_dfi[i]->GetcdmFieldFileNameFormat();
    cdm_DFI_FUB *dfi_fub = dynamic_cast<cdm_DFI_FUB*>(m_in_dfi[i]);
    cdm_Array * buf_xyz = NULL;
    if( dfi_fub || m_param->Get_OutputFormat() == CDM::E_CDM_FMT_FUB ) {
      buf_xyz = cdm_Array::instanceArray
      ( m_in_dfi[i]->GetDataType(),
        m_in_dfi[i]->GetArrayShape(),
        sz,
        outGc,
        3);
    }
//20160425.fub.e

    //出力タイプのセット
    if( m_param->Get_OutputDataType() == CDM::E_CDM_DTYPE_UNKNOWN )
    {
      d_type = m_in_dfi[i]->GetDataType();
    } else {
      d_type = m_param->Get_OutputDataType();
    }

    //出力バッファのインスタンス
    cdm_Array* src = cdm_Array::instanceArray
    ( d_type,
      //m_in_dfi[i]->GetArrayShape(),
      m_param->Get_OutputArrayShape(),
      szS,
      //0,
      outGc,
      //m_in_dfi[i]->GetNumVariables());
      nVari);

//20160425.fub.s
    cdm_Array* src_xyz = NULL;
//20160425.fub.e
   
    //DFI_FInfoクラスの取得
//20160425.fub.s
  //const cdm_FileInfo* DFI_FInfo = m_in_dfi[i]->GetcdmFileInfo();
    cdm_FileInfo* DFI_FInfo = (cdm_FileInfo *)m_in_dfi[i]->GetcdmFileInfo();
    const cdm_Domain* DFI_Domain = m_in_dfi[i]->GetcdmDomain();
//20160425.fub.e
    prefix=DFI_FInfo->Prefix; 

    //TimeSliceクラスの取得
    const cdm_TimeSlice* TSlice = m_in_dfi[i]->GetcdmTimeSlice();

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
      if( ret != CDM::E_CDM_SUCCESS ) {
        printf("ReadData Error\n");
        return false;
      }

//20160425.fub.s
      if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_FUB )
      {

        if( dfi_fub ) {
          CDM::E_CDM_FORMAT t_fmt = DFI_FInfo->FileFormat;
          DFI_FInfo->FileFormat = CDM::E_CDM_FMT_FUB_COD;
          int t_val = DFI_FInfo->NumVariables;
          DFI_FInfo->NumVariables = 3;

          ret = m_in_dfi[i]->ReadData(buf_xyz,
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

          DFI_FInfo->FileFormat = t_fmt;
          DFI_FInfo->NumVariables = t_val;

          src_xyz = NULL;

          if( ret == CDM::E_CDM_SUCCESS ) {
            //座標値出力バッファのインスタンス
            src_xyz = cdm_Array::instanceArray
            ( d_type,
            m_param->Get_OutputArrayShape(),
            szS,
            outGc,
            3);

          }
        }
 
        if( src_xyz == NULL && j==0 ) {
          //座標値出力バッファのインスタンス
          src_xyz = cdm_Array::instanceArray
          ( d_type,
          m_param->Get_OutputArrayShape(),
          szS,
          outGc,
          3);

          if( d_type == CDM::E_CDM_FLOAT64 ) {
            double *buf_p = (double *)buf_xyz->getData();
            //座標値を計算する
            for(int k=0-outGc, kk=m_Head[2]; kk<=m_Tail[2]; k++, kk++ ) {
            for(int j=0-outGc, jj=m_Head[1]; jj<=m_Tail[1]; j++, jj++ ) {
            for(int i=0-outGc, ii=m_Head[0]; ii<=m_Tail[0]; i++, ii++ ) {
              buf_p[_CDM_IDX_IJKN(i,j,k,0,sz[0],sz[1],sz[2],outGc)] = DFI_Domain->CellX(ii-1);
              buf_p[_CDM_IDX_IJKN(i,j,k,1,sz[0],sz[1],sz[2],outGc)] = DFI_Domain->CellY(jj-1);
              buf_p[_CDM_IDX_IJKN(i,j,k,2,sz[0],sz[1],sz[2],outGc)] = DFI_Domain->CellZ(kk-1);
            }}}
          } else if( d_type == CDM::E_CDM_FLOAT32 ) {
            float *buf_p = (float *)buf_xyz->getData();
            //座標値を計算する
            for(int k=0-outGc, kk=m_Head[2]; kk<=m_Tail[2]; k++, kk++ ) {
            for(int j=0-outGc, jj=m_Head[1]; jj<=m_Tail[1]; j++, jj++ ) {
            for(int i=0-outGc, ii=m_Head[0]; ii<=m_Tail[0]; i++, ii++ ) {
              buf_p[_CDM_IDX_IJKN(i,j,k,0,sz[0],sz[1],sz[2],outGc)] = DFI_Domain->CellX(ii-1);
              buf_p[_CDM_IDX_IJKN(i,j,k,1,sz[0],sz[1],sz[2],outGc)] = DFI_Domain->CellY(jj-1);
              buf_p[_CDM_IDX_IJKN(i,j,k,2,sz[0],sz[1],sz[2],outGc)] = DFI_Domain->CellZ(kk-1);
            }}}
          }
        }
      }
//20160425.fub.e

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

        for(int n=0; n<nVari; n++) convertXY(buf,src,headS,tailS,n);
//20160425.s
        if( src_xyz ) {
          for(int n=0; n<3; n++) convertXY(buf_xyz,src_xyz,headS,tailS,n);
        }
//20160425.e
      }

      CDM::E_CDM_OUTPUT_FNAME output_fname = m_param->Get_OutputFilenameFormat();
      m_out_dfi[i]->set_output_fname(output_fname);

      //minmaxの初期化
      int nsize = nVari;
      if( nVari > 1 ) nsize++;
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
  
      m_out_dfi[i]->SetcdmTimeSlice(*TSlice);

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


//20160425.fub.s
      if( src_xyz ) {
      //out DFI_FInfoクラスの取得
//20160425.fub.s
        cdm_FileInfo* out_DFI_FInfo = (cdm_FileInfo *)m_out_dfi[i]->GetcdmFileInfo();
//20160425.fub.e
        CDM::E_CDM_FORMAT t_fmt = out_DFI_FInfo->FileFormat;
        out_DFI_FInfo->FileFormat = CDM::E_CDM_FMT_FUB_COD;
        int t_val = out_DFI_FInfo->NumVariables;
        out_DFI_FInfo->NumVariables = 3;

        ret = m_out_dfi[i]->WriteData(
                                      (unsigned)TSlice->SliceList[j].step,
                                      outGc,
                                      rtime,
                                      src_xyz,
                                      tmp_minmax,
                                      true,
                                      idummy,
                                      ddummy);

        out_DFI_FInfo->FileFormat = t_fmt;
        out_DFI_FInfo->NumVariables = t_val;
      }
//20160425.fub.e

    } 
    delete src;
//20160425.fub.s
    if( src_xyz ) delete src_xyz;
//20160425.fub.e
  }

  return true;

}

