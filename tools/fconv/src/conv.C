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
 * @file   conv.C
 * @brief  CONV Class
 * @author aics
 */

#include "conv.h"
#include "convMx1.h"
#include "convMxM.h"
#include "convMxN.h"

// #################################################################
// コンストラクタ
CONV::CONV()
{
  m_procGrp = 0;
  m_myRank  = -1;
  m_numProc = 0;
  
  m_pflag=0;
  m_pflagv=0;
  m_lflag=0;
  m_lflagv=0;
 
  m_bgrid_interp_flag = false;
 
  m_in_dfi.clear();
  
  m_staging=0;
  
}


// #################################################################
// デストラクタ
CONV::~CONV()
{
  for(int i=0; i<m_in_dfi.size(); i++ ) if( !m_in_dfi[i] ) delete m_in_dfi[i];
  LOG_OUT_ fclose(m_fplog);
}

// #################################################################
// CONVのインスタンス
CONV* CONV::ConvInit(InputParam* param)
{
  CONV* conv = NULL;

  if(param->Get_ConvType() == E_CONV_OUTPUT_Mx1 ) {
    conv = new convMx1();
  } else if( param->Get_ConvType() == E_CONV_OUTPUT_MxM ) {
    conv = new convMxM();
  } else if( param->Get_ConvType() == E_CONV_OUTPUT_MxN ) {
    conv = new convMxN();
  }

  //InputParamのポインタをセット
  conv->m_param = param;

  //格子点補間フラグのセット
  if( param->Get_OutputFormat() == CDM::E_CDM_FMT_PLOT3D ||
      param->Get_OutputFormat() == CDM::E_CDM_FMT_AVS || 
      param->Get_OutputFormat() == CDM::E_CDM_FMT_VTK ) {
    conv->m_bgrid_interp_flag = true;
  } else {
    conv->m_bgrid_interp_flag = false;
  }
 
  return conv;
}

// #################################################################
//
CDM::E_CDM_ERRORCODE CONV::ReadDfiFiles()
{
  
  // dfiファイルの読込み
  int tempg[3];
  int tempd[3];
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;
  for(int i=0; i<m_param->m_dfiList.size(); i++) {
    cdm_DFI* dfi_in = cdm_DFI::ReadInit(MPI_COMM_WORLD,
                                        m_param->m_dfiList[i].in_dfi_name,
                                        tempg,
                                        tempd,
                                        ret);
    if( dfi_in == NULL ) return ret;
    if( ret != CDM::E_CDM_SUCCESS && ret != CDM::E_CDM_ERROR_INVALID_DIVNUM ) return ret;
    m_in_dfi.push_back(dfi_in);
    m_param->m_dfiList[i].in_dfi = dfi_in;
  }
 
  LOG_OUTV_ {
    PrintDFI(m_fplog);
    m_param->PrintParam(m_fplog);
  }
  
  STD_OUTV_ {
    if( m_myRank == 0 ) {
      PrintDFI(stdout);
      m_param->PrintParam(stdout);
    }
  }
  
  // dfi毎の成分数、出力ガイドセルのチェックと更新、等
  if( !CheckDFIdata() ) return CDM::E_CDM_ERROR;

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// dfi毎の成分数、出力ガイドセルのチェックと更新、等
bool CONV::CheckDFIdata()
{

  bool ierr=true;
  bool upsta=true;
  bool upend=true;
  int *sta;
  int *end;

  for( int i=0; i<m_in_dfi.size(); i++) {
    //コンバート成分数のチェック
    if( m_param->Get_OutputFormat() == CDM::E_CDM_FMT_SPH ) {
      if( m_in_dfi[i]->GetNumComponent() > 3 ) {
        printf("\tCan't Converter OutputFormat.\n");
        ierr=false;
      }
    }

    //ガイドセル数のチェックと更新
    if( m_param->Get_OutputGuideCell() > 0 ) {
      if( m_in_dfi[i]->GetNumGuideCell() < m_param->Get_OutputGuideCell() ) {
        m_param->Set_OutputGuideCell(m_in_dfi[i]->GetNumGuideCell());
        printf("\tupdate OutputGudeCell %d\n",m_param->Get_OutputGuideCell() );
      }
    }

    //入力指示範囲のチェック
    if( m_param->Get_CropIndexStart_on() ) {
      sta = m_param->Get_CropIndexStart();
      for(int j=0; j<3; j++) {
        if( sta[j]<1 ) {
          sta[j]=1;
          upsta=false;
        }
      }
    }
    if( m_param->Get_CropIndexEnd_on() ) {
      const cdm_Domain* DFI_Domian = m_in_dfi[i]->GetcdmDomain();
      end = m_param->Get_CropIndexEnd();
      for(int j=0; j<3; j++) {
        if( end[j]>DFI_Domian->GlobalVoxel[j] ) {
          end[j]=DFI_Domian->GlobalVoxel[j];
          upend=false;
        }
      }
    }

    //Prefixの重複チェック
    const cdm_FileInfo* Finfo1 = m_in_dfi[i]->GetcdmFileInfo();
    for(int j=i+1; j<m_in_dfi.size(); j++) {
      const cdm_FileInfo* Finfo2 = m_in_dfi[j]->GetcdmFileInfo();
      if( Finfo1->Prefix == Finfo2->Prefix ) {
        printf("\tCan't duplicate Prefix \"%s\" dfi : %s dfi : %s\n",
                Finfo1->Prefix.c_str(),m_in_dfi[i]->get_dfi_fname().c_str(),
                m_in_dfi[j]->get_dfi_fname().c_str());
        ierr=false;
      } 
    }

    //voxel size のチェック
    if( m_param->Get_ConvType() == E_CONV_OUTPUT_MxN ) {
      const cdm_Domain* domain1 = m_in_dfi[i]->GetcdmDomain();
      for(int j=i+1; j<m_in_dfi.size(); j++) {
        const cdm_Domain* domain2 = m_in_dfi[j]->GetcdmDomain();
        if( domain1->GlobalVoxel[0] != domain2->GlobalVoxel[0] ||
            domain1->GlobalVoxel[1] != domain2->GlobalVoxel[1] ||
            domain1->GlobalVoxel[2] != domain2->GlobalVoxel[2]   ) {
          printf("\tCan't mismatch ClobalVoxel size %s (%d,%d,%d) \n"
                 "\t                                %s (%d,%d,%d) \n",
            m_in_dfi[i]->get_dfi_fname().c_str(),
            domain1->GlobalVoxel[0],domain1->GlobalVoxel[1],domain1->GlobalVoxel[2],
            m_in_dfi[j]->get_dfi_fname().c_str(),
            domain2->GlobalVoxel[0],domain2->GlobalVoxel[1],domain2->GlobalVoxel[2]);

          ierr=false;
        }
      }
    }
  }

  //入力指示範囲の更新
  if( !upsta ) {
    m_param->Set_CropIndexStart(sta);
    printf("\tupdate input CropIndexStart : %d %d %d\n",sta[0],sta[1],sta[2]);
  }
  if( !upend ) {
    m_param->Set_CropIndexEnd(end);
    printf("\tupdate input CropIndexEnd   : %d %d %d\n",end[0],end[1],end[2]);
  }

  return ierr;

}

// #################################################################
//
void CONV::CheckDir(string dirstr)
{
  Hostonly_
  {
    
#ifndef _WIN32
    
    if( dirstr.size() == 0 ) {
      //printf("\toutput current directory\n");
      return;
    }
    
    DIR* dir;
    if( !(dir = opendir(dirstr.c_str())) ) {
      if( errno == ENOENT ) {
        mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
        if ( cdm_DFI::MakeDirectorySub(dirstr) != 0 )
        {
          printf("\tCan't generate directory(%s).\n", dirstr.c_str());
          Exit(0);
        }
      }
      else {
        printf("Directory open error.(%s)", dirstr.c_str());
        Exit(0);
      }
    }
    else {
      if( closedir(dir) == -1 ) {
        printf("Directory close error.(%s)", dirstr.c_str());
        Exit(0);
      }
    }
    
#else // for windows
    
    if( dirstr.size() == 0 ) {
      printf("\toutput current directory\n");
      return;
    }
    
    // check to exist directory
    if (IsDirExsist(dirstr)) {
      // exist directory
      return;
    }
    
    // make directory
    if(!CreateDirectory(dirstr.c_str(), NULL)){
      printf("\tCan't generate directory(%s).\n", dirstr.c_str());
      Exit(0);
    }
    
#endif  // _WIN32
    
  }
  
  return;
}

// #################################################################
//
void CONV::OpenLogFile()
{
  //log file open
  string prefix="log_comb_id";
  int len = prefix.size()+10;//+6+4
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  sprintf(tmp, "%s_%06d.%s", prefix.c_str(), m_myRank, "txt");
  
  std::string logname(tmp);
  if ( tmp ) delete [] tmp;
  
  if ( !(m_fplog = fopen(logname.c_str(), "w")) ){
    printf("\tFile Open Error : '%s'\n",logname.c_str());
    Exit(0);
  }
  fprintf(m_fplog,"####################\n",logname.c_str());
  fprintf(m_fplog,"### log_comb.txt ###\n",logname.c_str());
  fprintf(m_fplog,"####################\n",logname.c_str());
  fprintf(m_fplog,"\n");
  
  fprintf(m_fplog,"procGrp  = %d\n", m_procGrp);
  fprintf(m_fplog,"myRank   = %d\n", m_myRank);
  fprintf(m_fplog,"numProc  = %d\n", m_numProc);
  fprintf(m_fplog,"HostName = %s\n", m_HostName.c_str());
  fprintf(m_fplog,"\n");
  
}

// #################################################################
//
/*
void CONV::CloseLogFile()
{
  fclose(m_fplog);
}
*/

// #################################################################
// dfiをログ出力
void CONV::PrintDFI(FILE* fp)
{

  fprintf(fp,"\n");
  fprintf(fp,"*** dfi file info ***\n");
  fprintf(fp,"\n");
  for(int i=0; i<m_param->m_dfiList.size(); i++) {
    fprintf(fp,"\tDFI File Name : %s\n",m_param->m_dfiList[i].in_dfi_name.c_str());
    cdm_DFI *dfi = m_param->m_dfiList[i].in_dfi;
    const cdm_FileInfo *DFI_Info = dfi->GetcdmFileInfo();
    if( DFI_Info->DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {
      fprintf(fp,"\tDFI_Info->DFIType                  = \"Cartesian\"\n");
    } else {
      fprintf(fp,"\tDFI_Info->DFIType                  = \"\"\n");
    }
    fprintf(fp,"\tDFI_Info->DirectoryPath            = \"%s\"\n",DFI_Info->DirectoryPath.c_str());
    if( DFI_Info->TimeSliceDirFlag == CDM::E_CDM_ON ) {
      fprintf(fp,"\tDFI_Info->TimeSliceDirFlag         = \"on\"\n");
    } else {
      fprintf(fp,"\tDFI_Info->TimeSliceDirFlag         = \"off\"\n");
    }
    fprintf(fp,"\tDFI_Info->Prefix                   = \"%s\"\n",DFI_Info->Prefix.c_str());
    fprintf(fp,"\tDFI_Info->FileFormat               = \"%s\"\n",
            dfi->GetFileFormatString().c_str());
    if( DFI_Info->FieldFilenameFormat == CDM::E_CDM_FNAME_RANK_STEP ) {
      fprintf(fp,"\tDFI_Info->FieldFilenameFormat      = \"rank_step\"\n");
    }else {
      fprintf(fp,"\tDFI_Info->FieldFilenameFormat      = \"step_rank\"\n");
    }
    fprintf(fp,"\tDFI_Info->GuideCell                = %d\n",DFI_Info->GuideCell);
    fprintf(fp,"\tDFI_Info->DataType                 = \"%s\"\n",
            dfi->GetDataTypeString().c_str());
    if( DFI_Info->Endian == CDM::E_CDM_LITTLE ) {
      fprintf(fp,"\tDFI_Info->Endian                   = \"little\"\n");
    }else if( DFI_Info->Endian == CDM::E_CDM_BIG ) {
      fprintf(fp,"\tDFI_Info->Endian                   = \"big\"\n");
    }else {
      fprintf(fp,"\tDFI_Info->Endian                   = \"\"\n");
    }
    fprintf(fp,"\tDFI_Info->ArrayShape               = \"%s\"\n",
            dfi->GetArrayShapeString().c_str());
    fprintf(fp,"\tDFI_Info->Component                = %d\n",DFI_Info->Component);
    for(int j=0; j<DFI_Info->ComponentVariable.size(); j++ ) {
      fprintf(fp,"\t  DFI_Info->ComponentVariable[%d]        = %s\n",j,
              DFI_Info->ComponentVariable[j].c_str());
    }

    const cdm_MPI *DFI_MPI = dfi->GetcdmMPI();
    fprintf(fp,"\tDFI_MPI->NumberOfRank              = %d\n",DFI_MPI->NumberOfRank); 
    fprintf(fp,"\tDFI_MPI->NumberOfGroup             = %d\n",DFI_MPI->NumberOfGroup); 
   
    const cdm_Domain *DFI_Domain = dfi->GetcdmDomain();
    fprintf(fp,"\tDFI_Domain->GlobalVoxel[0]         = %d\n",DFI_Domain->GlobalVoxel[0]); 
    fprintf(fp,"\tDFI_Domain->GlobalVoxel[1]         = %d\n",DFI_Domain->GlobalVoxel[1]); 
    fprintf(fp,"\tDFI_Domain->GlobalVoxel[2]         = %d\n",DFI_Domain->GlobalVoxel[2]); 
    fprintf(fp,"\tDFI_Domain->GlobalDivision[0]      = %d\n",DFI_Domain->GlobalDivision[0]);
    fprintf(fp,"\tDFI_Domain->GlobalDivision[1]      = %d\n",DFI_Domain->GlobalDivision[1]);
    fprintf(fp,"\tDFI_Domain->GlobalDivision[2]      = %d\n",DFI_Domain->GlobalDivision[2]);

    const cdm_Process *DFI_Process = dfi->GetcdmProcess();
    fprintf(fp,"\tDFI_Process->RankList.size()       = %d\n",DFI_Process->RankList.size());
    for(int j=0; j<DFI_Process->RankList.size(); j++) {
      fprintf(fp,"\t  DFI_Process->RankList[%d].RankID       = %d\n",j,
              DFI_Process->RankList[j].RankID);
      fprintf(fp,"\t  DFI_Process->RankList[%d].HostName     = %s\n",j,
              DFI_Process->RankList[j].HostName.c_str());
      fprintf(fp,"\t  DFI_Process->RankList[%d].VoxelSize[0] = %d\n",j,
              DFI_Process->RankList[j].VoxelSize[0]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].VoxelSize[1] = %d\n",j,
              DFI_Process->RankList[j].VoxelSize[1]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].VoxelSize[2] = %d\n",j,
              DFI_Process->RankList[j].VoxelSize[2]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].HeadIndex[0] = %d\n",j,
              DFI_Process->RankList[j].HeadIndex[0]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].HeadIndex[1] = %d\n",j,
              DFI_Process->RankList[j].HeadIndex[1]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].HeadIndex[2] = %d\n",j,
              DFI_Process->RankList[j].HeadIndex[2]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].TailIndex[0] = %d\n",j,
              DFI_Process->RankList[j].TailIndex[0]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].TailIndex[1] = %d\n",j,
              DFI_Process->RankList[j].TailIndex[1]);
      fprintf(fp,"\t  DFI_Process->RankList[%d].TailIndex[2] = %d\n",j,
              DFI_Process->RankList[j].TailIndex[2]);
    }
    
    const cdm_TimeSlice* DFI_TSlice = dfi->GetcdmTimeSlice();
    fprintf(fp,"\tDFI_TSlice->SliceList.size()       = %d\n",DFI_TSlice->SliceList.size());
    for(int j=0; j<DFI_TSlice->SliceList.size(); j++ ) {
      fprintf(fp,"\t  DFI_TSlice->SliceList[%d].step         = %d\n",j,
              DFI_TSlice->SliceList[j].step);
      fprintf(fp,"\t  DFI_TSlice->SliceList[%d}.time         = %f\n",j,
              DFI_TSlice->SliceList[j].time);
      if( !DFI_TSlice->SliceList[j].avr_mode ) {
        fprintf(fp,"\t  DFI_TSlice->SliceList[%d}.AveragedStep = %d\n",j,
              DFI_TSlice->SliceList[j].AveragedStep);
        fprintf(fp,"\t  DFI_TSlice->SliceList[%d}.AveragedTime = %f\n",j,
              DFI_TSlice->SliceList[j].AveragedTime);
      }
      if( DFI_Info->Component > 1 ) {
        fprintf(fp,"\t  DFI_TSlice->SliceList[%d}.VectorMin    = %f\n",j,
              DFI_TSlice->SliceList[j].VectorMin);
        fprintf(fp,"\t  DFI_TSlice->SliceList[%d}.VectorMax    = %f\n",j,
              DFI_TSlice->SliceList[j].VectorMax);
        for(int k=0; k<DFI_TSlice->SliceList[j].Min.size(); k++) {
          fprintf(fp,"\t    DFI_TSlice->SliceList[%d}.Min[%d]      = %f\n",j,k,
              DFI_TSlice->SliceList[j].Min[k]);
        }
        for(int k=0; k<DFI_TSlice->SliceList[j].Max.size(); k++) {
          fprintf(fp,"\t    DFI_TSlice->SliceList[%d}.Max[%d]      = %f\n",j,k,
              DFI_TSlice->SliceList[j].Max[k]);
        }
      }
    }

    fprintf(fp,"\n");
  }

}

// #################################################################
//
void CONV::WriteTime(double* tt)
{
  fprintf(m_fplog,"\n\n");
  fprintf(m_fplog,"TIME : ReadInit      %10.3f sec.\n", tt[0]);
  fprintf(m_fplog,"TIME : ReadDfiFiles  %10.3f sec.\n", tt[1]);
  fprintf(m_fplog,"TIME : Converter     %10.3f sec.\n", tt[2]);
  fprintf(m_fplog,"TIME : Total Time    %10.3f sec.\n", tt[3]);
}

// #################################################################
// メモリ使用量を表示する
void CONV::MemoryRequirement(const double Memory, FILE* fp)
{
  const double mem = Memory;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional
  
  // Global memory
  fprintf (fp," MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)\n", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)\n", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)\n", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)\n", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)\n", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)\n", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  fflush(fp);
}

// #################################################################
// メモリ使用量を表示する
void CONV::MemoryRequirement(const double TotalMemory, const double sphMemory, const double plot3dMemory, const double thinMemory, FILE* fp)
{
  double mem;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional
  
  fprintf (fp,"*** Required MemorySize ***");
  fprintf (fp,"\n");
  
  mem = sphMemory;
  fprintf (fp,"  read SPH MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  mem = plot3dMemory;
  fprintf (fp,"  write PLOT3D MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  mem = thinMemory;
  fprintf (fp,"  write thin out MemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  mem = TotalMemory;
  fprintf (fp,"  TotalMemorySize = ");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  fprintf (fp,"\n");
  fprintf (fp,"\n");
  
  fflush(fp);
}

// #################################################################
// step番号からtimeを取得
double CONV::GetSliceTime(cdm_DFI* dfi, int step)
{

  const cdm_TimeSlice* Tslice = dfi->GetcdmTimeSlice();
  for(int i=0; i<Tslice->SliceList.size(); i++) {
    if( Tslice->SliceList[i].step == step ) return Tslice->SliceList[i].time;
  }

  return 0.0;

} 

// #################################################################
// X-Y面のコンバイン
bool CONV::convertXY(
                     cdm_Array* buf,
                     cdm_Array* &src,
                     int headS[3],
                     int tailS[3],
                     int n)
{
 
  //debug
  const int *tmp = src->getHeadIndex();
 
  //copy
  int gcB          = buf->getGcInt();
  const int *headB = buf->getHeadIndex();
  const int *tailB = buf->getTailIndex();
  int       gcS    = src->getGcInt();
  int sta[3],end[3];
  for( int i=0;i<3;i++ )
  {
    sta[i] = (headB[i]-gcB>=headS[i]-gcS) ? headB[i]-gcB : headS[i]-gcS;
    end[i] = (tailB[i]+gcB<=tailS[i]+gcS) ? tailB[i]+gcB : tailS[i]+gcS;
  }

  CDM::E_CDM_DTYPE buf_dtype =  buf->getDataType();

  //uint8
  if( buf_dtype == CDM::E_CDM_UINT8 ) {
    cdm_TypeArray<unsigned char> *B = dynamic_cast<cdm_TypeArray<unsigned char>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //int8
  else if( buf_dtype == CDM::E_CDM_INT8 ) {
    cdm_TypeArray<char> *B = dynamic_cast<cdm_TypeArray<char>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //uint16
  else if( buf_dtype == CDM::E_CDM_UINT16 ) {
    cdm_TypeArray<unsigned short> *B = dynamic_cast<cdm_TypeArray<unsigned short>*>(buf); 
    return copyArray(B, src, sta, end, n);
  }
  //int16
  else if( buf_dtype == CDM::E_CDM_INT16 ) {
    cdm_TypeArray<short> *B = dynamic_cast<cdm_TypeArray<short>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //uint32
  else if( buf_dtype == CDM::E_CDM_UINT32 ) {
    cdm_TypeArray<unsigned int> *B = dynamic_cast<cdm_TypeArray<unsigned int>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //int32
  else if( buf_dtype == CDM::E_CDM_INT32 ) {
    cdm_TypeArray<int> *B = dynamic_cast<cdm_TypeArray<int>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //uint64
  else if( buf_dtype == CDM::E_CDM_UINT64 ) {
    cdm_TypeArray<unsigned long long> *B = dynamic_cast<cdm_TypeArray<unsigned long long>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //int64
  else if( buf_dtype == CDM::E_CDM_INT64 ) {
    cdm_TypeArray<long long> *B = dynamic_cast<cdm_TypeArray<long long>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //float32
  else if( buf_dtype == CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *B = dynamic_cast<cdm_TypeArray<float>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  //float64
  else if( buf_dtype == CDM::E_CDM_FLOAT64 ) {
    cdm_TypeArray<double> *B = dynamic_cast<cdm_TypeArray<double>*>(buf);
    return copyArray(B, src, sta, end, n);
  }
  
  return false;
}

// #################################################################
// 出力ファイル形式から拡張子を求める
std::string CONV::GetFilenameExt(int file_format_type)
{

  if     ( file_format_type == CDM::E_CDM_FMT_SPH ) return D_CDM_EXT_SPH;
  else if( file_format_type == CDM::E_CDM_FMT_BOV ) return D_CDM_EXT_BOV;
  else if( file_format_type == CDM::E_CDM_FMT_AVS ) return D_CDM_EXT_SPH;
  else if( file_format_type == CDM::E_CDM_FMT_PLOT3D ) return D_CDM_EXT_FUNC;
  else if( file_format_type == CDM::E_CDM_FMT_VTK ) return D_CDM_EXT_VTK;

  return "";
}

// #################################################################
// step基準のリスト生成
void CONV::makeStepList(vector<step_rank_info> &StepRankList)
{

  //総ステップ数を求める
  int Total_step = 0;
  for(int i=0; i<m_in_dfi.size(); i++){
    const cdm_TimeSlice* TSlice = m_in_dfi[i]->GetcdmTimeSlice();
    Total_step+=TSlice->SliceList.size();
  }

  //自ランクで担当するステップ数を求める
  int nStep = Total_step/m_numProc;
  if( Total_step%m_numProc != 0 ) {
    for(int i=0; i<Total_step%m_numProc; i++) {
      if( m_myRank == i ) nStep++;
    }
  }

  //自ランクが担当するステップのスタートとエンドを求める
  int sta,end;
  sta = m_myRank * nStep;
  if( Total_step%m_numProc != 0 ) {
    if( m_myRank >= Total_step%m_numProc ) sta = sta+Total_step%m_numProc;
  }
  end = sta+nStep-1;

  //処理ステップリストの生成
  int cnt=0;
  for(int i=0; i<m_in_dfi.size(); i++){
    step_rank_info info;
    info.stepStart = -1;
    const cdm_TimeSlice* TSlice = m_in_dfi[i]->GetcdmTimeSlice();
    for( int j=0; j<TSlice->SliceList.size(); j++) {

      if( sta > cnt ) { cnt++; continue; }
      if( info.stepStart == -1 ) {
        info.dfi=m_in_dfi[i];
        info.stepStart=j;
      }
      info.stepEnd = j;
      cnt++;
      if( end < cnt ) break;
    }
    if( info.stepStart > -1 ) StepRankList.push_back(info);
    if( end < cnt ) break;
  }

  //rantStart,rankEndのセット
  for(int i=0; i<StepRankList.size(); i++) {
    const cdm_Process* DFI_Process = StepRankList[i].dfi->GetcdmProcess();
    StepRankList[i].rankStart=0;
    StepRankList[i].rankEnd=DFI_Process->RankList.size()-1;
  }
}

//#################################################################
// rank基準のリスト生成
void CONV::makeRankList(vector<step_rank_info> &StepRankList)
{

  //総rank数を求める
  int Total_rank = 0;
  for(int i=0; i<m_in_dfi.size(); i++){
    const cdm_Process* DFI_Process = m_in_dfi[i]->GetcdmProcess();
    Total_rank+=DFI_Process->RankList.size();
  }

  //自ランクで担当するrank数を求める
  int nRank = Total_rank/m_numProc;
  if( Total_rank%m_numProc != 0 ) {
    for(int i=0; i<Total_rank%m_numProc; i++) {
      if( m_myRank == i ) nRank++;
    }
  }

  //自ランクが担当するランクのスタートとエンドを求める
  int sta,end;
  sta = m_myRank * nRank;
  if( Total_rank%m_numProc != 0 ) {
    if( m_myRank >= Total_rank%m_numProc ) sta = sta+Total_rank%m_numProc;
  }
  end = sta+nRank-1;

  //処理rankリストの生成
  int cnt=0;
  for(int i=0; i<m_in_dfi.size(); i++) {
    step_rank_info info;
    info.rankStart = -1;
    const cdm_Process* DFI_Process = m_in_dfi[i]->GetcdmProcess();
    for(int j=0; j<DFI_Process->RankList.size(); j++) {
      if( sta > cnt ) { cnt++; continue; }
      if( info.rankStart == -1 ) {
        info.dfi = m_in_dfi[i];
        info.rankStart = j;
      }
      info.rankEnd = j;
      cnt++; 
      if( end < cnt ) break;
    }    
    if( info.rankStart > -1 ) StepRankList.push_back(info);
    if( end < cnt ) break;
  }

  //stepStart,stepEndのセット
  for(int i=0; i<StepRankList.size(); i++) {
    const cdm_TimeSlice* TSlice = StepRankList[i].dfi->GetcdmTimeSlice();
    StepRankList[i].stepStart=0;
    StepRankList[i].stepEnd=TSlice->SliceList.size()-1;
  }

}

//#################################################################
// データタイプ毎にminmaxを求める
bool CONV::DtypeMinMax(cdm_Array* src,
                       double *min,
                       double *max)
{
  CDM::E_CDM_DTYPE d_type = src->getDataType();
  if( d_type == CDM::E_CDM_UINT8 ) {
    cdm_TypeArray<unsigned char> *data = dynamic_cast<cdm_TypeArray<unsigned char>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_INT8 ) {
    cdm_TypeArray<char> *data = dynamic_cast<cdm_TypeArray<char>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_UINT16 ) {
    cdm_TypeArray<unsigned short> *data = dynamic_cast<cdm_TypeArray<unsigned short>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_INT16 ) {
    cdm_TypeArray<short> *data = dynamic_cast<cdm_TypeArray<short>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_UINT32 ) {
    cdm_TypeArray<unsigned int> *data = dynamic_cast<cdm_TypeArray<unsigned int>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_INT32 ) {
    cdm_TypeArray<int> *data = dynamic_cast<cdm_TypeArray<int>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_UINT64 ) {
    cdm_TypeArray<unsigned long long> *data = dynamic_cast<cdm_TypeArray<unsigned long long>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_INT64 ) {
    cdm_TypeArray<long long> *data = dynamic_cast<cdm_TypeArray<long long>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_FLOAT32 ) {
    cdm_TypeArray<float> *data = dynamic_cast<cdm_TypeArray<float>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  } else if( d_type == CDM::E_CDM_FLOAT64 ) {
    cdm_TypeArray<double> *data = dynamic_cast<cdm_TypeArray<double>*>(src);
    if( !calcMinMax(data,min,max) ) return false;
  }
 
  return true;

}

//#################################################################
// index.dfiファイル出力
bool CONV::WriteIndexDfiFile(vector<dfi_MinMax*> minmaxList)
{

  if( !m_param->Get_Outputdfi_on() ) return false;

  //出力dfiファイル名の取得
  /*
  vector<std::string> out_dfi_name  = m_InputCntl->Get_OutdfiNameList();
  vector<std::string> out_proc_name = m_InputCntl->Get_OutprocNameList();

  if( minmaxList.size() != out_dfi_name.size() && 
      minmaxList.size() != out_proc_name.size() ) return false;
  */
  for(int i=0; i<minmaxList.size(); i++) {

    cdm_DFI* dfi = minmaxList[i]->dfi;

    std::string out_dfi_name  = m_param->m_dfiList[i].out_dfi_name;
    std::string out_proc_name = m_param->m_dfiList[i].out_proc_name;

    FILE* fp=NULL;
    if( !(fp = fopen(out_dfi_name.c_str(), "w")) )
    {
      printf("Can't open file.(%s)\n", out_dfi_name.c_str());
      return false;
    }

    //成分数の取得
    int nComp = dfi->GetNumComponent();

    cdm_FileInfo *dfi_Finfo = (cdm_FileInfo *)dfi->GetcdmFileInfo();

    CDM::E_CDM_ARRAYSHAPE shape = dfi_Finfo->ArrayShape;
    //if( dfi_Finfo->FileFormat == (CDM::E_CDM_FMT_BOV) ) {
    if( m_param->Get_OutputFormat() == (CDM::E_CDM_FMT_BOV) ) {
      shape = m_param->Get_OutputArrayShape();
    }

    int outGc = m_param->Get_OutputGuideCell();
    if( outGc > 0 ) {
      if( outGc > dfi_Finfo->GuideCell ) outGc=dfi_Finfo->GuideCell;
    }

    CDM::E_CDM_DTYPE out_d_type = m_param->Get_OutputDataType();
    if( out_d_type == CDM::E_CDM_DTYPE_UNKNOWN ) out_d_type=dfi->GetDataType();

    CDM::E_CDM_OUTPUT_FNAME FieldFilenameFormat;
    FieldFilenameFormat=m_param->Get_OutputFilenameFormat();

    //FileInfoの出力
    cdm_FileInfo *Finfo = new cdm_FileInfo(CDM::E_CDM_DFITYPE_CARTESIAN,
                              FieldFilenameFormat,
                              m_param->Get_OutputDir(),
                              dfi_Finfo->TimeSliceDirFlag,
                              dfi_Finfo->Prefix,
                              m_param->Get_OutputFormat(),
                              outGc,
                              //(CDM::E_CDM_DTYPE)m_InputCntl->Get_OutputDataType(),
                              out_d_type,
                              dfi_Finfo->Endian,
                              //(CDM::E_CDM_ARRAYSHAPE)m_InputCntl->Get_OutputArrayShape(),
                              shape,
                              nComp);

    for(int n=0; n<nComp; n++) {
       std::string variable = dfi->getComponentVariable(n);
       if( variable != "" ) Finfo->setComponentVariable(n,variable);
    }

    if( Finfo->Write(fp, 0) != CDM::E_CDM_SUCCESS ) {
      fclose(fp);
      return false;
    }
    delete Finfo;

    //FilePathの出力
    cdm_FilePath *dfi_Fpath = (cdm_FilePath *)dfi->GetcdmFilePath();
    cdm_FilePath *Fpath = new cdm_FilePath(out_proc_name);
    if( Fpath->Write(fp, 1) != CDM::E_CDM_SUCCESS )
    {
      fclose(fp);
      return false;
    }
    delete Fpath;

    //Unitの出力
    cdm_Unit *dfi_Unit = (cdm_Unit *)dfi->GetcdmUnit();
    if( dfi_Unit->Write(fp, 0) != CDM::E_CDM_SUCCESS )
    {
      fclose(fp);
      return false;
    }

    //TimeSliceの出力
    const cdm_TimeSlice *dfi_TSlice = dfi->GetcdmTimeSlice();
    cdm_TimeSlice *TSlice = new cdm_TimeSlice();
    int nsize = nComp;
    if( nComp > 1 ) nsize++;
    double* minmax = new double[nsize*2];
    for(int j=0; j<dfi_TSlice->SliceList.size(); j++) {
      for(int n=0; n<nsize; n++) {
        minmax[n*2+0] = minmaxList[i]->Min[j*nsize+n];
        minmax[n*2+1] = minmaxList[i]->Max[j*nsize+n];
      }
      TSlice->AddSlice(dfi_TSlice->SliceList[j].step,
                       dfi_TSlice->SliceList[j].time,
                       minmax,
                       nComp,
                       true,
                       0,
                       0.0);

    }

    if( TSlice->Write(fp, 1) != CDM::E_CDM_SUCCESS )
    {
      fclose(fp);
      return false;
    }

    //fclose(fp);
    
  }

  return true;
}

//#################################################################
bool CONV::makeProcInfo(cdm_DFI* dfi, 
                        cdm_Domain* &out_domain,
                        cdm_MPI* &out_mpi, 
                        cdm_Process* &out_process,
                        int numProc)
{
  
  //間引き数の取得
  int thin_count = m_param->Get_ThinOut();

  //MPI情報の生成
  out_mpi = new cdm_MPI(numProc,0);

  //Domain 情報の生成
  cdm_Domain* dfi_domain = (cdm_Domain *)dfi->GetcdmDomain();
  double Gorigin[3];
  double Gregion[3];
  int Gvoxel[3];
  int Gdiv[3];
  int IndexStart[3];
  int IndexEnd[3];
  for(int i=0; i<3; i++) {
    Gorigin[i] = dfi_domain->GlobalOrigin[i];
    Gregion[i] = dfi_domain->GlobalRegion[i];
    Gvoxel[i]  = dfi_domain->GlobalVoxel[i];
    Gdiv[i]    = dfi_domain->GlobalDivision[i];
    IndexStart[i]=1;
    IndexEnd[i]=Gvoxel[i];
  }

  //pitを計算
  double pit[3];
  for(int i=0; i<3; i++) pit[i]=Gregion[i]/(double)Gvoxel[i];
 
  //入力領域指示ありのときボクセルサイズを更新
  if( m_param->Get_CropIndexStart_on() ) {
    const int* Start=m_param->Get_CropIndexStart();
    for(int i=0; i<3; i++) IndexStart[i]=Start[i];
    Gvoxel[0]=Gvoxel[0]-IndexStart[0]+1;
    Gvoxel[1]=Gvoxel[1]-IndexStart[1]+1;
    Gvoxel[2]=Gvoxel[2]-IndexStart[2]+1;
  }
  if( m_param->Get_CropIndexEnd_on() ) {
    const int* End=m_param->Get_CropIndexEnd();
    for(int i=0; i<3; i++) IndexEnd[i]=End[i];
    Gvoxel[0]=Gvoxel[0]-(dfi_domain->GlobalVoxel[0]-IndexEnd[0]);
    Gvoxel[1]=Gvoxel[1]-(dfi_domain->GlobalVoxel[1]-IndexEnd[1]);
    Gvoxel[2]=Gvoxel[2]-(dfi_domain->GlobalVoxel[2]-IndexEnd[2]);
  }

  //Gregionの更新
  for(int i=0; i<3; i++) Gregion[i]=(double)Gvoxel[i]*pit[i];

  //間引きありのときボクセルサイズを更新
  if( thin_count > 1 ) {
    for(int i=0; i<3; i++) {
      if( Gvoxel[i]%thin_count != 0 ) Gvoxel[i]=Gvoxel[i]/thin_count+1;
      else                            Gvoxel[i]=Gvoxel[i]/thin_count;
    }
  }
  //numProcが1のとき（Mx1のとき）GlobalDivisionを１にする
  if( numProc == 1 ) for(int i=0; i<3; i++) Gdiv[i]=1;

  //out_domainの生成 
  int *iblank = NULL; //cdm_Domainコンストラクタのiblank追加による一時的な対処(AS田中)
  out_domain = new cdm_Domain(Gorigin,Gregion,Gvoxel,Gdiv,iblank);

  //Process 情報の生成
  const cdm_Process* dfi_Process = dfi->GetcdmProcess();
  out_process = new cdm_Process();
  cdm_Rank rank;
  if( numProc == dfi_Process->RankList.size() ) {
    for(int i=0; i<numProc; i++) {
      rank.RankID = dfi_Process->RankList[i].RankID;
      for(int j=0; j<3; j++) {
        rank.VoxelSize[j]=dfi_Process->RankList[i].VoxelSize[j];
        rank.HeadIndex[j]=dfi_Process->RankList[i].HeadIndex[j];
        rank.TailIndex[j]=dfi_Process->RankList[i].TailIndex[j];
      }
      if( thin_count > 1 ) {
        for(int j=0; j<3; j++) {
          if( rank.VoxelSize[j]%thin_count != 0 ) rank.VoxelSize[j]=rank.VoxelSize[j]/thin_count+1;
          else                                    rank.VoxelSize[j]=rank.VoxelSize[j]/thin_count;
          if( (rank.HeadIndex[j]-1)%thin_count != 0 ) {
             rank.HeadIndex[j]=(rank.HeadIndex[j]-1)/thin_count+1;
          } else {
             rank.HeadIndex[j]=(rank.HeadIndex[j]-1)/thin_count;
          }
          rank.HeadIndex[j]=rank.HeadIndex[j]+1;
          rank.TailIndex[j]=rank.HeadIndex[j]+rank.VoxelSize[j]-1;
        }
      }  
      out_process->RankList.push_back(rank);
    }
  } else if( numProc == 1 ) {
    rank.RankID=0;
    for(int i=0; i<3; i++) {
      rank.VoxelSize[i]=Gvoxel[i];
      //rank.HeadIndex[i]=IndexStart[i];
      rank.HeadIndex[i]=1;
      rank.TailIndex[i]=rank.HeadIndex[i]+Gvoxel[i]-1;
    }
    out_process->RankList.push_back(rank);
  }

  return true;
}


//#################################################################
// index.dfiファイル出力
bool CONV::WriteProcDfiFile(std::string proc_name,
                            cdm_Domain* out_domain,
                            cdm_MPI* out_mpi,
                            cdm_Process* out_process)
{

  FILE* fp = NULL;

  if( out_domain == NULL || out_mpi == NULL || out_process == NULL ) return false;

  //proc.dfiファイルオープン
  if( !(fp = fopen(proc_name.c_str(), "w")) )
  {
     printf("Can't open file.(%s)\n", proc_name.c_str());
     return false;
  }

  //Domain {} の出力
  if( out_domain->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return false;
  }

  //MPI {} の出力
  if( out_mpi->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return false;
  }

  //Process {} の出力
  if( out_process->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return false;
  }

  fclose(fp);

  return true;
}
