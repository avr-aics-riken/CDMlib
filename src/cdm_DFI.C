/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI.C
 * @brief  cdm_DFI Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K
#include <typeinfo>
#include "cdm_DFI_SPH.h"
#include "cdm_DFI_BOV.h"
//FCONV 20131122.s
#include "cdm_DFI_AVS.h"
#include "cdm_DFI_PLOT3D.h"
#include "cdm_DFI_VTK.h"
//FCONV 20131122.e
#include "cdm_NonUniformDomain.h"

// #################################################################
// コンストラクタ
cdm_DFI::cdm_DFI()
{

 m_read_type = CDM::E_CDM_READTYPE_UNKNOWN;
 m_RankID = 0;

 m_input_type = CDM::E_CDM_FILE_TYPE_DEFAULT;
 m_output_type = CDM::E_CDM_FILE_TYPE_DEFAULT;
 m_output_fname = CDM::E_CDM_FNAME_DEFAULT;

}


// #################################################################
// デストラクタ
cdm_DFI::~cdm_DFI()
{
  if( DFI_Domain != NULL ){ delete DFI_Domain; }
}

// #################################################################
// DFI read インスタンス
cdm_DFI* cdm_DFI::ReadInit(const MPI_Comm comm, 
                           const std::string DfiName,
                           const int G_Voxel[3],
                           const int G_Div[3],
                           CDM::E_CDM_ERRORCODE &ret)
{

  /** DFIのディレクトリパスの取得 */
  std::string dirName = CDM::cdmPath_DirName(DfiName);

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  cdm_TextParser tpCntl;

  /** index.dfi read */
  /** TPインスタンス */
  tpCntl.getTPinstance();

  FILE*fp = NULL;
  if( !(fp=fopen(DfiName.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_READ_INDEXFILE_OPENERROR;
    return NULL;
  }
  fclose(fp);

  /** 入力ファイル index.dfi をセット */
  int ierror = 0;
  ierror = tpCntl.readTPfile(DfiName);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_TEXTPARSER;
    return NULL;
  }

  /** Fileinfoの読込み */
  cdm_FileInfo F_info;
  if( F_info.Read(tpCntl) != CDM::E_CDM_SUCCESS ) 
  {
    printf("\tFileInfo Data Read error %s\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_READ_FILEINFO;
    return NULL;
  }

  /** FilePathの読込み */
  cdm_FilePath F_path;
  if( F_path.Read(tpCntl) != CDM::E_CDM_SUCCESS )
  {
    printf("\tFilePath Data Read error %s\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_READ_FILEPATH;
    return NULL;
  }
  
  /** VisItオプションの読込み */
  cdm_VisIt visit;
  if( visit.Read(tpCntl) != CDM::E_CDM_SUCCESS )
  {
    printf("\tVisIt Data Read error %s\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_READ_FILEPATH;
    return NULL;
  }

  /** Unitの読込み */
  cdm_Unit unit;
  if( unit.Read(tpCntl) != CDM::E_CDM_SUCCESS )
  {
    printf("\tUnit Data Read error %s\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_READ_UNIT;
    return NULL;
  }

  /** TimeSliceの読込み */
  cdm_TimeSlice TimeSlice;
  if( TimeSlice.Read(tpCntl, F_info.FileFormat) != CDM::E_CDM_SUCCESS )
  {
    printf("\tTimeSlice Data Read error %s\n",DfiName.c_str());
    ret = CDM::E_CDM_ERROR_READ_TIMESLICE;
    return NULL;
  }

  /** TextParserの破棄 */
  tpCntl.remove();

  /** proc.dfi file name の取得 */
  std::string dfiname  = CDM::cdmPath_FileName(F_path.ProcDFIFile,".dfi");
  std::string procfile = CDM::cdmPath_ConnectPath(dirName,dfiname);

  /** proc.dfi read */
  /** TPインスタンス */
  tpCntl.getTPinstance();

  fp = NULL;
  if( !(fp=fopen(procfile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",procfile.c_str());
    ret = CDM::E_CDM_ERROR_READ_PROCFILE_OPENERROR;
    return NULL;
  }
  fclose(fp);

  /** 入力ファイル proc.dfi をセット */
  ierror = tpCntl.readTPfile(procfile);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",procfile.c_str());
    ret = CDM::E_CDM_ERROR_TEXTPARSER;
    return NULL;
  }

  /** Domainの読込み */
  cdm_Domain* domain = NULL;
  if( F_info.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {
    domain = new cdm_Domain;
  } else if(F_info.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN) {
    domain = new cdm_NonUniformDomain<double>;
  }
  if( domain->Read(tpCntl, dirName) != CDM::E_CDM_SUCCESS )
  {
    printf("\tDomain Data Read error %s\n",procfile.c_str());
    ret = CDM::E_CDM_ERROR_READ_DOMAIN;
    return NULL;
  }

  /** MPIの読込み */
  cdm_MPI mpi;
  if( mpi.Read(tpCntl,domain) != CDM::E_CDM_SUCCESS )
  {
    printf("\tMPI Data Read error %s\n",procfile.c_str());
    ret = CDM::E_CDM_ERROR_READ_MPI;
    return NULL;
  }

  /** Processの読込み */
  cdm_Process process;
  if( process.Read(tpCntl) != CDM::E_CDM_SUCCESS )
  {
    printf("\tProcess Data Read error %s\n",procfile.c_str());
    ret = CDM::E_CDM_ERROR_READ_PROCESS;
    return NULL;
  }

  /** TextParserの破棄 */
  tpCntl.remove();

  /** dfiのインスタンス **/
  cdm_DFI *dfi = NULL;
  // CIOlibではSPHをサポートしていたが、CDMlibではSPHはサポートしない。
  // ↓SPHを復活させた(2014.10.17)
  if( F_info.FileFormat == CDM::E_CDM_FMT_SPH )
  {
    dfi = new cdm_DFI_SPH(F_info, F_path, visit, unit, domain, mpi, TimeSlice, process);
  }
  else if( F_info.FileFormat == CDM::E_CDM_FMT_BOV )
  {
    dfi = new cdm_DFI_BOV(F_info, F_path, visit, unit, domain, mpi, TimeSlice, process);
  }
  else if( F_info.FileFormat == CDM::E_CDM_FMT_PLOT3D )
  {
    dfi = new cdm_DFI_PLOT3D(F_info, F_path, visit, unit, domain, mpi, TimeSlice, process);
  }
  else
  {
    return NULL;
  }

  //読込みタイプのチェック
  dfi->m_read_type = dfi->CheckReadType(G_Voxel,dfi->DFI_Domain->GlobalVoxel,
                                         G_Div,dfi->DFI_Domain->GlobalDivision);
  if( dfi->m_read_type == CDM::E_CDM_READTYPE_UNKNOWN ) {
    //printf("\tDimension size error (%d %d %d)\n", 
    //       G_Voxel[0], G_Voxel[1], G_Voxel[2]);
    ret = CDM::E_CDM_ERROR_INVALID_DIVNUM;
    dfi->m_comm = comm;
    dfi->m_indexDfiName = DfiName;
    dfi->m_RankID = RankID;
    return dfi;
  }

#if 0
  if( dfi->m_start_type == E_CDM_SAMEDIV_SAMERES ) {
    printf("***** SAMEDIV_SAMERES\n");
  } else if( dfi->m_start_type == E_CDM_SAMEDIV_REFINEMENT ) {
    printf("***** SAMEDIV_REFINEMENT\n");
  } else if( dfi->m_start_type == E_CDM_DIFFDIV_SAMERES ) {
    printf("***** DIFFDIV_SAMERES\n");
  } else if( dfi->m_start_type == E_CDM_DIFFDIV_REFINEMENT ) {
    printf("***** DIFFDIV_REFINEMENT\n");
  }
#endif

  dfi->m_comm = comm;
  dfi->m_indexDfiName = DfiName;
  dfi->m_RankID = RankID;

  ret = CDM::E_CDM_SUCCESS;

  return dfi;

}

// #################################################################
// cdm_FileInfoクラスのポインタ取得
const
cdm_FileInfo* cdm_DFI::GetcdmFileInfo()
{
  return &DFI_Finfo;
}

// #################################################################
void 
cdm_DFI::SetcdmFilePath(cdm_FilePath FPath)
{
  DFI_Fpath = FPath;
}

// #################################################################
// cdm_FilePathクラスのポインタ取得
const
cdm_FilePath* cdm_DFI::GetcdmFilePath()
{
  return &DFI_Fpath;
}

// #################################################################
void
cdm_DFI::SetcdmVisIt(cdm_VisIt Visit)
{
  DFI_VisIt = Visit;
}

// #################################################################
const
cdm_VisIt* cdm_DFI::GetcdmVisIt()
{
  return &DFI_VisIt;
}


// #################################################################
// cdm_Unitクラスのポインタ取得
const
cdm_Unit* cdm_DFI::GetcdmUnit()
{
  return &DFI_Unit;
}

// #################################################################
//
void cdm_DFI::SetcdmUnit(cdm_Unit unit)
{
  DFI_Unit = unit;
}

// #################################################################
// cdm_Domainクラスのポインタ取得
const
cdm_Domain* cdm_DFI::GetcdmDomain()
{
  return DFI_Domain;
}

// #################################################################
//
void 
cdm_DFI::SetcdmDomain(cdm_Domain* domain)
{
  DFI_Domain = domain;
}


// #################################################################
// cdm_MPIクラスのポインタ取得
const
cdm_MPI* cdm_DFI::GetcdmMPI()
{
  return &DFI_MPI;
}

// #################################################################
//
void 
cdm_DFI::SetcdmMPI(cdm_MPI mpi)
{
  DFI_MPI = mpi;
}

// #################################################################
// cdm_TimeSliceクラスのポインタ取得
const
cdm_TimeSlice* cdm_DFI::GetcdmTimeSlice()
{
  return &DFI_TimeSlice;
}

// #################################################################
void
cdm_DFI::SetcdmTimeSlice(cdm_TimeSlice TSlice)
{
  DFI_TimeSlice = TSlice;
}

// #################################################################
// cdm_Processクラスのポインタ取得
const
cdm_Process* cdm_DFI::GetcdmProcess()
{
  return &DFI_Process;
}

// #################################################################
//
void
cdm_DFI::SetcdmProcess(cdm_Process Process)
{
  DFI_Process = Process;
}

// #################################################################
// 配列形状を文字列で返す
std::string cdm_DFI::GetArrayShapeString()
{
  if( DFI_Finfo.ArrayShape == CDM::E_CDM_IJKN ) return D_CDM_IJNK;
  if( DFI_Finfo.ArrayShape == CDM::E_CDM_NIJK ) return D_CDM_NIJK;
  return " ";
}

// #################################################################
// 配列形状を返す(e_num番号)
CDM::E_CDM_ARRAYSHAPE cdm_DFI::GetArrayShape()
{
  return (CDM::E_CDM_ARRAYSHAPE)DFI_Finfo.ArrayShape;
}

// #################################################################
// データタイプの取り出し(文字列)
std::string cdm_DFI::GetDataTypeString()
{
  return ConvDatatypeE2S((CDM::E_CDM_DTYPE)DFI_Finfo.DataType);
}

// #################################################################
// データタイプの取り出し(e_num)
CDM::E_CDM_DTYPE cdm_DFI::GetDataType()
{
  return (CDM::E_CDM_DTYPE)DFI_Finfo.DataType;
}

// #################################################################
// dfi種別の取り出し(文字列)
std::string cdm_DFI::GetDFITypeString()
{
  if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) return D_CDM_DFITYPE_CARTESIAN;
  if( DFI_Finfo.DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) return D_CDM_DFITYPE_NON_UNIFORM_CARTESIAN;
  return "";
}

// #################################################################
// dfi種別の取り出し(e_num)
CDM::E_CDM_DFITYPE cdm_DFI::GetDFIType()
{
  return (CDM::E_CDM_DFITYPE)DFI_Finfo.DFIType;
}

//FCONV 20140123.s
// #################################################################
// FileFormatの取り出し(文字列)
std::string cdm_DFI::GetFileFormatString()
{
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_UNKNOWN ) return "";
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_SPH ) return "sph";
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_BOV ) return "bov";
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_AVS ) return "avs";
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_PLOT3D ) return "plot3d";
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_VTK ) return "vtk";
  return "";
}

// #################################################################
// FileFormatの取り出し(e_num)
CDM::E_CDM_FORMAT cdm_DFI::GetFileFormat()
{
  return (CDM::E_CDM_FORMAT)DFI_Finfo.FileFormat;
}


// #################################################################
// 変数の個数の取り出し
int cdm_DFI::GetNumVariables()
{
  return DFI_Finfo.NumVariables;
}

// #################################################################
// 仮想セル数の取り出し
int cdm_DFI::GetNumGuideCell()
{
  return DFI_Finfo.GuideCell;
}

// #################################################################
// データタイプを文字列からe_num番号へ変換
CDM::E_CDM_DTYPE cdm_DFI::ConvDatatypeS2E(const std::string datatype)
{

  if     ( !strcasecmp(datatype.c_str(),"Int8"   ) ) return CDM::E_CDM_INT8;
  else if( !strcasecmp(datatype.c_str(),"Int16"  ) ) return CDM::E_CDM_INT16;
  else if( !strcasecmp(datatype.c_str(),"Int32"  ) ) return CDM::E_CDM_INT32;
  else if( !strcasecmp(datatype.c_str(),"Int64"  ) ) return CDM::E_CDM_INT64;
  else if( !strcasecmp(datatype.c_str(),"UInt8"  ) ) return CDM::E_CDM_UINT8;
  else if( !strcasecmp(datatype.c_str(),"UInt16" ) ) return CDM::E_CDM_UINT16;
  else if( !strcasecmp(datatype.c_str(),"UInt32" ) ) return CDM::E_CDM_UINT32;
  else if( !strcasecmp(datatype.c_str(),"UInt64" ) ) return CDM::E_CDM_UINT64;
  else if( !strcasecmp(datatype.c_str(),"Float32") ) return CDM::E_CDM_FLOAT32;
  else if( !strcasecmp(datatype.c_str(),"Float64") ) return CDM::E_CDM_FLOAT64;

  return CDM::E_CDM_DTYPE_UNKNOWN;
}

// #################################################################
// データタイプをe_num番号から文字列に変換
std::string cdm_DFI::ConvDatatypeE2S(const CDM::E_CDM_DTYPE Dtype)
{
  if     ( Dtype == CDM::E_CDM_INT8    ) return D_CDM_INT8;
  else if( Dtype == CDM::E_CDM_INT16   ) return D_CDM_INT16;
  else if( Dtype == CDM::E_CDM_INT32   ) return D_CDM_INT32;
  else if( Dtype == CDM::E_CDM_INT64   ) return D_CDM_INT64;
  else if( Dtype == CDM::E_CDM_UINT8   ) return D_CDM_UINT8;
  else if( Dtype == CDM::E_CDM_UINT16  ) return D_CDM_UINT16;
  else if( Dtype == CDM::E_CDM_UINT32  ) return D_CDM_UINT32;
  else if( Dtype == CDM::E_CDM_UINT64  ) return D_CDM_UINT64;
  else if( Dtype == CDM::E_CDM_FLOAT32 ) return D_CDM_FLOAT32;
  else if( Dtype == CDM::E_CDM_FLOAT64 ) return D_CDM_FLOAT64;
  else return "dummy";

}
// #################################################################
// データサイズの取り出し
int cdm_DFI::get_cdm_Datasize(CDM::E_CDM_DTYPE Dtype)
{

  if     ( Dtype == CDM::E_CDM_INT8    ) return sizeof(char);
  else if( Dtype == CDM::E_CDM_INT16   ) return sizeof(short);
  else if( Dtype == CDM::E_CDM_INT32   ) return sizeof(int);
  else if( Dtype == CDM::E_CDM_INT64   ) return sizeof(long long);
  else if( Dtype == CDM::E_CDM_UINT8   ) return sizeof(unsigned char);
  else if( Dtype == CDM::E_CDM_UINT16  ) return sizeof(unsigned short);
  else if( Dtype == CDM::E_CDM_UINT32  ) return sizeof(unsigned int);
  else if( Dtype == CDM::E_CDM_UINT64  ) return sizeof(unsigned long long);
  else if( Dtype == CDM::E_CDM_FLOAT32 ) return sizeof(float);
  else if( Dtype == CDM::E_CDM_FLOAT64 ) return sizeof(double);
  else return 0;

}

// #################################################################
// DFI DomainのGlobalVoxelの取り出し
const int* cdm_DFI::GetDFIGlobalVoxel()
{
  return DFI_Domain->GlobalVoxel;
}

// #################################################################
// DFI DomainのGlobalDivisionの取り出し
const int* cdm_DFI::GetDFIGlobalDivision()
{
  return DFI_Domain->GlobalDivision;
}
// #################################################################
// Create Domain & Process  
void cdm_DFI::cdm_Create_dfiProcessInfo(const MPI_Comm comm, 
                                        cdm_Process &G_Process)
{

  cdm_Rank G_Rank;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  if( nrank > 1 ) {
    int *headtail = NULL;
    if( RankID == 0 ) {
      headtail = new int[6*nrank];
    }

    int sbuff[6];
    for(int i=0; i<3; i++) {
      sbuff[i]   = DFI_Process.RankList[RankID].HeadIndex[i];
      sbuff[i+3] = DFI_Process.RankList[RankID].TailIndex[i];
    }

    MPI_Gather(sbuff,6,MPI_INT,headtail,6,MPI_INT,0,comm);

    if( RankID == 0 ) {
      for(int i=0; i<nrank; i++) {
        G_Rank.RankID=i;
        for(int j=0; j<3; j++) {
          G_Rank.HeadIndex[j]=headtail[i*6+j];
          G_Rank.TailIndex[j]=headtail[i*6+j+3];
          G_Rank.VoxelSize[j]=G_Rank.TailIndex[j]-G_Rank.HeadIndex[j]+1;
        }
        G_Process.RankList.push_back(G_Rank);
      }
    }

    if ( RankID == 0 ) delete [] headtail;

  } else {
    G_Rank.RankID=0;
    for(int i=0; i<3; i++) {
      G_Rank.HeadIndex[i]=DFI_Process.RankList[0].HeadIndex[i]; 
      G_Rank.TailIndex[i]=DFI_Process.RankList[0].TailIndex[i];
      G_Rank.VoxelSize[i]=G_Rank.TailIndex[i]-G_Rank.HeadIndex[i]+1;
    }
    G_Process.RankList.push_back(G_Rank);
  }
}

// #################################################################
// 読込み判定
CDM::E_CDM_READTYPE cdm_DFI::CheckReadType(const int G_voxel[3], 
                                           const int DFI_GlobalVoxel[3],
                                           const int G_Div[3],
                                           const int DFI_GlobalDivision[3])
{

  bool isSameDiv=true;

  //分割数チェック
  for(int i=0; i<3; i++ ) {
    if( DFI_GlobalDivision[i] != G_Div[i] ) {
      isSameDiv = false;
    }
  }

  if( isSameDiv ) {
    if( G_voxel[0] == DFI_GlobalVoxel[0]   &&
        G_voxel[1] == DFI_GlobalVoxel[1]   &&
        G_voxel[2] == DFI_GlobalVoxel[2]   ) return CDM::E_CDM_SAMEDIV_SAMERES;

    if( G_voxel[0] == DFI_GlobalVoxel[0]*2 &&
        G_voxel[1] == DFI_GlobalVoxel[1]*2 &&
        G_voxel[2] == DFI_GlobalVoxel[2]*2 ) return CDM::E_CDM_SAMEDIV_REFINEMENT;
  } else {
    if( G_voxel[0] == DFI_GlobalVoxel[0]   &&
        G_voxel[1] == DFI_GlobalVoxel[1]   &&
        G_voxel[2] == DFI_GlobalVoxel[2]   ) return CDM::E_CDM_DIFFDIV_SAMERES;

    if( G_voxel[0] == DFI_GlobalVoxel[0]*2 &&
        G_voxel[1] == DFI_GlobalVoxel[1]*2 &&
        G_voxel[2] == DFI_GlobalVoxel[2]*2 ) return CDM::E_CDM_DIFFDIV_REFINEMENT;
  }
 
  return CDM::E_CDM_READTYPE_UNKNOWN;
}

// #################################################################
// 読込み範囲を求める
void cdm_DFI::CreateReadStartEnd(bool isSame,
                            const int head[3], 
                            const int tail[3], 
                            const int gc, 
                            const int DFI_head[3],
                            const int DFI_tail[3], 
                            const int DFI_gc, 
                            const CDM::E_CDM_READTYPE readflag, 
                            int copy_sta[3], 
                            int copy_end[3],
                            int read_sta[3],
                            int read_end[3])
{

  int src_head[3],src_tail[3],src_gc;
  if( !isSame ) {
  // 粗密のとき密に変換、ガイドセルは倍にする
    src_gc = DFI_gc*2;
    for(int i=0; i<3; i++) {
      src_head[i]=DFI_head[i]*2-1;
      src_tail[i]=DFI_tail[i]*2;
    }
  } else {
  // 粗密でない時各値をコピー
    src_gc = DFI_gc;
    for(int i=0; i<3; i++) {
      src_head[i]=DFI_head[i];
      src_tail[i]=DFI_tail[i];
    }
  }

//スタート、エンドをセット
  for(int i=0; i<3; i++) {
    copy_sta[i] = max(head[i],src_head[i]);
    copy_end[i] = min(tail[i],src_tail[i]);
    
    //仮想セルが読込みの実セル内のときの処理（スタート）
    if( copy_sta[i] == 1 ) {
      copy_sta[i] -= min(gc,src_gc);
    } else if( head[i]>src_head[i] ) {
      copy_sta[i] = max(head[i]-gc,src_head[i]);
    }
 
    //仮想セルが読込みの実セル内のときの処理
    if( ( isSame  && copy_end[i] == DFI_Domain->GlobalVoxel[i] )  ||
        ( !isSame && copy_end[i] == DFI_Domain->GlobalVoxel[i]*2 ) ) {
      copy_end[i] += min(gc,src_gc);
    } else if( tail[i]<src_tail[i] ) {
      copy_end[i] = min(tail[i]+gc,src_tail[i]);
    }

    //read satrt/end のセット
    if( !isSame ) {
      if( copy_sta[i]>0 ) read_sta[i] = (copy_sta[i]+1)/2;
      else                read_sta[i] =  copy_sta[i]/2;

      if( copy_end[i]>0 ) read_end[i] = (copy_end[i]+1)/2;
      else                read_end[i] =  copy_end[i]/2;

    } else {
      read_sta[i] = copy_sta[i];
      read_end[i] = copy_end[i];
    }

  }
}

// #################################################################
// ファイル名を作成(ディレクトリパス付加）
std::string cdm_DFI::Generate_FieldFileName(int RankID, 
                                       int step, 
                                       const bool mio)
{

  if( DFI_Finfo.DirectoryPath.empty() ) return NULL;
  if( DFI_Finfo.Prefix.empty() ) return NULL;

  std::string fmt;
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_SPH ) {
    fmt=D_CDM_EXT_SPH;
  } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_BOV ) {
    fmt=D_CDM_EXT_BOV_DATAFILE;
//FCONV 20131122.s
  } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_AVS ) {
    //fmt=D_CDM_EXT_SPH;
    fmt=D_CDM_EXT_BOV_DATAFILE;
  } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_VTK ) {
    fmt=D_CDM_EXT_VTK;
  } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_PLOT3D ) {
    fmt=D_CDM_EXT_FUNC;
//FCONV 20131122.e
  }

  int len = DFI_Finfo.DirectoryPath.size() + DFI_Finfo.Prefix.size() + fmt.size() + 25; 
  // id(6) + step(10) + 1(\0) + "_"(2) + "."(1)+"id"(2)
  if( DFI_Finfo.TimeSliceDirFlag == CDM::E_CDM_ON ) len += 11;

  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);

  if( mio ) {
    if( DFI_Finfo.TimeSliceDirFlag == CDM::E_CDM_ON ) {
      sprintf(tmp, "%s/%010d/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),step,DFI_Finfo.Prefix.c_str(), 
            step,RankID,fmt.c_str());
    } else {
      sprintf(tmp, "%s/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,RankID,fmt.c_str());
    }
  } else {
    if( DFI_Finfo.TimeSliceDirFlag == CDM::E_CDM_ON ) {
      sprintf(tmp, "%s/%010d/%s_%010d.%s",DFI_Finfo.DirectoryPath.c_str(),step,DFI_Finfo.Prefix.c_str(), 
            step,fmt.c_str());
    } else {
      sprintf(tmp, "%s/%s_%010d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,fmt.c_str());
    }
  }
  
  std::string fname(tmp);
  if( tmp ) delete [] tmp;

  return fname;
}

//FCONV 20131128.s
// #################################################################
// ファイル名のみの生成（static 関数）
std::string cdm_DFI::Generate_FileName(std::string prefix,
                                       int RankID,
                                       int step,
                                       std::string ext,
                                       CDM::E_CDM_OUTPUT_FNAME output_fname,
                                       bool mio,
                                       CDM::E_CDM_ONOFF TimeSliceDirFlag)
{

  int len = prefix.size()+ext.size()+100;
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);

  //step出力なしのファイル名生成
  if( step < 0 ) 
  {
    if( mio ) {
      sprintf(tmp,"%s_id%06d.%s",prefix.c_str(),RankID,ext.c_str());
    } else {
      sprintf(tmp,"%s.%s",prefix.c_str(),ext.c_str());
    }
    std::string fname(tmp);
    if( tmp ) delete [] tmp;
    return fname;
  }

  //RankID出力なしのファイル名生成
  if( !mio ) {
    sprintf(tmp,"%s_%010d.%s",prefix.c_str(),step,ext.c_str());
    std::string fname(tmp);
    if( tmp ) delete [] tmp;
    return fname;
  }

  //step_rank
  if( output_fname != CDM::E_CDM_FNAME_RANK_STEP ) 
  {
    if( TimeSliceDirFlag == CDM::E_CDM_ON ) {
      sprintf(tmp,"%010d/%s_%010d_id%06d.%s",step,prefix.c_str(),step,RankID,ext.c_str());
    } else {
      sprintf(tmp,"%s_%010d_id%06d.%s",prefix.c_str(),step,RankID,ext.c_str());
    }
  } else if( output_fname == CDM::E_CDM_FNAME_RANK_STEP ) 
  {
  //rank_step
    if( TimeSliceDirFlag == CDM::E_CDM_ON ) {
      sprintf(tmp,"%010d/%s_id%06d_%010d.%s",step,prefix.c_str(),RankID,step,ext.c_str());
    } else {
      sprintf(tmp,"%s_id%06d_%010d.%s",prefix.c_str(),RankID,step,ext.c_str());
    }
  }

  std::string fname(tmp);
  if( tmp ) delete [] tmp;
  return fname;

}
//FCONV 20131128.s


// #################################################################
// ディレクトリがなければ作成、既存なら何もしない
int cdm_DFI::MakeDirectory(const std::string path)
{
  int ret = MakeDirectorySub(path);
  if( ret != 0 )
  {
     // 既存以外のエラー
     if ( EEXIST != errno )
     {
        printf( "\tError(errno)=[%s]\n", strerror(errno) );
        return 0;
     }
  }

  // failed
  return 1;       
}

// #################################################################
// ディレクトリがなければ作成、既存なら何もしない
int cdm_DFI::MakeDirectoryPath()
{
  // DirectoryPath with TimeSlice
  std::string path = Generate_Directory_Path();

  return MakeDirectory(path);
}

// #################################################################
int cdm_DFI::MakeDirectorySub( std::string path )
{

  umask(022);

  int ret = mkdir(path.c_str(), 0777);
  if( ret != 0 )
  {
    if( errno == EEXIST ) return 0;

    std::string parent = CDM::cdmPath_DirName(path);
    int ret2 = MakeDirectorySub( parent );
    if( ret2 != 0 )
    {
      return ret2;
    }
    ret = MakeDirectorySub( path );
  }

  return ret;

}

// #################################################################
// Directoryパスを生成する関数
std::string cdm_DFI::Generate_Directory_Path()
{

  // dfiのパスとDirectoryPathを連結する関数
  // ただし、絶対パスのときはdfiのパスは無視
  // CDM::cdmPath_isAbsoluteがtrueのとき絶対パス
  // DirectoryPath + TimeSliceDir
  std::string path = m_directoryPath;
  if( DFI_Finfo.TimeSliceDirFlag == CDM::E_CDM_ON )
  {
    //path = CDM::cdmPath_ConnectPath(path, m_timeSliceDir);
    path = CDM::cdmPath_ConnectPath(path, "");
  }

  // absolute path
  if( CDM::cdmPath_isAbsolute(path) )
  {
    return path;
  }

  // relative path
  std::string dfidir = CDM::cdmPath_DirName(m_indexDfiName);
  path = CDM::cdmPath_ConnectPath(dfidir, path);
  return path;

}

// #################################################################
// 出力DFIファイル名を作成する
std::string cdm_DFI::Generate_DFI_Name(const std::string prefix)
{

  // directory path
  std::string dirName  = CDM::cdmPath_DirName(prefix);

  // file extension
  std::string dfiname = CDM::cdmPath_FileName(prefix,".dfi");

  // filename
  std::string fname = CDM::cdmPath_ConnectPath( dirName, dfiname );

#if 0 // for debug
  printf("prefix    =%s\n", prefix.c_str() );
  printf("  dirName =%s\n", dirName.c_str() );
  printf("  dfiname =%s\n", dfiname.c_str() );
  printf("  fname   =%s\n", fname.c_str() );
  printf("\n");
#endif

  return fname;
}

// #################################################################
// Unitを追加する
void cdm_DFI::AddUnit(const std::string Name,
                      const std::string Unit,
                      const double reference,
                      const double difference,
                      const bool BsetDiff)
{

  /** UnitElemの生成 */
  cdm_UnitElem unit =  cdm_UnitElem(Name,Unit,reference,difference,BsetDiff);
  
  /** UnilListへのセット */
  DFI_Unit.UnitList.insert(map<std::string,cdm_UnitElem>::value_type(Name,unit));

}

// #################################################################
// UuitElemを取得する
CDM::E_CDM_ERRORCODE cdm_DFI::GetUnitElem(const std::string Name,
                                          cdm_UnitElem &unit)
{
  return DFI_Unit.GetUnitElem(Name, unit);
}

// #################################################################
// UnitElemのメンバ変数毎に取得する
CDM::E_CDM_ERRORCODE cdm_DFI::GetUnit(const std::string Name,
                                      std::string &unit,
                                      double &ref,
                                      double &diff,
                                      bool &bSetDiff)
{
  return DFI_Unit.GetUnit(Name, unit, ref, diff, bSetDiff);
}


// #################################################################
// TimeSlice OnOff フラグをセットする
void cdm_DFI::SetTimeSliceFlag(const CDM::E_CDM_ONOFF ONOFF)
{
  DFI_Finfo.TimeSliceDirFlag = ONOFF;
}

// #################################################################
// FileInfoの変数名を登録する
void cdm_DFI::setVariableName(int pvari, std::string variName)
{

  DFI_Finfo.setVariableName(pvari, variName);

}

// #################################################################
// FileInfoの変数名を取得する
std::string cdm_DFI::getVariableName(int pvari)
{

  return DFI_Finfo.getVariableName(pvari);

}

// #################################################################
// DFIに出力されているminmaxの合成値を取得
CDM::E_CDM_ERRORCODE cdm_DFI::getVectorMinMax(const unsigned step,
                                              double &vec_min,
                                              double &vec_max)
{

  return DFI_TimeSlice.getVectorMinMax(step,vec_min,vec_max);

}

// #################################################################
// DFIに出力されているminmaxの合成値を取得
CDM::E_CDM_ERRORCODE cdm_DFI::getMinMax(const unsigned step,
                                        const int variNo,
                                        double &min_value,
                                        double &max_value)
{

  return DFI_TimeSlice.getMinMax(step,variNo,min_value,max_value);

}

// #################################################################
// 読込みランクリストの作成
CDM::E_CDM_ERRORCODE
cdm_DFI::CheckReadRank(const cdm_Domain* dfi_domain,
                       const int head[3],
                       const int tail[3],
                       CDM::E_CDM_READTYPE readflag,
                       vector<int> &readRankList)
{

  return DFI_Process.CheckReadRank(dfi_domain,head,tail,readflag,readRankList);

}
