/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI.C
 * @brief  cio_DFI Class
 * @author aics    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K
#include "cio_DFI_SPH.h"
#include "cio_DFI_BOV.h"
//FCONV 20131122.s
#include "cio_DFI_AVS.h"
#include "cio_DFI_PLOT3D.h"
#include "cio_DFI_VTK.h"
//FCONV 20131122.e

// #################################################################
// コンストラクタ
cio_DFI::cio_DFI()
{

 m_read_type = CIO::E_CIO_READTYPE_UNKNOWN;
 m_RankID = 0;

 m_output_type = CIO::E_CIO_OUTPUT_TYPE_DEFAULT;
 m_output_fname = CIO::E_CIO_FNAME_DEFAULT;

}


// #################################################################
// デストラクタ
cio_DFI::~cio_DFI()
{

}

// #################################################################
// DFI read インスタンス
cio_DFI* cio_DFI::ReadInit(const MPI_Comm comm, 
                           const std::string DfiName,
                           const int G_Voxel[3],
                           const int G_Div[3],
                           CIO::E_CIO_ERRORCODE &ret)
{

  /** DFIのディレクトリパスの取得 */
  std::string dirName = CIO::cioPath_DirName(DfiName);

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  cio_TextParser tpCntl;

  /** index.dfi read */
  /** TPインスタンス */
  tpCntl.getTPinstance();

  FILE*fp = NULL;
  if( !(fp=fopen(DfiName.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",DfiName.c_str());
    ret = CIO::E_CIO_ERROR_READ_INDEXFILE_OPENERROR;
    return NULL;
  }
  fclose(fp);

  /** 入力ファイル index.dfi をセット */
  int ierror = 0;
  ierror = tpCntl.readTPfile(DfiName);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",DfiName.c_str());
    ret = CIO::E_CIO_ERROR_TEXTPARSER;
    return NULL;
  }

  /** Fileinfoの読込み */
  cio_FileInfo F_info;
  if( F_info.Read(tpCntl) != CIO::E_CIO_SUCCESS ) 
  {
    printf("\tFileInfo Data Read error %s\n",DfiName.c_str());
    ret = CIO::E_CIO_ERROR_READ_FILEINFO;
    return NULL;
  }

  /** FilePathの読込み */
  cio_FilePath F_path;
  if( F_path.Read(tpCntl) != CIO::E_CIO_SUCCESS )
  {
    printf("\tFilePath Data Read error %s\n",DfiName.c_str());
    ret = CIO::E_CIO_ERROR_READ_FILEPATH;
    return NULL;
  }

  /** Unitの読込み */
  cio_Unit unit;
  if( unit.Read(tpCntl) != CIO::E_CIO_SUCCESS )
  {
    printf("\tUnit Data Read error %s\n",DfiName.c_str());
    ret = CIO::E_CIO_ERROR_READ_UNIT;
    return NULL;
  }

  /** TimeSliceの読込み */
  cio_TimeSlice TimeSlice;
  if( TimeSlice.Read(tpCntl) != CIO::E_CIO_SUCCESS )
  {
    printf("\tTimeSlice Data Read error %s\n",DfiName.c_str());
    ret = CIO::E_CIO_ERROR_READ_TIMESLICE;
    return NULL;
  }

  /** TextParserの破棄 */
  tpCntl.remove();

  /** proc.dfi file name の取得 */
  std::string dfiname  = CIO::cioPath_FileName(F_path.ProcDFIFile,".dfi");
  std::string procfile = CIO::cioPath_ConnectPath(dirName,dfiname);

  /** proc.dfi read */
  /** TPインスタンス */
  tpCntl.getTPinstance();

  fp = NULL;
  if( !(fp=fopen(procfile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",procfile.c_str());
    ret = CIO::E_CIO_ERROR_READ_PROCFILE_OPENERROR;
    return NULL;
  }
  fclose(fp);

  /** 入力ファイル proc.dfi をセット */
  ierror = tpCntl.readTPfile(procfile);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",procfile.c_str());
    ret = CIO::E_CIO_ERROR_TEXTPARSER;
    return NULL;
  }

  /** Domainの読込み */
  cio_Domain domain;
  if( domain.Read(tpCntl) != CIO::E_CIO_SUCCESS ) 
  {
    printf("\tDomain Data Read error %s\n",procfile.c_str());
    ret = CIO::E_CIO_ERROR_READ_DOMAIN;
    return NULL;
  }

  /** MPIの読込み */
  cio_MPI mpi;
  if( mpi.Read(tpCntl,domain) != CIO::E_CIO_SUCCESS )
  {
    printf("\tMPI Data Read error %s\n",procfile.c_str());
    ret = CIO::E_CIO_ERROR_READ_MPI;
    return NULL;
  }

  /** Processの読込み */
  cio_Process process;
  if( process.Read(tpCntl) != CIO::E_CIO_SUCCESS )
  {
    printf("\tProcess Data Read error %s\n",procfile.c_str());
    ret = CIO::E_CIO_ERROR_READ_PROCESS;
    return NULL;
  }

  /** TextParserの破棄 */
  tpCntl.remove();

  /** dfiのインスタンス **/
  cio_DFI *dfi = NULL;
  if( F_info.FileFormat == CIO::E_CIO_FMT_SPH ) {
    dfi = new cio_DFI_SPH(F_info, F_path, unit, domain, mpi, TimeSlice, process);
  } else if( F_info.FileFormat == CIO::E_CIO_FMT_BOV ) {
    dfi = new cio_DFI_BOV(F_info, F_path, unit, domain, mpi, TimeSlice, process);
  } else {
    return NULL;
  }

  //読込みタイプのチェック
  dfi->m_read_type = dfi->CheckReadType(G_Voxel,dfi->DFI_Domain.GlobalVoxel,
                                         G_Div,dfi->DFI_Domain.GlobalDivision);
  if( dfi->m_read_type == CIO::E_CIO_READTYPE_UNKNOWN ) {
    //printf("\tDimension size error (%d %d %d)\n", 
    //       G_Voxel[0], G_Voxel[1], G_Voxel[2]);
    ret = CIO::E_CIO_ERROR_INVALID_DIVNUM;
    dfi->m_comm = comm;
    dfi->m_indexDfiName = DfiName;
    dfi->m_RankID = RankID;
    return dfi;
  }

#if 0
  if( dfi->m_start_type == E_CIO_SAMEDIV_SAMERES ) {
    printf("***** SAMEDIV_SAMERES\n");
  } else if( dfi->m_start_type == E_CIO_SAMEDIV_REFINEMENT ) {
    printf("***** SAMEDIV_REFINEMENT\n");
  } else if( dfi->m_start_type == E_CIO_DIFFDIV_SAMERES ) {
    printf("***** DIFFDIV_SAMERES\n");
  } else if( dfi->m_start_type == E_CIO_DIFFDIV_REFINEMENT ) {
    printf("***** DIFFDIV_REFINEMENT\n");
  }
#endif

  dfi->m_comm = comm;
  dfi->m_indexDfiName = DfiName;
  dfi->m_RankID = RankID;

  ret = CIO::E_CIO_SUCCESS;

  return dfi;

}

// #################################################################
// cio_FileInfoクラスのポインタ取得
const
cio_FileInfo* cio_DFI::GetcioFileInfo()
{
  return &DFI_Finfo;
}

// #################################################################
void 
cio_DFI::SetcioFilePath(cio_FilePath FPath)
{
  DFI_Fpath = FPath;
}

// #################################################################
// cio_FilePathクラスのポインタ取得
const
cio_FilePath* cio_DFI::GetcioFilePath()
{
  return &DFI_Fpath;
}

// #################################################################
// cio_Unitクラスのポインタ取得
const
cio_Unit* cio_DFI::GetcioUnit()
{
  return &DFI_Unit;
}

// #################################################################
//
void cio_DFI::SetcioUnit(cio_Unit unit)
{
  DFI_Unit = unit;
}

// #################################################################
// cio_Domainクラスのポインタ取得
const
cio_Domain* cio_DFI::GetcioDomain()
{
  return &DFI_Domain;
}

// #################################################################
//
void 
cio_DFI::SetcioDomain(cio_Domain domain)
{
  DFI_Domain = domain;
}


// #################################################################
// cio_MPIクラスのポインタ取得
const
cio_MPI* cio_DFI::GetcioMPI()
{
  return &DFI_MPI;
}

// #################################################################
//
void 
cio_DFI::SetcioMPI(cio_MPI mpi)
{
  DFI_MPI = mpi;
}

// #################################################################
// cio_TimeSliceクラスのポインタ取得
const
cio_TimeSlice* cio_DFI::GetcioTimeSlice()
{
  return &DFI_TimeSlice;
}

// #################################################################
void
cio_DFI::SetcioTimeSlice(cio_TimeSlice TSlice)
{
  DFI_TimeSlice = TSlice;
}

// #################################################################
// cio_Processクラスのポインタ取得
const
cio_Process* cio_DFI::GetcioProcess()
{
  return &DFI_Process;
}

// #################################################################
//
void
cio_DFI::SetcioProcess(cio_Process Process)
{
  DFI_Process = Process;
}

// #################################################################
// DFI Write インスタンス float 型
cio_DFI* cio_DFI::WriteInit(const MPI_Comm comm,
                            const std::string DfiName,
                            const std::string Path,
                            const std::string prefix,
                            const CIO::E_CIO_FORMAT format,
                            const int GCell,
                            const CIO::E_CIO_DTYPE DataType,
                            const CIO::E_CIO_ARRAYSHAPE ArrayShape, 
                            const int nComp,
                            const std::string proc_fname,
                            const int G_size[3],
                            const float pitch[3],
                            const float G_origin[3],
                            const int division[3],
                            const int head[3],
                            const int tail[3],
                            const std::string hostname,
                            const CIO::E_CIO_ONOFF TSliceOnOff)
{

  // float型をdouble型に変換してdouble版WriteInit関数を呼ぶ
 
  double d_pch[3],d_org[3];
  for(int i=0; i<3; i++) {
    d_pch[i]=(double)pitch[i];
    d_org[i]=(double)G_origin[i];
  }

  return WriteInit(comm, 
                   DfiName, 
                   Path, 
                   prefix, 
                   format, 
                   GCell, 
                   DataType,
                   ArrayShape, 
                   nComp, 
                   proc_fname,
                   G_size, 
                   d_pch,
                   d_org, 
                   division, 
                   head, 
                   tail, 
                   hostname,
                   TSliceOnOff);

}

// #################################################################
// DFI Write インスタンス double 型
cio_DFI* cio_DFI::WriteInit(const MPI_Comm comm,
                            const std::string DfiName,
                            const std::string Path,
                            const std::string prefix,
                            const CIO::E_CIO_FORMAT format,
                            const int GCell,
                            const CIO::E_CIO_DTYPE DataType,
                            const CIO::E_CIO_ARRAYSHAPE ArrayShape,
                            const int nComp,
                            const std::string proc_fname,
                            const int G_size[3],
                            const double pitch[3],
                            const double G_origin[3],
                            const int division[3],
                            const int head[3],
                            const int tail[3],
                            const std::string hostname,
                            const CIO::E_CIO_ONOFF TSliceOnOff)
{

//FCONV 20140131.s
  if( format == CIO::E_CIO_FMT_SPH ) {
    if( nComp > 1 && ArrayShape == CIO::E_CIO_IJKN ) {
      printf("\tCIO error sph file undefined ijkn component>1.\n");
      return NULL;
    }
  }
//FCONV 20140131.e

  cio_DFI *dfi = NULL;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  cio_FileInfo out_F_info;
  out_F_info.DirectoryPath    = Path;
  out_F_info.TimeSliceDirFlag = TSliceOnOff;
  out_F_info.Prefix           = prefix;
  out_F_info.FileFormat       = format;
  out_F_info.GuideCell        = GCell;
  out_F_info.DataType         = DataType;
  out_F_info.ArrayShape       = ArrayShape;
  out_F_info.Component        = nComp;

  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  if( cdumy[0] == 0x01 ) out_F_info.Endian = CIO::E_CIO_LITTLE;
  if( cdumy[0] == 0x00 ) out_F_info.Endian = CIO::E_CIO_BIG;

  cio_FilePath out_F_path;
  out_F_path.ProcDFIFile = proc_fname;

  cio_Unit out_unit;

  cio_MPI out_mpi;
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  cio_Domain out_domain;
  cio_Process out_Process;
  cio_Rank out_Rank;

  for(int i=0; i<nrank; i++ ) {
     out_Process.RankList.push_back(out_Rank);
  }

  out_Process.RankList[RankID].RankID=RankID;
  out_Process.RankList[RankID].HostName=hostname;
  for(int i=0; i<3; i++) {
    out_Process.RankList[RankID].HeadIndex[i]=head[i];
    out_Process.RankList[RankID].TailIndex[i]=tail[i];
    out_Process.RankList[RankID].VoxelSize[i]=tail[i]-head[i]+1;
  }

  for(int i=0; i<3; i++) {
    out_domain.GlobalVoxel[i]  = G_size[i];
    out_domain.GlobalDivision[i] = division[i];
    out_domain.GlobalOrigin[i] = G_origin[i];
    out_domain.GlobalRegion[i] = pitch[i]*G_size[i];
  }

  cio_TimeSlice out_TSlice;

  char tmpname[512];
  memset(tmpname,0x00,sizeof(char)*512);
  if( gethostname(tmpname, 512) != 0 ) printf("*** error gethostname() \n");

  if( out_F_info.FileFormat == CIO::E_CIO_FMT_SPH ) {
    dfi = new cio_DFI_SPH(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_Process);
  } else if( out_F_info.FileFormat == CIO::E_CIO_FMT_BOV ) {
    dfi = new cio_DFI_BOV(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_Process);
//FCONV 20131122.s
  } else if( out_F_info.FileFormat == CIO::E_CIO_FMT_AVS ) {
    dfi = new cio_DFI_AVS(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_Process);
  } else if( out_F_info.FileFormat == CIO::E_CIO_FMT_PLOT3D ) {
    dfi = new cio_DFI_PLOT3D(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_Process);
  } else if( out_F_info.FileFormat == CIO::E_CIO_FMT_VTK ) {
    dfi = new cio_DFI_VTK(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_Process);
//FCONV 20131122.e
  } else return NULL;


  dfi->m_indexDfiName = DfiName;
  dfi->m_directoryPath = CIO::cioPath_DirName(DfiName);
  dfi->m_comm = comm;
  dfi->m_RankID = RankID;

  return dfi;

}

// #################################################################
// 配列形状を文字列で返す
std::string cio_DFI::GetArrayShapeString()
{
  if( DFI_Finfo.ArrayShape == CIO::E_CIO_IJKN ) return D_CIO_IJNK;
  if( DFI_Finfo.ArrayShape == CIO::E_CIO_NIJK ) return D_CIO_NIJK;
  return " ";
}

// #################################################################
// 配列形状を返す(e_num番号)
CIO::E_CIO_ARRAYSHAPE cio_DFI::GetArrayShape()
{
  return (CIO::E_CIO_ARRAYSHAPE)DFI_Finfo.ArrayShape;
}

// #################################################################
// データタイプの取り出し(文字列)
std::string cio_DFI::GetDataTypeString()
{
  return ConvDatatypeE2S((CIO::E_CIO_DTYPE)DFI_Finfo.DataType);
}

// #################################################################
// データタイプの取り出し(e_num)
CIO::E_CIO_DTYPE cio_DFI::GetDataType()
{
  return (CIO::E_CIO_DTYPE)DFI_Finfo.DataType;
}

//FCONV 20140123.s
// #################################################################
// FileFormatの取り出し(文字列)
std::string cio_DFI::GetFileFormatString()
{
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_UNKNOWN ) return "";
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_SPH ) return "sph";
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_BOV ) return "bov";
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_AVS ) return "avs";
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_PLOT3D ) return "plot3d";
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_VTK ) return "vtk";
  return "";
}

// #################################################################
// FileFormatの取り出し(e_num)
CIO::E_CIO_FORMAT cio_DFI::GetFileFormat()
{
  return (CIO::E_CIO_FORMAT)DFI_Finfo.FileFormat;
}


// #################################################################
// 成分数の取り出し
int cio_DFI::GetNumComponent()
{
  return DFI_Finfo.Component;
}

// #################################################################
// 仮想セル数の取り出し
int cio_DFI::GetNumGuideCell()
{
  return DFI_Finfo.GuideCell;
}

// #################################################################
// データタイプを文字列からe_num番号へ変換
CIO::E_CIO_DTYPE cio_DFI::ConvDatatypeS2E(const std::string datatype)
{

  if     ( !strcasecmp(datatype.c_str(),"Int8"   ) ) return CIO::E_CIO_INT8;
  else if( !strcasecmp(datatype.c_str(),"Int16"  ) ) return CIO::E_CIO_INT16;
  else if( !strcasecmp(datatype.c_str(),"Int32"  ) ) return CIO::E_CIO_INT32;
  else if( !strcasecmp(datatype.c_str(),"Int64"  ) ) return CIO::E_CIO_INT64;
  else if( !strcasecmp(datatype.c_str(),"UInt8"  ) ) return CIO::E_CIO_UINT8;
  else if( !strcasecmp(datatype.c_str(),"UInt16" ) ) return CIO::E_CIO_UINT16;
  else if( !strcasecmp(datatype.c_str(),"UInt32" ) ) return CIO::E_CIO_UINT32;
  else if( !strcasecmp(datatype.c_str(),"UInt64" ) ) return CIO::E_CIO_UINT64;
  else if( !strcasecmp(datatype.c_str(),"Float32") ) return CIO::E_CIO_FLOAT32;
  else if( !strcasecmp(datatype.c_str(),"Float64") ) return CIO::E_CIO_FLOAT64;

  return CIO::E_CIO_DTYPE_UNKNOWN;
}

// #################################################################
// データタイプをe_num番号から文字列に変換
std::string cio_DFI::ConvDatatypeE2S(const CIO::E_CIO_DTYPE Dtype)
{
  if     ( Dtype == CIO::E_CIO_INT8    ) return D_CIO_INT8;
  else if( Dtype == CIO::E_CIO_INT16   ) return D_CIO_INT16;
  else if( Dtype == CIO::E_CIO_INT32   ) return D_CIO_INT32;
  else if( Dtype == CIO::E_CIO_INT64   ) return D_CIO_INT64;
  else if( Dtype == CIO::E_CIO_UINT8   ) return D_CIO_UINT8;
  else if( Dtype == CIO::E_CIO_UINT16  ) return D_CIO_UINT16;
  else if( Dtype == CIO::E_CIO_UINT32  ) return D_CIO_UINT32;
  else if( Dtype == CIO::E_CIO_UINT64  ) return D_CIO_UINT64;
  else if( Dtype == CIO::E_CIO_FLOAT32 ) return D_CIO_FLOAT32;
  else if( Dtype == CIO::E_CIO_FLOAT64 ) return D_CIO_FLOAT64;
  else return "dummy";

}
// #################################################################
// データサイズの取り出し
int cio_DFI::get_cio_Datasize(CIO::E_CIO_DTYPE Dtype)
{

  if     ( Dtype == CIO::E_CIO_INT8    ) return sizeof(char);
  else if( Dtype == CIO::E_CIO_INT16   ) return sizeof(short);
  else if( Dtype == CIO::E_CIO_INT32   ) return sizeof(int);
  else if( Dtype == CIO::E_CIO_INT64   ) return sizeof(long long);
  else if( Dtype == CIO::E_CIO_UINT8   ) return sizeof(unsigned char);
  else if( Dtype == CIO::E_CIO_UINT16  ) return sizeof(unsigned short);
  else if( Dtype == CIO::E_CIO_UINT32  ) return sizeof(unsigned int);
  else if( Dtype == CIO::E_CIO_UINT64  ) return sizeof(unsigned long long);
  else if( Dtype == CIO::E_CIO_FLOAT32 ) return sizeof(float);
  else if( Dtype == CIO::E_CIO_FLOAT64 ) return sizeof(double);
  else return 0;

}

// #################################################################
// DFI DomainのGlobalVoxelの取り出し
int* cio_DFI::GetDFIGlobalVoxel()
{
  return DFI_Domain.GlobalVoxel;
}

// #################################################################
// DFI DomainのGlobalDivisionの取り出し
int* cio_DFI::GetDFIGlobalDivision()
{
  return DFI_Domain.GlobalDivision;
}
// #################################################################
// Create Domain & Process  
void cio_DFI::cio_Create_dfiProcessInfo(const MPI_Comm comm, 
                                        cio_Process &G_Process)
{

  cio_Rank G_Rank;

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
CIO::E_CIO_READTYPE cio_DFI::CheckReadType(const int G_voxel[3], 
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
        G_voxel[2] == DFI_GlobalVoxel[2]   ) return CIO::E_CIO_SAMEDIV_SAMERES;

    if( G_voxel[0] == DFI_GlobalVoxel[0]*2 &&
        G_voxel[1] == DFI_GlobalVoxel[1]*2 &&
        G_voxel[2] == DFI_GlobalVoxel[2]*2 ) return CIO::E_CIO_SAMEDIV_REFINEMENT;
  } else {
    if( G_voxel[0] == DFI_GlobalVoxel[0]   &&
        G_voxel[1] == DFI_GlobalVoxel[1]   &&
        G_voxel[2] == DFI_GlobalVoxel[2]   ) return CIO::E_CIO_DIFFDIV_SAMERES;

    if( G_voxel[0] == DFI_GlobalVoxel[0]*2 &&
        G_voxel[1] == DFI_GlobalVoxel[1]*2 &&
        G_voxel[2] == DFI_GlobalVoxel[2]*2 ) return CIO::E_CIO_DIFFDIV_REFINEMENT;
  }
 
  return CIO::E_CIO_READTYPE_UNKNOWN;
}

// #################################################################
// 読込み範囲を求める
void cio_DFI::CreateReadStartEnd(bool isSame,
                            const int head[3], 
                            const int tail[3], 
                            const int gc, 
                            const int DFI_head[3],
                            const int DFI_tail[3], 
                            const int DFI_gc, 
                            const CIO::E_CIO_READTYPE readflag, 
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
    if( ( isSame  && copy_end[i] == DFI_Domain.GlobalVoxel[i] )  ||
        ( !isSame && copy_end[i] == DFI_Domain.GlobalVoxel[i]*2 ) ) {
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
std::string cio_DFI::Generate_FieldFileName(int RankID, 
                                       int step, 
                                       const bool mio)
{

  if( DFI_Finfo.DirectoryPath.empty() ) return NULL;
  if( DFI_Finfo.Prefix.empty() ) return NULL;

  std::string fmt;
  if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_SPH ) {
    fmt=D_CIO_EXT_SPH;
  } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_BOV ) {
    fmt=D_CIO_EXT_BOV;
//FCONV 20131122.s
  } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_AVS ) {
    //fmt=D_CIO_EXT_SPH;
    fmt=D_CIO_EXT_BOV;
  } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_VTK ) {
    fmt=D_CIO_EXT_VTK;
  } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_PLOT3D ) {
    fmt=D_CIO_EXT_FUNC;
//FCONV 20131122.e
  }

  int len = DFI_Finfo.DirectoryPath.size() + DFI_Finfo.Prefix.size() + fmt.size() + 25; 
  // id(6) + step(10) + 1(\0) + "_"(2) + "."(1)+"id"(2)
  if( DFI_Finfo.TimeSliceDirFlag == CIO::E_CIO_ON ) len += 11;

  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);

  if( mio ) {
    if( DFI_Finfo.TimeSliceDirFlag == CIO::E_CIO_ON ) {
      sprintf(tmp, "%s/%010d/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),step,DFI_Finfo.Prefix.c_str(), 
            step,RankID,fmt.c_str());
    } else {
      sprintf(tmp, "%s/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,RankID,fmt.c_str());
    }
  } else {
    if( DFI_Finfo.TimeSliceDirFlag == CIO::E_CIO_ON ) {
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
std::string cio_DFI::Generate_FileName(std::string prefix,
                                       int RankID,
                                       int step,
                                       std::string ext,
                                       CIO::E_CIO_OUTPUT_FNAME output_fname,
                                       bool mio,
                                       CIO::E_CIO_ONOFF TimeSliceDirFlag)
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
  if( output_fname != CIO::E_CIO_FNAME_RANK_STEP ) 
  {
    if( TimeSliceDirFlag == CIO::E_CIO_ON ) {
      sprintf(tmp,"%010d/%s_%010d_id%06d.%s",step,prefix.c_str(),step,RankID,ext.c_str());
    } else {
      sprintf(tmp,"%s_%010d_id%06d.%s",prefix.c_str(),step,RankID,ext.c_str());
    }
  } else if( output_fname == CIO::E_CIO_FNAME_RANK_STEP ) 
  {
  //rank_step
    if( TimeSliceDirFlag == CIO::E_CIO_ON ) {
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
int cio_DFI::MakeDirectory(const std::string path)
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
int cio_DFI::MakeDirectoryPath()
{
  // DirectoryPath with TimeSlice
  std::string path = Generate_Directory_Path();

  return MakeDirectory(path);
}

// #################################################################
int cio_DFI::MakeDirectorySub( std::string path )
{

  umask(022);

  int ret = mkdir(path.c_str(), 0777);
  if( ret != 0 )
  {
    if( errno == EEXIST ) return 0;

    std::string parent = CIO::cioPath_DirName(path);
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
std::string cio_DFI::Generate_Directory_Path()
{

  // dfiのパスとDirectoryPathを連結する関数
  // ただし、絶対パスのときはdfiのパスは無視
  // CIO::cioPath_isAbsoluteがtrueのとき絶対パス
  // DirectoryPath + TimeSliceDir
  std::string path = m_directoryPath;
  if( DFI_Finfo.TimeSliceDirFlag == CIO::E_CIO_ON )
  {
    //path = CIO::cioPath_ConnectPath(path, m_timeSliceDir);
    path = CIO::cioPath_ConnectPath(path, "");
  }

  // absolute path
  if( CIO::cioPath_isAbsolute(path) )
  {
    return path;
  }

  // relative path
  std::string dfidir = CIO::cioPath_DirName(m_indexDfiName);
  path = CIO::cioPath_ConnectPath(dfidir, path);
  return path;

}

// #################################################################
// 出力DFIファイル名を作成する
std::string cio_DFI::Generate_DFI_Name(const std::string prefix)
{

  // directory path
  std::string dirName  = CIO::cioPath_DirName(prefix);

  // file extension
  std::string dfiname = CIO::cioPath_FileName(prefix,".dfi");

  // filename
  std::string fname = CIO::cioPath_ConnectPath( dirName, dfiname );

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
void cio_DFI::AddUnit(const std::string Name,
                      const std::string Unit,
                      const double reference,
                      const double difference,
                      const bool BsetDiff)
{

  /** UnitElemの生成 */
  cio_UnitElem unit =  cio_UnitElem(Name,Unit,reference,difference,BsetDiff);
  
  /** UnilListへのセット */
  DFI_Unit.UnitList.insert(map<std::string,cio_UnitElem>::value_type(Name,unit));

}

// #################################################################
// UuitElemを取得する
CIO::E_CIO_ERRORCODE cio_DFI::GetUnitElem(const std::string Name,
                                          cio_UnitElem &unit)
{
  return DFI_Unit.GetUnitElem(Name, unit);
}

// #################################################################
// UnitElemのメンバ変数毎に取得する
CIO::E_CIO_ERRORCODE cio_DFI::GetUnit(const std::string Name,
                                      std::string &unit,
                                      double &ref,
                                      double &diff,
                                      bool &bSetDiff)
{
  return DFI_Unit.GetUnit(Name, unit, ref, diff, bSetDiff);
}


// #################################################################
// TimeSlice OnOff フラグをセットする
void cio_DFI::SetTimeSliceFlag(const CIO::E_CIO_ONOFF ONOFF)
{
  DFI_Finfo.TimeSliceDirFlag = ONOFF;
}

// #################################################################
// FileInfoの成分名を登録する
void cio_DFI::setComponentVariable(int pcomp, std::string compName)
{

  DFI_Finfo.setComponentVariable(pcomp, compName);

}

// #################################################################
// FileInfoの成分名を取得する
std::string cio_DFI::getComponentVariable(int pcomp)
{

  return DFI_Finfo.getComponentVariable(pcomp);

}

// #################################################################
// DFIに出力されているminmaxの合成値を取得
CIO::E_CIO_ERRORCODE cio_DFI::getVectorMinMax(const unsigned step,
                                              double &vec_min,
                                              double &vec_max)
{

  return DFI_TimeSlice.getVectorMinMax(step,vec_min,vec_max);

}

// #################################################################
// DFIに出力されているminmaxの合成値を取得
CIO::E_CIO_ERRORCODE cio_DFI::getMinMax(const unsigned step,
                                        const int compNo,
                                        double &min_value,
                                        double &max_value)
{

  return DFI_TimeSlice.getMinMax(step,compNo,min_value,max_value);

}

// #################################################################
// 読込みランクリストの作成
CIO::E_CIO_ERRORCODE
cio_DFI::CheckReadRank(cio_Domain dfi_domain,
                       const int head[3],
                       const int tail[3],
                       CIO::E_CIO_READTYPE readflag,
                       vector<int> &readRankList)
{

  return DFI_Process.CheckReadRank(dfi_domain,head,tail,readflag,readRankList);

}
