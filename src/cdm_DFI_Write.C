/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_Write.C
 * @brief  cdm_DFI Class
 * @author aics     
 */

#include "cdm_DFI.h"

// #################################################################
// Index DFIファイルの出力 (API関数)
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteIndexDfiFile()
{

  CDM::E_CDM_ERRORCODE err = CDM::E_CDM_SUCCESS;

  if( m_indexDfiName != "" ) {
    //index dfi ファイルのディレクトリ作成
    cdm_DFI::MakeDirectory(m_directoryPath);
    std::string dfiname = CDM::cdmPath_FileName(m_indexDfiName,".dfi");
    std::string fname = CDM::cdmPath_ConnectPath( m_directoryPath, dfiname );

    //index dfi のファイル出力
    if( m_RankID == 0 ) {
      err = WriteIndexDfiFile(fname);
    }
  } else {
    printf("\tError : dfi file name is not set\n");
    return CDM::E_CDM_ERROR_WRITE_INDEXFILENAME_EMPTY;
  }

  return err;

}

// #################################################################
// Index DFIファイルの出力
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteIndexDfiFile(const std::string dfi_name)
{

  if ( dfi_name.empty() ) return CDM::E_CDM_ERROR_WRITE_INDEXFILENAME_EMPTY;
  if ( DFI_Finfo.Prefix.empty() ) return CDM::E_CDM_ERROR_WRITE_PREFIX_EMPTY;

  FILE* fp = NULL;

  // File exist ?
  bool flag = false;
  if ( fp = fopen(dfi_name.c_str(), "r") )
  {
    flag = true;
    fclose(fp);
  }

  if( !(fp = fopen(dfi_name.c_str(), "w")) )
  {
    fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
    return CDM::E_CDM_ERROR_WRITE_INDEXFILE_OPENERROR;
  }

  //FileInfo {} の出力
  CDM::E_CDM_OUTPUT_FNAME orgFname = DFI_Finfo.FieldFilenameFormat;
  DFI_Finfo.FieldFilenameFormat = m_output_fname;
  if( DFI_Finfo.Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_FILEINFO;
  }
  DFI_Finfo.FieldFilenameFormat = orgFname;

  //FilePath {} の出力
  if( DFI_Fpath.Write(fp, 1) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_FILEPATH;
  }

  // VisIt用のオプション出力
  if( DFI_VisIt.Write(fp, 1) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_VISIT;
  }

  //Unit {} の出力
  if( DFI_Unit.Write(fp, 0) != CDM::E_CDM_SUCCESS ) 
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_UNIT;
  }

  //TimeSlice {} の出力
  if ( DFI_TimeSlice.Write(fp, 1, DFI_Finfo.FileFormat) != CDM::E_CDM_SUCCESS ) 
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_TIMESLICE;
  }

  //追加情報の出力
//20150918.NetCDF.s
#ifdef _WITH_NETCDF4_
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_NETCDF4 )
  {
    cdm_DFI_NETCDF *dfi_nc = (cdm_DFI_NETCDF*)this;
    if( dfi_nc->WriteAdditionalTP(fp, 0) != CDM::E_CDM_SUCCESS )
    {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_DFI_NETCDF;
    }
  }
#endif
//20150918.NetCDF.e


  fclose(fp);
  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// proc DFIファイルの出力コントロール (float 版)
/*
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteProcDfiFile(const MPI_Comm comm,
                          bool out_host,
                          float* org)
{

  //origin の再設定
  double d_org[3];
  if( org != NULL ) {
    for(int i=0; i<3; i++) {
      d_org[i]=(double)org[i];
    }
  } else {
    for(int i=0; i<3; i++) {
      d_org[i]=DFI_Domain.GlobalOrigin[i];
    }
  }

  return WriteProcDfiFile(comm, out_host, d_org);

}
*/
// #################################################################
// proc DFIファイルの出力コントロール
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteProcDfiFile(const MPI_Comm comm,
                          const bool out_host,
                          const int cell_id,
                          const int bcf_id)
                          //double* org)
{

  //procファイル名の生成
  std::string procFileName = CDM::cdmPath_DirName(m_indexDfiName)+"/"+CDM::cdmPath_FileName(DFI_Fpath.ProcDFIFile,".dfi");

  if( procFileName.empty() ) return CDM::E_CDM_ERROR_WRITE_PROCFILENAME_EMPTY;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  cdm_MPI out_mpi;
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  //CellID,境界IDをセット
  DFI_Process.RankList[RankID].c_id = cell_id;
  DFI_Process.RankList[RankID].bc_id = bcf_id;

  cdm_Process out_Process;

  //出力するProcess情報の生成
  cdm_Create_dfiProcessInfo(comm, out_Process);

  //ホスト名出力指示ありの時、各ランクのホスト名を集める
  if( out_host ) {
    const int LEN=256;
    char *recbuf = new char[out_Process.RankList.size()*LEN];
    char  sedbuf[LEN];
    //sprintf(sedbuf,"%s",hostname.c_str());
    sprintf(sedbuf,"%s",DFI_Process.RankList[RankID].HostName.c_str());
    MPI_Gather(sedbuf,LEN,MPI_CHAR,recbuf,LEN,MPI_CHAR,0,comm);

    for( int i=0; i<out_Process.RankList.size(); i++ ) {
     char* hn =&(recbuf[i*LEN]);
     out_Process.RankList[i].HostName=(std::string(hn));
    }

    if( recbuf ) delete [] recbuf;
  }

  //proc.dfの出力
  if( RankID == 0 ) {

    FILE* fp = NULL;
    if( !(fp = fopen(procFileName.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", procFileName.c_str());
      return CDM::E_CDM_ERROR_WRITE_PROCFILE_OPENERROR;
    }

    //Domain {} の出力
    if( DFI_Domain->Write(fp, 0) != CDM::E_CDM_SUCCESS )
    {
      if (fp) fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_DOMAIN;
    }

    //MPI {} の出力
    if( out_mpi.Write(fp, 0) != CDM::E_CDM_SUCCESS )
    {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_MPI;
    }

    //Process {} の出力
    if( out_Process.Write(fp, 0) != CDM::E_CDM_SUCCESS )
    {
      fclose(fp);
      return CDM::E_CDM_ERROR_WRITE_PROCESS;
    }

    fclose(fp);
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// grid ファイル出力
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteGridFile(const int* iblank)
{

  bool flag;
  //ファイルフォーマットチェック。gridファイルがあるのは、PLOT3DとAVSのみ。
  if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_PLOT3D )
  {
    flag = write_GridData(iblank);
    if ( !flag ) return CDM::E_CDM_ERROR_WRITE_GRIDFILE;
  }
  else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_AVS )
  {
    std::cout << "Output cod file (Not supported yet)" << std::endl;
    return CDM::E_CDM_ERROR_WRITE_GRIDFILE;
  }
  else
  {
    printf("\tError : File format \"%s\" has no grid file.\n", GetFileFormatString().c_str());
    return CDM::E_CDM_ERROR_WRITE_GRIDFILE;
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// fileld data 出力
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteData(const unsigned step,
                   const int gc,
                   double time,
                   cdm_Array* val,
                   double* minmax,
                   const bool avr_mode,
                   const unsigned step_avr,
                   double time_avr)
{

  //printf("WriteData RankID : %d\n",m_RankID);

  bool mio=false;
  if( DFI_MPI.NumberOfRank > 1 ) mio=true;

  std::string outFile,tmp;
//FCONV 20131128.s
//if( m_output_fname != CDM::E_CDM_FNAME_RANK_STEP ) {
  if( m_output_fname != CDM::E_CDM_FNAME_RANK_STEP && m_output_fname != CDM::E_CDM_FNAME_RANK ) {
    tmp = Generate_FieldFileName(m_RankID,step,mio);
    if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      outFile = tmp;
    } else {
      outFile = m_directoryPath + "/"+ tmp;
    }
  } else {
    std::string ext;
    if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_SPH ) {
      ext = D_CDM_EXT_SPH;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_BOV ) {
      ext = D_CDM_EXT_BOV_DATAFILE;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_AVS ) {
      //ext = D_CDM_EXT_SPH;
      ext = D_CDM_EXT_BOV_DATAFILE;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_VTK ) {
      ext = D_CDM_EXT_VTK;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_PLOT3D ) {
      ext = D_CDM_EXT_FUNC;
//20150918.NetCDF.s
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_NETCDF4 ) {
      ext = D_CDM_EXT_NC;
//20150918.NetCDF.e
//20160331.fub.s
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_FUB ) {
      ext = D_CDM_EXT_FUB;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_FUB_COD ) {
      ext = D_CDM_EXT_XYZ;
//20160331.fub.e
    }
    tmp = Generate_FileName(DFI_Finfo.Prefix,
                            m_RankID,
                            step,ext,
                            m_output_fname,
                            mio,
                            DFI_Finfo.TimeSliceDirFlag,
                            DFI_Finfo.RankNoPrefix);
    if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      outFile = DFI_Finfo.DirectoryPath +"/"+ tmp;
    } else {
      outFile = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
    }
  }
//FCONV 20131128.e

  std::string dir = CDM::cdmPath_DirName(outFile);

  if( MakeDirectory(dir) != 1 ) return CDM::E_CDM_ERROR_MAKEDIRECTORY;

  cdm_Array *outArray = val;
  if( gc > DFI_Finfo.GuideCell ) {
    //出力用バッファのインスタンス
    outArray = cdm_Array::instanceArray
               ( DFI_Finfo.DataType
               , DFI_Finfo.ArrayShape
               , DFI_Process.RankList[m_RankID].VoxelSize
               , DFI_Finfo.GuideCell
               , DFI_Finfo.NumVariables); 
    //配列のコピー val -> outArray
    int ret = val->copyArray(outArray);
  }
  else if( gc < DFI_Finfo.GuideCell )
  {
    //出力用に用意したデータのガイドセル値より出力するガイドセル値の方が大きい場合はエラー
    printf("\tError : Number of guide cells %d %d\n", gc, DFI_Finfo.GuideCell);
    return CDM::E_CDM_ERROR_NUM_OF_GUIDECELLS;
  }

  // フィールドデータの出力
  CDM::E_CDM_ERRORCODE err = CDM::E_CDM_SUCCESS;
  err = WriteFieldData(outFile, step, time, outArray, avr_mode, step_avr, time_avr);

  //出力バッファのメモリ解放
  if( val != outArray ) {
    delete outArray;
  }

  if( err != CDM::E_CDM_SUCCESS ) return err;

//FCONV 20131218.s
  if( m_indexDfiName != "" ) {
    //index dfi ファイルのディレクトリ作成
    cdm_DFI::MakeDirectory(m_directoryPath);
    std::string dfiname = CDM::cdmPath_FileName(m_indexDfiName,".dfi");
    std::string fname = CDM::cdmPath_ConnectPath( m_directoryPath, dfiname );

    //Slice へのセット
    bool bExist = false;
    for( int i=0;i<DFI_TimeSlice.SliceList.size();i++ )
    {
      if( DFI_TimeSlice.SliceList[i].step == step )
      {
        bExist = true;
        break;
      }
    }
    if( !bExist )
    {
      DFI_TimeSlice.AddSlice(step, time, minmax, DFI_Finfo.NumVariables, DFI_Finfo.FileFormat,
                             avr_mode, step_avr, time_avr);
    }

    //index dfi のファイル出力
    if( m_RankID == 0 ) {
      err = WriteIndexDfiFile(fname);
    }
  }
//FCONV 20131218.e
//FCONV 20131125.s
  if( !write_ascii_header(step,time) ) return CDM::E_CDM_ERROR;
//FCONV 20131125.e

  return err;
}

// #################################################################
// fileld data 出力(dfi fileの出力なし)
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteFieldDataFile(const unsigned step,
                            const int gc,
                            double time,
                            cdm_Array* val,
                            const bool avr_mode,
                            const unsigned step_avr,
                            double time_avr)
{

  bool mio=false;
  if( DFI_MPI.NumberOfRank > 1 ) mio=true;

  std::string outFile,tmp;
//FCONV 20131128.s
  if( m_output_fname != CDM::E_CDM_FNAME_RANK_STEP ) {
    tmp = Generate_FieldFileName(m_RankID,step,mio);
    if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      outFile = tmp;
    } else {
      outFile = m_directoryPath + "/"+ tmp;
    }
  } else {
    std::string ext;
    if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_SPH ) {
      ext = D_CDM_EXT_SPH;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_BOV ) {
      ext = D_CDM_EXT_BOV_DATAFILE;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_AVS ) {
      //ext = D_CDM_EXT_SPH;
      ext = D_CDM_EXT_BOV_DATAFILE;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_VTK ) {
      ext = D_CDM_EXT_VTK;
    } else if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_PLOT3D ) {
      ext = D_CDM_EXT_FUNC;
    }
    tmp = Generate_FileName(DFI_Finfo.Prefix,
                            m_RankID,
                            step,ext,
                            m_output_fname,
                            mio,
                            DFI_Finfo.TimeSliceDirFlag,
                            DFI_Finfo.RankNoPrefix);
    if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      outFile = DFI_Finfo.DirectoryPath +"/"+ tmp;
    } else {
      outFile = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
    }
  }
//FCONV 20131128.e

  std::string dir = CDM::cdmPath_DirName(outFile);

  if( MakeDirectory(dir) != 1 ) return CDM::E_CDM_ERROR_MAKEDIRECTORY;

  cdm_Array *outArray = val;
  if( gc > DFI_Finfo.GuideCell ) {
    //出力用バッファのインスタンス
    outArray = cdm_Array::instanceArray
               ( DFI_Finfo.DataType
               , DFI_Finfo.ArrayShape
               , DFI_Process.RankList[m_RankID].VoxelSize
               , DFI_Finfo.GuideCell
               , DFI_Finfo.NumVariables); 
    //配列のコピー val -> outArray
    int ret = val->copyArray(outArray);
  }
  else if( gc < DFI_Finfo.GuideCell )
  {
    //出力用に用意したデータのガイドセル値より出力するガイドセル値の方が大きい場合はエラー
    printf("\tError : Number of guide cells %d %d\n", gc, DFI_Finfo.GuideCell);
    return CDM::E_CDM_ERROR_NUM_OF_GUIDECELLS;
  }

  // フィールドデータの出力
  CDM::E_CDM_ERRORCODE err = CDM::E_CDM_SUCCESS;
  err = WriteFieldData(outFile, step, time, outArray, avr_mode, step_avr, time_avr);

  //出力バッファのメモリ解放
  if( val != outArray ) {
    delete outArray;
  }

  if( err != CDM::E_CDM_SUCCESS ) return err;

//FCONV 20131125.s
  if( !write_ascii_header(step,time) ) return CDM::E_CDM_ERROR;
//FCONV 20131125.e

  return err;
}

// #################################################################
// フィールドデータ出力
CDM::E_CDM_ERRORCODE
cdm_DFI::WriteFieldData(std::string fname,
                        const unsigned step,
                        double time,
                        cdm_Array *val,
                        const bool avr_mode,
                        const unsigned step_avr,
                        const double time_avr)
{
  //追記モードにするかどうかをチェック(NetCDF対応)
  bool addMode = CheckAddWriteMode();

  // ファイルオープン
//FILE* fp;
//if( (fp = fopen(fname.c_str(),"wb")) == NULL ) {
  cdm_FILE* fp;
  if( (fp = cdm_FILE::OpenWriteBinary(fname,DFI_Finfo.FileFormat, addMode)) == NULL ) {
    fprintf(stderr,"Can't open file.(%s)\n",fname.c_str());
    return CDM::E_CDM_ERROR_OPEN_FIELDDATA;
  }

  //printf("field file name : %s\n",fname.c_str());

  //ヘッダー出力
  if( write_HeaderRecord(fp, step, time, m_RankID) != CDM::E_CDM_SUCCESS ) {
//  fclose(fp);
    cdm_FILE::CloseFile(fp);
    return CDM::E_CDM_ERROR_WRITE_FIELD_HEADER_RECORD;
  }

  cdm_Array *outArray = val;

  //格子点補間処理ありの場合、図心データから格子点への補間を行う
  if( m_bgrid_interp_flag ) {
    //配列サイズの取得
    const int *szVal = val->getArraySizeInt();
    //配列における変数の個数の取得
    int nVari = val->getNvari();
    //格子点データ配列サイズのセット
    int szOut[3];
    for(int i=0; i<3; i++) szOut[i]=szVal[i]+1;
    //出力バッファのインスタンス
    outArray =  cdm_Array::instanceArray
                           (val->getDataType(),
                            val->getArrayShape(),
                            szOut,
                            0,
                            nVari);
    //unsigned char
    if( val->getDataType() == CDM::E_CDM_UINT8 ) {
      cdm_TypeArray<unsigned char> *V = dynamic_cast<cdm_TypeArray<unsigned char>*>(val);
      cdm_TypeArray<unsigned char> *P = dynamic_cast<cdm_TypeArray<unsigned char>*>(outArray);
      setGridData(P,V);

    //char
    }else if( val->getDataType() == CDM::E_CDM_INT8 ) {
      cdm_TypeArray<char> *V = dynamic_cast<cdm_TypeArray<char>*>(val);
      cdm_TypeArray<char> *P = dynamic_cast<cdm_TypeArray<char>*>(outArray);
      setGridData(P,V);

    //unsigned short
    } else if( val->getDataType() == CDM::E_CDM_UINT16 ) {
      cdm_TypeArray<unsigned short> *V = dynamic_cast<cdm_TypeArray<unsigned short>*>(val);
      cdm_TypeArray<unsigned short> *P = dynamic_cast<cdm_TypeArray<unsigned short>*>(outArray);
      setGridData(P,V);

    //short
    } else if( val->getDataType() == CDM::E_CDM_INT16 ) {
      cdm_TypeArray<short> *V = dynamic_cast<cdm_TypeArray<short>*>(val);
      cdm_TypeArray<short> *P = dynamic_cast<cdm_TypeArray<short>*>(outArray);
      setGridData(P,V);

    //uint
    } else if( val->getDataType() == CDM::E_CDM_UINT32 ) {
      cdm_TypeArray<unsigned int> *V = dynamic_cast<cdm_TypeArray<unsigned int>*>(val);
      cdm_TypeArray<unsigned int> *P = dynamic_cast<cdm_TypeArray<unsigned int>*>(outArray);
      setGridData(P,V);

    //int
    } else if( val->getDataType() == CDM::E_CDM_INT32 ) {
      cdm_TypeArray<int> *V = dynamic_cast<cdm_TypeArray<int>*>(val);
      cdm_TypeArray<int> *P = dynamic_cast<cdm_TypeArray<int>*>(outArray);
      setGridData(P,V);

    //ulong
    } else if( val->getDataType() == CDM::E_CDM_UINT64 ) {
      cdm_TypeArray<unsigned long long> *V = dynamic_cast<cdm_TypeArray<unsigned long long>*>(val);
      cdm_TypeArray<unsigned long long> *P = dynamic_cast<cdm_TypeArray<unsigned long long>*>(outArray);
      setGridData(P,V);

    //long
    } else if( val->getDataType() == CDM::E_CDM_INT64 ) {
      cdm_TypeArray<long long> *V = dynamic_cast<cdm_TypeArray<long long>*>(val);
      cdm_TypeArray<long long> *P = dynamic_cast<cdm_TypeArray<long long>*>(outArray);
      setGridData(P,V);

    //float
    } else if( val->getDataType() == CDM::E_CDM_FLOAT32 ) {
      cdm_TypeArray<float> *V = dynamic_cast<cdm_TypeArray<float>*>(val);
      cdm_TypeArray<float> *P = dynamic_cast<cdm_TypeArray<float>*>(outArray);
      setGridData(P,V);

    //double
    } else if( val->getDataType() == CDM::E_CDM_FLOAT64 ) {
      cdm_TypeArray<double> *V = dynamic_cast<cdm_TypeArray<double>*>(val);
      cdm_TypeArray<double> *P = dynamic_cast<cdm_TypeArray<double>*>(outArray);
      setGridData(P,V);

    }

  } 

  //データ出力
  //if( write_DataRecord(fp, val, DFI_Finfo.GuideCell, m_RankID) != CDM::E_CDM_SUCCESS) {
  if( write_DataRecord(fp, outArray, DFI_Finfo.GuideCell, m_RankID) != CDM::E_CDM_SUCCESS) {
//  fclose(fp);
    cdm_FILE::CloseFile(fp);
    return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;
  }

  //average 出力
  if( !avr_mode ) {
    if( write_averaged(fp, step_avr, time_avr) != CDM::E_CDM_SUCCESS ) {
//    fclose(fp);
      cdm_FILE::CloseFile(fp);
      return CDM::E_CDM_ERROR_WRITE_FIELD_AVERAGED_RECORD;
    }
  }

//fclose(fp);
  cdm_FILE::CloseFile(fp);

  return CDM::E_CDM_SUCCESS;

}

