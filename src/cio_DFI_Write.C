/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_Write.C
 * @brief  cio_DFI Class
 * @author aics     
 */

#include "cio_DFI.h"

// #################################################################
// Index DFIファイルの出力
CIO::E_CIO_ERRORCODE
cio_DFI::WriteIndexDfiFile(const std::string dfi_name)
{

  if ( dfi_name.empty() ) return CIO::E_CIO_ERROR_WRITE_INDEXFILENAME_EMPTY;
  if ( DFI_Finfo.Prefix.empty() ) return CIO::E_CIO_ERROR_WRITE_PREFIX_EMPTY;

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
    return CIO::E_CIO_ERROR_WRITE_INDEXFILE_OPENERROR;
  }

  //FileInfo {} の出力
  if( DFI_Finfo.Write(fp, 0) != CIO::E_CIO_SUCCESS )
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FILEINFO;
  }

  //FilePath {} の出力
  if( DFI_Fpath.Write(fp, 1) != CIO::E_CIO_SUCCESS )
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FILEPATH;
  }


  //Unit {} の出力
  if( DFI_Unit.Write(fp, 0) != CIO::E_CIO_SUCCESS ) 
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_UNIT;
  }

  //TimeSlice {} の出力
  if ( DFI_TimeSlice.Write(fp, 1) != CIO::E_CIO_SUCCESS ) 
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_TIMESLICE;
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// proc DFIファイルの出力コントロール (float 版)
/*
CIO::E_CIO_ERRORCODE
cio_DFI::WriteProcDfiFile(const MPI_Comm comm,
                          bool out_host,
                          float* org)
{

  //orign の再設定
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
CIO::E_CIO_ERRORCODE
cio_DFI::WriteProcDfiFile(const MPI_Comm comm,
                          bool out_host)
                          //double* org)
{

  //procファイル名の生成
  std::string procFileName = CIO::cioPath_DirName(m_indexDfiName)+"/"+CIO::cioPath_FileName(DFI_Fpath.ProcDFIFile,".dfi");

  if( procFileName.empty() ) return CIO::E_CIO_ERROR_WRITE_PROCFILENAME_EMPTY;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  cio_MPI out_mpi;
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  cio_Domain out_domain;
  cio_Process out_Process;

  //出力するProcess情報の生成
  cio_Create_dfiProcessInfo(comm, out_Process);

  //orign の設定
  /*
  if( org!=NULL ) {
    for(int i=0; i<3; i++) {
      out_domain.GlobalOrigin[i] = org[i];
    }
  } else {
  */
    for(int i=0; i<3; i++) {
      out_domain.GlobalOrigin[i] = DFI_Domain.GlobalOrigin[i];
    }
  //}

  //Domain の設定
  for(int i=0; i<3; i++) {
    out_domain.GlobalVoxel[i]    = DFI_Domain.GlobalVoxel[i];
    out_domain.GlobalDivision[i] = DFI_Domain.GlobalDivision[i];
    out_domain.GlobalRegion[i]   = DFI_Domain.GlobalRegion[i];
  }

  //ホスト名出力指示ありの時、各ランクのホスト名を集める
  if( out_host ) {
    const int LEN=256;
    char *recbuf = new char[out_Process.RankList.size()*LEN];
    char  sedbuf[LEN];
    //sprintf(sedbuf,"%s",hostname.c_str());
    sprintf(sedbuf,"%s",DFI_Process.RankList[RankID].HostName.c_str());
    MPI_Gather(sedbuf,LEN,MPI_CHAR,recbuf,LEN,MPI_CHAR,0,MPI_COMM_WORLD);

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
      return CIO::E_CIO_ERROR_WRITE_PROCFILE_OPENERROR;
    }

    //Domain {} の出力
    if( out_domain.Write(fp, 0) != CIO::E_CIO_SUCCESS )
    {
      if (fp) fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_DOMAIN;
    }

    //MPI {} の出力
    if( out_mpi.Write(fp, 0) != CIO::E_CIO_SUCCESS )
    {
      fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_MPI;
    }

    //Process {} の出力
    if( out_Process.Write(fp, 0) != CIO::E_CIO_SUCCESS )
    {
      fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_PROCESS;
    }

    fclose(fp);
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// fileld data 出力
CIO::E_CIO_ERRORCODE
cio_DFI::WriteData(const unsigned step,
                   const int gc,
                   double time,
                   cio_Array* val,
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
  if( m_output_fname != CIO::E_CIO_FNAME_RANK_STEP ) {
    tmp = Generate_FieldFileName(m_RankID,step,mio);
    if( CIO::cioPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      outFile = tmp;
    } else {
      outFile = m_directoryPath + "/"+ tmp;
    }
  } else {
    std::string ext;
    if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_SPH ) {
      ext = D_CIO_EXT_SPH;
    } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_BOV ) {
      ext = D_CIO_EXT_BOV;
    } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_AVS ) {
      //ext = D_CIO_EXT_SPH;
      ext = D_CIO_EXT_BOV;
    } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_VTK ) {
      ext = D_CIO_EXT_VTK;
    } else if( DFI_Finfo.FileFormat == CIO::E_CIO_FMT_PLOT3D ) {
      ext = D_CIO_EXT_FUNC;
    }
    tmp = Generate_FileName(DFI_Finfo.Prefix,
                            m_RankID,
                            step,ext,
                            m_output_fname,
                            mio,
                            DFI_Finfo.TimeSliceDirFlag);
    if( CIO::cioPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      outFile = DFI_Finfo.DirectoryPath +"/"+ tmp;
    } else {
      outFile = m_directoryPath + "/" + DFI_Finfo.DirectoryPath +"/"+ tmp;
    }
  }
//FCONV 20131128.e

  std::string dir = CIO::cioPath_DirName(outFile);

  if( MakeDirectory(dir) != 1 ) return CIO::E_CIO_ERROR_MAKEDIRECTORY;

  cio_Array *outArray = val;
  if( gc != DFI_Finfo.GuideCell ) {
    //出力用バッファのインスタンス
    outArray = cio_Array::instanceArray
               ( DFI_Finfo.DataType
               , DFI_Finfo.ArrayShape
               , DFI_Process.RankList[m_RankID].VoxelSize
               , DFI_Finfo.GuideCell
               , DFI_Finfo.Component); 
    //配列のコピー val -> outArray
    int ret = val->copyArray(outArray);
  }

  // フィールドデータの出力
  CIO::E_CIO_ERRORCODE err = CIO::E_CIO_SUCCESS;
  err = WriteFieldData(outFile, step, time, outArray, avr_mode, step_avr, time_avr);

  //出力バッファのメモリ解放
  if( val != outArray ) {
    delete outArray;
  }

  if( err != CIO::E_CIO_SUCCESS ) return err;

//FCONV 20131218.s
  if( m_indexDfiName != "" ) {
    //index dfi ファイルのディレクトリ作成
    cio_DFI::MakeDirectory(m_directoryPath);
    std::string dfiname = CIO::cioPath_FileName(m_indexDfiName,".dfi");
    std::string fname = CIO::cioPath_ConnectPath( m_directoryPath, dfiname );

    //Slice へのセット
    DFI_TimeSlice.AddSlice(step, time, minmax, DFI_Finfo.Component, avr_mode,
                           step_avr, time_avr);

    //index dfi のファイル出力
    if( m_RankID == 0 ) {
      err = WriteIndexDfiFile(fname);
    }
  }
//FCONV 20131218.e
//FCONV 20131125.s
  if( !write_ascii_header(step,time) ) return CIO::E_CIO_ERROR;
//FCONV 20131125.e

  return err;
}

// #################################################################
// フィールドデータ出力
CIO::E_CIO_ERRORCODE
cio_DFI::WriteFieldData(std::string fname,
                        const unsigned step,
                        double time,
                        cio_Array *val,
                        const bool avr_mode,
                        const unsigned step_avr,
                        const double time_avr)
{

  FILE* fp;
  if( (fp = fopen(fname.c_str(),"wb")) == NULL ) {
    fprintf(stderr,"Can't open file.(%s)\n",fname.c_str());
    return CIO::E_CIO_ERROR_OPEN_FIELDDATA;
  }

  //printf("field file name : %s\n",fname.c_str());

  //ヘッダー出力
  if( write_HeaderRecord(fp, step, time, m_RankID) != CIO::E_CIO_SUCCESS ) {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FIELD_HEADER_RECORD;
  }

  cio_Array *outArray = val;

  //格子点補間処理ありの場合、図心データから格子点への補間を行う
  if( m_bgrid_interp_flag ) {
    //配列サイズの取得
    const int *szVal = val->getArraySizeInt();
    //配列成分の取得
    int nComp = val->getNcomp();
    //格子点データ配列サイズのセット
    int szOut[3];
    for(int i=0; i<3; i++) szOut[i]=szVal[i]+1;
    //出力バッファのインスタンス
    outArray =  cio_Array::instanceArray
                           (val->getDataType(),
                            val->getArrayShape(),
                            szOut,
                            0,
                            nComp);
    //unsigned char
    if( val->getDataType() == CIO::E_CIO_UINT8 ) {
      cio_TypeArray<unsigned char> *V = dynamic_cast<cio_TypeArray<unsigned char>*>(val);
      cio_TypeArray<unsigned char> *P = dynamic_cast<cio_TypeArray<unsigned char>*>(outArray);
      setGridData(P,V);

    //char
    }else if( val->getDataType() == CIO::E_CIO_INT8 ) {
      cio_TypeArray<char> *V = dynamic_cast<cio_TypeArray<char>*>(val);
      cio_TypeArray<char> *P = dynamic_cast<cio_TypeArray<char>*>(outArray);
      setGridData(P,V);

    //unsigned short
    } else if( val->getDataType() == CIO::E_CIO_UINT16 ) {
      cio_TypeArray<unsigned short> *V = dynamic_cast<cio_TypeArray<unsigned short>*>(val);
      cio_TypeArray<unsigned short> *P = dynamic_cast<cio_TypeArray<unsigned short>*>(outArray);
      setGridData(P,V);

    //short
    } else if( val->getDataType() == CIO::E_CIO_INT16 ) {
      cio_TypeArray<short> *V = dynamic_cast<cio_TypeArray<short>*>(val);
      cio_TypeArray<short> *P = dynamic_cast<cio_TypeArray<short>*>(outArray);
      setGridData(P,V);

    //uint
    } else if( val->getDataType() == CIO::E_CIO_UINT32 ) {
      cio_TypeArray<unsigned int> *V = dynamic_cast<cio_TypeArray<unsigned int>*>(val);
      cio_TypeArray<unsigned int> *P = dynamic_cast<cio_TypeArray<unsigned int>*>(outArray);
      setGridData(P,V);

    //int
    } else if( val->getDataType() == CIO::E_CIO_INT32 ) {
      cio_TypeArray<int> *V = dynamic_cast<cio_TypeArray<int>*>(val);
      cio_TypeArray<int> *P = dynamic_cast<cio_TypeArray<int>*>(outArray);
      setGridData(P,V);

    //ulong
    } else if( val->getDataType() == CIO::E_CIO_UINT64 ) {
      cio_TypeArray<unsigned long long> *V = dynamic_cast<cio_TypeArray<unsigned long long>*>(val);
      cio_TypeArray<unsigned long long> *P = dynamic_cast<cio_TypeArray<unsigned long long>*>(outArray);
      setGridData(P,V);

    //long
    } else if( val->getDataType() == CIO::E_CIO_INT64 ) {
      cio_TypeArray<long long> *V = dynamic_cast<cio_TypeArray<long long>*>(val);
      cio_TypeArray<long long> *P = dynamic_cast<cio_TypeArray<long long>*>(outArray);
      setGridData(P,V);

    //float
    } else if( val->getDataType() == CIO::E_CIO_FLOAT32 ) {
      cio_TypeArray<float> *V = dynamic_cast<cio_TypeArray<float>*>(val);
      cio_TypeArray<float> *P = dynamic_cast<cio_TypeArray<float>*>(outArray);
      setGridData(P,V);

    //double
    } else if( val->getDataType() == CIO::E_CIO_FLOAT64 ) {
      cio_TypeArray<double> *V = dynamic_cast<cio_TypeArray<double>*>(val);
      cio_TypeArray<double> *P = dynamic_cast<cio_TypeArray<double>*>(outArray);
      setGridData(P,V);

    }

  } 

  //データ出力
  //if( write_DataRecord(fp, val, DFI_Finfo.GuideCell, m_RankID) != CIO::E_CIO_SUCCESS) {
  if( write_DataRecord(fp, outArray, DFI_Finfo.GuideCell, m_RankID) != CIO::E_CIO_SUCCESS) {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FIELD_DATA_RECORD;
  }

  //average 出力
  if( !avr_mode ) {
    if( write_averaged(fp, step_avr, time_avr) != CIO::E_CIO_SUCCESS ) {
      fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_FIELD_AVERAGED_RECORD;
    }
  }

  fclose(fp);

  return CIO::E_CIO_SUCCESS;

}

