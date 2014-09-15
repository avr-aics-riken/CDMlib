/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_DFI_Read.C
 * @brief  cdm_DFI Class
 * @author aics    
 */

#include "cdm_DFI.h"

// #################################################################
// フィールドデータの読込み（引数で渡された配列にデータを読込み返す）
CDM::E_CDM_ERRORCODE
cdm_DFI::ReadData(cdm_Array *dst,
                  const unsigned step, 
                  const int gc, 
                  const int Gvoxel[3], 
                  const int Gdivision[3], 
                  const int head[3], 
                  const int tail[3],
                  double &time,
                  const bool avr_mode, 
                  unsigned &avr_step, 
                  double &avr_time)
{

  CDM::E_CDM_ERRORCODE ret;

  /** dts にHead/Tailをセット */

  int Shead[3];
  for(int i=0; i<3; i++) Shead[i] = head[i]; 
  dst->setHeadIndex(Shead);

  /** index DFIファイルの ディレクトリパスを取得 */
  std::string dir = CDM::cdmPath_DirName(m_indexDfiName);

  bool mio = false;             ///< DFIファイルの並列フラグ
  bool isSame =true;            ///< 粗密フラグ true:密 false:粗
  CDM::E_CDM_READTYPE readflag; ///<読込み判定フラグ

  /** 読込みフラグ取得 */
  readflag = CheckReadType(Gvoxel, DFI_Domain.GlobalVoxel,
                           Gdivision, DFI_Domain.GlobalDivision);

  /** 粗密フラグセット */
  if( readflag == CDM::E_CDM_SAMEDIV_REFINEMENT || readflag == CDM::E_CDM_DIFFDIV_REFINEMENT ) isSame = false; 

  /**読込みランクリストの生成 */
  ret = DFI_Process.CheckReadRank(DFI_Domain, head, tail, readflag, m_readRankList);
  if( ret != CDM::E_CDM_SUCCESS ) {
    printf("error code : %d\n",(int)ret);
    return ret;
  }

  if( DFI_Process.RankList.size() > 1 ) mio = true; ///<Processが１より大きい時並列

  for(int i=0; i<m_readRankList.size(); i++) {
    int n = m_readRankList[i];
    int ID= DFI_Process.RankList[n].RankID;

    /** ファイル名の生成 */
    std::string fname;
    if( CDM::cdmPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
      fname = Generate_FieldFileName(ID,step,mio);
    } else {
      std::string tmp = Generate_FieldFileName(ID,step,mio);
      fname = CDM::cdmPath_ConnectPath( dir, tmp );
    }

    int copy_sta[3],copy_end[3],read_sta[3],read_end[3];

    /** 読込み領域 start end の取得 */
    CreateReadStartEnd(isSame,head, tail, gc, DFI_Process.RankList[n].HeadIndex, 
                       DFI_Process.RankList[n].TailIndex,
                       DFI_Finfo.GuideCell, readflag, 
                       copy_sta, copy_end, read_sta, read_end);

    /** 読込み方法の取得 */

    /** フィールドデータの読込み */

    cdm_Array* src = ReadFieldData(fname, step, time, read_sta, read_end,
                                   DFI_Process.RankList[n].HeadIndex, 
                                   DFI_Process.RankList[n].TailIndex, 
                                   avr_mode, avr_step, avr_time, ret);
    if( ret != CDM::E_CDM_SUCCESS ) {
      delete src;
      return ret;
    }

    /** 読込めたファイル名の出力（ランク0のみ) */
    if( m_RankID == 0 ) {
      printf("\t[%s] has read :\tstep=%d  time=%e ]\n",fname.c_str(), step, time);
    }


    /** src にHead/Tailをセット */
    src->setHeadIndex(read_sta);

    /** 粗密処理 */
    if( !isSame ) {
      cdm_Array *temp = src;
      int err;
      src = cdm_Array::interp_coarse(temp,err,false);
      delete temp;
    }

    src->copyArray(copy_sta,copy_end,dst);
    delete src;

  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// フィールドデータの読込み
cdm_Array* cdm_DFI::ReadFieldData(std::string fname,
                                  const unsigned step,
                                  double &time,
                                  const int sta[3],
                                  const int end[3],
                                  const int DFI_head[3],
                                  const int DFI_tail[3],
                                  bool avr_mode,
                                  unsigned &avr_step,
                                  double &avr_time,
                                  CDM::E_CDM_ERRORCODE &ret )
{

  ret = CDM::E_CDM_SUCCESS;


  if( !fname.c_str() || !DFI_Finfo.Component ) {
    ret = CDM::E_CDM_ERROR_READ_FIELDDATA_FILE;
    return NULL;
  }

  /** ファイルオープン */
  FILE* fp;
  if( !(fp=fopen(fname.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",fname.c_str());
    ret = CDM::E_CDM_ERROR_OPEN_FIELDDATA;
    return NULL;
  }

  /** Endian セット */
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  CDM::E_CDM_ENDIANTYPE Endian=CDM::E_CDM_ENDIANTYPE_UNKNOWN;
  if( cdumy[0] == 0x01 ) Endian = CDM::E_CDM_LITTLE;
  if( cdumy[0] == 0x00 ) Endian = CDM::E_CDM_BIG;

  bool matchEndian = true;
  if( Endian != DFI_Finfo.Endian ) matchEndian = false;

  //RealType real_type;
  int voxsize[3];
  /** ヘッダーレコードの読込み */
  ret = read_HeaderRecord(fp, matchEndian, step, DFI_head, DFI_tail, 
                         DFI_Finfo.GuideCell, voxsize, time);
  if( ret != CDM::E_CDM_SUCCESS )
  {
    ret = CDM::E_CDM_ERROR_READ_FIELD_HEADER_RECORD;
    printf("**** read error\n");
    fclose(fp);
    return NULL;
  }

  int sz[3]; ///< voxsize - 2*gc : 実セル数
  for(int i=0; i<3; i++) sz[i]=voxsize[i]-2*DFI_Finfo.GuideCell;

  int szB[3],headB[3];
  for(int i=0; i<3; i++) {
    szB[i] = voxsize[i];
    if( m_extend_arraysize_flag ) { //格子点の場合、配列サイズを1大きくする
      szB[i]+=1;
    }
    headB[i] = DFI_head[i] - DFI_Finfo.GuideCell;
  }
  //１層ずつ読み込むので、バッファのZサイズは１にしておく
  szB[2]=1;

//FCONV 20121216.s
  //読み込みバッファ
  cdm_Array* buf=NULL;
  //配列形状がIJKNのときは成分数を１にしてインスタンスする
  if( DFI_Finfo.ArrayShape == CDM::E_CDM_NIJK ) {
    buf = cdm_Array::instanceArray
                   ( DFI_Finfo.DataType
                   , DFI_Finfo.ArrayShape
                   , szB
                   , 0 
                   , DFI_Finfo.Component );
  } else if( DFI_Finfo.ArrayShape == CDM::E_CDM_IJKN ) {
    buf = cdm_Array::instanceArray
                   ( DFI_Finfo.DataType
                   , DFI_Finfo.ArrayShape
                   , szB
                   , 0
                   , 1 );
  }
//FCONV 20121216.e

  int szS[3];
  int headS[3];
  for(int i=0; i<3; i++) {
    szS[i]=end[i]-sta[i]+1;
    if( m_extend_arraysize_flag ) { //格子点の場合、配列サイズを1大きくする
      szS[i]+=1;
    }
    headS[i]=sta[i];
  }

  cdm_Array* src = cdm_Array::instanceArray
                   ( DFI_Finfo.DataType
                   , DFI_Finfo.ArrayShape
                   , szS
                   , 0
                   , DFI_Finfo.Component );
  src->setHeadIndex( headS );


  //data 読込み
  //if( !read_Datarecord(fp, matchEndian, buf, headB, voxsize[2], src ) ) {
  ret = read_Datarecord(fp, matchEndian, buf, headB, voxsize[2], src );
  if( ret != CDM::E_CDM_SUCCESS) {
    ret = CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
    fclose(fp);
    printf("ERROR Data Record Read error!!!!\n");
    delete buf;
    return NULL;
  }

  //read average
  if( !avr_mode ) {
    //if( !read_averaged(fp, matchEndian, step, avr_step, avr_time) )
    ret = read_averaged(fp, matchEndian, step, avr_step, avr_time);
    if( ret !=CDM::E_CDM_SUCCESS )
    {
      ret = CDM::E_CDM_ERROR_READ_FIELD_AVERAGED_RECORD;
      delete buf;
      return src;
    }
  }

  fclose(fp);
  delete buf;

  return src;

}
