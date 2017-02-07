/*
 *  CDMlib - Cartesian Data Management library
 *
 *  Copyright (c) 2013-2016 Advanced Institute for Computational Science, RIKEN.
 *  All rights reserved.
 *
 */

/**
 * @file   cdm_DFI_FUB.C
 * @brief  cdm_DFI_FUB Class
 * @author aics
 */

#include "cdm_DFI.h"
#include "cdm_DFI_FUB.h"

// #################################################################
// コンストラクタ
cdm_DFI_FUB::cdm_DFI_FUB()
{

  //fubFlist.clear();

}

// #################################################################
// デストラクタ
cdm_DFI_FUB::~cdm_DFI_FUB()
{

}

// #################################################################
// ファイルのヘッダーレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::read_HeaderRecord( cdm_FILE* pFile,
                                bool matchEndian,
                                unsigned step,
                                const int head[3],
                                const int tail[3],
                                int gc,
                                int voxsize[3],
                                double &time)
{

  unsigned int dmy,type_dmy;

  FILE *fp = pFile->m_fp;

  //REC1
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_FUB_REC1;
  }
  if( !matchEndian ) BSWAP32(dmy);

  //REC2
  int size[4];
  if( (fread(size, sizeof(int), 4 , fp)) == 4 ) {
    if( !matchEndian ) {
      BSWAP32(size[0]);
      BSWAP32(size[1]);
      BSWAP32(size[2]);
      BSWAP32(size[3]);
    }
  }else{
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_FUB_REC1;
  }

  voxsize[0]=size[0];
  voxsize[1]=size[1];
  voxsize[2]=size[2];

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_FUB_REC1;
  }
  if( !matchEndian ) BSWAP32(dmy);

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// フィールドデータファイルのデータレコード読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::read_Datarecord(cdm_FILE* pFile,
                             bool matchEndian,
                             unsigned step,
                             cdm_Array* buf,
                             int head[3],
                             int nz,
                             cdm_Array* &src)
{

  FILE *fp = pFile->m_fp;

  // fortran record の読み込み
  unsigned int idmy;
  if( fread(&idmy,sizeof(int),1,fp) != 1 ) {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_FUB_REC2;
  }
  if( !matchEndian ) BSWAP32(idmy);

  // 1層ずつ読み込み
  int hzB = head[2];

  CDM::E_CDM_ARRAYSHAPE shape = buf->getArrayShape();

#ifdef CDM_BUFFER_MB_SIZE
  size_t ndata = src->getArrayLength();
  if( src->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
#else

  //NIJKの読込み
  if( shape == CDM::E_CDM_NIJK ) {
    for( int k=0; k<nz; k++ ) {
      //headインデクスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      //１層読み込み
      size_t ndata = buf->getArrayLength();
      if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;

      // コピー
      buf->copyArray(src);
    }
  }
  //IJKNの読込み
  else if( shape == CDM::E_CDM_IJKN ) {
    for(int n=0; n<src->getNvari(); n++) {
    for(int k=0; k<nz; k++) {
      //headインデックスをずらす
      head[2]=hzB+k;
      buf->setHeadIndex(head);

      //１層読み込み
      size_t ndata = buf->getArrayLength();
      if( buf->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;

      //コピー
      buf->copyArrayNvari(src,n);
    }}
  }
#endif

  // fortran record の読み込み
  if( fread(&idmy,sizeof(int),1,fp) != 1 ) {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_FUB_REC2;
  }
  if( !matchEndian ) BSWAP32(idmy);

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// averageデータレコードの読込み
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::read_averaged(cdm_FILE* pFile,
                           bool matchEndian,
                           unsigned step,
                           unsigned &avr_step,
                           double &avr_time)
{

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ヘッダファイルの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::write_HeaderRecord(cdm_FILE* pFile,
                                const unsigned step,
                                const double time,
                                const int RankID)
{

  FILE *fp = pFile->m_fp;

  int size[4];
  for(int i=0; i<3; i++) {
    size[i] = (int)DFI_Process.RankList[RankID].VoxelSize[i]+(int)(2*DFI_Finfo.GuideCell);
  }
  size[3] = DFI_Finfo.NumVariables;

  unsigned int dmy = sizeof(int)*4;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_FIELD_HEADER_RECORD;
  if( fwrite(size, sizeof(int), 4, fp) != 4 ) return CDM::E_CDM_ERROR_WRITE_FIELD_HEADER_RECORD;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_FIELD_HEADER_RECORD;

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
//  データレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::write_DataRecord(cdm_FILE* pFile,
                              cdm_Array* val,
                              const int gc,
                              const int RankID)
{

  FILE *fp =  pFile->m_fp;

  CDM::E_CDM_DTYPE Dtype = (CDM::E_CDM_DTYPE)DFI_Finfo.DataType;
  int Real_size = get_cdm_Datasize(Dtype);

  int size[3];
  for(int i=0; i<3; i++ ) {
   size[i] = (int)DFI_Process.RankList[RankID].VoxelSize[i]+(int)(2*gc);
  }

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.NumVariables > 1 ) dLen *= (size_t)DFI_Finfo.NumVariables;

  unsigned int dmy = dLen * Real_size;

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;
  if( val->writeBinary(fp) != dLen ) return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return CDM::E_CDM_ERROR_WRITE_FIELD_DATA_RECORD;

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// Averageレコードの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::write_averaged(cdm_FILE* pFile,
                            const unsigned step_avr,
                            const double time_avr)
{

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// fub特有のFileList要素を読み込む
/*
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::readFileListTP( cdm_TextParser tpCntl )
{

  std::string label_base,label, str;
  double ct;

  // TP
  TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    return CDM::E_CDM_ERROR_TEXTPARSER;
  }

  //FileList要素の存在チェック
  label_base = "/FileList";
  if( !tpCntl.chkNode(label_base) )
  {
    printf("\tCDM Parsing infomation : No Elem name [%s]\n", label_base.c_str());
    return CDM::E_CDM_ERROR_READ_PROCESS;
  }

  // /FileListに移動
  tp->changeNode(label_base);

  // 子のラベルを取得
  vector<std::string> labels;
  tp->getNodes(labels,1);

  // 子のRank読み込み
  for(size_t i=0; i<labels.size(); i++ )
  {
    // Rank要素か確認
    label = labels[i];
    if( strcasecmp(label.substr(0,4).c_str(), "Rank") ) continue;

    // Rank要素の読み込み
    tp->changeNode("/FileList/"+label);

    //ID
    label = "ID";
    if( !(tpCntl.GetValue(label, &ct, false )) ) {
      printf("\tCDM Parsing error : fail to get \"FileList/Rank/ID\"\n");
      return CDM::E_CDM_ERROR_READ_FILELIST_ID;
    }
    fubFname fub;
    fub.id = (int)ct;

    //CoordinateFileName
    label = "CoordinateFileName";
    if( !(tpCntl.GetValue(label, &str, false )) ) {
      printf("\tCDM Parsing error : fail to get \"FileList/Rank/CoordinateFileName\"\n");
      return CDM::E_CDM_ERROR_READ_FILELIST_CCORDINATEFILENAME;
    }
    fub.cFname = str;

    //FieldDataFileName
    label = "FieldDataFileName";
    if( !(tpCntl.GetValue(label, &str, false )) ) {
      printf("\tCDM Parsing error : fail to get \"FileList/Rank/FieldDataFileName\"\n");
      return CDM::E_CDM_ERROR_READ_FILELIST_FIELDDATAFILENAME;
    }
    fub.fFname = str;

    fubFlist.push_back(fub);

  }

  return CDM::E_CDM_SUCCESS;
}
*/
// #################################################################
// field data file nameの取得
/*
std::string
cdm_DFI_FUB::getFileNameFromFileList( const int ID )
{

  for(int i=0; fubFlist.size(); i++)
  {
    if(fubFlist[i].id == ID ) {
      if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_FUB ) return fubFlist[i].fFname;
      if( DFI_Finfo.FileFormat == CDM::E_CDM_FMT_FUB_COD ) return fubFlist[i].cFname;
    }
  }

  return "";

}
*/
// #################################################################
// coordinate file nameの取得
std::string
cdm_DFI_FUB::getCoordinateFileName(std::string FieldFileName)
{

  if( FieldFileName.empty() ) return "";

/*
  //FileListのフィールドファイル名を検索して座標値ファイル名を取得
  for(int i=0; fubFlist.size(); i++)
  {
    if( fubFlist[i].fFname == FieldFileName ) return fubFlist[i].cFname;
  }
*/

  //FileListがないとき拡張子をxyzに変更して座標値ファイル名にする
  std::string fub_ext = D_CDM_EXT_FUB;
  std::string::size_type pos = FieldFileName.find("."+fub_ext);
  
  std::string fname,xyz_fub;
  xyz_fub = D_CDM_EXT_XYZ;
  fname = FieldFileName;
  fname.replace(pos,4,"."+xyz_fub);

  return fname;

}



// #################################################################
// 座標値データの読み込み
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::ReadCoordinateData(cdm_Array *dst,
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

  CDM::E_CDM_FORMAT xfmt;
  xfmt = DFI_Finfo.FileFormat;
  int xnum;
  xnum = DFI_Finfo.NumVariables;

  DFI_Finfo.FileFormat = CDM::E_CDM_FMT_FUB_COD;
  DFI_Finfo.NumVariables = 3;
 
  CDM::E_CDM_ERRORCODE ret;
  ret = ReadData(dst, step, gc, Gvoxel, Gdivision, head, tail, time,
                 avr_mode, avr_step, avr_time);

  DFI_Finfo.FileFormat = xfmt;
  DFI_Finfo.NumVariables = xnum; 

  return ret;
}

// #################################################################
// 座標値データの出力
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::WriteCoordinateData(const unsigned step,
                                 const int gc,
                                 double time,
                                 cdm_Array* val)
{

  CDM::E_CDM_FORMAT xfmt;
  xfmt = DFI_Finfo.FileFormat;
  DFI_Finfo.FileFormat = CDM::E_CDM_FMT_FUB_COD;
  //std::string tname = DFI_Finfo.Prefix;
  //DFI_Finfo.Prefix = "nod";
  std::string tdfiname = m_indexDfiName;
  m_indexDfiName="";
  int xnum = DFI_Finfo.NumVariables;
  DFI_Finfo.NumVariables = 3;

  CDM::E_CDM_ERRORCODE ret;
  ret = WriteData(step,
                  gc,
                  time,
                  val,
                  (double*)NULL,
                  true,
                  0,
                  0.0);

  DFI_Finfo.FileFormat = xfmt;
  //DFI_Finfo.Prefix = tname;
  m_indexDfiName = tdfiname;
  DFI_Finfo.NumVariables = xnum;

  return ret;

}

// #################################################################
// cdm_FieldFileNameFormatクラスのポインタ取得
const cdm_FieldFileNameFormat*
cdm_DFI_FUB::GetcdmFieldFileNameFormat()
{

  if( DFI_FieldFileNameFormat.LabelList.size() <= 0 ) return NULL;

  return &DFI_FieldFileNameFormat;

}

// #################################################################
// FileListの出力
/*
CDM::E_CDM_ERRORCODE
cdm_DFI_FUB::WriteFileList(FILE* fp,
                           const unsigned tab)
{

  if( fubFlist.size() <= 0 ) CDM::E_CDM_SUCCESS;

  fprintf(fp, "FileList {\n");
  fprintf(fp, "\n");

  for(int i=0; i<fubFlist.size(); i++) {
    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "Rank[@] {\n");
    fprintf(fp, "\n");
    
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "ID        = %d\n", fubFlist[i].id);

    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "CoordinateFileName = \"%s\"\n",fubFlist[i].cFname.c_str());

    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "FieldDataFileName  = \"%s\"\n",fubFlist[i].fFname.c_str());

    fprintf(fp, "\n");
    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "}\n");
  }
  fprintf(fp, "}\n");

  return CDM::E_CDM_SUCCESS;
}
*/
