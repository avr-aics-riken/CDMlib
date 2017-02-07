/*
 * frm (File Rank Mapper)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   Staging.C
 * @brief  Staging class 関数
 * @author aics
 */

#include <stdlib.h>
#include "Staging_Utility.h"

//20160407.fub.s
//////////////////////////////////////////////////////////////////
//引数のエラーチェック
stg_ErrorCode
Staging::ArgErrorCheck()
{

  if( m_step < 0 && m_CoordinateStep >= 0 ) {
    printf("ERROR input option -c : %d\n",m_CoordinateStep);
    return STG_ERROR_ILLEGAL_OPTION_C;
  }

  return STG_SUCCESS;
}

//////////////////////////////////////////////////////////////////
// fubフィルコピー（フィールドデータファイルと座標値データファイルコピー
void Staging::FubFileCopy(int readRankID, int step, const bool mio,
                          std::string path2, cdm_DFI *dfi, char *cmd)
{

  cdm_DFI_FUB *dfi_fub = dynamic_cast<cdm_DFI_FUB*>(dfi);
  
/*
  //FileListからのファイル名取得
  std::string fname = dfi_fub->getFileNameFromFileList(readRankID);
*/
  cdm_FieldFileNameFormat* 
  Ffformat = (cdm_FieldFileNameFormat *)dfi->GetcdmFieldFileNameFormat();

  std::string fname = "";
  if( Ffformat ) {
     fname = Ffformat->GenerateFileName("FieldFile","",step,readRankID);
  }

  //FileListから取得したファイル名がディレクトリ付のとき、ディレクトリ名取得
  std::string fname_Path="";
  if( !fname.empty() ) {
    fname_Path = CDM::cdmPath_DirName(fname);
    if( !fname_Path.empty() ) {
      std::string::size_type pos = fname_Path.find("/");
      if( pos < 2 ) fname_Path.erase(0,pos+1); 
      pos = fname_Path.find("/");
      if( (int)pos < 0 || pos < fname_Path.size()-1 ) fname_Path += "/";
      MakeDirectory(path2+fname_Path);
    }
  }

  std::string Ffname;
  //FileListがないときはファイル名を生成
  if( fname.empty() ) {
    Ffname = CDM::cdmPath_ConnectPath(m_inPath,Generate_FileName(readRankID,step,mio));
  } else {
    if( dfi_Finfo->DirectoryPath.empty() ) {
      Ffname = fname;
    } else {
      Ffname = dfi_Finfo->DirectoryPath + "/" + fname;
    }
  }
   
  //フィールドデータファイルのコピーコマンド生成と実行
  memset(cmd, 0, sizeof(char)*512 );
  sprintf(cmd,"cp %s %s%s\n",Ffname.c_str(),path2.c_str(),fname_Path.c_str());
  system(cmd); 

  //Coordenate data file名の取得
  std::string Cfname, OutCfname;
  OutCfname="";
  //if( m_step < 0 || dfi_Finfo->FieldFilenameFormat == CDM::E_CDM_FNAME_CUSTOM ) {
  if( m_step < 0 ) {
    if( fname.empty() ) {
      Cfname = dfi_fub->getCoordinateFileName(Ffname);
    } else {
      Cfname = Ffformat->GenerateFileName("CoordinateFile",dfi_Finfo->DirectoryPath,
                                          step,readRankID);
    }
  } else if( m_CoordinateStep >= 0 && fname.empty() ) {
    Ffname = CDM::cdmPath_ConnectPath(m_inPath,Generate_FileName(readRankID,m_CoordinateStep,mio));
    Cfname = dfi_fub->getCoordinateFileName(Ffname);
    Ffname = CDM::cdmPath_FileName(Generate_FileName(readRankID,step,mio));
    OutCfname = dfi_fub->getCoordinateFileName(Ffname);
  } else {
    if( fname.empty() ) {
      Cfname = CDM::cdmPath_ConnectPath(m_inPath,Generate_FileName(readRankID,step,mio));
    } else {
      if( m_CoordinateStep < 0 ) {
        Cfname = Ffformat->GenerateFileName("CoordinateFile",dfi_Finfo->DirectoryPath,
                                            step,readRankID);
      } else {
        Cfname = Ffformat->GenerateFileName("CoordinateFile",dfi_Finfo->DirectoryPath,
                                            m_CoordinateStep,readRankID);
        OutCfname = Ffformat->GenerateFileName("CoordinateFile","",
                                            step,readRankID);
      }
    }
  }

  //座標値データファイルのコピーコマンド生成と実行
  FILE *tmp_fp;
  if( (tmp_fp=fopen(Cfname.c_str(), "rb")) ) {
    fclose(tmp_fp);
    memset(cmd, 0, sizeof(char)*512 );
    sprintf(cmd,"cp %s %s%s%s\n",Cfname.c_str(),path2.c_str(), fname_Path.c_str()
           ,OutCfname.c_str());
    system(cmd);
  } else if( m_step>=0 ) {
    printf("Error undefined Coordinate data file : %s\n",Cfname.c_str());
  }

  return;

}  
//20160407.fub.e

///////////////////////////////////////////////////////////////////
// 初期化、ファイルの読込み
//bool Staging::Initial( string infofile, string dfiname )
bool Staging::Initial( string infofile )
{

  //string fconvfile; //FCONV 入力ファイル名

  //infoファイル読込み
  //if( !ReadInfo(infofile, fconvfile) ) return false;
  if( !ReadInfo() ) return false;

  //FCONV入力ファイルの読込み
  //if( !ReadFconvInputFile(fconvfile) ) return false;
  if( !ReadFconvInputFile() ) return false;

  if( !FconvInputCheck() ) return false;

  DFI.clear();

  //dfiファイル読込み
  CDM::E_CDM_ERRORCODE ret;
  for(int i=0; i<m_dfi_fname.size(); i++) {
    cdm_DFI* dfi_in = cdm_DFI::ReadInit(MPI_COMM_WORLD, m_dfi_fname[i], m_GVoxel, m_Gdiv, ret);

    if( dfi_in == NULL ) return ret;
    if( ret != CDM::E_CDM_SUCCESS && ret != CDM::E_CDM_ERROR_INVALID_DIVNUM ) return ret;

    DFI.push_back(dfi_in);
  }

  if( m_GVoxel[0] == 0 || m_GVoxel[1] == 0 || m_GVoxel[2] == 0 ) {
    dfi_Domain = DFI[0]->GetcdmDomain();
    for(int i=0; i<3; i++) m_GVoxel[i]=dfi_Domain->GlobalVoxel[i]; 
  }
  if( m_Gdiv[0] == 0 || m_Gdiv[1] == 0 || m_Gdiv[2] == 0 ) {
    dfi_Domain = DFI[0]->GetcdmDomain();
    for(int i=0; i<3; i++) m_Gdiv[i]=dfi_Domain->GlobalDivision[i];
  }

  //dfi情報が格納されたクラスポインタの取得
  /*
  dfi_Finfo = DFI->GetcdmFileInfo();
  dfi_Fpath = DFI->GetcdmFilePath();
  dfi_Visit = DFI->GetcdmVisIt();
  dfi_Unit  = DFI->GetcdmUnit();
  dfi_Domain= DFI->GetcdmDomain();
  dfi_MPI   = DFI->GetcdmMPI();
  dfi_TSlice= DFI->GetcdmTimeSlice();
  dfi_Process=(cdm_Process *)DFI->GetcdmProcess();
  */

  if( !m_ActiveSubdomain.empty() ) {
    int divSudomain[3] = {0,0,0};
    //ActiveSubdomainファイルの読込み
    stg_ErrorCode ret = ReadActiveSubdomainFile( m_ActiveSubdomain, m_subDomainInfo, divSudomain);
    if( ret != STG_SUCCESS ) return false;
  } else {
    int nRank = m_GRankInfo.size();
    if( nRank == 0 ) nRank = m_Gdiv[0]*m_Gdiv[1]*m_Gdiv[2];
    //ActiveSubdomainファイルがない場合、全ランク活性ファイルとする
    if( CheckData(nRank) != STG_SUCCESS ) return false; 
  }

  return true;

}

///////////////////////////////////////////////////////////////////
// infoファイルの読込み
//bool Staging::ReadInfo(string infofile, string &fconvfile)
bool Staging::ReadInfo()
{

  if( m_infofile.empty() ) return true;

  //実行プログラム情報ローダのインスタンス生成
  cdm_TextParser tp_stg;
  tp_stg.getTPinstance();

  FILE* fp = NULL;
  if( !(fp=fopen(m_infofile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",m_infofile.c_str());
    return false;
  }
  fclose(fp);

  int ierr = tp_stg.readTPfile(m_infofile);
  if( ierr ) {
    printf("\tinput file not found '%s'\n",m_infofile.c_str());
    return false;
  }

  string label,label_base,label_leaf;
  double v[3];
  string str;
  int ct;

  //GlobalVoxel
  for (int n=0; n<3; n++) v[n]=0.0;
  label = "/Domain/GlobalVoxel";
  if( !(tp_stg.GetVector(label, v, 3 )) ) {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return false;
  }
  for(int i=0;i<3;i++ ) m_GVoxel[i]=v[i];

  //GlobalDivision
  for (int n=0; n<3; n++) v[n]=0.0;
  label = "/Domain/GlobalDivision";
  if( !(tp_stg.GetVector(label, v, 3 )) ) {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return false;
  }
  for(int i=0;i<3;i++ ) m_Gdiv[i]=v[i];

  //ActiveSubdomain
  label = "/Domain/ActiveSubdomainFile";
  if( !(tp_stg.GetValue(label, &str )) ) {
    m_ActiveSubdomain="";
  } else {
    m_ActiveSubdomain=str;
  }

//FCONV 20140116.s
/*
  //InputFile
  label = "/FCONVInfo/InputFile";
  if( !(tp_stg.GetValue(label, &str )) ) {
    fconvfile="";
  } else {
    fconvfile=str;
  }

  //NumberOfProcess
  label = "/FCONVInfo/NumberOfProcess";
  if( !(tp_stg.GetValue(label, &ct )) ) {
    m_fconv_numproc=0;
  } else {
    m_fconv_numproc=ct;
  }
*/
  m_ConvType = STG_E_OUTPUT_MxN;

  for(int i=0; i<3; i++) {
    m_CropStart[i]=1;
    m_CropEnd[i]=m_GVoxel[i];
  }

//FCONV 20140116.e

  //NumberOfRank
  label = "/MPI/NumberOfRank";
  if( !(tp_stg.GetValue(label, &ct )) ) {
    m_NumberOfRank=0;
  } else {
    m_NumberOfRank=ct;
  }

  m_GRankInfo.clear();

  Rank rank;
  //Rank
  int nnode=0;
  label_base = "/Process";
  if( tp_stg.chkNode(label_base) ) {
    nnode = tp_stg.countLabels(label_base);
    if( m_NumberOfRank != nnode ) return false;
  } else {
    if( m_NumberOfRank > 0 ) return false;
    return true;
  }

  TextParser *tp = tp_stg.getTPPtr();
  if( !tp ) 
  {
    return false;
  }

  // /Processに移動
  tp->changeNode(label_base);

  // 子のラベルを取得
  vector<std::string> labels;
  tp->getNodes(labels,1);

/*
  for (int i=0; i<nnode; i++) {
    if(!tp_stg.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      return false;
    }
    if( strcasecmp(str.substr(0,4).c_str(), "Rank") ) continue;
    label_leaf=label_base+"/"+str;
*/

  //子のRankを読み込み
  for( size_t i=0; i<labels.size(); i++) {
    //Rank要素かどうか確認
    label=labels[i];
    if( strcasecmp(label.substr(0,4).c_str(), "Rank") ) continue;

    //Rankに移動
    label_leaf=label_base + "/" + label;
    tp->changeNode(label_leaf);

    //ID
    label = "ID";
    if ( !(tp_stg.GetValue(label, &ct, false )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return false;
    }
    else {
      rank.RankID= ct;
    }

    //VoxelSize
    label = "VoxelSize";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tp_stg.GetVector(label, v, 3, false )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return false;
    }
    for (int n=0; n<3; n++) rank.VoxelSize[n]=v[n];

    //HeadIndex&TailIndexの初期化
    for(int n=0; n<3; n++) {
      rank.HeadIndex[n]=0;
      rank.TailIndex[n]=0;
    }

    //HeadIndex
    label = "HeadIndex";
    if( tp_stg.GetVector(label, v, 3, false ) ) {
      for (int n=0; n<3; n++) rank.HeadIndex[n]=v[n];
      //TailIndex
      label = "TailIndex";
      if( tp_stg.GetVector(label, v, 3, false ) ) {
        for (int n=0; n<3; n++) rank.TailIndex[n]=v[n];
      } else {
        //TailIndexがないとき自動セット
        for (int n=0; n<3; n++) {
          rank.TailIndex[n]=rank.HeadIndex[n]+rank.VoxelSize[n]-1;
        }
      }
    //HeadIndexがなしでTailIndexがあるとき
    } else {
      label = "TailIndex";
      if( tp_stg.GetVector(label, v, 3, false ) ) {
        for (int n=0; n<3; n++) {
          rank.TailIndex[n]=v[n];
          rank.HeadIndex[n]=rank.TailIndex[n]-rank.VoxelSize[n]+1;
        }
      }
    }

    //RankPosition
    label = "RankPosition";
    for (int n=0; n<3; n++) v[n]=-1.0;
    tp_stg.GetVector(label, v, 3, false ); 
    for (int n=0; n<3; n++) rank.RankPosition[n]=v[n]; 

    m_GRankInfo.push_back(rank);

  }

  if( m_GRankInfo.size() != m_NumberOfRank ) return false;

  tp_stg.remove();

  return true;
}

///////////////////////////////////////////////////////////////////
// FCONV 入力ファイルの読込み
//bool Staging::ReadFconvInputFile(string fconvfile)
bool Staging::ReadFconvInputFile()
{

  if( m_fconvfile.empty() ) {
    m_fconv_inputfile = false;
    return true;
  }

  string label;
  string str;
  int v[3];
  int ct;

  //実行プログラム情報ローダのインスタンス生成
  //cdm_TextParser tp_stg;
  TextParser tp_stg;
  //tp_stg.getTPinstance();

  /*
  FILE* fp = NULL;
  if( !(fp=fopen(m_fconvfile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",m_fconvfile.c_str());
    return false;
  }
  fclose(fp);
  */

  //int ierr = tp_stg.readTPfile(m_fconvfile);
  int ierr = tp_stg.read(m_fconvfile);

  int nnode=0;
  string label_base = "/ConvData/InputDFI";
  if( tp_stg.chkNode(label_base) ) {
    nnode = tp_stg.countLabels(label_base);
  } else return false;

  int ncnt=0;
  label_base = "/ConvData";
  ncnt++;
  //dfiのファイル名の読込み
  for (int i=0; i<nnode; i++) {
    if(!tp_stg.getNodeStr(label_base, ncnt, str)) {
      printf("\tParsing error : No Elem name\n");
      return false;
    }
    if( !strcasecmp(str.substr(0,8).c_str(), "InputDFI") ) {
      label = label_base+"/"+str;
      if( !(tp_stg.getInspectedValue(label, str )) ) {
        printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      m_dfi_fname.push_back(str);
      ncnt++;
    }
  }

  //コンバートタイプの読込み
  label = "/ConvData/ConvType";
  //if( (tp_stg.GetValue(label,&str)) ) {
  if( (tp_stg.getInspectedValue(label,str)) ) {
    if     ( !strcasecmp(str.c_str(), "Mx1") )  m_ConvType = STG_E_OUTPUT_Mx1;
    else if( !strcasecmp(str.c_str(), "MxN") )  m_ConvType = STG_E_OUTPUT_MxN;
    else if( !strcasecmp(str.c_str(), "MxM") )  m_ConvType = STG_E_OUTPUT_MxM;
    else
    {
      printf("\tInvalid keyword is described for '%s'\n", label.c_str());
      return false;
    }
  }

  //MxN用出力分割情報の取得
  label = "/ConvData/OutputDivision";
  if( (tp_stg.getInspectedVector(label, v, 3 )) ) {
    for(int i=0; i<3; i++) m_Gdiv[i]=v[i];
  }

  //出力ガイドセル数
  label = "/ConvData/OutputGuideCell";
  if( !(tp_stg.getInspectedValue(label,ct)) ) {
    m_outGc=0;
  } else {
    m_outGc=ct;
  }

  //ファイル割振り方法の取得
  label = "/ConvData/MultiFileCasting";
  //if( (tp_stg.GetValue(label,&str)) ) {
  if( (tp_stg.getInspectedValue(label,str)) ) {
    if     ( !strcasecmp(str.c_str(), "step") ) m_outList = STG_E_OUTPUT_TYPE_STEP;
    else if( !strcasecmp(str.c_str(), "rank") ) m_outList = STG_E_OUTPUT_TYPE_RANK;
  } else m_outList = STG_E_OUTPUT_TYPE_STEP;


  //入力指示の読込み
  label = "/ConvData/CropIndexStart";
  if( (tp_stg.getInspectedVector(label, v, 3 )) ) {
    for(int i=0; i<3; i++) m_CropStart[i]=v[i];
    m_cropIndexStart_on=true; 
  }
  label = "/ConvData/CropIndexEnd";
  if( (tp_stg.getInspectedVector(label, v, 3 )) ) {
    for(int i=0; i<3; i++) m_CropEnd[i]=v[i];
    m_cropIndexEnd_on=true;
  }

  tp_stg.remove();

  m_fconv_inputfile = true;

  return true;

}

///////////////////////////////////////////////////////////////////
// FCONVファイルの入力チェック
//
bool Staging::FconvInputCheck()
{

  if( m_ConvType == STG_E_OUTPUT_Mx1 && m_outList == STG_E_OUTPUT_TYPE_RANK ) {
    printf("\tInvalid MultiFileCasting \"rank\"\n");
    return false;
  }
  
  if( m_ConvType == STG_E_OUTPUT_MxM ) {
    if( m_cropIndexStart_on || m_cropIndexEnd_on ) {
      printf("\tInvalid CropIndexStart or CropIndexEnd\n");
      return false;
    }
  }

  if( m_ConvType != STG_E_OUTPUT_MxN ) {
    for(int i=0; i<3; i++) m_Gdiv[i]=0;
  }

  return true; 
}

///////////////////////////////////////////////////////////////////
// step基準のリスト生成
void Staging::makeStepList(int myID)
{

  m_StepRankList.clear();

  //総ステップ数を求める
  int Total_step = 0;
  for(int i=0; i<DFI.size(); i++) {
    const cdm_TimeSlice* TSlice = DFI[i]->GetcdmTimeSlice();
    for(int j=0; j<TSlice->SliceList.size(); j++) {
      if( m_step > 0 && TSlice->SliceList[j].step != m_step ) continue;
      Total_step++;
    }
  }

  //各ランクで担当するステップ数を求める
  int nStep = Total_step/m_fconv_numproc;
  if( Total_step%m_fconv_numproc != 0 ) {
    for(int i=0; i<Total_step%m_fconv_numproc; i++) {
      if( myID==i ) nStep++;
    }
  }

  //各ランクが担当するステップのスタートとエンドを求める
  int sta,end;
  sta = myID * nStep;
  if( Total_step%m_fconv_numproc != 0 ) {
    if( myID >= Total_step%m_fconv_numproc ) sta = sta+Total_step%m_fconv_numproc;
  }
  end = sta+nStep-1;
  if( nStep == 0 ) return;

  //処理ステップリストの生成
  int cnt=0;
  for(int i=0; i<DFI.size(); i++) {
    step_rank_info info;
    info.stepStart = -1;
    const cdm_TimeSlice* TSlice = DFI[i]->GetcdmTimeSlice();
    for(int j=0; j<TSlice->SliceList.size(); j++) {
//20160510.fub.s
      if( m_step >= 0 && TSlice->SliceList[j].step != m_step ) continue;
//20160510.fub.e
      if( sta > cnt ) { cnt++; continue; }
      if( info.stepStart == -1 ) {
        info.dfi = DFI[i];
        info.stepStart = j;
      }
      info.stepEnd = j;
      cnt++;
      if( end < cnt ) break;
    }
    if( info.stepStart > -1 ) m_StepRankList.push_back(info);
    if( end < cnt ) break;
  }    

  //rantStart,rankEndのセット
  for(int i=0; i<m_StepRankList.size(); i++) {
    const cdm_Process* DFI_Process = m_StepRankList[i].dfi->GetcdmProcess();
    m_StepRankList[i].rankStart=0;
    m_StepRankList[i].rankEnd=DFI_Process->RankList.size()-1;
  }

}

///////////////////////////////////////////////////////////////////
// rank基準のリスト生成
void Staging::makeRankList(int myID)
{

  m_StepRankList.clear();

  //総ランク数を求める
  int Total_rank = 0;
  for(int i=0; i<DFI.size(); i++) {
    const cdm_Process* DFI_Process = DFI[i]->GetcdmProcess();
    Total_rank += DFI_Process->RankList.size();
  }

  //各ランクで担当するrank数を求める
  int nRank = Total_rank/m_fconv_numproc;
  if( Total_rank%m_fconv_numproc != 0 ) {
    for(int i=0; i<Total_rank%m_fconv_numproc; i++) {
      if( myID == i ) nRank++;
    }
  }

  //自ランクが担当するランクのスタートとエンドを求める
  int sta,end;
  sta = myID * nRank;
  if( Total_rank%m_fconv_numproc != 0 ) {
    if( myID >= Total_rank%m_fconv_numproc ) sta = sta+Total_rank%m_fconv_numproc;
  }
  end = sta+nRank-1;
  if( nRank == 0 ) return;

  //処理rankリストの生成
  int cnt=0;
  for(int i=0; i<DFI.size(); i++) {
    step_rank_info info;
    info.rankStart = -1;
    const cdm_Process* DFI_Process = DFI[i]->GetcdmProcess();
    for(int j=0; j<DFI_Process->RankList.size(); j++) {
      if( sta > cnt ) { cnt++; continue; }
      if( info.rankStart == -1 ) {
        info.dfi = DFI[i];
        info.rankStart = j;
      }
      info.rankEnd = j;
      cnt++;
      if( end < cnt ) break;
    } 
    if( info.rankStart > -1 ) m_StepRankList.push_back(info);
    if( end < cnt ) break;
  }

  //tepStart,stepEndのセット
  for(int i=0; i<m_StepRankList.size(); i++) {
    const cdm_TimeSlice* TSlice = m_StepRankList[i].dfi->GetcdmTimeSlice();
    m_StepRankList[i].stepStart=0;
    m_StepRankList[i].stepEnd=TSlice->SliceList.size()-1;
  }
  //printf("myID : %d nRank : %d sta : %d end : %d m_StepRankList.size : %d\n",
  //        myID,nRank,sta,end,(int)m_StepRankList.size());

}

///////////////////////////////////////////////////////////////////
// 領域分割数の取得
const int* Staging::GetDivNum() const
{
  return m_Gdiv;
}

///////////////////////////////////////////////////////////////////
// 活性サブドメイン情報の存在チェック
bool Staging::IsExistSubdomain( ActiveSubDomain subDomain)
{
  for( size_t i=0;i<m_subDomainInfo.size();i++ )
  {
    ActiveSubDomain dom = m_subDomainInfo[i];
    if( dom == subDomain ) return true;
  }
  return false;
}

///////////////////////////////////////////////////////////////////
// 活性サブドメイン情報の追加
bool Staging::AddSubdomain(ActiveSubDomain subDomain )
{

  //既存チェック
  if( IsExistSubdomain(subDomain) ) return false;

  //追加
  m_subDomainInfo.push_back(subDomain);
  return true;
}

///////////////////////////////////////////////////////////////////
// 活性サブドメインの数を取得
int Staging::GetSubdomainNum() const
{
  if( m_subDomainInfo.size() > 0 )
  {
    return (int)m_subDomainInfo.size();
  }
  if( m_Gdiv[0] <= 0 || m_Gdiv[1] <= 0 || m_Gdiv[2] <= 0 )
  {
    return 0;
  }
  return m_Gdiv[0] * m_Gdiv[1] * m_Gdiv[2];
}

///////////////////////////////////////////////////////////////////
// 活性サブドメインの数を取得
const ActiveSubDomain* Staging::GetSubdomainInfo( size_t idx ) const
{
  if( int(idx) >= GetSubdomainNum() ) return NULL;
  return &(m_subDomainInfo[idx]);
}

///////////////////////////////////////////////////////////////////
// ランクマップを生成
//bool Staging::CreateRankMap()
int* Staging::CreateRankMap()
{

  //m_rankMap = NULL;

  // 領域分割数を取得
  const int* div = GetDivNum();
  if( !div ) return NULL;

  // マップ領域を確保(初期値NULL)
  size_t ndiv = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  int *rankMap = new int[ndiv];
  if( !rankMap )
  {
    return NULL;
  }
  for( size_t i=0;i<ndiv;i++ ) rankMap[i] = -1;

  // 活性サブドメイン情報配置位置に0をセット
  for( int i=0;i<GetSubdomainNum();i++ )
  {
    //サブドメイン情報
    const ActiveSubDomain* dom = GetSubdomainInfo(i);
    if( !dom )
    {
      delete [] rankMap;
      return NULL;
    }

    // 位置を取得
    const int *pos = dom->GetPos();
    if( !pos )
    {
      delete [] rankMap;
      return NULL;
    }

    // 0をセット
    rankMap[_IDX_S3D(pos[0],pos[1],pos[2],div[0],div[1],div[2],0)] = 0;
  }

  bool flg = true;
  for(int i=0; i<m_GRankInfo.size(); i++)
  {
    if( m_GRankInfo[i].RankPosition[0]<0 || m_GRankInfo[i].RankPosition[1]<0 || 
        m_GRankInfo[i].RankPosition[2]<0 ) {
      flg = false;
      break;
    }
  }

  if( !flg || m_GRankInfo.size()<1 ) {
    // i->j->kの優先順で活性サブドメインにランク番号をセット
    int rankCount = 0;
    for( int k=0;k<div[2];k++ ){
    for( int j=0;j<div[1];j++ ){
    for( int i=0;i<div[0];i++ ){
      if( rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] == 0 )
      {
        rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] = rankCount;
        rankCount++;
      }
    }}}
  } else {
    //printf("RankPosition set\n");
    // i->j->kの優先順で活性サブドメインにランク番号をセット
    for( int n=0; n<m_GRankInfo.size(); n++ ) {
      int i=m_GRankInfo[n].RankPosition[0];
      int j=m_GRankInfo[n].RankPosition[1];
      int k=m_GRankInfo[n].RankPosition[2];
      if( rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] == 0 )
      {
        rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] = n;
      }
    }
  }

  // ランクマップをセット
  //if( m_rankMap ) delete [] m_rankMap;
  //m_rankMap = rankMap;

  //for(int i=0;i<ndiv;i++) printf("rankMap[%d] : %d\n",i,rankMap[i]);

  return rankMap;
}

///////////////////////////////////////////////////////////////////
// ランクマップを生成
int* Staging::CreateActiveRankMap()
{

  int i,j,k;

  int div[3];
  div[0]=m_mapX.size();
  div[1]=m_mapY.size();
  div[2]=m_mapZ.size();

  size_t ndiv = div[0]*div[1]*div[2];

  int *rankMap = new int[ndiv];
  for(int i=0; i<ndiv; i++) rankMap[i]=-1;

  headT::iterator it;

  for(int n=0; n<m_GRankInfo.size(); n++)
  {
    it=m_mapX.find(m_GRankInfo[n].HeadIndex[0]);
    i=it->second;

    it=m_mapY.find(m_GRankInfo[n].HeadIndex[1]);
    j=it->second;

    it=m_mapZ.find(m_GRankInfo[n].HeadIndex[2]);
    k=it->second;

    rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)] = n;
  }

  return rankMap; 

}

///////////////////////////////////////////////////////////////////
// ActiveSubdomainファイルのエンディアンチェック 
stg_EMatchType Staging::isMatchEndianSbdmMagick( int ident )
{
  char magick_c[] = "SBDM";
  int  magick_i=0;

  //cheak match
  magick_i = (magick_c[3]<<24) + (magick_c[2]<<16) + (magick_c[1]<<8) + magick_c[0];
  if( magick_i == ident )
  {
    return STG_Match;
  }

  //chack unmatch
  magick_i = (magick_c[0]<<24) + (magick_c[1]<<16) + (magick_c[2]<<8) + magick_c[3];
  if( magick_i == ident )
  {
    return STG_UnMatch;
  }

  //unknown format
  return STG_UnKnown;
}


///////////////////////////////////////////////////////////////////
// ActiveSubdomainファイルの読み込み
stg_ErrorCode Staging::ReadActiveSubdomainFile( std::string subDomainFile )
{

  //読込み
  int div[3];
  stg_ErrorCode ret = ReadActiveSubdomainFile( subDomainFile, m_subDomainInfo, div );
  if( ret != STG_SUCCESS ) return ret;

  // 領域分割数のチェック
  if( m_Gdiv[0] > 0 && m_Gdiv[1] > 0 && m_Gdiv[2] > 0 )
  {
    if( div[0] != m_Gdiv[0] || div[0] != m_Gdiv[0] ||div[0] != m_Gdiv[0] )
    {
      return STG_ERROR_MISMATCH_DIV_SUBDOMAIN;
    }
  } 

  return STG_SUCCESS;
}


///////////////////////////////////////////////////////////////////
// ActiveSubdomainファイルの読み込み(static関数)
stg_ErrorCode Staging::ReadActiveSubdomainFile( std::string subDomainFile,
                                                std::vector<ActiveSubDomain>& subDomainInfo,
                                                int div[3] )
{
  if( subDomainFile.empty() ) return STG_ERROR_OPEN_SBDM;

  // ファイルオープン
  FILE*fp = fopen( subDomainFile.c_str(), "rb" );
  if( !fp ) return STG_ERROR_OPEN_SBDM; 

  //エンディアン識別子
  int ident;
    if( fread( &ident, sizeof(int), 1, fp ) != 1 )
  {
    fclose(fp);
    return STG_ERROR_READ_SBDM_HEADER;
  }

  stg_EMatchType endian = isMatchEndianSbdmMagick( ident );
  if( endian == STG_UnKnown )
  {
    fclose(fp);
    return STG_ERROR_READ_SBDM_FORMAT;
  }

  // 領域分割数
  if( fread( div, sizeof(int), 3, fp ) != 3 )
  {
    fclose(fp);
    return STG_ERROR_READ_SBDM_DIV;
  }
  if( endian == STG_UnMatch ) BSWAPVEC(div,3);

  // contents
  size_t nc = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  unsigned char *contents = new unsigned char[nc];
  if( fread( contents, sizeof(unsigned char), nc, fp ) != nc )
  {
    delete [] contents;
    fclose(fp);
    return STG_ERROR_READ_SBDM_CONTENTS;
  }

  // ファイルクローズ
  fclose(fp);

  size_t ptr = 0;
  // 活性ドメイン情報の生成
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    if( contents[ptr] == 0x01 )
    {
      int pos[3] = {i,j,k};
      ActiveSubDomain dom( pos );
      subDomainInfo.push_back(dom);
    }
    ptr++;
  }}}

  // contentsのdelete
  delete [] contents;

  // 活性ドメインの数をチェック
  if( subDomainInfo.size() == 0 )
  {
    return STG_ERROR_SBDM_NUMDOMAIN_ZERO;
  }
 
  return STG_SUCCESS;
}
 
///////////////////////////////////////////////////////////////////
// 領域情報のチェック
stg_ErrorCode Staging::CheckData( int nRank )
{

  // 領域分割数
  if( m_Gdiv[0] <= 0 || m_Gdiv[1] <= 0 || m_Gdiv[2] <= 0 ) 
    return STG_ERROR_INVALID_DIVNUM;

  //活性サブドメイン情報
  int ndom = m_subDomainInfo.size();

  bool flg=true;
  for( int i=0; i<m_GRankInfo.size(); i++ ) {
    if( m_GRankInfo[i].RankPosition[0]<0 || m_GRankInfo[i].RankPosition[1]<0 ||
        m_GRankInfo[i].RankPosition[2]<0 ) {
      flg=false;
      break;
    }
  }

  if( ndom==0 && m_GRankInfo.size()==0 ) flg=false;

  if( ndom == 0 ) {
    //活性サブドメイン情報が空のとき、全領域を活性サブドメインとする
    //if( nRank != m_Gdiv[0]*m_Gdiv[1]*m_Gdiv[2] )
    //{
    //  return STG_ERROR_MISMATCH_NP_SUBDOMAIN;
    //}
    if( flg ) {
      for( int n=0; n<m_GRankInfo.size(); n++ ) {
        int i=m_GRankInfo[n].RankPosition[0];
        int j=m_GRankInfo[n].RankPosition[1];
        int k=m_GRankInfo[n].RankPosition[2];
        int pos[3] = {i,j,k};
        ActiveSubDomain dom( pos );
        AddSubdomain( dom );
      }
    } else {
      for( int k=0;k<m_Gdiv[2];k++ ){
      for( int j=0;j<m_Gdiv[1];j++ ){
      for( int i=0;i<m_Gdiv[0];i++ ){
        int pos[3] = {i,j,k};
        ActiveSubDomain dom( pos );
        AddSubdomain( dom );
      }}}
    }
  } else {
    //if( nRank != ndom )
    //{
    //  return STG_ERROR_MISMATCH_NP_SUBDOMAIN;
    //}
  }

  return STG_SUCCESS;
}

///////////////////////////////////////////////////////////////////
// VOXEL数の取得
const int* Staging::GetVoxNum() const
{
  return m_GVoxel;
}

///////////////////////////////////////////////////////////////////
// head&tail の生成
bool Staging::CreateHeadTail(int* rankMap, vector<Rank>& RankInfo)
{

  if( !rankMap ) return false;
 
  //領域分割数
  const int *div = GetDivNum();
  if( !div ) return false;
 
  int numrank=0;
  for(int i=0; i<div[0]*div[1]*div[2]; i++ ) {
    if( rankMap[i] >= 0 ) numrank++;
  }

  if( numrank != RankInfo.size() ) {
    printf("ERROR MISMATCH NumberOfRank\n");
    return false;
  }
 
  bool flg=true;
  for( int i=0; i<RankInfo.size(); i++ ) {
    if( RankInfo[i].HeadIndex[0]<1 || RankInfo[i].HeadIndex[1]<1 || RankInfo[i].HeadIndex[2]<1 )
    {
      flg=false;
      break;
    }
  }

  if( flg ) return true;


  //全体ボクセル数
  const int *gvox = GetVoxNum();
  if( !gvox ) return false;

  // ローカルのVOXEL数
  int *nvX = new int[div[0]];
  int *nvY = new int[div[1]];
  int *nvZ = new int[div[2]];
  int *nv[3] = {nvX,nvY,nvZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    //基準のボクセル数
    int nbase = gvox[n] / div[n];

    //余り
    int amari = gvox[n] % div[n];

    //ボクセル数をセット
    for( int i=0;i<div[n];i++ )
    {
      nvd[i] = nbase;
      if( i<amari ) nvd[i]++;
    }
  }

  //head
  int *headX = new int[div[0]];
  int *headY = new int[div[1]];
  int *headZ = new int[div[2]];
  int *head[3] = {headX,headY,headZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    int *hd = head[n];
    hd[0] = 1;

    for( int i=1;i<div[n];i++ )
    {
      hd[i] = hd[i-1]+nvd[i-1];
    }
  }

  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    int rankNo = rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)];

    if( rankNo < 0 ) continue;

    m_HeadTail[rankNo][0] = headX[i];
    m_HeadTail[rankNo][1] = headY[j];
    m_HeadTail[rankNo][2] = headZ[k];
    m_HeadTail[rankNo][3] = headX[i]+nvX[i]-1;
    m_HeadTail[rankNo][4] = headY[j]+nvY[j]-1;
    m_HeadTail[rankNo][5] = headZ[k]+nvZ[k]-1;

    RankInfo[rankNo].RankID=rankNo;
    RankInfo[rankNo].VoxelSize[0]=m_HeadTail[rankNo][3]-m_HeadTail[rankNo][0]+1;
    RankInfo[rankNo].VoxelSize[1]=m_HeadTail[rankNo][4]-m_HeadTail[rankNo][1]+1;
    RankInfo[rankNo].VoxelSize[2]=m_HeadTail[rankNo][5]-m_HeadTail[rankNo][2]+1;
    RankInfo[rankNo].HeadIndex[0]=m_HeadTail[rankNo][0];
    RankInfo[rankNo].HeadIndex[1]=m_HeadTail[rankNo][1];
    RankInfo[rankNo].HeadIndex[2]=m_HeadTail[rankNo][2];
    RankInfo[rankNo].TailIndex[0]=m_HeadTail[rankNo][3];
    RankInfo[rankNo].TailIndex[1]=m_HeadTail[rankNo][4];
    RankInfo[rankNo].TailIndex[2]=m_HeadTail[rankNo][5];

  }}}

  delete [] nvX;
  delete [] nvY;
  delete [] nvZ;
  delete [] headX;
  delete [] headY;
  delete [] headZ;

  return true;
}
///////////////////////////////////////////////////////////////////
// head&tail の生成
bool Staging::CreateHeadTail(int* rankMap)
{

  if( !rankMap ) return false;

  //領域分割数
  const int *div = GetDivNum();
  if( !div ) return false;

  //全体ボクセル数
  const int *tvox = GetVoxNum();
  if( !tvox ) return false;
  int gvox[3];
  for(int i=0; i<3; i++) gvox[i]=tvox[i];

/*
  int outGc = m_outGc;
  if( outGc > 0 ) {
    if( dfi_Finfo->GuideCell < outGc ) outGc=dfi_Finfo->GuideCell;
  }
*/

  int CropStart[3];
  int CropEnd[3];

  if( !m_cropIndexStart_on ) {
    for(int i=0; i<3; i++ ) CropStart[i]=1;
  } else {
     for(int i=0; i<3; i++ ) CropStart[i]=m_CropStart[i];
  }

  if( !m_cropIndexEnd_on ) {
    for(int i=0; i<3; i++ ) CropEnd[i]=dfi_Domain->GlobalVoxel[i];
  } else {
    for(int i=0; i<3; i++ ) CropEnd[i]=m_CropEnd[i];
  }

/*
  if( outGc > 0 ) {
    for(int i=0; i<3; i++) {
      if( CropStart[i]>1 ) CropStart[i]=CropStart[i]-outGc;
      if( CropEnd[i]<dfi_Domain->GlobalVoxel[i] ) CropEnd[i]=CropEnd[i]+outGc;
    }
  }
*/

  //MxNの入力指示あれている場合のボクセル数を更新
  //if( m_cropIndexStart_on || m_cropIndexEnd_on ) {
    gvox[0] = CropEnd[0]-CropStart[0]+1;
    gvox[1] = CropEnd[1]-CropStart[1]+1;
    gvox[2] = CropEnd[2]-CropStart[2]+1;
  //} else {
  //  for(int i=0; i<3; i++) m_CropStart[i]=1;
  //  for(int i=0; i<3; i++) m_CropEnd[i]=gvox[i];
  //}


  // ローカルのVOXEL数
  int *nvX = new int[div[0]];
  int *nvY = new int[div[1]];
  int *nvZ = new int[div[2]];
  int *nv[3] = {nvX,nvY,nvZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    //基準のボクセル数
    int nbase = gvox[n] / div[n];

    //余り
    int amari = gvox[n] % div[n];

    //ボクセル数をセット
    for( int i=0;i<div[n];i++ )
    {
      nvd[i] = nbase;
      if( i<amari ) nvd[i]++;
    }
  }

  //head
  int *headX = new int[div[0]];
  int *headY = new int[div[1]];
  int *headZ = new int[div[2]];
  int *head[3] = {headX,headY,headZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    int *hd = head[n];
    //hd[0] = 1;
    hd[0] = CropStart[n];

    for( int i=1;i<div[n];i++ )
    {
      hd[i] = hd[i-1]+nvd[i-1];
    }
  }

  Rank rank;

  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    int rankNo = rankMap[_IDX_S3D(i,j,k,div[0],div[1],div[2],0)];

    if( rankNo < 0 ) continue;

    m_HeadTail[rankNo][0] = headX[i];
    m_HeadTail[rankNo][1] = headY[j];
    m_HeadTail[rankNo][2] = headZ[k];
    m_HeadTail[rankNo][3] = headX[i]+nvX[i]-1;
    m_HeadTail[rankNo][4] = headY[j]+nvY[j]-1;
    m_HeadTail[rankNo][5] = headZ[k]+nvZ[k]-1;

    rank.RankID=rankNo;
    rank.VoxelSize[0]=m_HeadTail[rankNo][3]-m_HeadTail[rankNo][0]+1;
    rank.VoxelSize[1]=m_HeadTail[rankNo][4]-m_HeadTail[rankNo][1]+1;
    rank.VoxelSize[2]=m_HeadTail[rankNo][5]-m_HeadTail[rankNo][2]+1;
    rank.HeadIndex[0]=m_HeadTail[rankNo][0];
    rank.HeadIndex[1]=m_HeadTail[rankNo][1];
    rank.HeadIndex[2]=m_HeadTail[rankNo][2];
    rank.TailIndex[0]=m_HeadTail[rankNo][3];
    rank.TailIndex[1]=m_HeadTail[rankNo][4];
    rank.TailIndex[2]=m_HeadTail[rankNo][5];
    m_GRankInfo.push_back(rank);

  }}}

  delete [] nvX;
  delete [] nvY;
  delete [] nvZ;
  delete [] headX;
  delete [] headY;
  delete [] headZ;

  return true;
}

///////////////////////////////////////////////////////////////////
// head&tailをガイドセルで更新
void Staging::UpdateHeadTail(const int* mst_head, const int* mst_tail,
                             int* head,int* tail)
{

  //コピー元のヘッドとテイルをコピー
  for(int i=0; i<3; i++) {
    head[i]=mst_head[i];
    tail[i]=mst_tail[i];
  }

  //出力ガイドセル数がDFIより大きいときはDFIのガイドセルに更新
  int outGc = m_outGc;
  if( dfi_Finfo->GuideCell < outGc ) outGc=dfi_Finfo->GuideCell;
  
  //出力ガイドセルが無いときはリターン
  if( outGc < 1 ) return;

  for(int i=0; i<3; i++) {
    if( head[i]>1 ) head[i]=head[i]-outGc;
    if( tail[i]<dfi_Domain->GlobalVoxel[i] ) tail[i]=tail[i]+outGc;
  }

}


///////////////////////////////////////////////////////////////////
// head&tail の生成
void Staging::SetHeadTail()
{
  for(int i=0; i<m_GRankInfo.size(); i++) {
    m_HeadTail[i][0]=m_GRankInfo[i].HeadIndex[0];
    m_HeadTail[i][1]=m_GRankInfo[i].HeadIndex[1];
    m_HeadTail[i][2]=m_GRankInfo[i].HeadIndex[2];
    m_HeadTail[i][3]=m_GRankInfo[i].TailIndex[0];
    m_HeadTail[i][4]=m_GRankInfo[i].TailIndex[1];
    m_HeadTail[i][5]=m_GRankInfo[i].TailIndex[2];

    m_headX.insert(m_HeadTail[i][0]);
    m_headY.insert(m_HeadTail[i][1]);
    m_headZ.insert(m_HeadTail[i][2]);
  }
 
}

///////////////////////////////////////////////////////////////////
// head map の生成
void Staging::CreateHeadMap(std::set<int>head, headT &map)
{

  int cnt=0;
  for(std::set<int>::iterator it=head.begin(); it!=head.end(); it++)
  {
    int key=*it;
    map.insert(headT::value_type(key,cnt));
    cnt++;
  }

}
///////////////////////////////////////////////////////////////////
// 粗密判定
stg_EGlobalVoxel Staging::CheckGlobalVoxel(int Gvoxel[3], int DFI_Gvoxel[3])
{

  if( Gvoxel[0] == DFI_Gvoxel[0]   &&
      Gvoxel[1] == DFI_Gvoxel[1]   &&
      Gvoxel[2] == DFI_Gvoxel[2]   ) return STG_E_GV_SAME;

  if( Gvoxel[0] == DFI_Gvoxel[0]*2 &&
      Gvoxel[1] == DFI_Gvoxel[1]*2 &&
      Gvoxel[2] == DFI_Gvoxel[2]*2 ) return STG_E_GVX2_SAME;

  return STG_E_OTHER;
}

///////////////////////////////////////////////////////////////////
// ファイルコピー
bool Staging::FileCopy(vector<int>readRankList, int myRank, cdm_DFI *dfi)
{

  bool mio;
  if( dfi_Process->RankList.size()>1 ) mio=true;
  else mio=false;

  char*cmd = new char[512];
  char tmp[20];

  if( m_outPath == "" ) m_outPath=".";
  int len = m_outPath.size()+7;

  string fname_dfi;

  sprintf(tmp,"%06d", myRank);
  string path = m_outPath + string("/") + string(tmp);
  MakeDirectory(path);
  path += string("/");

  for(int i=0; i<readRankList.size(); i++) {

    int ostep = m_step;
    if( dfi_Finfo->FieldFilenameFormat == CDM::E_CDM_FNAME_RANK ) ostep = 0;

    //出力時刻指定あり
    if( ostep >= 0 ) {
      string path2;
      if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON ) {
        sprintf(tmp,"%010d",ostep);
        path2 = path+tmp;
        MakeDirectory(path2);
      } else {
        path2 = path;
      }

//20160407.fub.s
      if( dfi_Finfo->FileFormat != CDM::E_CDM_FMT_FUB ) {
//20160407.fub.e    
        fname_dfi = CDM::cdmPath_ConnectPath(m_inPath,Generate_FileName(readRankList[i],ostep,mio)); 
        memset(cmd, 0, sizeof(char)*512 );
        sprintf(cmd,"cp %s %s\n",fname_dfi.c_str(),path2.c_str());
        system(cmd);
//20160406.fub.s
      } else if ( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_FUB ) {
      //FubFileCopy(readRankList[i],ostep,mio,path2,dfi,cmd);
        FubFileCopy(dfi_Process->RankList[readRankList[i]].RankID,
                    ostep,mio,path2,dfi,cmd);
      }
//20160406.fub.e

    //出力指定なし
    }else if( ostep<0 ) {
      for(int j=0; j<dfi_TSlice->SliceList.size(); j++) {
        int step=dfi_TSlice->SliceList[j].step;
        string path2;
        if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON ) {
          sprintf(tmp,"%010d",step);
          path2 = path+tmp;
          MakeDirectory(path2);
        } else {
          path2 = path;
        }

//20160407.fub.s
        if( dfi_Finfo->FileFormat != CDM::E_CDM_FMT_FUB ) {
//20160407.fub.e
          fname_dfi = CDM::cdmPath_ConnectPath(m_inPath,Generate_FileName(readRankList[i],step,mio)); 
          memset(cmd, 0, sizeof(char)*512 );
          sprintf(cmd,"cp %s %s\n",fname_dfi.c_str(),path2.c_str());
          system(cmd);
//20160406.fub.s
        } else if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_FUB ) {
        //FubFileCopy(readRankList[i],step,mio,path2,dfi,cmd);
          FubFileCopy(dfi_Process->RankList[readRankList[i]].RankID,
                      step,mio,path2,dfi,cmd);
        }
//20160406.fub.e
      }
    }
  }
  memset(cmd, 0, sizeof(char)*512 );
  std::string procfname = CDM::cdmPath_ConnectPath(m_inPath,dfi_Fpath->ProcDFIFile);
  sprintf(cmd,"cp %s %s/%s_proc.dfi\n",procfname.c_str(),
          path.c_str(),dfi_Finfo->Prefix.c_str());
  system(cmd);

  if( cmd ) delete [] cmd;

  return true;
}

///////////////////////////////////////////////////////////////////
// ファイルコピー(FCONV用)
//20160509.fub.s
//bool Staging::FileCopy(step_rank_info info, int myRank)
bool Staging::FileCopy(step_rank_info info, int myRank, int ndfi)
//20160509.fub.e
{

  dfi_Process = (cdm_Process *)info.dfi->GetcdmProcess();
//20160408.fub.s
//dfi_Finfo   = info.dfi->GetcdmFileInfo();
  dfi_Finfo   = (cdm_FileInfo *)info.dfi->GetcdmFileInfo();
//20160408.fub.e
  dfi_TSlice  = info.dfi->GetcdmTimeSlice();
  dfi_Fpath   = info.dfi->GetcdmFilePath();
  dfi_Visit   = info.dfi->GetcdmVisIt();
  dfi_Unit    = info.dfi->GetcdmUnit();
  dfi_Domain  = info.dfi->GetcdmDomain();
  dfi_MPI     = info.dfi->GetcdmMPI();

  string dfi_fname = info.dfi->get_dfi_fname();
  m_inPath = CDM::cdmPath_DirName(dfi_fname);

  int outGc=m_outGc;
  if( outGc > 0 ) {
    if( dfi_Finfo->GuideCell < outGc ) outGc=dfi_Finfo->GuideCell;
  }

  int CropStart[3];
  int CropEnd[3];

  if( !m_cropIndexStart_on ) {
    for(int i=0; i<3; i++ ) CropStart[i]=1;
  } else {
    for(int i=0; i<3; i++ ) CropStart[i]=m_CropStart[i];
  }

  if( !m_cropIndexEnd_on ) {
    for(int i=0; i<3; i++ ) CropEnd[i]=dfi_Domain->GlobalVoxel[i];
  } else {
    for(int i=0; i<3; i++ ) CropEnd[i]=m_CropEnd[i];
  }

  if( outGc > 0 ) {
    for(int i=0; i<3; i++) {
      if( CropStart[i]>1 ) CropStart[i]=CropStart[i]-outGc;
      if( CropEnd[i]<dfi_Domain->GlobalVoxel[i] ) CropEnd[i]=CropEnd[i]+outGc;
    }
  }


  //printf("m_CropStart : %d %d %d\n",m_CropStart[0],m_CropStart[1],m_CropStart[2]);
  //printf("m_CropEnd   : %d %d %d\n",m_CropEnd[0],m_CropEnd[1],m_CropEnd[2]);

  bool mio;
  if( dfi_Process->RankList.size()>1 ) mio = true;
  else mio=false;

  char* cmd = new char[512];
  char tmp[20];
  string fname;

  if( m_outPath == "" ) m_outPath=".";
  int len = m_outPath.size()+7;

  //コピーするディレクトリの生成
  sprintf(tmp,"%06d", myRank);
  string path = m_outPath + string("/") + string(tmp);
  MakeDirectory(path);
  path += string("/"); 
 
  string path2;
  //ファイル名を生成しディレクトリのファイルをコピーする
  for(int i=info.stepStart; i<=info.stepEnd; i++) {
    int step = dfi_TSlice->SliceList[i].step;
    //TimeSlice 出力の場合のディレクトリの生成
    if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON ) {
      sprintf(tmp,"%010d",step);
      path2 = path+tmp;
      MakeDirectory(path2);
    } else {
      path2 = path;
    }
    for(int j=info.rankStart; j<=info.rankEnd; j++) {

      if( CropStart[0] > dfi_Process->RankList[j].TailIndex[0] ||
          CropEnd[0]   < dfi_Process->RankList[j].HeadIndex[0] ) continue;
      if( CropStart[1] > dfi_Process->RankList[j].TailIndex[1] ||
          CropEnd[1]   < dfi_Process->RankList[j].HeadIndex[1] ) continue;
      if( CropStart[2] > dfi_Process->RankList[j].TailIndex[2] ||
          CropEnd[2]   < dfi_Process->RankList[j].HeadIndex[2] ) continue;

//20160509.fub.s
      if( dfi_Finfo->FileFormat != CDM::E_CDM_FMT_FUB ) {
//20160509.fub.e
        fname=CDM::cdmPath_ConnectPath(m_inPath,Generate_FileName(j,step,mio));
        memset(cmd, 0, sizeof(char)*512 );
        sprintf(cmd,"cp %s %s\n",fname.c_str(),path2.c_str());
        system(cmd);
//20160509.fub.s
      } else {
        FubFileCopy(dfi_Process->RankList[j].RankID,step,mio,path2,info.dfi,
                    cmd); 
      }
//20160509.fub.e
    }
  }

  //proc.dfiのコピー
//20160509.fub.s
  if( dfi_Finfo->Prefix.empty() ) {
    char prefix[128];
    sprintf(prefix,"dfi%d",ndfi);
    dfi_Finfo->Prefix = prefix;
  }
//20160509.fub.e
  memset(cmd, 0, sizeof(char)*512 );
  string procfname = CDM::cdmPath_ConnectPath(m_inPath,dfi_Fpath->ProcDFIFile);
  sprintf(cmd,"cp %s %s/%s_proc.dfi\n",procfname.c_str(),path.c_str(),
          dfi_Finfo->Prefix.c_str());
  system(cmd);
 
  if( cmd ) delete [] cmd;

  //index.dfiの生成、出力
  string indexfname = path+CDM::cdmPath_FileName(dfi_fname,".dfi");
  WriteIndexDfiFile(indexfname, info);

  return true;
}

///////////////////////////////////////////////////////////////////
// dfiファイル出力
bool Staging::OutputDFI(string fname, int* rankMap, cdm_DFI *dfi)
{
  if( m_outPath == "" ) m_outPath="./";
  int len = m_outPath.size()+7;
  char tmp[20];

  string fname_dfi;
  int sta_x,sta_y,sta_z;
  int end_x,end_y,end_z;
  int div[3];
  for(int i=0; i<3; i++ ) div[i]=dfi_Domain->GlobalDivision[i];

  for(int k=0; k<m_Gdiv[2]; k++) {
  for(int j=0; j<m_Gdiv[1]; j++) {
  for(int i=0; i<m_Gdiv[0]; i++) {
    int myRank = rankMap[_IDX_S3D(i,j,k,m_Gdiv[0],m_Gdiv[1],m_Gdiv[2],0)];
    if( myRank < 0 ) continue;

    sprintf(tmp,"%06d",myRank);
    string path =m_outPath + string("/") + string(tmp);

    MakeDirectory(path);
    path += string("/");

    string DfiName = path+CDM::cdmPath_FileName(fname,".dfi");

    WriteIndexDfiFile(DfiName, dfi);

  }}}

  return true;
}
///////////////////////////////////////////////////////////////////
// ファイル名を作成
std::string Staging::Generate_FileName(int RankID, int step, const bool mio)
{

  if( dfi_Finfo->DirectoryPath.empty() ) return NULL;
  if( dfi_Finfo->Prefix.empty() ) return NULL;

  string fmt;
  if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_SPH ) {
    fmt=D_CDM_EXT_SPH;
  } else if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_BOV ) {
    fmt=D_CDM_EXT_BOV_DATAFILE;
//20150918.NetCDF.s
  } else if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_NETCDF4 ) {
    fmt=D_CDM_EXT_NC;
//20150918.NetCDF.e
//20160406.fub.s
  } else if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_FUB ) {
    fmt=D_CDM_EXT_FUB;
//20160406.fub.e
  }

#if 0
  int len = dfi_Finfo->DirectoryPath.size() + dfi_Finfo->Prefix.size() +
            fmt.size() + 25;

  if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON ) len += 11;

  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);

  if( mio ) {
    if( dfi_Finfo->FieldFilenameFormat == CDM::E_CDM_FNAME_RANK_STEP ) {
    //ファイル命名順がrank_step
      if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON ) {
        sprintf(tmp, "%s/%010d/%s_id%06d_%010d.%s",dfi_Finfo->DirectoryPath.c_str(),
                                         step,dfi_Finfo->Prefix.c_str(),
                                         RankID,step,fmt.c_str());
      } else {
        sprintf(tmp, "%s/%s_id%06d_%010d.%s",dfi_Finfo->DirectoryPath.c_str(),
                                         dfi_Finfo->Prefix.c_str(),
                                         RankID,step,fmt.c_str());
      }
    } else {
    //ファイル命名順がstep_rank または、未指定
      if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON ) {
        sprintf(tmp, "%s/%010d/%s_%010d_id%06d.%s",dfi_Finfo->DirectoryPath.c_str(),
                                         step,dfi_Finfo->Prefix.c_str(),
                                         step,RankID,fmt.c_str());
      } else {
        sprintf(tmp, "%s/%s_%010d_id%06d.%s",dfi_Finfo->DirectoryPath.c_str(),
                                         dfi_Finfo->Prefix.c_str(),
                                         step,RankID,fmt.c_str());
      }
    }
  } else {
    sprintf(tmp, "%s/%s_%010d.%s",dfi_Finfo->DirectoryPath.c_str(),
                                  dfi_Finfo->Prefix.c_str(),
                                  step,fmt.c_str());
  }
  std::string fname(tmp);
  if( tmp ) delete [] tmp;
#else
  std::string tmp = cdm_DFI::Generate_FileName(dfi_Finfo->Prefix, RankID, step, fmt, dfi_Finfo->FieldFilenameFormat,
                                               mio, dfi_Finfo->TimeSliceDirFlag, dfi_Finfo->RankNoPrefix);
  std::string fname = dfi_Finfo->DirectoryPath;
  fname += "/";
  fname += tmp;
#endif

  return fname;
}
///////////////////////////////////////////////////////////////////
// ディレクトリがなければ作成、既存なら何もしない
int Staging::MakeDirectory(string path)
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

  return 1;
}

///////////////////////////////////////////////////////////////////
// ディレクトリがなければ作成、既存なら何もしない
int Staging::MakeDirectoryPath()
{

  std::string path = Generate_Directory_Path();

  return MakeDirectory(path);
}

///////////////////////////////////////////////////////////////////
int Staging::MakeDirectorySub( std::string path )
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
std::string Staging::Generate_Directory_Path()
{

  // dfiのパスとDirectoryPathを連結する関数
  // ただし、絶対パスのときはdfiのパスは無視
  // CDM::cdmPath_isAbsoluteがtrueのとき絶対パス
  // DirectoryPath + TimeSliceDir
  std::string path = m_inPath;
  if( dfi_Finfo->TimeSliceDirFlag == CDM::E_CDM_ON )
  {
    path = CDM::cdmPath_ConnectPath(path, "");
  }

  if( CDM::cdmPath_isAbsolute(path) )
  {
    return path;
  }

  path = CDM::cdmPath_ConnectPath(m_inPath, path);
  return path;

}

// #################################################################
// Index DFIファイルの出力
CDM::E_CDM_ERRORCODE
Staging::WriteIndexDfiFile(const std::string dfi_name, cdm_DFI *dfi)
{
  if ( dfi_name.empty() ) return CDM::E_CDM_ERROR_WRITE_INDEXFILENAME_EMPTY;
  if ( dfi_Finfo->Prefix.empty() ) return CDM::E_CDM_ERROR_WRITE_PREFIX_EMPTY;

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
  cdm_FileInfo *t_Finfo = (cdm_FileInfo *)dfi_Finfo;
//FEAST 20140203.s
  t_Finfo->DirectoryPath="./";
//FEAST 20140203.s
  if( t_Finfo->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_FILEINFO;
  }

  //FilePath {} の出力
  cdm_FilePath t_Fpath;
  t_Fpath.ProcDFIFile = dfi_Finfo->Prefix+"_proc.dfi";
  if( t_Fpath.Write(fp, 1) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_FILEPATH;
  }
  
  //VisIt {} の出力
  cdm_VisIt *t_Visit = (cdm_VisIt *)dfi_Visit;
  if( t_Visit->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_VISIT;
  }

  //Unit {} の出力
  cdm_Unit *t_Unit = (cdm_Unit *)dfi_Unit;
  if( t_Unit->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_UNIT;
  }

  //TimeSlice {} の出力
  cdm_TimeSlice *t_Slice = new cdm_TimeSlice();
  int nsize = dfi_Finfo->NumVariables;
  if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_SPH && dfi_Finfo->NumVariables > 1 )
  {
    nsize++;
  }

  double* minmax = new double[nsize*2];
  for(int i=0; i<dfi_TSlice->SliceList.size(); i++)
  {
    int step = dfi_TSlice->SliceList[i].step;

    if( m_step >= 0 && m_step != step ) continue;

    double *ominmax = minmax;
    if( dfi_TSlice->SliceList[i].Min.size()!=nsize ||
        dfi_TSlice->SliceList[i].Max.size()!=nsize )
    {
      ominmax = NULL;
    }
    else
    {
      for(int n=0; n<nsize; n++)
      {
        minmax[n*2+0] = dfi_TSlice->SliceList[i].Min[n];
        minmax[n*2+1] = dfi_TSlice->SliceList[i].Max[n];
      }
    }

    t_Slice->AddSlice(step,
                      dfi_TSlice->SliceList[i].time,
                      ominmax,
                      dfi_Finfo->NumVariables,
                      dfi_Finfo->FileFormat,
                      dfi_TSlice->SliceList[i].avr_mode,
                      dfi_TSlice->SliceList[i].AveragedStep,
                      dfi_TSlice->SliceList[i].AveragedTime);
  }
  delete [] minmax;

  if ( t_Slice->Write(fp, 1, dfi_Finfo->FileFormat) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_TIMESLICE;
  }

/*
  //FileListの出力
  if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_FUB ) {
    cdm_DFI_FUB *dfi_fub = dynamic_cast<cdm_DFI_FUB*>(dfi);
    return dfi_fub->WriteFileList(fp, 1);
  }
*/
  //FieldFileNameFormatの出力
  if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_FUB ) { 
    cdm_FieldFileNameFormat* Ffformat = 
        (cdm_FieldFileNameFormat *)dfi->GetcdmFieldFileNameFormat();
    if( Ffformat ) Ffformat->Write(fp, 1);
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// Index DFIファイルの出力(FCONV用)
CDM::E_CDM_ERRORCODE
Staging::WriteIndexDfiFile(const std::string dfi_name, const step_rank_info info)
{

  if( dfi_name.empty() ) return CDM::E_CDM_ERROR_WRITE_INDEXFILENAME_EMPTY;
  if( dfi_Finfo->Prefix.empty() ) return CDM::E_CDM_ERROR_WRITE_PREFIX_EMPTY;

  FILE* fp = NULL;

  //File exist ?
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
  cdm_FileInfo *t_Finfo = (cdm_FileInfo *)dfi_Finfo;
  std::string tmp = t_Finfo->DirectoryPath;
  t_Finfo->DirectoryPath="./";
  if( t_Finfo->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_FILEINFO;
  }
  t_Finfo->DirectoryPath=tmp;

  //FilePath {} の出力
  cdm_FilePath t_Fpath;
  t_Fpath.ProcDFIFile = dfi_Finfo->Prefix+"_proc.dfi";
  if( t_Fpath.Write(fp, 1) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_FILEPATH;
  }

  //VisIt {} の出力
  cdm_VisIt *t_Visit = (cdm_VisIt *)dfi_Visit;
  if( t_Visit->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_VISIT;
  }
  
  //Unit {} の出力
 cdm_Unit *t_Unit = (cdm_Unit *)dfi_Unit;
  if( t_Unit->Write(fp, 0) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_UNIT;
  }

  //TimeSlice {} の出力
  cdm_TimeSlice *t_Slice = new cdm_TimeSlice();
  int nsize = dfi_Finfo->NumVariables;

//20160510.fub.s
/*
  if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_SPH && dfi_Finfo->NumVariables > 1 )
  {
    nsize++;
  }
*/
//20160510.fub.e

  double* minmax = new double[nsize*2];
  for(int i=info.stepStart; i<=info.stepEnd; i++) {
    int step = dfi_TSlice->SliceList[i].step;

//20160510.fub.s
    //if( m_step >= 0 && m_step != step ) continue;

    double* ominmax = minmax;
    if( dfi_TSlice->SliceList[i].Min.size() != nsize ||
        dfi_TSlice->SliceList[i].Max.size() != nsize ) {
      ominmax = NULL;  
    } else {
//20160510.fub.e

      for(int n=0; n<nsize; n++) {
        minmax[n*2+0] = dfi_TSlice->SliceList[i].Min[n];
        minmax[n*2+1] = dfi_TSlice->SliceList[i].Max[n];
      }

//20160510.fub.s
    }
//20160510.fub.e

    t_Slice->AddSlice(step,
                      dfi_TSlice->SliceList[i].time,
//20160510.fub.s
                    //minmax,
                      ominmax,
//20160510.fub.e
                      dfi_Finfo->NumVariables,
                      dfi_Finfo->FileFormat,
                      dfi_TSlice->SliceList[i].avr_mode,
                      dfi_TSlice->SliceList[i].AveragedStep,
                      dfi_TSlice->SliceList[i].AveragedTime);
  }
  if ( t_Slice->Write(fp, 1, dfi_Finfo->FileFormat) != CDM::E_CDM_SUCCESS )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_WRITE_TIMESLICE;
  }

//20160510.fub.s
//FieldFileNameFormatの出力
  if( dfi_Finfo->FileFormat == CDM::E_CDM_FMT_FUB )
  {
    cdm_FieldFileNameFormat* Ffformat =
       (cdm_FieldFileNameFormat *)info.dfi->GetcdmFieldFileNameFormat();
    if( Ffformat ) Ffformat->Write(fp,1);
  }
//20160510.fub.e

  return CDM::E_CDM_SUCCESS;
}

