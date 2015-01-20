/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_Process.C
 * @brief  cdm_Rank & cdm_Process Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cdm_Rank::cdm_Rank()
{

 RankID = 0;
 HostName = "";
 for(int i=0; i<3; i++) {
   VoxelSize[i]=0;
   HeadIndex[i]=0;
   TailIndex[i]=0;
 }
 c_id = -1;
 bc_id = -1;

}


// #################################################################
// デストラクタ
cdm_Rank::~cdm_Rank()
{

}

// #################################################################
// DFIファイル:Rank要素を読み込む
CDM::E_CDM_ERRORCODE
cdm_Rank::Read(cdm_TextParser tpCntl,
               std::string label_leaf)
{
#if 0
  std::string str;
  std::string label;
  int ct;
  int iv[3];

  //ID
  label = label_leaf + "/ID";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_ID;
  }
  else {
    RankID= ct;
  }

  //HostName
  label = label_leaf + "/HostName";
  if ( !(tpCntl.GetValue(label, &str )) ) {
    printf("\tCDM Parsing error : fail to get '%s'\n", label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_HOSTNAME;
  }
  HostName= str;

  //VoxelSize
  label = label_leaf + "/VoxelSize";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.GetVector(label, iv, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_VOXELSIZE;
  }
  VoxelSize[0]=iv[0];
  VoxelSize[1]=iv[1];
  VoxelSize[2]=iv[2];

  //HeadIndex
  label = label_leaf + "/HeadIndex";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.GetVector(label, iv, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_HEADINDEX;
  }
  HeadIndex[0]=iv[0];
  HeadIndex[1]=iv[1];
  HeadIndex[2]=iv[2];

  //TailIndex
  label = label_leaf + "/TailIndex";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.GetVector(label, iv, 3 )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_TAILINDEX;
  }
  TailIndex[0]=iv[0];
  TailIndex[1]=iv[1];
  TailIndex[2]=iv[2];

#else
  std::string str;
  std::string label;
  int ct;
  int iv[3];

  // TextParser
  TextParser *tp = tpCntl.getTPPtr();
  tp->changeNode(label_leaf);

  //ID
  label = "ID";
  if ( !(tpCntl.GetValue(label, &ct, false )) ) {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_ID;
  }
  else {
    RankID= ct;
  }

  //HostName
  label = "HostName";
  if ( !(tpCntl.GetValue(label, &str, false )) ) {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_HOSTNAME;
  }
  HostName= str;

  //VoxelSize
  label = "VoxelSize";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.GetVector(label, iv, 3, false )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_VOXELSIZE;
  }
  VoxelSize[0]=iv[0];
  VoxelSize[1]=iv[1];
  VoxelSize[2]=iv[2];

  //HeadIndex
  label = "HeadIndex";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.GetVector(label, iv, 3, false )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_HEADINDEX;
  }
  HeadIndex[0]=iv[0];
  HeadIndex[1]=iv[1];
  HeadIndex[2]=iv[2];

  //TailIndex
  label = "TailIndex";
  for (int n=0; n<3; n++) iv[n]=0.0;
  if ( !(tpCntl.GetVector(label, iv, 3, false )) ) 
  {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_TAILINDEX;
  }
  TailIndex[0]=iv[0];
  TailIndex[1]=iv[1];
  TailIndex[2]=iv[2];

  //CellID
  label = "CellID";
  if ( !(tpCntl.GetValue(label, &ct, false )) ) {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_CELLID;
  }
  else {
    c_id= ct;
  }

  //BCflagID
  label = "BCflagID";
  if ( !(tpCntl.GetValue(label, &ct, false )) ) {
    printf("\tCDM Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_BCFLAGID;
  }
  else {
    bc_id= ct;
  }
#endif

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// DFIファイル:Rank要素を出力する
CDM::E_CDM_ERRORCODE
cdm_Rank::Write(FILE* fp, 
                const unsigned tab) 
{

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "ID        = %d\n", RankID);

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "HostName  = \"%s\"\n", HostName.c_str());

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "VoxelSize = (%d, %d, %d)\n", VoxelSize[0],
                                              VoxelSize[1],
                                              VoxelSize[2]);

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "HeadIndex = (%d, %d, %d)\n", HeadIndex[0],
                                              HeadIndex[1],
                                              HeadIndex[2]);

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "TailIndex = (%d, %d, %d)\n", TailIndex[0],
                                              TailIndex[1],
                                              TailIndex[2]);

    fprintf(fp, "\n");

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "CellID    = %d\n", c_id);

    _CDM_WRITE_TAB(fp, tab);
    fprintf(fp, "BCflagID  = %d\n", bc_id);

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// コンストラクタ
cdm_Process::cdm_Process()
{

  m_rankMap=NULL;

}

// #################################################################
// デストラクタ
cdm_Process::~cdm_Process()
{
  if( m_rankMap ) delete m_rankMap;
}

// #################################################################
// DFIファイル:Process要素を読み込む
CDM::E_CDM_ERRORCODE
cdm_Process::Read(cdm_TextParser tpCntl)
{
#if 0
  std::string str;
  std::string label_base,label_leaf;
  int nnode=0;
  CDM::E_CDM_ERRORCODE iret;  ///<リターンコード

  cdm_Rank rank;

  //Process 
  nnode=0;
  label_base = "/Process";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }

  for (int i=0; i<nnode; i++) {

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tCDM Parsing error : No Elem name\n");
      return CDM::E_CDM_ERROR_READ_DFI_NO_RANK;
    }
    if( strcasecmp(str.substr(0,4).c_str(), "Rank") ) continue;
    label_leaf=label_base+"/"+str;

    /** Rankの読込み */
    iret = rank.Read(tpCntl, label_leaf);
    if( iret == CDM::E_CDM_SUCCESS ) {
      RankList.push_back(rank); 
    } else  return iret;

  }
#else
  CDM::E_CDM_ERRORCODE iret;

  // TP
  TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    return CDM::E_CDM_ERROR_TEXTPARSER;
  }

  // Process要素の存在チェック
  std::string label_base = "/Process";
  if( !tpCntl.chkNode(label_base) )
  {
    printf("\tCDM Parsing error : No Elem name [%s]\n", label_base.c_str());
    return CDM::E_CDM_ERROR_READ_PROCESS;
  }

  // /Processに移動
  tp->changeNode(label_base);

  // 子のラベルを取得
  vector<std::string> labels;
  tp->getNodes(labels,1);

  // 子のRankを読み込み
  for( size_t i=0;i<labels.size();i++ )
  {
    // Rank要素かどうか確認
    std::string label = labels[i];
    if( strcasecmp(label.substr(0,4).c_str(), "Rank") ) continue;

    // Rank要素の読込み
    cdm_Rank rank;
    std::string leaf = label_base + "/" + label;
    if( (iret = rank.Read(tpCntl,leaf)) == CDM::E_CDM_SUCCESS )
    {
      RankList.push_back(rank);
    }
    else
    {
      return iret;
    }
  }
#endif

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// 読込みランクリストの作成
CDM::E_CDM_ERRORCODE 
cdm_Process::CheckReadRank(const cdm_Domain* dfi_domain, 
                           const int head[3],
                           const int tail[3],
                           CDM::E_CDM_READTYPE readflag,
                           vector<int> &ReadRankList)
{

  headT mapHeadX,mapHeadY,mapHeadZ;

  //DFIにProcess/Rank[@]がない処理
  if( RankList.empty() ) {
    CDM::E_CDM_ERRORCODE ret = CreateRankList(dfi_domain, mapHeadX, mapHeadY, mapHeadZ);   
    if( ret != CDM::E_CDM_SUCCESS ) return ret;
  }

  //rankMapが未定義（DFIにProcess/Rank[@]がある場合）
  if( m_rankMap == NULL ) {
    m_rankMap = CreateRankMap(dfi_domain->GlobalDivision, mapHeadX, mapHeadY, mapHeadZ); 
  }

  //mapHeadX,mapHeadY,mapHeadZが未定義
  if( mapHeadX.empty() || mapHeadY.empty() || mapHeadZ.empty() ) {
    std::set<int>headx,heady,headz;
    for(int i=0; i<RankList.size(); i++ ) {
      headx.insert(RankList[i].HeadIndex[0]);
      heady.insert(RankList[i].HeadIndex[1]);
      headz.insert(RankList[i].HeadIndex[2]);
    }
    CreateHeadMap(headx,mapHeadX);
    CreateHeadMap(heady,mapHeadY);
    CreateHeadMap(headz,mapHeadZ);
  }
  
  return CheckStartEnd(dfi_domain, head, tail, readflag, mapHeadX, mapHeadY, mapHeadZ, ReadRankList);

}

// #################################################################
// DFIのProcessにHeadIndex,TailIndex指定が無い場合RankListを生成
CDM::E_CDM_ERRORCODE 
cdm_Process::CreateRankList(const cdm_Domain* dfi_domain,
                            map<int,int> &mapHeadX,
                            map<int,int> &mapHeadY,
                            map<int,int> &mapHeadZ)
{

  vector<cdm_ActiveSubDomain> subDomainInfo;

  CDM::E_CDM_ERRORCODE ret;
 
  ret = CreateSubDomainInfo(dfi_domain,subDomainInfo);
  if( ret != CDM::E_CDM_SUCCESS ) return ret;
  
  m_rankMap = CreateRankMap(dfi_domain->GlobalDivision,subDomainInfo);

  ret = CreateRankList(dfi_domain->GlobalDivision, dfi_domain->GlobalVoxel,
                       mapHeadX, mapHeadY, mapHeadZ);

  return ret;

}

// #################################################################
// ActiveSubDomain情報を作成
CDM::E_CDM_ERRORCODE 
cdm_Process::CreateSubDomainInfo(const cdm_Domain* domain,
                                 vector<cdm_ActiveSubDomain> &subDomainInfo)
{
  if( !domain->ActiveSubdomainFile.empty() ) {
    int divSudomain[3] = {0,0,0};
    CDM::E_CDM_ERRORCODE ret = ReadActiveSubdomainFile( domain->ActiveSubdomainFile,
                                                   subDomainInfo, divSudomain);
    if( ret != CDM::E_CDM_SUCCESS ) return ret;
  } else {
    //活性サブドメイン情報
    for( int k=0;k<domain->GlobalDivision[2];k++ ){
    for( int j=0;j<domain->GlobalDivision[1];j++ ){
    for( int i=0;i<domain->GlobalDivision[0];i++ ){
      int pos[3] = {i,j,k};
      cdm_ActiveSubDomain dom( pos );
      subDomainInfo.push_back(dom);
    }}}

  }

  return CDM::E_CDM_SUCCESS;
}
// #################################################################
// subDomainをもとにCPM同様の分割方法でRankListを生成する
CDM::E_CDM_ERRORCODE 
cdm_Process::CreateRankList(const int div[3],
                            const int gvox[3],
                            map<int,int> &mapHeadX,
                            map<int,int> &mapHeadY,
                            map<int,int> &mapHeadZ )
{
  if( !m_rankMap ) return CDM::E_CDM_ERROR;
  int ndiv = div[0]*div[1]*div[2];

  cdm_Rank rank;

  //ローカルのVOXEL数
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

  CreateHeadMap(headX,div[0],mapHeadX);
  //CreateHeadMap(headY,div[0],mapHeadY);
  CreateHeadMap(headY,div[1],mapHeadY);
  //CreateHeadMap(headZ,div[0],mapHeadZ);
  CreateHeadMap(headZ,div[2],mapHeadZ);

  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    int rankNo = m_rankMap[_CDM_IDX_IJK(i,j,k,div[0],div[1],div[2],0)];
    if( rankNo < 0 ) continue;
    rank.RankID = rankNo;
    rank.VoxelSize[0]=(headX[i]+nvX[i]-1)-headX[i]+1;
    rank.VoxelSize[1]=(headY[j]+nvY[j]-1)-headY[j]+1;
    rank.VoxelSize[2]=(headZ[k]+nvZ[k]-1)-headZ[k]+1;
    rank.HeadIndex[0]=headX[i];
    rank.HeadIndex[1]=headY[j];
    rank.HeadIndex[2]=headZ[k];
    rank.TailIndex[0]=headX[i]+nvX[i]-1;
    rank.TailIndex[1]=headY[j]+nvY[j]-1;
    rank.TailIndex[2]=headZ[k]+nvZ[k]-1;
    RankList.push_back(rank);
  }}}

  delete [] nvX;
  delete [] nvY;
  delete [] nvZ;
  delete [] headX;
  delete [] headY;
  delete [] headZ;

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ActiveSubdomainファイルの読み込み(static関数)
CDM::E_CDM_ERRORCODE 
cdm_Process::ReadActiveSubdomainFile(
                     std::string subDomainFile,
                     std::vector<cdm_ActiveSubDomain>& subDomainInfo,
                     int div[3] )
{
  if( subDomainFile.empty() ) return CDM::E_CDM_ERROR_OPEN_SBDM;

  // ファイルオープン
  FILE*fp = fopen( subDomainFile.c_str(), "rb" );
  if( !fp ) return CDM::E_CDM_ERROR_OPEN_SBDM;

  //エンディアン識別子
  int ident;
    if( fread( &ident, sizeof(int), 1, fp ) != 1 )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_SBDM_HEADER;
  }

 //エンディアンチェック
  if( isMatchEndianSbdmMagick( ident ) < 0 )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_SBDM_FORMAT;
  }

  // 領域分割数
  if( fread( div, sizeof(int), 3, fp ) != 3 )
  {
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_SBDM_DIV;
  }

  if( isMatchEndianSbdmMagick( ident ) == 0 ) BSWAPVEC(div,3);

  // contents
  size_t nc = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  unsigned char *contents = new unsigned char[nc];
  if( fread( contents, sizeof(unsigned char), nc, fp ) != nc )
  {
    delete [] contents;
    fclose(fp);
    return CDM::E_CDM_ERROR_READ_SBDM_CONTENTS;
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
      cdm_ActiveSubDomain dom( pos );
      subDomainInfo.push_back(dom);
    }
    ptr++;
  }}}

  // contentsのdelete
  delete [] contents;

  // 活性ドメインの数をチェック
  if( subDomainInfo.size() == 0 )
  {
    return CDM::E_CDM_ERROR_SBDM_NUMDOMAIN_ZERO;
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// ActiveSubdomainファイルのエンディアンチェック
int cdm_Process::isMatchEndianSbdmMagick( int ident )
{
  char magick_c[] = "SBDM";
  int  magick_i=0;

  //cheak match
  magick_i = (magick_c[3]<<24) + (magick_c[2]<<16) + (magick_c[1]<<8) + magick_c[0];
  if( magick_i == ident )
  {
    return 1;
  }

  //chack unmatch
  magick_i = (magick_c[0]<<24) + (magick_c[1]<<16) + (magick_c[2]<<8) + magick_c[3];
  if( magick_i == ident )
  {
    return 0;
  }

  //unknown format
  return -1;
}

// #################################################################
// ランクマップを生成 （非活性を含む）
int* cdm_Process::CreateRankMap(
                                const int div[3],
                                std::vector<cdm_ActiveSubDomain> &subDomainInfo)
{

  size_t ndiv = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  int *rankMap = new int[ndiv];
  if( !rankMap ) return NULL;

  for( size_t i=0;i<ndiv;i++ ) rankMap[i] = -1;

  // 活性サブドメイン情報配置位置に0をセット
  for( int i=0; i<subDomainInfo.size(); i++ )
  {
    //サブドメイン情報
    const cdm_ActiveSubDomain dom = subDomainInfo[i];

    //位置を取得
    const int *pos = dom.GetPos();
    if( !pos )
    {
      delete [] rankMap;
      return NULL;
    }

    //0をセット
    rankMap[_CDM_IDX_IJK(pos[0],pos[1],pos[2],div[0],div[1],div[2],0)] = 0;
  }

 //i->j->kの優先順で活性サブドメインにランク番号をセット
  int rankCount = 0;
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    if( rankMap[_CDM_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] == 0 )
    {
      rankMap[_CDM_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] = rankCount;
      rankCount++;
    }
  }}}

  return rankMap;
}

// #################################################################
//head map の生成 set版
void cdm_Process::CreateHeadMap(std::set<int>head,
                                headT &map)
{

  map.clear();

  int cnt=0;
  for(std::set<int>::iterator it=head.begin();it!=head.end();it++)
  {
    int key=*it;
    map.insert(headT::value_type(key,cnt));
    cnt++;
  }
}

// #################################################################
//head map の生成 int配列版
void cdm_Process::CreateHeadMap(int* head,
                                int ndiv,
                                headT &map)
{
  
  map.clear();

  for(int i=0; i<ndiv; i++)
  {
    map.insert(headT::value_type(head[i],i));
  }
}


// #################################################################
//DFI用ランクマップを生成
int* cdm_Process::CreateRankMap(const int div[3],
                                headT &mapHeadX,
                                headT &mapHeadY,
                                headT &mapHeadZ)
{

  int i,j,k;

  std::set<int>headx,heady,headz;
  for(int i=0; i<RankList.size(); i++ ) {
    headx.insert(RankList[i].HeadIndex[0]);
    heady.insert(RankList[i].HeadIndex[1]);
    headz.insert(RankList[i].HeadIndex[2]);
  }
  CreateHeadMap(headx,mapHeadX);
  CreateHeadMap(heady,mapHeadY);
  CreateHeadMap(headz,mapHeadZ);

  size_t ndiv = div[0]*div[1]*div[2];

  int *rankMap = new int[ndiv];
  for(int i=0; i<ndiv; i++) rankMap[i]=-1;

  headT::iterator it;

  for(int n=0; n<RankList.size(); n++)
  {
    it=mapHeadX.find(RankList[n].HeadIndex[0]);
    i=it->second;

    it=mapHeadY.find(RankList[n].HeadIndex[1]);
    j=it->second;

    it=mapHeadZ.find(RankList[n].HeadIndex[2]);
    k=it->second;

    int rnkPos=_CDM_IDX_IJK(i,j,k,div[0],div[1],div[2],0);

    rankMap[_CDM_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] = n;
  }

  return rankMap;

}

// #################################################################
//  読込みランクファイルリストの作成
CDM::E_CDM_ERRORCODE 
cdm_Process::CheckStartEnd(const cdm_Domain* dfi_domain,
                           const int head[3],
                           const int tail[3],
                           CDM::E_CDM_READTYPE readflag,
                           headT mapHeadX,
                           headT mapHeadY,
                           headT mapHeadZ,
                           vector<int> &readRankList)
{

  int StartEnd[6];

  int ndiv = dfi_domain->GlobalDivision[0]*
             dfi_domain->GlobalDivision[1]*
             dfi_domain->GlobalDivision[2];
  int head2[3],tail2[3];
  if( readflag == CDM::E_CDM_SAMEDIV_SAMERES || readflag == CDM::E_CDM_DIFFDIV_SAMERES ) {
    for(int i=0; i<3; i++) {
      head2[i]=head[i];
      tail2[i]=tail[i];
    }
  } else {
    for(int i=0; i<3; i++) {
      if( head[i] < 0 ) head2[i]=head[i]/2;
      else              head2[i]=(head[i]+1)/2;
      if( tail[i] < 0 ) tail2[i]=tail[i]/2;
      else              tail2[i]=(tail[i]+1)/2;
    }
  }

//x方向の絞り込み
  for( headT::iterator it=mapHeadX.begin();it!=mapHeadX.end();it++ )
  {
    if( head2[0] >= (*it).first ) StartEnd[0] = (*it).second;
    if( tail2[0] >= (*it).first ) StartEnd[3] = (*it).second;
    else break;
  }

//y方向の絞り込み
  for( headT::iterator it=mapHeadY.begin();it!=mapHeadY.end();it++ )
  {
    if( head2[1] >= (*it).first ) StartEnd[1] = (*it).second;
    if( tail2[1] >= (*it).first ) StartEnd[4] = (*it).second;
    else break;
  }

//z方向の絞り込み
  for( headT::iterator it=mapHeadZ.begin();it!=mapHeadZ.end();it++ )
  {
    if( head2[2] >= (*it).first ) StartEnd[2] = (*it).second;
    if( tail2[2] >= (*it).first ) StartEnd[5] = (*it).second;
    else break;
  }

  readRankList.clear();

  for(int k=StartEnd[2]; k<=StartEnd[5]; k++) {
  for(int j=StartEnd[1]; j<=StartEnd[4]; j++) {
  for(int i=StartEnd[0]; i<=StartEnd[3]; i++) {
    int rank = m_rankMap[_CDM_IDX_IJK(i,j,k,dfi_domain->GlobalDivision[0],
                                          dfi_domain->GlobalDivision[1],
                                          dfi_domain->GlobalDivision[2],0)];
    if( rank<0 ) continue;

    readRankList.push_back(rank);

  }}}

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// DFIファイル:Process要素を出力する
CDM::E_CDM_ERRORCODE
cdm_Process::Write(FILE* fp, 
                   const unsigned tab) 
{

  fprintf(fp, "Process {\n");
  fprintf(fp, "\n");

  cdm_Rank rank;

  for(int i=0; i<RankList.size(); i++) {
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Rank[@] {\n");
    fprintf(fp, "\n");

    rank = RankList[i];

    //Rank要素の出力
    if( rank.Write(fp,tab+2) != CDM::E_CDM_SUCCESS ) return CDM::E_CDM_ERROR;

    fprintf(fp, "\n");
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "}\n");
  }

  fprintf(fp, "\n");
  fprintf(fp, "}\n");

  return CDM::E_CDM_SUCCESS;
}


