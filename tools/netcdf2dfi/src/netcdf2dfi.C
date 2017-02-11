/*
###################################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

/**
 * @file   netcdf2dfi.C
 * @brief  netcdf2dfi class
 * @author aics
 */

#include "netcdf2dfi.h"

// コンストラクタ
NetCDF2DFI::NetCDF2DFI()
{
  m_num_ncfiles = 0;
  m_DirectoryPath = "";
  m_Prefix = "";
  m_RankNoPrefix = CDM::C_CDM_RANKNOPREFIX;
  m_GuideCell = 0;
  m_num_steps = 0;
  m_OutputDirectoryPath = "";
  m_index_fname = "index.dfi";
  m_proc_fname = "proc.dfi";
  m_vecVariable.clear();
  m_FieldFilenameFormat = CDM::E_CDM_FNAME_DEFAULT;

  m_nameX = "x";
  m_nameY = "y";
  m_nameZ = "z";
  m_nameT = "time";

  m_nrank = 1;
  m_myrank = 0;

  m_indexFsta = 0;
  m_indexFend = 0;

  DFI_Domain = new cdm_Domain();
}

// デストラクタ
NetCDF2DFI::~NetCDF2DFI()
{
}

// インスタンス
NetCDF2DFI*
NetCDF2DFI::get_instance( int argc, char **argv, int nrank, int myrank )
{
  if( argc != 2 )
  {
    Usage(argv[0]);
    return NULL;
  }

  // 入力パラメータファイル名
  const char *ifname = argv[1];
  if( myrank == 0 )
  {
    printf( "input parameter file : %s\n", ifname );
  }

  // インスタンス
  NetCDF2DFI *ncPtr = new NetCDF2DFI;
  if( !ncPtr )
  {
    return NULL;
  }

  // MPI
  ncPtr->m_nrank  = nrank;
  ncPtr->m_myrank = myrank;

  // 入力パラメータファイルの読み込み
  if( !ncPtr->ReadParam(ifname) )
  {
    stmpd_printf("[%d] read param : %s\n", myrank, ifname );
    delete ncPtr;
    return NULL;
  }

  // 担当ファイルの開始終了インデクスのセット
  vector<int> vecFiles; //各ランクが担当するファイル数
  for( int i=0;i<nrank;i++ )
  {
    int nf = ncPtr->m_num_ncfiles / nrank;
    if( ncPtr->m_num_ncfiles % nrank > i )
    {
      nf++;
    }
    vecFiles.push_back(nf);
  }
  int indexFsta = 0;
  int indexFend = ncPtr->m_num_ncfiles - 1;
  for( int i=1;i<=myrank;i++ )
  {
    indexFsta += vecFiles[i-1];
  }
  indexFend = indexFsta + vecFiles[myrank] - 1;
  ncPtr->m_indexFsta = indexFsta;
  ncPtr->m_indexFend = indexFend;
//stmpd_printf( "[%d] start/end file index = %d - %d\n", myrank, indexFsta, indexFend );

  return ncPtr;
}

// usage
void
NetCDF2DFI::Usage( const char *prog )
{
  printf( "Usage :: %s ifname\n"
          "  ifname : input parameter filename\n"
        , prog);
}

// 入力パラメータファイルの読み込み
bool
NetCDF2DFI::ReadParam( const char *ifname )
{
  int ret;
  std::string leaf = "/NetCDF2DFI";
  std::string label, sdmy;
  int idmy;

  // ファイルの存在チェック
  FILE*fp = NULL;
  if( !(fp=fopen(ifname,"rb")) ) {
    printf("[%d] Can't open file. (%s)\n",m_myrank, ifname);
    return false;
  }
  fclose(fp);

  // オープン
  cdm_TextParser tpCntl;
  tpCntl.getTPinstance();
  if( (ret = tpCntl.readTPfile(std::string(ifname))) != TP_NO_ERROR )
  {
    printf("[%d] read error (%s). tp error code=%d\n", m_myrank, ifname, ret);
    return false;
  }
  TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    printf("[%d] read error (%s). tp pointer is NULL\n", m_myrank, ifname);
    return false;
  }

  // num_ncfiles
  label = leaf + "/num_ncfiles";
  if( !(tpCntl.GetValue(label, &idmy)) )
  {
    printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
    return false;
  }
  m_num_ncfiles = idmy;

  // DirectoryPath
  label = leaf + "/DirectoryPath";
  if( !(tpCntl.GetValue(label, &sdmy)) )
  {
    printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
    return false;
  }
  m_DirectoryPath = sdmy;

  // Prefix
  label = leaf + "/Prefix";
  if( !(tpCntl.GetValue(label, &sdmy)) )
  {
    printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
    return false;
  }
  m_Prefix = sdmy;

  // FieldFilenameFormat
  label = leaf + "/FieldFilenameFormat";
  if( tpCntl.GetValue(label, &sdmy) )
  {
    if( !strcasecmp(sdmy.c_str(),"step_rank" ) )
    {
      m_FieldFilenameFormat = CDM::E_CDM_FNAME_STEP_RANK;
    }
    else if( !strcasecmp(sdmy.c_str(),"rank_step" ) )
    {
      m_FieldFilenameFormat = CDM::E_CDM_FNAME_RANK_STEP;
    }
    else if( !strcasecmp(sdmy.c_str(),"rank" ) )
    {
      m_FieldFilenameFormat = CDM::E_CDM_FNAME_RANK;
    }
    else
    {
      printf("[%d] read error (%s). label=%s\n",m_myrank, ifname, label.c_str());
      return false;
    }
  }

  // RankNoPrefix
  label = leaf + "/RankNoPrefix";
  if( tpCntl.GetValue(label, &sdmy) )
  {
    m_RankNoPrefix = sdmy;
  }

  // GuideCell
  label = leaf + "/GuideCell";
  if( !(tpCntl.GetValue(label, &idmy)) )
  {
    printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
    return false;
  }
  m_GuideCell = idmy;

  // Variable
  int nvar=0;
  label = leaf + "/Variable";
  if( !tpCntl.chkNode(label) )
  {
    printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
    return false;
  }
  nvar = tpCntl.countLabels(label);
  if( nvar>0 )
  {
    // 子ノードのラベルを取得
    vector<std::string> labels;
    tp->changeNode(leaf);
    tp->getNodes(labels,1);
    for( int i=0;i<labels.size();i++ )
    {
      label = labels[i];
      if( strcasecmp(label.substr(8,1).c_str(), "[") ) continue;
      if( !strcasecmp(label.substr(0,8).c_str(), "variable") )
      {
        label = leaf+"/" + label + "/name";
        if ( !(tpCntl.GetValue(label, &sdmy )) )
        {
          printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
          return false;
        }
        else
        {
          m_vecVariable.push_back(sdmy);
        }
      }
    }
  }

  // num_steps
  if( m_FieldFilenameFormat != CDM::E_CDM_FNAME_RANK )
  {
    label = leaf + "/num_steps";
    if( !(tpCntl.GetValue(label, &idmy)) )
    {
      printf("[%d] read error (%s). label=%s\n", m_myrank, ifname, label.c_str());
      return false;
    }
    m_num_steps = idmy;
  }

  // OutputDirectoryPath
  label = leaf + "/OutputDirectoryPath";
  if( tpCntl.GetValue(label, &sdmy) )
  {
    m_OutputDirectoryPath = sdmy;
  }

  // index_fname
  label = leaf + "/index_fname";
  if( tpCntl.GetValue(label, &sdmy) )
  {
    m_index_fname = sdmy;
  }

  // proc_fname
  label = leaf + "/proc_fname";
  if( tpCntl.GetValue(label, &sdmy) )
  {
    m_proc_fname = sdmy;
  }

  // VariableName
  std::string leaf_leaf = leaf + "/VariableName";
  if( tpCntl.chkNode(leaf_leaf) )
  {
    // xを取得
    label = leaf_leaf + "/x";
    if ( tpCntl.GetValue(label, &sdmy) )
    {
      m_nameX = sdmy;
    }

    // yを取得
    label = leaf_leaf + "/y";
    if ( tpCntl.GetValue(label, &sdmy) )
    {
      m_nameY = sdmy;
    }

    // zを取得
    label = leaf_leaf + "/z";
    if ( tpCntl.GetValue(label, &sdmy) )
    {
      m_nameZ = sdmy;
    }

    // timeを取得
    label = leaf_leaf + "/time";
    if ( tpCntl.GetValue(label, &sdmy) )
    {
      m_nameT = sdmy;
    }
  }

  return true;
}

// 入力パラメータファイルの印刷
bool
NetCDF2DFI::PrintInputParam()
{
  printf( "\n---- input parmeters ----\n" );

  std::string leaf = "  /NetCDF2DFI";
  std::string label, sdmy;
  int idmy;

  // num_ncfiles
  label = leaf + "/num_ncfiles";
  printf( "%s = %d\n", label.c_str(), m_num_ncfiles);

  // DirectoryPath
  label = leaf + "/DirectoryPath";
  printf( "%s = %s\n", label.c_str(), m_DirectoryPath.c_str());

  // Prefix
  label = leaf + "/Prefix";
  printf( "%s = %s\n", label.c_str(), m_Prefix.c_str());

  // FieldFilenameFormat
  label = leaf + "/FieldFilenameFormat";
  const char *fdmy[4] = {"default", "step_rank", "rank_step", "rank"};
  printf( "%s = %s\n", label.c_str(), fdmy[(int)m_FieldFilenameFormat+1]);

  // RankNoPrefix
  label = leaf + "/RankNoPrefix";
  printf( "%s = %s\n", label.c_str(), m_RankNoPrefix.c_str());

  // GuideCell
  label = leaf + "/GuideCell";
  printf( "%s = %d\n", label.c_str(), m_GuideCell);

  // Variable
  label = leaf + "/Variable";
  for( int i=0;i<m_vecVariable.size();i++ )
  {
    printf( "%s[%d]/name = %s\n", label.c_str(), i, m_vecVariable[i].c_str());
  }

  // num_steps
  if( m_FieldFilenameFormat != CDM::E_CDM_FNAME_RANK )
  {
    label = leaf + "/num_steps";
    printf( "%s = %d\n", label.c_str(), m_num_steps);
  }

  // OutputDirectoryPath
  label = leaf + "/OutputDirectoryPath";
  printf( "%s = %s\n", label.c_str(), m_OutputDirectoryPath.c_str());

  // index_fname
  label = leaf + "/index_fname";
  printf( "%s = %s\n", label.c_str(), m_index_fname.c_str());

  // proc_fname
  label = leaf + "/proc_fname";
  printf( "%s = %s\n", label.c_str(), m_proc_fname.c_str());

  // VariableName
  std::string leaf_leaf = leaf + "/VariableName";
  printf( "%s/x = %s\n", leaf_leaf.c_str(), m_nameX.c_str());
  printf( "%s/y = %s\n", leaf_leaf.c_str(), m_nameY.c_str());
  printf( "%s/z = %s\n", leaf_leaf.c_str(), m_nameZ.c_str());
  printf( "%s/time = %s\n", leaf_leaf.c_str(), m_nameT.c_str());

  return true;
}

// ヘッダーレコードの読み込み
int
NetCDF2DFI::ReadHeader()
{
  CDM::E_CDM_ERRORCODE ret;

  if( m_myrank==0 )
  {
    printf( "\n---- read header record ----\n" );
  }

  // ファイル命名用変数のセット
  bool mio = false; //分散ファイルフラグ
  if( m_num_ncfiles > 1 ) mio = true;
  std::string ext = D_CDM_EXT_NC; //拡張子

  // 先頭座標値格納領域
  double *bufCoord = new double[m_num_ncfiles*3];
  for( int i=0;i<m_num_ncfiles*3;i++ )
  {
    bufCoord[i] = -DBL_MAX;
  }

  // 自ランクが担当するファイルを読み込み
  for( int ID=m_indexFsta;ID<=m_indexFend;ID++ )
  {
    //ファイル名(ステップ0)
    std::string fname = cdm_DFI::Generate_FileName(m_Prefix, ID, 0, ext, m_FieldFilenameFormat,
                                                   mio, CDM::E_CDM_OFF, m_RankNoPrefix);
    fname = m_DirectoryPath + "/" + fname;
//stmpd_printf( "[%d] fname = %s\n", m_myrank, fname.c_str() );

    // ファイルオープン
    cdm_FILE *pFile = cdm_FILE::OpenReadBinary( fname, CDM::E_CDM_FMT_NETCDF4 );
    if( !pFile )
    {
      printf("[%d] Can't open file. (%s)\n", m_myrank, fname.c_str());
      return 9;
    }

    // ヘッダー読み込み
    cdm_DFI_NETCDF::stVarInfo varInfoX, varInfoY, varInfoZ, varInfoT;
    ret = cdm_DFI_NETCDF::read_HeaderRecord( pFile->m_ncid, m_nameX, m_nameY, m_nameZ, m_nameT
                                           , varInfoX, varInfoY, varInfoZ, varInfoT );
    if( ret != CDM::E_CDM_SUCCESS )
    {
      printf("[%d] ERROR : read header record. [%s]\n", m_myrank, fname.c_str());
      printf("             CDM error code = %d\n", ret);
      cdm_FILE::CloseFile(pFile);
      return 9;
    }

    // 先頭の座標値をセット
    bufCoord[ID*3+0] = varInfoX.data[0+m_GuideCell];
    bufCoord[ID*3+1] = varInfoY.data[0+m_GuideCell];
    bufCoord[ID*3+2] = varInfoZ.data[0+m_GuideCell];

    // ファイルクローズ
    cdm_FILE::CloseFile(pFile);

    // 時刻数をセット(rankのとき)
    if( m_FieldFilenameFormat == CDM::E_CDM_FNAME_RANK )
    {
      m_num_steps = varInfoT.dims[0].len;
//stmpd_printf( "[%d] rank step = %d\n", m_myrank, m_num_steps );
    }
  }

  // 先頭座標値のマージ
  double *topCoord = new double[m_num_ncfiles*3];
  for( int i=0;i<m_num_ncfiles*3;i++ )
  {
    topCoord[i] = bufCoord[i];
  }
  m_setX.clear();
  m_setY.clear();
  m_setZ.clear();
  if( m_nrank > 1 )
  {
    MPI_Allreduce(bufCoord, topCoord, m_num_ncfiles*3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  for( int i=0;i<m_num_ncfiles;i++ )
  {
    m_setX.insert(topCoord[i*3+0]);
    m_setY.insert(topCoord[i*3+1]);
    m_setZ.insert(topCoord[i*3+2]);
  }

  if( m_myrank == 0 )
  {
#if 0
    stmpd_printf("X div = %d\n", m_setX.size());
    for( set<double>::iterator it=m_setX.begin();it!=m_setX.end();it++ )
    {
      stmpd_printf( "  %e\n", *it );
    }
    stmpd_printf("Y div = %d\n", m_setY.size());
    for( set<double>::iterator it=m_setY.begin();it!=m_setY.end();it++ )
    {
      stmpd_printf( "  %e\n", *it );
    }
    stmpd_printf("Z div = %d\n", m_setZ.size());
    for( set<double>::iterator it=m_setZ.begin();it!=m_setZ.end();it++ )
    {
      stmpd_printf( "  %e\n", *it );
    }
#endif
    printf( "  Global Division = %d %d %d\n", m_setX.size(), m_setY.size(), m_setZ.size() );
    if( m_FieldFilenameFormat == CDM::E_CDM_FNAME_RANK )
    {
      printf( "  num_steps in NC file = %d\n", m_num_steps );
    }
  }

  delete [] topCoord;
  delete [] bufCoord;

  // 正常終了
  return 0;
}

// データレコードの読み込み
int
NetCDF2DFI::ReadData()
{
  CDM::E_CDM_ERRORCODE ret;

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if( m_myrank==0 )
  {
    printf( "\n---- read data record ----\n" );
  }
  fflush(stdout);

  // ファイル命名用変数のセット
  bool mio = false; //分散ファイルフラグ
  if( m_num_ncfiles > 1 ) mio = true;
  std::string ext = D_CDM_EXT_NC; //拡張子

  // ファイルごとにため込む情報の宣言
  CDM::E_CDM_DTYPE dType = CDM::E_CDM_DTYPE_UNKNOWN;
  vector<cdm_TimeSlice> vecFileTimeSlice; //ファイルごとのTimeSlice情報
  vector<cdm_Rank> vecFileRank;           //ファイルごとのRank情報
  double coord1[3] = {0.0, 0.0, 0.0};     //下側の座標
  double coord2[3] = {0.0, 0.0, 0.0};     //上側の座標

  // 自ランクが担当するファイルを読み込み
  for( int ID=m_indexFsta;ID<=m_indexFend;ID++ )
  {
    printf( "  [%d] read file rank = %d (%d/%d)\n", m_myrank, ID, ID+1, m_num_ncfiles );

    // ファイルのTimeSliceを追加
    cdm_TimeSlice FileTimeSlice;

    // 時刻ループ
    for( int step=0;step<m_num_steps;step++ )
    {
      //ファイル名
      std::string fname = cdm_DFI::Generate_FileName(m_Prefix, ID, step, ext, m_FieldFilenameFormat,
                                                     mio, CDM::E_CDM_OFF, m_RankNoPrefix);
      fname = m_DirectoryPath + "/" + fname;
//stmpd_printf( "[%d] fname = %s\n", m_myrank, fname.c_str() );

      // ファイルオープン
      cdm_FILE *pFile = cdm_FILE::OpenReadBinary( fname, CDM::E_CDM_FMT_NETCDF4 );
      if( !pFile )
      {
        printf("[%d] Can't open file. (%s)\n", m_myrank, fname.c_str());
        return 9;
      }

      // ヘッダー読み込み
      cdm_DFI_NETCDF::stVarInfo varInfoX, varInfoY, varInfoZ, varInfoT;
      ret = cdm_DFI_NETCDF::read_HeaderRecord( pFile->m_ncid, m_nameX, m_nameY, m_nameZ, m_nameT
                                             , varInfoX, varInfoY, varInfoZ, varInfoT );
      if( ret != CDM::E_CDM_SUCCESS )
      {
        printf("[%d] ERROR : read header record. [%s]\n", m_myrank, fname.c_str());
        printf("             CDM error code = %d\n", ret);
        cdm_FILE::CloseFile(pFile);
        return 9;
      }

      // 読み込むステップインデクス
      unsigned istep = 0;
      if( m_FieldFilenameFormat == CDM::E_CDM_FNAME_RANK )
      {
        istep = step;
      }

      // データの読み込み
      std::vector<cdm_DFI_NETCDF::stVarInfo> vecVarInfo;
      cdm_Array* pArray = cdm_DFI_NETCDF::read_Datarecord
                            ( pFile->m_ncid, istep, m_vecVariable, CDM::E_CDM_DTYPE_UNKNOWN
                            , varInfoX, varInfoY, varInfoZ, varInfoT, vecVarInfo, ret );
      if( ret != CDM::E_CDM_SUCCESS )
      {
        printf("[%d] ERROR : read data record. [%s] step=%d\n", m_myrank, fname.c_str(), step);
        printf("             CDM error code = %d\n", ret);
        cdm_FILE::CloseFile(pFile);
        return 9;
      }

      // Rank情報
      cdm_Rank rank;
      rank.RankID = ID;
      rank.HostName = "unknown";
      rank.VoxelSize[0] = varInfoX.dims[0].len - 2*m_GuideCell;
      rank.VoxelSize[1] = varInfoY.dims[0].len - 2*m_GuideCell;
      rank.VoxelSize[2] = varInfoZ.dims[0].len - 2*m_GuideCell;

      // DFI情報のセット
      if( step == 0 )
      {
        // データタイプのセット
        dType = pArray->getDataType();

        // 座標値、時刻のUnit追加
        DFI_Unit.AddUnit(varInfoX.name, varInfoX.unit, 0.0);
        DFI_Unit.AddUnit(varInfoY.name, varInfoY.unit, 0.0);
        DFI_Unit.AddUnit(varInfoZ.name, varInfoZ.unit, 0.0);
        DFI_Unit.AddUnit(varInfoT.name, varInfoT.unit, 0.0);

        // 成分ごとのUnit追加
        for( int n=0;n<vecVarInfo.size();n++ )
        {
          DFI_Unit.AddUnit(vecVarInfo[n].name, vecVarInfo[n].unit, 0.0);
        }

        //headにはポジションをセットしておく
        set<double>::iterator it;
        int pos;
        for( it=m_setX.begin(), pos=0; it!=m_setX.end(); it++, pos++ )
        {
          if( *it == varInfoX.data[0+m_GuideCell] )
          {
            rank.HeadIndex[0] = pos;
            break;
          }
        }
        for( it=m_setY.begin(), pos=0; it!=m_setY.end(); it++, pos++ )
        {
          if( *it == varInfoY.data[0+m_GuideCell] )
          {
            rank.HeadIndex[1] = pos;
            break;
          }
        }
        for( it=m_setZ.begin(), pos=0; it!=m_setZ.end(); it++, pos++ )
        {
          if( *it == varInfoZ.data[0+m_GuideCell] )
          {
            rank.HeadIndex[2] = pos;
            break;
          }
        }
        vecFileRank.push_back(rank);

        // pitch
        if( rank.HeadIndex[0] == 0 )
        {
          double pitch = varInfoX.data[0+m_GuideCell+1] - varInfoX.data[0+m_GuideCell];
          coord1[0] = varInfoX.data[0+m_GuideCell] - pitch*0.5;
        }
        if( rank.HeadIndex[1] == 0 )
        {
          double pitch = varInfoY.data[0+m_GuideCell+1] - varInfoY.data[0+m_GuideCell];
          coord1[1] = varInfoY.data[0+m_GuideCell] - pitch*0.5;
        }
        if( rank.HeadIndex[2] == 0 )
        {
          double pitch = varInfoZ.data[0+m_GuideCell+1] - varInfoZ.data[0+m_GuideCell];
          coord1[2] = varInfoZ.data[0+m_GuideCell] - pitch*0.5;
        }
        if( rank.HeadIndex[0] == m_setX.size()-1 )
        {
          double pitch = varInfoX.data[varInfoX.dims[0].len-1-m_GuideCell] - varInfoX.data[varInfoX.dims[0].len-1-m_GuideCell-1];
          coord2[0] = varInfoX.data[varInfoX.dims[0].len-1-m_GuideCell] + pitch*0.5;
        }
        if( rank.HeadIndex[1] == m_setY.size()-1 )
        {
          double pitch = varInfoY.data[varInfoY.dims[0].len-1-m_GuideCell] - varInfoY.data[varInfoY.dims[0].len-1-m_GuideCell-1];
          coord2[1] = varInfoY.data[varInfoY.dims[0].len-1-m_GuideCell] + pitch*0.5;
        }
        if( rank.HeadIndex[2] == m_setZ.size()-1 )
        {
          double pitch = varInfoZ.data[varInfoZ.dims[0].len-1-m_GuideCell] - varInfoZ.data[varInfoZ.dims[0].len-1-m_GuideCell-1];
          coord2[2] = varInfoZ.data[varInfoZ.dims[0].len-1-m_GuideCell] + pitch*0.5;
        }
      }

      // 時刻
      double time;
      if( m_FieldFilenameFormat == CDM::E_CDM_FNAME_RANK )
      {
        time = varInfoT.data[step];
      }
      else
      {
        time = varInfoT.data[0];
      }

      // Slice情報を取得
      cdm_Slice slice = GetSliceInfo(step, time, pArray);
      FileTimeSlice.SliceList.push_back(slice);

      // ファイルクローズ
      cdm_FILE::CloseFile(pFile);

      // 変数抽出ファイルの出力
      if( m_OutputDirectoryPath.length() > 0 )
      {
        WriteFieldData(ID, step, time, rank.VoxelSize,
                       varInfoX, varInfoY, varInfoZ,
                       varInfoT,
                       vecVarInfo, pArray);
      }

      // 配列削除
      delete pArray;

#if 0
{
  string ff = cdm_DFI::Generate_FileName(m_Prefix, ID, step, "txt", CDM::E_CDM_FNAME_STEP_RANK,
                                         mio, CDM::E_CDM_OFF, m_RankNoPrefix);
  FILE *fp = fopen(ff.c_str(), "wt");
  slice.Write(fp, 0, CDM::E_CDM_FMT_NETCDF4);
  fclose(fp);
}
#endif
    }// for step

    // TimeSliceを追加
    vecFileTimeSlice.push_back(FileTimeSlice);

  } // for file


  // 出力DFI情報の生成
  {
    // FileInfo
    DFI_Finfo.DFIType             = CDM::E_CDM_DFITYPE_CARTESIAN;
    DFI_Finfo.FieldFilenameFormat = m_FieldFilenameFormat;
    DFI_Finfo.DirectoryPath       = m_DirectoryPath;
    if( m_OutputDirectoryPath.length() > 0 )
    {
      DFI_Finfo.DirectoryPath     = m_OutputDirectoryPath;
    }
    DFI_Finfo.TimeSliceDirFlag    = CDM::E_CDM_OFF;
    DFI_Finfo.Prefix              = m_Prefix;
    DFI_Finfo.FileFormat          = CDM::E_CDM_FMT_NETCDF4;
    DFI_Finfo.GuideCell           = m_GuideCell;
    DFI_Finfo.DataType            = dType;
//    DFI_Finfo.Endian              = CDM::E_CDM_ENDIANTYPE_UNKNOWN;
    {
      // endianはnaitiveをセット
      int idumy = 1;
      char* cdumy = (char*)(&idumy);
      if( cdumy[0] == 0x01 ) DFI_Finfo.Endian = CDM::E_CDM_LITTLE;
      if( cdumy[0] == 0x00 ) DFI_Finfo.Endian = CDM::E_CDM_BIG;
    }
    DFI_Finfo.ArrayShape          = CDM::E_CDM_IJKN;
    DFI_Finfo.NumVariables        = (int)m_vecVariable.size();
    DFI_Finfo.VariableName.clear();
    for( int n=0;n<m_vecVariable.size();n++ )
    {
      DFI_Finfo.VariableName.push_back(m_vecVariable[n]);
    }
    DFI_Finfo.RankNoPrefix        = m_RankNoPrefix;

    // NetCDF4
    // 出力用のcdm_DFI_NETCDFインスタンス時にセット

    // Visit
    // 何もなし

    // FilePath
    DFI_Fpath.ProcDFIFile = m_proc_fname;

    // UnitList
    // データ読み込み時に収集済み

    // TimeSlice
    // ファイル毎のTimeSlice情報をマージ
    mergeTimeSlice(vecFileTimeSlice);

    // Domain情報, Process情報
    mergeRank(vecFileRank, coord1, coord2);

    // MPI情報
    DFI_MPI.NumberOfRank  = m_num_ncfiles;
    DFI_MPI.NumberOfGroup = 1;
  }

  // 正常終了
  return 0;
}

// Slice情報を取得
cdm_Slice
NetCDF2DFI::GetSliceInfo(int step, double time, cdm_Array *pArray)
{
  // 基本情報のセット
  cdm_Slice slice;
  slice.step = step;
  slice.time = time;
  slice.avr_mode = true;
  slice.AveragedStep=0;
  slice.AveragedTime=0.0;

  // min/maxの取得(double or floatの場合のみ)
  if( pArray->getDataType() == CDM::E_CDM_FLOAT32 )
  {
    cdm_TypeArray<float> *pTypeArray = (cdm_TypeArray<float>*)pArray;
    GetMinMax(pTypeArray, m_GuideCell, slice.Min, slice.Max);
  }
  else if( pArray->getDataType() == CDM::E_CDM_FLOAT64 )
  {
    cdm_TypeArray<double> *pTypeArray = (cdm_TypeArray<double>*)pArray;
    GetMinMax(pTypeArray, m_GuideCell, slice.Min, slice.Max);
  }

  return slice;
}

// ファイル毎のTimeSlice情報をマージ
void
NetCDF2DFI::mergeTimeSlice(vector<cdm_TimeSlice> &vecFileTimeSlice)
{
  // 成分ごとにmin/maxをマージする
  int nbuf = m_vecVariable.size()*m_num_steps;
  double *sbuf = new double[nbuf];
  double *rbuf = new double[nbuf];

  // 1番目のファイルの情報をコピー
  DFI_TimeSlice.SliceList = vecFileTimeSlice[0].SliceList;

  //*** min ***//
  // 成分ループ
  for( int n=0;n<m_vecVariable.size();n++ )
  {
    // 時間ループ
    for( int step=0;step<m_num_steps;step++ )
    {
      double minv = DBL_MAX;
      for( int i=0;i<vecFileTimeSlice.size();i++ )
      {
        minv = min( minv, vecFileTimeSlice[i].SliceList[step].Min[n] );
      }
      sbuf[n*m_num_steps+step] = minv;
    }
  }

  // プロセス間でマージ
  if( m_nrank == 1 )
  {
    for( int n=0;n<nbuf;n++ )
    {
      rbuf[n] = sbuf[n];
    }
  }
  else
  {
    MPI_Allreduce(sbuf, rbuf, nbuf, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }

  //*** max ***//
  // 成分ループ
  for( int n=0;n<m_vecVariable.size();n++ )
  {
    // 時間ループ
    for( int step=0;step<m_num_steps;step++ )
    {
      double maxv = -DBL_MAX;
      for( int i=0;i<vecFileTimeSlice.size();i++ )
      {
        maxv = max( maxv, vecFileTimeSlice[i].SliceList[step].Max[n] );
      }
      sbuf[n*m_num_steps+step] = maxv;
    }
  }

  // プロセス間でマージ
  if( m_nrank == 1 )
  {
    for( int n=0;n<nbuf;n++ )
    {
      rbuf[n] = sbuf[n];
    }
  }
  else
  {
    MPI_Allreduce(sbuf, rbuf, nbuf, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }

  // TimeSliceにセット
  for( int n=0;n<m_vecVariable.size();n++ )
  {
    for( int step=0;step<m_num_steps;step++ )
    {
      DFI_TimeSlice.SliceList[step].Max[n] = rbuf[n*m_num_steps+step];
    }
  }


  delete [] sbuf;
  delete [] rbuf;
}

// Domain情報, Process情報のマージ
void
NetCDF2DFI::mergeRank(vector<cdm_Rank> &vecFileRank, double coord1[3], double coord2[3])
{
  // 領域情報をマージする
  int nbuf = m_num_ncfiles*6;
  int *sbuf = new int[nbuf];
  int *rbuf = new int[nbuf];
  for( int n=0;n<nbuf;n++ )
  {
    sbuf[n] = 0;
  }

  // ファイル毎にID位置にVoxelSizeと位置情報を格納
  for( int n=0;n<vecFileRank.size();n++ )
  {
    int ID = vecFileRank[n].RankID;
    sbuf[ID*6+0] = vecFileRank[n].VoxelSize[0];
    sbuf[ID*6+1] = vecFileRank[n].VoxelSize[1];
    sbuf[ID*6+2] = vecFileRank[n].VoxelSize[2];
    sbuf[ID*6+3] = vecFileRank[n].HeadIndex[0];
    sbuf[ID*6+4] = vecFileRank[n].HeadIndex[1];
    sbuf[ID*6+5] = vecFileRank[n].HeadIndex[2];
  }

  // プロセス間でマージ
  if( m_nrank == 1 )
  {
    for( int n=0;n<nbuf;n++ )
    {
      rbuf[n] = sbuf[n];
    }
  }
  else
  {
    MPI_Allreduce(sbuf, rbuf, nbuf, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  }

  // 領域分割位置のVoxelSize格納領域
  vector<int> vecVoxelSizeX(m_setX.size());
  vector<int> vecVoxelSizeY(m_setY.size());
  vector<int> vecVoxelSizeZ(m_setZ.size());

  // ランク情報を生成
  DFI_Process.RankList.clear();
  for( int ID=0;ID<m_num_ncfiles;ID++ )
  {
    //ランク情報を格納
    cdm_Rank rank;
    rank.RankID = ID;
    rank.HostName = "unknown";
    rank.VoxelSize[0] = rbuf[ID*6+0];
    rank.VoxelSize[1] = rbuf[ID*6+1];
    rank.VoxelSize[2] = rbuf[ID*6+2];
    rank.HeadIndex[0] = rbuf[ID*6+3]; //この時点では領域分割位置
    rank.HeadIndex[1] = rbuf[ID*6+4];
    rank.HeadIndex[2] = rbuf[ID*6+5];
    DFI_Process.RankList.push_back(rank);

    //  領域分割位置のVoxelSizeをセット
    vecVoxelSizeX[rank.HeadIndex[0]] = rank.VoxelSize[0];
    vecVoxelSizeY[rank.HeadIndex[1]] = rank.VoxelSize[1];
    vecVoxelSizeZ[rank.HeadIndex[2]] = rank.VoxelSize[2];
  }

  // HeadIndexとGlobalVoxel
  int head;
  DFI_Domain->GlobalVoxel[0] = DFI_Domain->GlobalVoxel[1] = DFI_Domain->GlobalVoxel[2] = 0;
  vector<int> vecHeadIndexX;
  vector<int> vecHeadIndexY;
  vector<int> vecHeadIndexZ;
  head = 1;
  for( int i=0;i<vecVoxelSizeX.size();i++ )
  {
    vecHeadIndexX.push_back(head);
    head += vecVoxelSizeX[i];
    DFI_Domain->GlobalVoxel[0] += vecVoxelSizeX[i];
  }
  head = 1;
  for( int i=0;i<vecVoxelSizeY.size();i++ )
  {
    vecHeadIndexY.push_back(head);
    head += vecVoxelSizeY[i];
    DFI_Domain->GlobalVoxel[1] += vecVoxelSizeY[i];
  }
  head = 1;
  for( int i=0;i<vecVoxelSizeZ.size();i++ )
  {
    vecHeadIndexZ.push_back(head);
    head += vecVoxelSizeZ[i];
    DFI_Domain->GlobalVoxel[2] += vecVoxelSizeZ[i];
  }

  // head/tailをセット
  for( int n=0;n<DFI_Process.RankList.size();n++ )
  {
    int posX = DFI_Process.RankList[n].HeadIndex[0];
    int posY = DFI_Process.RankList[n].HeadIndex[1];
    int posZ = DFI_Process.RankList[n].HeadIndex[2];

    DFI_Process.RankList[n].HeadIndex[0] = vecHeadIndexX[posX];
    DFI_Process.RankList[n].HeadIndex[1] = vecHeadIndexY[posY];
    DFI_Process.RankList[n].HeadIndex[2] = vecHeadIndexZ[posZ];

    DFI_Process.RankList[n].TailIndex[0] = DFI_Process.RankList[n].HeadIndex[0] + vecVoxelSizeX[posX] - 1;
    DFI_Process.RankList[n].TailIndex[1] = DFI_Process.RankList[n].HeadIndex[1] + vecVoxelSizeY[posY] - 1;
    DFI_Process.RankList[n].TailIndex[2] = DFI_Process.RankList[n].HeadIndex[2] + vecVoxelSizeZ[posZ] - 1;
  }

  // GlobalDivision
  DFI_Domain->GlobalDivision[0] = m_setX.size();
  DFI_Domain->GlobalDivision[1] = m_setY.size();
  DFI_Domain->GlobalDivision[2] = m_setZ.size();

  // coord
  double scoord[6], rcoord[6];
  for( int i=0;i<3;i++ )
  {
    scoord[i]   = rcoord[i]   = coord1[i];
    scoord[i+3] = rcoord[i+3] = coord2[i];
  }
  if( m_nrank > 1 )
  {
    MPI_Allreduce(scoord, rcoord, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  }
  DFI_Domain->GlobalOrigin[0] = rcoord[0];
  DFI_Domain->GlobalOrigin[1] = rcoord[1];
  DFI_Domain->GlobalOrigin[2] = rcoord[2];
  DFI_Domain->GlobalRegion[0] = rcoord[3] - rcoord[0];
  DFI_Domain->GlobalRegion[1] = rcoord[4] - rcoord[1];
  DFI_Domain->GlobalRegion[2] = rcoord[5] - rcoord[2];

  delete [] sbuf;
  delete [] rbuf;
}

// DFIファイルの出力
int
NetCDF2DFI::WriteDFI()
{
  CDM::E_CDM_ERRORCODE ret;

  // DFIの出力
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if( m_myrank==0 )
  {
    printf( "\n---- write DFI file ----\n" );
  }
  fflush(stdout);

  // ランク0だけが出力
  if( m_myrank==0 )
  {
    // 出力用のdfiインスタンス
    cdm_DFI_NETCDF *oDFI = new cdm_DFI_NETCDF(DFI_Finfo, DFI_Fpath, DFI_VisIt, DFI_Unit, DFI_Domain, DFI_MPI, DFI_TimeSlice, DFI_Process);
    if( !oDFI )
    {
      printf( "output DFI instance error\n" );
      return 9;
    }

    // dimension名のセット
    oDFI->SetDimName( m_nameX, m_nameY, m_nameZ, m_nameT );

    FILE *fp = NULL;

    //*** index.dfi ***//
    fp = fopen(m_index_fname.c_str(), "wt");
    if( !fp )
    {
      printf( "output file open error. [%s]\n", m_index_fname.c_str() );
      delete oDFI;
      return 9;
    }

    // FileInfo
    DFI_Finfo.Write(fp, 0);

    // NetCDF4
    oDFI->WriteAdditionalTP(fp, 0);

    // FilePath
    DFI_Fpath.Write(fp, 1);

    // Visit
    DFI_VisIt.Write(fp, 1);

    // UnitList
    DFI_Unit.Write(fp, 0);

    // TimeSlice
    DFI_TimeSlice.Write(fp, 1, CDM::E_CDM_FMT_NETCDF4);

    // file close
    fclose(fp);
    printf( "  output %s\n", m_index_fname.c_str() );

    //*** proc.dfi ***//
    fp = fopen(m_proc_fname.c_str(), "wt");
    if( !fp )
    {
      printf( "output file open error. [%s]\n", m_proc_fname.c_str() );
      delete oDFI;
      return 9;
    }

    // Domain
    DFI_Domain->Write(fp, 0);

    // MPI
    DFI_MPI.Write(fp, 0);

    // Process
    DFI_Process.Write(fp, 0);

    // file close
    fclose(fp);
    printf( "  output %s\n", m_proc_fname.c_str() );

    delete oDFI;

    // extract message
    if( m_OutputDirectoryPath.length() > 0 )
    {
      printf( "  output extract files to [%s] directory\n", m_OutputDirectoryPath.c_str() );
    }
  }

  return 0;
}

// NetCDFファイル出力(変数抽出時のみ)
int
NetCDF2DFI::WriteFieldData(int ID, int step, double time, int VoxelSize[3],
                           cdm_DFI_NETCDF::stVarInfo &varInfoX, cdm_DFI_NETCDF::stVarInfo &varInfoY, cdm_DFI_NETCDF::stVarInfo &varInfoZ,
                           cdm_DFI_NETCDF::stVarInfo &varInfoT,
                           vector<cdm_DFI_NETCDF::stVarInfo> &vecVarInfo, cdm_Array *pArray)
{
  // 出力(変数抽出時のみ)
  if( m_OutputDirectoryPath.length() == 0 )
  {
    return 0;
  }

  // ファイル命名用変数のセット
  bool mio = false; //分散ファイルフラグ
  if( m_num_ncfiles > 1 ) mio = true;
  std::string ext = D_CDM_EXT_NC; //拡張子

  CDM::E_CDM_OUTPUT_FNAME FieldFilenameFormat = m_FieldFilenameFormat;

  //ファイル名
  std::string fname = cdm_DFI::Generate_FileName(m_Prefix, ID, step, ext, FieldFilenameFormat,
                                                 mio, CDM::E_CDM_OFF, m_RankNoPrefix);
  fname = m_OutputDirectoryPath + "/" + fname;
//  printf( "  [%d] write file [%s]. rank = %d (%d/%d). step=%d\n", m_myrank, fname.c_str(), ID, ID+1, m_num_ncfiles, step );

  // ディレクトリの生成
  std::string dir = CDM::cdmPath_DirName(fname);
  if( cdm_DFI::MakeDirectorySub(dir) != 0 )
  {
    printf( "error : make directory [%s]\n", dir.c_str() );
    return 9;
  }

  // ファイルオープン
  cdm_FILE* fp = NULL;
  bool addMode = false; //追記モードフラグ
  if( step>0 && FieldFilenameFormat==CDM::E_CDM_FNAME_RANK )
  {
    addMode = true;
  }
  if( (fp = cdm_FILE::OpenWriteBinary(fname,CDM::E_CDM_FMT_NETCDF4, addMode)) == NULL ) {
    printf("Can't open file.(%s)\n",fname.c_str());
    return CDM::E_CDM_ERROR_OPEN_FIELDDATA;
  }

  // ヘッダーの出力
  CDM::E_CDM_ERRORCODE eRet;
  eRet = cdm_DFI_NETCDF::write_HeaderRecord(fp, (unsigned)step, time, VoxelSize, m_GuideCell,
                                            pArray->getDataType(), m_vecVariable, DFI_Unit,
                                            varInfoX, varInfoY, varInfoZ, varInfoT, vecVarInfo);
  if( eRet != CDM::E_CDM_SUCCESS )
  {
    printf( "[%d] error : write header record [%s]\n  CDM error code = %d\n", m_myrank, fname.c_str(), (int)eRet );
    cdm_FILE::CloseFile(fp);
    return 9;
  }



  // データの出力
  eRet = cdm_DFI_NETCDF::write_DataRecord(fp, pArray, m_GuideCell, VoxelSize, varInfoT, vecVarInfo);
  if( eRet != CDM::E_CDM_SUCCESS )
  {
    printf( "[%d] error : write data record [%s]\n  CDM error code = %d\n", m_myrank, fname.c_str(), (int)eRet );
    cdm_FILE::CloseFile(fp);
    return 9;
  }



  // ファイルクローズ
  cdm_FILE::CloseFile(fp);

  return 0;
}
