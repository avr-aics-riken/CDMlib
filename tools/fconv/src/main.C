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
 * @file main.C
 * @brief convのmain関数
 * @author aics
 */

#include "InputParam.h"
#include "conv.h"

void usage(const char *progname)
{
  std::cerr
    << "Usage: " << progname << " <option> filename <options>\n"
    << " filename: file name when -f specified\n"
    << " Options:\n"
    << "  -f filename  : input file name for combine (ex : comb.tp)\n"
    << "  -d dirname   : output directory name for combine (this option is given priority over input file) \n"
    << "  -v verbose   : print more info\n"
    << "  -l log out   : print out logfile\n"
    << "  -h           : Show usage and exit\n"
    << std::endl;
}

int main( int argc, char **argv )
{
  char *progname = argv[0];
  bool out_comb = false;
  bool out_log  = false;
  bool thin_out = false;
  int thin_count=1;
  string fname;
  string dname="";
  //int pflag=0;//出力しない
  int pflag=1;//出力する
  int pflagv=0;
  int lflag;
  int lflagv;

  // タイミング用変数
  double t0, t1, t2, t3, t4, t5;

  // ##################################################################
  // 初期処理

  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う
  cpm_ParaManager* paraMngr = cpm_ParaManager::get_instance(argc, argv);
  if( !paraMngr ) {
    return CPM_ERROR_PM_INSTANCE;
  }

  // ##################################################################
  // 入力オプション処理

  int opt;
  while ((opt = getopt(argc, argv, "avlf:d:hs:")) != -1) {
    switch (opt) {
    case 'f':
      fname = optarg;
      out_comb = true;
      break;
/*
    case 'd':
      dname = optarg;
      break;
*/
    case 'v':
      pflagv = 1;
      break;
    case 'l':
      out_log = true;
      break;
/*
    case 's':
      thin_out = true;
      thin_count = atoi(optarg);
      break;
*/
    case 'h': // Show usage and exit
      usage(progname);
      return 0;
      break;
    case ':': // not find argument
      usage(progname);
      return 0;
      break;
    case '?': // unknown option
      usage(progname);
      return 0;
      break;
    }
  }

  // 入力ファイルが存在するかどうか
  if( !(out_comb) ){
    usage(progname);
    return 0;
  }

  // 画面出力、ログ出力の整理
  if(pflagv==1) pflag =1; 
  if(pflag ==0) pflagv=0; 
  lflag=0;
  lflagv=0;
  if(out_log){
    lflag=pflag;
    lflagv=pflagv;
  }

  // ##################################################################
  // InputParamインスタンス
  InputParam param(paraMngr);

  // ##################################################################
  // 入力ファイルの読込みとチェック
  cout << endl;
  cout << "Input Parameter File Read" << endl;
  t0 = cpm_Base::GetWTime();
  if( !param.Read(fname) ) {
    cout << "Input Parameter File Read Error" << endl;
    return 0;
  };

  // ##################################################################
  // conv のインスタンス
  CONV* conv = CONV::ConvInit(&param);
  if( conv == NULL ) return 0;
  if ( !conv->importCPM(paraMngr) )
  {
    return CPM_ERROR_PM_INSTANCE;
  }

  // ##################################################################
  // ログファイルのオープン
  int m_lflag=lflag;
  LOG_OUT_ conv->OpenLogFile();

  // ##################################################################
  // 出力指定ディレクトリのチェック
  if( dname.size() > 0 ) conv->CheckDir(dname);
  conv->CheckDir(param.Get_OutputDir());

  // ##################################################################
  // 引数のセット
  conv->m_pflag=pflag;
  conv->m_pflagv=pflagv;
  conv->m_lflag=lflag;
  conv->m_lflagv=lflagv;

  // ##################################################################
  // dfiファイルの読み込み
  cout << endl;
  cout << "ReadDfiFiles" << endl;
  t1 = cpm_Base::GetWTime();
  if( conv->ReadDfiFiles() != CDM::E_CDM_SUCCESS ) return 0;
 
  t2 = cpm_Base::GetWTime();
  // ##################################################################
  // VoxelInit
  conv->VoxelInit();

  // ##################################################################
  // 実行
  if( !conv->exec() ) return 0;

  // ##################################################################
  // 終了処理
  cout << endl;
  cout << "converter finish" << endl;
  t3 = cpm_Base::GetWTime();

  double tt[4];
  tt[0]=t1-t0;
  tt[1]=t2-t1;
  tt[2]=t3-t2;
  tt[3]=t3-t0;

  printf("\n\n");
  printf("TIME : ReadInit      %10.3f sec.\n", tt[0]);
  printf("TIME : ReadDfiFiles  %10.3f sec.\n", tt[1]);
  printf("TIME : ConvertFiles  %10.3f sec.\n", tt[2]);
  printf("TIME : Total Time    %10.3f sec.\n", tt[3]);
  LOG_OUT_ conv->WriteTime(tt);

  // ##################################################################
  // ログファイルのクローズ(デストラクタへ移動)
  //LOG_OUT_ conv->CloseLogFile();

  // ##################################################################
  // 並列環境の終了

  return 0;
}

