/*
 *
 *
 *
 *
 */

#define MPI_COMPILE

#include <mpi.h>
#include <set>
#include "fub_DATA.h"

/** 定義点タイプ */
enum DataType
{
  CELL = 0 ///< セル
, NODE = 1 ///< ノード
};

bool getParameter(TextParser *tp, int &vc, DataType &Dtype, double &tol,
                  int &StartFileNo, int &NumberOfBlock,
                  int &StartStepNo, int &EndStepNo,
                  string &DirPath);

void printParameter(int vc, DataType Dtype, double tol, 
                    int StartFileNo, int NumberOfBlock,
                    int StartStepNo, int EndStepNo,
                    string DirPath); 

void WriteDFI(FILE *fp, int vc, vector<string> valname,string endian,
              string dftype, string DirPath, set<int>StepNoList,
              unsigned int tab);

// #################################################################
// メイン

int main( int argc, char *argv[] ) {

#ifdef MPI_COMPILE 
  int nprocs,myRank;
  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if( nprocs > 1 ) 
  {
    if( myRank == 0 ) {
      printf("Error Number of Proccess is Greater Than 1 : %d\n",nprocs);
    }
    MPI_Finalize();  
    return 0;
  }  
#endif
  int ret=0;
  double tol = 0.0;
  int vc=0;

  //パラメータローダのインスタンス生成
  cdm_TextParser tpCntl;
  tpCntl.getTPinstance();

  //パラメータのロード
  if( argc>1 )
  {
    //入力ファイルの指定
    string input_file = argv[1];

    //file open check
    FILE *fp;
    if( !(fp = fopen(input_file.c_str(),"r")) ) 
    {
      printf("Error at file open : %s file\n",input_file.c_str());
      return 0;
    }

    if( (ret=tpCntl.readTPfile(input_file)) != TP_NO_ERROR )
    {
      printf("Error at reading %s file : %d\n",input_file.c_str(),ret);
      return 0;
    }
  } else {
    printf("Error undefined input file\n");
    return 0;
  }

  int      GuideCell;     ///< 仮想セル数
  DataType Dtype;         ///< DataType "Cell" "Node"
  double   Tolerance;     ///< 2点間距離判定用のトレランス値
  int      StartFileNo;   ///< スタートファイル（ランク）番号
  int      NumberOfBlock; ///< スタートファイル番号からの数
  int      StartStepNo;   ///< スタートステップ番号
  int      EndStepNo;     ///< エンドステップ番号
  string   DirPath;       ///< 入力ファイルのディレクトリパス

  //TP
    TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    return CDM::E_CDM_ERROR_TEXTPARSER;
  }

  if( !getParameter(tp, GuideCell, Dtype, Tolerance,
                    StartFileNo, NumberOfBlock,
                    StartStepNo, EndStepNo,
                     DirPath) ) {
    return CDM::E_CDM_ERROR_TEXTPARSER; 
  }

  vc=GuideCell;

  // FieldFileNameFormatクラスのインスタンス
  cdm_FieldFileNameFormat FFF;

  //FieldFileNameFormatクラスでのParse
  if( FFF.Read(tpCntl) != CDM::E_CDM_SUCCESS ) {
    //FFF.Print();
    return 0;
  }
  FFF.Print();

  vector<int>    fileNoList; ///< File No (Rank No) list
  vector<fub_DATA *> vecFubData;

  //ステップ番号リスト、ランク番号リストの生成
  set<int> CStepList;
  set<int> FStepList;
  set<int> FRankList;
  for(int n=0,i=StartFileNo; n<NumberOfBlock; n++,i++)
  {
    for(int j=StartStepNo; j<=EndStepNo; j++)
    {
      if( FFF.FileExist("CoordinateFile", DirPath, j, i) ) {
        CStepList.insert(j);
      }
      if( FFF.FileExist("FieldFile", DirPath, j, i) ) {
        FStepList.insert(j);
        FRankList.insert(i);
      }
    }
  }

  //取得ステップ番号、ランク番号の出力（デバッグ用）
  set<int>::iterator it = CStepList.begin();
  while( it != CStepList.end() )
  {
    printf("CoordinateFile step no : %d\n",*it);
    it++;
  }

  it = FStepList.begin();
  while( it != FStepList.end() )
  {
    printf("FieldFile step no      : %d\n",*it);
    it++;
  }

  it = FRankList.begin();
  while( it != FRankList.end() )
  {
    printf("FieldFile rank id      : %d\n",*it);
    it++;
  }

  FILE *fp=NULL;

  string endian;
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  if( cdumy[0] == 0x01 ) endian="little";
  if( cdumy[0] == 0x00 ) endian="big";

  //座標値データの読み込み
  string fubCname,fubFname;
  string dftype="";
  set<int>::iterator itr;
  set<int>::iterator its;

  itr = FRankList.begin();
  while( itr != FRankList.end() ) 
  {
    its = CStepList.begin();
    while( its != CStepList.end() ) 
    {
      fubCname = FFF.GenerateFileName("CoordinateFile", DirPath, *its, *itr);
      fubFname = FFF.GenerateFileName("FieldFile", DirPath, *its, *itr);
      fub_DATA *fub = new fub_DATA(fubCname, fubFname);
      if( !fub->readFub(GuideCell, dftype) ) {
        printf("Error fub file read error\n");
        delete fub;
        return 0;
      }
      fub->fileID=*itr;
      vecFubData.push_back(fub);
      endian = fub->endian;
      its++;
    }
    itr++;
  }

  //フィールドデータの成分数の取得
  int nVal = vecFubData[0]->getNval();

  //隣接ノードの検索
  for(int i=0; i<vecFubData.size(); i++) {
    fub_DATA *data1 = vecFubData[i];
    int cnt=0;
    for(int j=i+1; j<vecFubData.size(); j++) {
      fub_DATA *data2 = vecFubData[j];
      data1->checkNeighbor(data2, cnt, tol);
      if( cnt >= 6 ) break;
    }
  }

// 一番下のノードを検索
  fub_DATA *firstNode=NULL;
  for(int i=0; i<vecFubData.size(); i++) {
    if( vecFubData[i]->nID[X_MINUS] == NULL &&
        vecFubData[i]->nID[Y_MINUS] == NULL &&
        vecFubData[i]->nID[Z_MINUS] == NULL ) {
      firstNode = vecFubData[i];
      break;
    }
  }

//分割数の獲得(idiv)
  fub_DATA *checkNode=firstNode;
  int idiv=1;
  while( checkNode->nID[X_PLUS] ) 
  {
    idiv++;
    checkNode = checkNode->nID[X_PLUS];
  }

//分割数の獲得(jdiv)
  checkNode=firstNode;
  int jdiv=1;
  while( checkNode->nID[Y_PLUS] ) 
  {
    jdiv++;
    checkNode = checkNode->nID[Y_PLUS];
  }

//分割数の獲得(kdiv)
  checkNode=firstNode;
  int kdiv=1;
  while( checkNode->nID[Z_PLUS] ) 
  {
    kdiv++;
    checkNode = checkNode->nID[Z_PLUS];
  }

// ノードマップの作成
  fub_DATA **fubDataMap = new fub_DATA*[idiv*jdiv*kdiv];
  for(int k=0; k<kdiv; k++) {
  for(int j=0; j<jdiv; j++) {
  for(int i=0; i<idiv; i++) {
    fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)] = NULL;
  }}}

  fubDataMap[_IDX_S3D(0,0,0,idiv,jdiv,kdiv,0)]=firstNode;

  fub_DATA::setFubDataMap(0,0,0,idiv,jdiv,kdiv,fubDataMap);

  for(int k=0; k<kdiv; k++) {
  for(int j=0; j<jdiv; j++) {
  for(int i=0; i<idiv; i++) {
    fub_DATA* tmp = fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)];
    if( !tmp ) {
      printf("fubDataMap[%d,%d,%d] is NULL !!!!",i,j,k);
      delete [] fubDataMap;
      return 0;
    }
  }}}

// Head のセット
  int xHead[idiv+1],yHead[jdiv+1],zHead[kdiv+1];
  xHead[0]=yHead[0]=zHead[0]=1;
  // Data Type "Cell"
  if( Dtype == CELL ) {
    // X
    for(int i=0; i<idiv; i++) {
      xHead[i+1] = xHead[i] + fubDataMap[_IDX_S3D(i,0,0,idiv,jdiv,kdiv,0)]->size[0]-2*vc;
    }
    // Y
    for(int j=0; j<jdiv; j++) {
      yHead[j+1] = yHead[j] + fubDataMap[_IDX_S3D(0,j,0,idiv,jdiv,kdiv,0)]->size[1]-2*vc;
    }
    // Z
    for(int k=0; k<kdiv; k++) {
      zHead[k+1] = yHead[k] + fubDataMap[_IDX_S3D(0,0,k,idiv,jdiv,kdiv,0)]->size[2]-2*vc;
    }
  } else if( Dtype == NODE ) {
    // X
    for(int i=0; i<idiv; i++) {
      xHead[i+1] = xHead[i] + fubDataMap[_IDX_S3D(i,0,0,idiv,jdiv,kdiv,0)]->size[0]-2*vc-1;
    }
    xHead[idiv]++;
    // Y
    for(int j=0; j<jdiv; j++) {
      yHead[j+1] = yHead[j] + fubDataMap[_IDX_S3D(0,j,0,idiv,jdiv,kdiv,0)]->size[1]-2*vc-1;
    }
    yHead[jdiv]++;
    // Z
    for(int k=0; k<kdiv; k++) {
      zHead[k+1] = yHead[k] + fubDataMap[_IDX_S3D(0,0,k,idiv,jdiv,kdiv,0)]->size[2]-2*vc-1;
    }
    zHead[kdiv]++;
  }

// proc.dfi 出力
  fub_DATA::WriteProc(firstNode,idiv,jdiv,kdiv,xHead,yHead,zHead,vc,fubDataMap);
  printf("created proc.dfi file\n");

  vector<string>valname;

// field.dfi 出力
  char vname[5];
  for(int i=0; i<nVal; i++) {
    sprintf(vname,"val%d",i);
    valname.push_back(vname);
  }

  FILE *dfi_fp;
  if( !(dfi_fp=fopen("field.dfi", "w")) )
  {
    printf("Error open fnaem : %s\n","field.dfi");
    return 0;
  }

  WriteDFI(dfi_fp,vc,valname,endian,dftype,DirPath,FStepList,1);

  FFF.Write(dfi_fp,1);

  fclose(dfi_fp);

  printf("created field.dfi file\n");


  delete [] fubDataMap; 

#ifdef MPI_COMPILE
  MPI_Finalize(); 
#endif 
 
  return 0;

}

// #################################################################
// パラメータの取得

bool getParameter(TextParser *tp, int &vc, DataType &Dtype, double &tol,
                  int &StartFileNo, int &NumberOfBlock,
                  int &StartStepNo, int &EndStepNo, string &DirPath)
{

  string label, str;
  double ct;

  //仮想セル数
  label = "/Parameter/GuideCell";
  if( !(tp->getInspectedValue(label, ct)) )
  {
    printf("Error GuideCell Undefined\n");
    return false;
  } else {
    vc = (int)ct;
  }

  //データタイプ
  label = "/Parameter/DataType";
  if( !(tp->getInspectedValue(label, str)) )
  {
    printf("Error DataType Undefined\n");
    return false;
  } else {
    if( (strcasecmp(str.c_str(),"cell") == 0 ) ) {
      Dtype = CELL;
    } else if( (strcasecmp(str.c_str(),"node") == 0 ) ) {
      Dtype = NODE;
    } else {
      printf("Error Dtype unknown %s\n",str.c_str());
      return false;
    }
  }

  //2点間の距離比較用のトレランス値
  label = "/Parameter/Tolerance";
  if( !(tp->getInspectedValue(label, ct)) ) {
    tol=1.0e-8;
    printf("Warn Tolerance sets uncertain 1.0e-8.\n");
  } else {
    tol = (double)ct;
  }

  //スタートファイル番号
  label = "/FileInfo/StartFileNo";
  if( !(tp->getInspectedValue(label, ct)) ) {
    StartFileNo=1;
    printf("Warn StartFileNo sets uncertain 1.\n");
  } else {
    StartFileNo=(int)ct;
  }

  //ブロック（ファイル）数
  label = "/FileInfo/NumberOfBlock";
  if( !(tp->getInspectedValue(label, ct)) ) {
    printf("Error NumberOfBlock Undefined\n");
    return false;
  } else {
    NumberOfBlock = (int)ct;
  }

  //スタートステップ番号
  label = "/FileInfo/StartStepNo";
  if( !(tp->getInspectedValue(label, ct)) ) {
    StartStepNo=0;
  } else {
    StartStepNo=(int)ct;
  }

  //エンドステップ番号
  label = "/FileInfo/EndStepNo";
  if( !(tp->getInspectedValue(label, ct)) ) {
    EndStepNo=0;
  } else {
    EndStepNo = (int)ct;
  }

  //DirectoryPath
  label = "/FileInfo/DirectoryPath";
  if( !(tp->getInspectedValue(label, str)) ) {
    printf("Error CoordinateFileName Undefined\n");
    return false;
  } else {
    DirPath = str;
  }

  //取得しらパラメータの表示
  printParameter(vc, Dtype, tol, StartFileNo, NumberOfBlock,
                 StartStepNo, EndStepNo, DirPath);
  return true;

}

// #################################################################
// パラメータの表示

void printParameter(int vc, DataType Dtype, double tol, int StartFileNo,
                    int NumberOfBlock, int StartStepNo, int EndStepNo,
                    string DirPath) 
{

  printf("\n");
  printf("-----------------------------------------------\n");
  printf("\n");
  printf("Parameter\n");
  printf("{\n");
  printf("  GuideCell     = %d\n",vc);
  if( Dtype == CELL ) {
    printf("  DataType      = \"Cell\"\n");
  } else if( Dtype == NODE ) {
    printf("  DataType      = \"Node\"\n");
  }
  printf("  Tolerance     = %e\n",tol);
  printf("}\n");
  printf("\n");
  printf("FileInfo\n");
  printf("{\n");
  printf("  StartFileNo   = %d\n",StartFileNo);
  printf("  NumberOfBlock = %d\n",NumberOfBlock);
  printf("  StartStepNo   = %d\n",StartStepNo);
  printf("  EndStepNo     = %d\n",EndStepNo);
  printf("  DirectoryPath   = \"%s\"\n",DirPath.c_str());
  printf("}\n");
  printf("\n");
  printf("\n");

  return;
}

void WriteDFI(FILE *fp, int vc, vector<string> valname,string endian,
              string dftype, string DirPath, set<int>StepNoList,
              unsigned int tab)
{

//FileInfo
  fprintf(fp,"FileInfo {\n");
  fprintf(fp,"\n");

  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"DFIType            = \"Cartesian\"\n");
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"DirectoryPath      = \"%s\"\n",DirPath.c_str());
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"TimeSliceDirectory = \"off\"\n");
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"Prefix             = \"\"\n");
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"FileFormat         = \"fub\"\n");
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"FieldFilenameFormat= \"custom\"\n");
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"GuideCell          = %d\n",vc);
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"DataType           = \"%s\"\n",dftype.c_str());
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"Endian             = \"%s\"\n",endian.c_str());
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"NumVariables       = %d\n",valname.size());
  for(int i=0; i<valname.size(); i++)
  {
    _FUB_WRITE_TAB(fp,tab);
    fprintf(fp,"Variable[@]{ name  = \"%s\" }\n",valname[i].c_str());
  }
  
  fprintf(fp,"\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");

//FilePath
  fprintf(fp,"FilePath {\n");
  fprintf(fp,"\n");

  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"Process = \"proc.dfi\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");

//VisIt
  fprintf(fp,"VisIt {\n");
  fprintf(fp,"\n");

  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"PlotGuideCell = \"off\"\n");
  _FUB_WRITE_TAB(fp,tab);
  fprintf(fp,"ResultFormat  = \"FBinary\"\n");
  fprintf(fp,"\n");
  fprintf(fp,"}\n");  
  fprintf(fp,"\n");

//UnitList
  fprintf(fp,"UnitList {\n");
  fprintf(fp,"\n");
  fprintf(fp,"\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");

//TimeSlice
  fprintf(fp,"TimeSlice {\n");
  fprintf(fp,"\n");
  
  set<int>::iterator it = StepNoList.begin();
  while( it != StepNoList.end() )
  {
    _FUB_WRITE_TAB(fp,tab);
    fprintf(fp,"Slice[@] {\n");
    _FUB_WRITE_TAB(fp,tab+1);
    fprintf(fp,"Step = %u\n",*it);
    _FUB_WRITE_TAB(fp,tab+1);
    fprintf(fp,"Time = %e\n",0.0e0);
    _FUB_WRITE_TAB(fp,tab);
    fprintf(fp,"}\n");
    it++;
  }

  fprintf(fp,"\n");
  fprintf(fp,"}\n");
  fprintf(fp,"\n");

  return;
}
