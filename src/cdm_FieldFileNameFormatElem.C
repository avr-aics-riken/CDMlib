/*
 *
 *
 *
 */

#include "cdm_FieldFileNameFormatElem.h"

// #################################################################
// コンストラクタ
cdm_FieldFileNameFormatElem::cdm_FieldFileNameFormatElem(const string label)
{

  FnameLabel = label;

  FileName    = "";
  StepNoKey   = "";
  RankIdKey   = "";
  StepNoDisit = -1;
  RankIdDisit = -1;

  FileNameFormat = "";

  FnameFormat = E_FUB_NONE;

}

// #################################################################
// デストラクタ
cdm_FieldFileNameFormatElem::~cdm_FieldFileNameFormatElem()
{

}

// #################################################################
// TextParserでのParse
bool
cdm_FieldFileNameFormatElem::Read(TextParser *tp)
{

  const string root = "/FieldFileNameFormat";

  string node,leaf,str;
  double ct;

  leaf = root + "/" + FnameLabel;

  // FieldFilenameformat/label/FileNname
  node = leaf + "/FileName";
  if( !(tp->getInspectedValue(node, str)) )
  {
    printf("Error Undefined FileName\n");
    return false;
  } else {
    FileName = str;
  }

  // FieldFilenameformat/label/StepNoKey
  node = leaf + "/StepNoKey";
  if( !(tp->getInspectedValue(node, str)) )
  {
    StepNoKey = "";
  } else {
    StepNoKey = str;
  }

  // FieldFilenameformat/label/StepNoDisit
  node = leaf + "/StepNoDisit";
  int err = 0;
  StepNoDisit = GetDisitNo(tp,node,err);
  if( err != 0 ) return false;

  // FieldFilenameformat/label/RankIdKey
  node = leaf + "/RankIdKey";
  if( !(tp->getInspectedValue(node, str)) )
  {
    RankIdKey = "";
  } else {
    RankIdKey = str;
  }
  // FieldFilenameformat/label/RankIdDisit
  node = leaf + "/RankIdDisit";
  RankIdDisit = GetDisitNo(tp,node,err);

  //FnameFormat(step_rank, rank_step,,,,)のセット,ファイル名生成用のFormatの作成
  SetFnameFormat();

  return true;
}

// #################################################################
// Disitの取得（shorten=0, 整数値)
int 
cdm_FieldFileNameFormatElem::GetDisitNo(TextParser *tp, const string label, int &err)
{
  
  string str;
  double ct;

  if( (tp->getType(label, &err) ) == TP_STRING_VALUE ) {
    tp->getInspectedValue(label, str);
    if( (strcasecmp(str.c_str(),"shorten") == 0 ) ) {
      return 0;
    } else {
      printf("Error %s Parameter Illegal string : %s\n",label.c_str(),str.c_str());
      err = -1;
    }
  } else if( (tp->getType(label, &err) ) == TP_NUMERIC_VALUE ) {
    tp->getInspectedValue(label, ct);
    return (int)ct;
    if( (int)ct < 0 ) {
      printf("Error %s Parameter Illegal number : %d\n",label.c_str(),(int)ct);
      err = -1;
      return (int)ct;
    } else if( (int)ct == 0 ) {
      return -1;
    }
  } else return -1;
  return 0;
}

// #################################################################
// Parse の出力
void
cdm_FieldFileNameFormatElem::PrintParse()
{
  string node,str;
  double ct;

  printf("  %s\n", FnameLabel.c_str());
  printf("  {\n");
  printf("    FileName    = \"%s\"\n", FileName.c_str());
  printf("    StepNoKey   = \"%s\"\n", StepNoKey.c_str());
  if( StepNoDisit == 0 ) {
    printf("    StepNoDisit = \"Shorten\"\n");
  } else if( StepNoDisit < 0 ) {
    printf("    StepNoDisit = 0\n");
  } else if( StepNoDisit > 0 ) {
    printf("    StepNoDisit = %d\n",StepNoDisit);
  }
  printf("    RankIdKey   = \"%s\"\n", RankIdKey.c_str());
  if( RankIdDisit == 0 ) {
    printf("    RankIdDisit = \"Shorten\"\n");
  } else if( RankIdDisit < 0 ) {
    printf("    RankIdDisit = 0\n");
  } else if( RankIdDisit > 0 ) {
    printf("    RankIdDisit = %d\n",RankIdDisit);
  }
  printf("  }\n");

  return;
}

// #################################################################
// Set FileNameFormat
void 
cdm_FieldFileNameFormatElem::SetFnameFormat()
{

//FileNameのStepNoKeyの検索
  char *step = NULL;
  if( !StepNoKey.empty() ) {
    step = strstr((char *)FileName.c_str(), (char *)StepNoKey.c_str());
  }

//FileNameのRankIdKeyの検索
  char *rank = NULL;
  if( !RankIdKey.empty() ) {
    rank = strstr((char *)FileName.c_str(), (char *)RankIdKey.c_str());
  }

//FnameFormat(FielFilenameFormat)の決定
  if( step ) {
    if( rank ) {
      if( step < rank ) {
        FnameFormat = E_FUB_STEP_RANK;  ///< StepNoKeyが先に指示(step_rank)
      } else {
        FnameFormat = E_FUB_RANK_STEP;  ///< RankIdKeyが先に指示(rank_step)
      }
    } else {
      FnameFormat = E_FUB_STEP;         ///< FileNameにStepNoKeyのみ指示
    }
  } else {
    if( rank ) {
      FnameFormat = E_FUB_RANK;         ///<FileNameにRankIdKeyのみ指示
    } else {
      FnameFormat = E_FUB_NONE;         ///<FileNameにStepNoKeyもRankIdKeyも指示なし
    }
  }

//FileNameFormat(ファイル名形式）の作成
  FileNameFormat = FileName;
  //StepNoKeyをsprintfの出力形式に置換（ファイル名生成用）
  if( step ) {
    string StepDisit = SetDisitNoFormat(StepNoDisit);
    string::size_type pos = FileNameFormat.find(StepNoKey);
    if( pos )FileNameFormat.replace(pos,StepNoKey.size(),StepDisit);
  }
  //RankIdKeyをsprintfの出力形式に置換（ファイル名生成用）
  if( rank ) {
    string  RankDisit = SetDisitNoFormat(RankIdDisit);
    string::size_type pos = FileNameFormat.find(RankIdKey);
    FileNameFormat.replace(pos,RankIdKey.size(),RankDisit);
  }

  return;

}

// #################################################################
// 桁数を考慮した出力形式をもとめる
string
cdm_FieldFileNameFormatElem::SetDisitNoFormat(const int DisitNo)
{
 
  string str="";
  char *fmt = new char[10];

  if( DisitNo == 0 ) {
     sprintf(fmt,"%s%s","%","d");
     str = fmt;
  } else if( DisitNo > 0 ) {
     sprintf(fmt,"%s0%d%s","%",DisitNo,"d");
     str = fmt;
  }

  delete [] fmt;

  return str;
}

// #################################################################
// fubファイル名の生成
string
cdm_FieldFileNameFormatElem::GenerateFileName(const string DirPath,
                                          const int nStep, 
                                          const int nId)
{

  string fname = "";

  //DirPsth 
  if( DirPath.empty() ) {
    fname = "./"; 
  } else {
    int sz = DirPath.size();
    string::size_type po = DirPath.find("/");
    if( (sz-1) != (int)po ) {
      fname = DirPath + "/";
    } else {
      fname = DirPath;
    }
  }

  char *name = new char[512];

  //step_rank
  if(        FnameFormat == E_FUB_STEP_RANK ) {
    sprintf(name,FileNameFormat.c_str(),nStep,nId);
  } else if( FnameFormat == E_FUB_RANK_STEP ) {
    sprintf(name,FileNameFormat.c_str(),nId,nStep);
  } else if( FnameFormat == E_FUB_RANK ) {
    sprintf(name,FileNameFormat.c_str(),nId);
  } else if( FnameFormat == E_FUB_STEP ) {
    sprintf(name,FileNameFormat.c_str(),nStep);
  } else {
    sprintf(name,FileNameFormat.c_str());
  }
  
  fname += name;
  delete [] name;  

  return fname;

}

// #################################################################
// ファイル有無判定
bool
cdm_FieldFileNameFormatElem::FileExist(const string DirPath, const int nStep, 
                                       const int nId)
{

  //ファイル名の生成
  string fname = GenerateFileName(DirPath, nStep, nId);

  //ファイルオープン
  FILE *fp;
  if( !(fp = fopen(fname.c_str(),"rb"))  ) return false;
  fclose(fp); 
  return true;

}

// #################################################################
// FieldFileNmaeFormatの出力
bool cdm_FieldFileNameFormatElem::Write(FILE *fp, const unsigned tab)
{

   _FUB_WRITE_TAB(fp,tab);
   fprintf(fp,"%s\n",FnameLabel.c_str());
   _FUB_WRITE_TAB(fp,tab);
   fprintf(fp,"{\n");

   _FUB_WRITE_TAB(fp,tab+1);
   fprintf(fp,"FileName    = \"%s\"\n",FileName.c_str());
   _FUB_WRITE_TAB(fp,tab+1);
   fprintf(fp,"StepNoKey   = \"%s\"\n",StepNoKey.c_str());
   _FUB_WRITE_TAB(fp,tab+1);
   if( StepNoDisit == 0 ) {
      fprintf(fp,"StepNoDisit = \"Shorten\"\n");
   } else if ( StepNoDisit > 0 ) {
      fprintf(fp,"StepNoDisit = %d\n",StepNoDisit);
   } else {
      fprintf(fp,"StepNoDisit = %d\n",0);
   }
   _FUB_WRITE_TAB(fp,tab+1);
   fprintf(fp,"RankIdKey   = \"%s\"\n",RankIdKey.c_str());
   if( RankIdDisit == 0 ) {
      _FUB_WRITE_TAB(fp,tab+1);
      fprintf(fp,"RankIdDisit = \"Shorten\"\n");
   } else if( RankIdDisit > 0 ) {
      _FUB_WRITE_TAB(fp,tab+1);
      fprintf(fp,"RankIdDisit = %d\n",RankIdDisit);
   } else {
      _FUB_WRITE_TAB(fp,tab+1);
      fprintf(fp,"RankIdDisit = %d\n",0);
   }

   _FUB_WRITE_TAB(fp,tab);
   fprintf(fp,"}\n");

  return true;

}
