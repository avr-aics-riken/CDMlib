/*
 *
 *
 */

#include "cdm_DFI.h"
//#include "cdm_FieldFileNameFormat.h"

// #################################################################
// コンストラクタ
cdm_FieldFileNameFormat::cdm_FieldFileNameFormat()
{

  LabelList.clear(); 

}

// #################################################################
// デストラクタ
cdm_FieldFileNameFormat::~cdm_FieldFileNameFormat()
{

}

// #################################################################
// TextParserでのParse
CDM::E_CDM_ERRORCODE
//cdm_FieldFileNameFormat::Read(TextParser *tp)
cdm_FieldFileNameFormat::Read(cdm_TextParser tpCntl)
{

  //TP
  TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    return CDM::E_CDM_ERROR_TEXTPARSER;
  }  

  vector<string> top_label;  ///< カレントの子ノードのラベル格納
  tp->getNodes(top_label);
  //FieldFileNameFormatタグがある
  int cnt=0;
  for(int i=0; i< top_label.size(); i++)
  {
    if( (strcasecmp(top_label[i].c_str(),"FieldFileNameFormat") == 0 ) ) cnt ++;
  }

  if( cnt<1 ) return CDM::E_CDM_ERROR_UNDEFINED_FIELDFILENAMEFORMAT;

  vector<string> label;  ///< FieldFileNameFormat子ノードのラベル格納

  //FieldFileNameFormatタグへの移動
  tp->changeNode("FieldFileNameFormat");

  //FieldFileNameFormatの子ノードを取得
  tp->getNodes(label,2);

  //FieldFileNameFormatの子ノード数のループ
  for(int i=0; i<label.size(); i++)
  {
    //FieldFileNameFormatElemのインスタンス
    cdm_FieldFileNameFormatElem elem(label[i]);

    //FieldFileNameFormatElemのパース
    if( elem.Read(tp) ) 
    {
      //FieldFileNameFormatElemをmapに追加
      AddFieldFileNameFormatElem(elem);
    } else {
      return CDM::E_CDM_ERROR;
    }
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// パラメータの出力
void
cdm_FieldFileNameFormat::Print()
{

  for(int i=0; i<LabelList.size(); i++)
  {
    GetFieldFileNameFormatElem(LabelList[i])->PrintParse();
  }
  return;

}

// #################################################################
// FieldFileNameFormatElemクラスの取得
cdm_FieldFileNameFormatElem* 
cdm_FieldFileNameFormat::GetFieldFileNameFormatElem( const string label )
{

  map<string,cdm_FieldFileNameFormatElem>::iterator it;
  it = mapElem.find(label);
  if( it == mapElem.end() ) return NULL;
  return &it->second;

}

// #################################################################
// FieldFileNameFormatElemクラスの追加
bool
cdm_FieldFileNameFormat::AddFieldFileNameFormatElem(cdm_FieldFileNameFormatElem elem)
{
  //mapに既に登録されているときはメッセージ出力し、falseを返す
  if( GetFieldFileNameFormatElem(elem.FnameLabel) )
  {
    printf("Duplicate FieldFileNameFormatElem label %s\n",elem.FnameLabel.c_str());
    return false;
  //mapに未登録のときelemを登録し、trueを返す
  } else {
    LabelList.push_back(elem.FnameLabel);
    mapElem.insert( map<string,cdm_FieldFileNameFormatElem>::value_type(elem.FnameLabel,elem) );
    return true;
  }
}

// #################################################################
// label list の取得
vector<string> 
cdm_FieldFileNameFormat::GetLabelList()
{

  return LabelList;

}

// #################################################################
// ファイルの存在を判定
bool
cdm_FieldFileNameFormat::FileExist(string label, string DirPath, int nStep, int nId)
{
  return GetFieldFileNameFormatElem(label)->FileExist(DirPath,nStep,nId);
}

// #################################################################
// ファイル名の生成
string
cdm_FieldFileNameFormat::GenerateFileName(string label, string DirPath, 
                                      int nStep, int nId)
{
   cdm_FieldFileNameFormatElem* elm = GetFieldFileNameFormatElem(label);
   if( elm == NULL ) return "";
   return elm->GenerateFileName(DirPath,nStep, nId);
}

// #################################################################
// index.dfi FieldFileNameFormat出力
void
cdm_FieldFileNameFormat::Write(FILE *fp, const unsigned tab)
{

  fprintf(fp,"FieldFileNameFormat\n");
  fprintf(fp,"{\n");

  for(int i=0; i<LabelList.size(); i++) {
    GetFieldFileNameFormatElem(LabelList[i])->Write(fp, tab);
  }

  fprintf(fp,"}\n");
}
