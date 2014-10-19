/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   cdm_TextParser.C
 * @brief  TextParser Control class
 * @author aics
 */

#include "cdm_TextParser.h"



// #################################################################
// ラベルの有無をチェック
bool cdm_TextParser::chkLabel(const std::string label, bool fullPath)
{
  int ierror;
  std::string value;
  
  if( !tp ) return false;
  
  // ラベルがあるかチェック
  vector<std::string> labels;
  if( fullPath )
  {
    ierror=tp->getAllLabels(labels);
    if (ierror != 0)
    {
      cout <<  "ERROR in TextParser::getAllLabels file: "
      << " ERROR CODE "<< ierror << endl;
      return false;
    }
  }
  else
  {
    ierror=tp->getLabels(labels);
    if (ierror != 0)
    {
      cout <<  "ERROR in TextParser::getLabels file: "
      << " ERROR CODE "<< ierror << endl;
      return false;
    }
  }
  
  int flag=0;
  for (int i = 0; i < labels.size(); i++)
  {
	  if( !strcasecmp(label.c_str(), labels[i].c_str()) )
    {
		  flag=1;
		  break;
	  }
  }
  
  if (flag==0)
  {
	  return false;
  }
  
  return true;
}



// #################################################################
// ノードの有無をチェック
bool cdm_TextParser::chkNode(const std::string label)
{
  int ierror;
  std::string node;
  vector<std::string> labels;
  int len=label.length();
  
  if( !tp ) return false;
  
  // Nodeがあるかチェック
  ierror = tp->getAllLabels(labels);
  
  if (ierror != 0)
  {
    cout <<  "ERROR in TextParser::getAllLabels file: "
    << " ERROR CODE "<< ierror << endl;
    return false;
  }
  
  int flag=0;
  for (int i = 0; i < labels.size(); i++) {
	  node = labels[i].substr(0,len);
    
    if ( !strcasecmp(node.c_str(), label.c_str()) )
    {
		  flag=1;
		  break;
	  }
  }
  
  if (flag==0)
  {
	  return false;
  }
  
  return true;
}



// #################################################################
// ノード以下のラベルの数を数える
int cdm_TextParser::countLabels(const std::string label)
{
  int ierror;
  std::string node,str,chkstr="";
  vector<std::string> labels;
  int len=label.length();
  int flag=0;
  int inode=0;
  int next=0;
  
  if( !tp ) return -1;
  
  // Nodeがあるかチェック
  ierror=tp->getAllLabels(labels);
  
  if (ierror != 0){
    cout <<  "ERROR in TextParser::getAllLabels file: "
    << " ERROR CODE "<< ierror << endl;
    return -1;
  }
  
  for (int i = 0; i < labels.size(); i++) {
	  node=labels[i].substr(0,len);
    
	  if( !strcasecmp(node.c_str(), label.c_str()) ){
		  str=labels[i].substr(len+1);
		  next=str.find("/");
      
		  if(next==0) inode++;
		  else{
			  if(chkstr!=str.substr(0,next)){
				  chkstr=str.substr(0,next);
				  inode++;
			  }
		  }
      
	  }
  }
  
  return inode;
}



// #################################################################
// TextParserLibraryのインスタンス生成
void cdm_TextParser::getTPinstance()
{
  tp = new TextParser;
}



// #################################################################
// ノード以下のnnode番目の文字列を取得する
bool cdm_TextParser::GetNodeStr(const std::string label, const int nnode, std::string *ct)
{
  if ( !tp ) return -1;
  
  int ierror;
  int len=label.length();
  int flag=0;
  int inode=0;
  int next=0;
  
  std::string node;
  std::string str;
  std::string chkstr="";
  vector<std::string> labels;

  
  // Nodeがあるかチェック
  ierror = tp->getAllLabels(labels);
  
  if (ierror != 0)
  {
    cout <<  "ERROR in TextParser::getAllLabels file: " << " ERROR CODE "<< ierror << endl;
    return false;
  }
  
  for (int i = 0; i < labels.size(); i++) {
	  node = labels[i].substr(0, len);
    
	  if ( !strcasecmp(node.c_str(), label.c_str()) )
    {
		  str = labels[i].substr(len+1);
		  next = str.find("/");
      
		  if ( next == 0 )
      {
			  inode++;
		  }
		  else
      {
			  if ( chkstr != str.substr(0, next) )
        {
				  chkstr = str.substr(0, next);
				  inode++;
			  }
		  }
      
		  if ( inode == nnode )
      {
			  *ct = chkstr;
			  return true;
		  }
	  }
  }
  return false;
}



// #################################################################
// TextParser入力ファイルからベクトル値を取得する（整数型）
bool cdm_TextParser::GetVector(const std::string label, int *vec, const int nvec, bool checkPath)
{
  int ierr = TP_NO_ERROR;
  std::string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label,checkPath)){
	  return false;
  }

  // get value
  if( (ierr = tp->getValue(label, value)) != TP_NO_ERROR ) return false;

  // get type
  TextParserValueType type = tp->getType(label, &ierr);
  if( ierr != TP_NO_ERROR ) return false;
  if( type != TP_VECTOR_NUMERIC ) return false;

  // split
  vector<std::string> vec_value;
  if( (ierr = tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return false;

  // check number of vector element
  if( vec_value.size() != nvec ) return false;

  // string to real
  for(int i=0; i<vec_value.size(); i++ )
  {
    vec[i] = tp->convertInt(vec_value[i], &ierr);
    if( ierr != TP_NO_ERROR ) return false;
  }

  return true;
}


// #################################################################
// TextParser入力ファイルからベクトル値を取得する（実数型）
bool cdm_TextParser::GetVector(const std::string label, double *vec, const int nvec, bool checkPath)
{
  int ierr = TP_NO_ERROR;
  std::string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label,checkPath)){
	  return false;
  }

  // get value
  if( (ierr = tp->getValue(label, value)) != TP_NO_ERROR ){
    cout << " GetVector debug 333" << endl;
	  return false;
  }

  // get type
  TextParserValueType type = tp->getType(label, &ierr);
  if( ierr != TP_NO_ERROR ) return false;
  if( type != TP_VECTOR_NUMERIC ) return false;

  // split
  vector<std::string> vec_value;
  if( (ierr = tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return false;

  // check number of vector element
  if( vec_value.size() != nvec ) return false;

  // string to real
  for(int i=0; i<vec_value.size(); i++ )
  {
    vec[i] = tp->convertDouble(vec_value[i], &ierr);
    if( ierr != TP_NO_ERROR ) return false;
  }

  return true;
}


// #################################################################
// TextParser入力ファイルからベクトル値を取得する（文字列型）
bool cdm_TextParser::GetVector(std::string label, std::string *vec, const int nvec, bool checkPath)
{
  int ierr = TP_NO_ERROR;
  std::string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label,checkPath)){
	  return false;
  }

  // get value
  if( (ierr = tp->getValue(label, value)) != TP_NO_ERROR ) return false;

  // get type
  TextParserValueType type = tp->getType(label, &ierr);
  if( ierr != TP_NO_ERROR ) return false;
  if( type != TP_VECTOR_NUMERIC ) return false;

  // split
  vector<std::string> vec_value;
  if( (ierr = tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return false;

  // check number of vector element
  if( vec_value.size() != nvec ) return false;

  // string to string
  for(int i=0;i<vec_value.size();i++ )
  {
    vec[i] = vec_value[i];
  }

  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（整数型）
bool cdm_TextParser::GetValue(const std::string label, int *ct, bool checkPath)
{
  int ierror;
  std::string value;
  
  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label,checkPath)){
    return false;
  }

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << label << ierror << endl;
	  return false;
  }

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  return false;
  }
  if( type != TP_NUMERIC_VALUE ){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  return false;
  }

  // string to real
  int val = tp->convertInt(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR convertInt " << ierror << endl;
	  return false;
  }

  *ct=val;

  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（実数型）
bool cdm_TextParser::GetValue(const std::string label, double *ct, bool checkPath)
{
  int ierror;
  std::string value;
  std::string node;

  if( !tp ) return false;
  
  // ラベルがあるかチェック
  if( !chkLabel(label,checkPath)){
    return false;
  }

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << ierror << endl;
	  return false;
  }

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  return false;
  }
  if( type != TP_NUMERIC_VALUE ){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  return false;
  }

  // string to real
  //REAL_TYPE val = tp->convertFloat(value, &ierror);
  double val = tp->convertFloat(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR convertInt " << ierror << endl;
	  return false;
  }

  *ct=val;

  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（文字列型）
bool cdm_TextParser::GetValue(const std::string label, std::string *ct, bool checkPath)
{
  int ierror;
  std::string value;
  
  if ( !tp ) return false;
  
  // ラベルがあるかチェック
  if ( !chkLabel(label,checkPath) )
  {
    return false;
  }
  
  //値の取得
  ierror = tp->getValue(label, value); //labelは絶対パスを想定
  
  if (ierror != TP_NO_ERROR)
  {
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << label << endl;
	  return false;
  }
  
  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR)
  {
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  return false;
  }
  
  if( type != TP_STRING_VALUE )
  {
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  return false;
  }
  
  *ct=value;
  
  return true;
}



// #################################################################
// TextParserオブジェクトに入力ファイルをセットする
int cdm_TextParser::readTPfile(const std::string filename)
{
  int ierr = TP_NO_ERROR;
  if( !tp ) return TP_ERROR;
  
  // read
  if( (ierr = tp->read(filename)) != TP_NO_ERROR )
  {
    cout << "ERROR : in input file: " << filename << endl
    << "  ERROR CODE = "<< ierr << endl;
    return ierr;
  }
  return ierr;
}
