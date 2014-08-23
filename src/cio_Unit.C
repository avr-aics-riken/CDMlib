/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_Unit.C
 * @brief  cio_Unit Class
 * @author aics    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K



/** cio_UnitElem class **/

// #################################################################
// コンストラクタ
cio_UnitElem::cio_UnitElem()
{

  Name="";
  Unit="";
  reference=0.0;
  difference=0.0;
  BsetDiff=false;

}

// #################################################################
// コンストラクタ
cio_UnitElem::cio_UnitElem(const std::string _Name,
                           const std::string _Unit,
                           const double _reference,
                           const double _difference,
                           const bool _BsetDiff)
{
  Name       = _Name;
  Unit       = _Unit;
  reference  = _reference;
  difference = _difference;
  BsetDiff   = _BsetDiff;
}


// #################################################################
// デストラクタ
cio_UnitElem::~cio_UnitElem()
{


}

// #################################################################
// Unit要素の読込み
CIO::E_CIO_ERRORCODE 
cio_UnitElem::Read(cio_TextParser tpCntl,
                   const std::string label_leaf)
{

  std::string str,label;
  double dt;

  //単位系のの読込み
  label = label_leaf + "/Unit";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    return CIO::E_CIO_WARN_GETUNIT;
  }
  Unit=str;

  //値の読込み
  label = label_leaf + "/Reference";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    dt=0.0;
  }
  reference=dt;

  //diffの読込み
  label = label_leaf + "/Difference";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    difference=0.0;
    BsetDiff=false;
  } else {
    difference=dt;
    BsetDiff=true;
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// Unit要素の出力
CIO::E_CIO_ERRORCODE
cio_UnitElem::Write(FILE* fp, const unsigned tab)
{

  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "Unit       = \"%s\"\n",Unit.c_str());
  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "Reference  = %e\n",reference);
  if( BsetDiff ) {
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "Difference = %e\n",difference);
  }

  return CIO::E_CIO_SUCCESS;

}

/** cio_Uit class **/
// #################################################################
// コンストラクタ
cio_Unit::cio_Unit()
{

}

// #################################################################
// デストラクタ
cio_Unit::~cio_Unit()
{

  UnitList.clear();

}

// #################################################################
// Unitの読込み
CIO::E_CIO_ERRORCODE 
cio_Unit::Read(cio_TextParser tpCntl) 
{

  std::string str;
  std::string label_base,label_leaf;
  int nnode=0;
  CIO::E_CIO_ERRORCODE iret = CIO::E_CIO_SUCCESS;

  //UnitList
  label_base = "/UnitList";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }

  for(int i=0; i<nnode; i++) {
    /** UnitElemの読込み */
    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      //printf("\tCIO Parsing error : No Elem name\n");
      return iret;
    }
    label_leaf=label_base+"/"+str;
    cio_UnitElem unit;
    unit.Name = str;
    if( unit.Read(tpCntl,label_leaf) == CIO::E_CIO_SUCCESS ) {
      UnitList.insert(map<std::string,cio_UnitElem>::value_type(str,unit));
    }
  }

  return iret; 

}

// #################################################################
// 該当するUnitElemの取り出し
CIO::E_CIO_ERRORCODE 
cio_Unit::GetUnitElem(const std::string Name,
                      cio_UnitElem &unit)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合はNULLを返す
  if( it == UnitList.end() ) {
    return CIO::E_CIO_ERROR; 
  }

  //UnitElemを返す
  unit = (*it).second;

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// UnitElemのメンバ変数毎に取得する
CIO::E_CIO_ERRORCODE cio_Unit::GetUnit(const std::string Name,
                                       std::string &unit,
                                       double &ref,
                                       double &diff,
                                       bool &BsetDiff)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は空白を返す
  if( it == UnitList.end() ) {
    return CIO::E_CIO_WARN_GETUNIT;
  }

  //単位を返す
  unit=(*it).second.Unit;
  ref =(*it).second.reference;
  diff=(*it).second.difference;
  BsetDiff=(*it).second.BsetDiff;

  return CIO::E_CIO_SUCCESS;

}
// #################################################################
// ベース名の取り出し
/*
std::string cio_Unit::GetBaseName(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は空白を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return ""; 
  }

  //ベース名を返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.BaseName;

}
*/
// #################################################################
// ベース値の取り出し
/*
double cio_Unit::GetBaseValue(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は0.0を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return 0.0;
  }

  //ベース値を返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.BaseValue;

}
*/
// #################################################################
// Diff Name の取り出し
/*
std::string cio_Unit::GetDiffName(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は空白を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return "";
  }

  //DiffNameを返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.DiffName;
}
*/
// #################################################################
// Diff Valueの取り出し
/*
double cio_Unit::GetDiffValue(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は0.0を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return 0.0;
  }

  //Diff Valueを返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.DiffValue;

}
*/
// #################################################################
// Unitの出力
CIO::E_CIO_ERRORCODE
cio_Unit::Write(FILE* fp, 
                const unsigned tab) 
{

  fprintf(fp, "UnitList {\n");
  fprintf(fp, "\n");

  map<std::string,cio_UnitElem>::iterator it;
  for( it=UnitList.begin(); it!=UnitList.end(); it++ ) {

    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "%s {\n",(*it).second.Name.c_str());

    if( (*it).second.Write(fp,tab+2) != CIO::E_CIO_SUCCESS ) return CIO::E_CIO_ERROR;
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "}\n");
  }

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CIO::E_CIO_SUCCESS;

}

