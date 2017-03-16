#ifndef _CDM_FIELDFILENAMEFORMAT_H_
#define _CDM_FIELDFILENAMEFORMAT_H_

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

#include<map>

#include "cdm_FieldFileNameFormatElem.h"

class cdm_FieldFileNameFormat {

public:

  vector<string> LabelList;

  map<string,cdm_FieldFileNameFormatElem> mapElem;

public:

  /** コンストラクタ */
  cdm_FieldFileNameFormat();

  /** デストラクタ */
  ~cdm_FieldFileNameFormat();

public:

  /** TextParser */
  CDM::E_CDM_ERRORCODE
  //Read(TextParser *tp);
  Read(cdm_TextParser tpCntl);

  /** パラメータの出力 */
  void Print();

  /** FieldFileNameFormatElemクラスの追加 */
  bool AddFieldFileNameFormatElem(cdm_FieldFileNameFormatElem elem);

  /** FieldFileNameFormatElemクラスの取得 */
  cdm_FieldFileNameFormatElem* GetFieldFileNameFormatElem(const string label);

  /** label list の取得 */
  vector<string> GetLabelList();

  /** File有無判定 */
  bool FileExist(string label, string DirPath, int nStep, int nId);

  /** File名生成 */
  string GenerateFileName(string label, string DirPath, int nStep, int nId);

  /** index.dfi FieldFileNameFormat{} 出力 */
  void Write(FILE *fp, const unsigned tab);

};

#endif
