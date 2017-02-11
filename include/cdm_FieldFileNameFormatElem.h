#ifndef _CDM_FIELDFILENAMEFORMATELEM_H_
#define _CDM_FIELDFILENAMEFORMATELEM_H_

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

#include <stdio.h>
#include <string.h>
#include "TextParser.h"

#define _FUB_TAB_STR "  "
#define _FUB_WRITE_TAB(_FP,_NTAB) {\
  for(int _NTCNT=0; _NTCNT<_NTAB; _NTCNT++) fprintf(_FP,_FUB_TAB_STR); \
}
using namespace std;

enum FUB_FNAME_TYPE
{
  E_FUB_NONE      = -1, //stepもrankもなし
  E_FUB_STEP_RANK =  0, //step_rank
  E_FUB_RANK_STEP =  1, //rank_step
  E_FUB_RANK      =  2, //rank
  E_FUB_STEP      =  3
};

class cdm_FieldFileNameFormatElem {

protected:

  string FileName;  ///< 読み込まれるフィールドファイル名の規約名
                    ///< （例）"cellPhys_STEPNO_idFILENO.fub"

  string StepNoKey; ///< ステップ番号が付与されるキー文字　（例）"STEPNO"
  string RankIdKey; ///< ランク番号が付与されるキー文字    （例）"FILENO"

  int StepNoDisit;  ///< ステップ番号桁数（正:桁数, 0:shorten, 負:ステップ番号なし)
  int RankIdDisit;  ///< ランク番号桁数（正:桁数, 0:shorten, 負:ランク番号なし)

  string FileNameFormat;  ///< ファイル名生成用のFormat
                          ///<（例）"cellPhy_%06d_id%10d.fub"

  FUB_FNAME_TYPE FnameFormat;  ///< FieldFilenameFormat(enum FUB_FNAME_TYPE)

public:

  string FnameLabel;      /// FieldFileNameFormatのファイル名（FieldFile,CoordinateFile)

public:

  /** コンストラクタ */
  cdm_FieldFileNameFormatElem(const string label);

  /** デストラクタ */
  ~cdm_FieldFileNameFormatElem();

protected:

  /** FileNameFormatをセットする */
  void SetFnameFormat();

  /** 桁数を考慮した出力形式をセットする */
  string SetDisitNoFormat(const int DisitNo);

public:

  /** TextParserでのParsea */
  bool Read(TextParser *tp);

  /** Parseの出力 */
  void PrintParse();

  /** Disitの取得（shorten=0, 整数値) */
  int GetDisitNo(TextParser *tp, const string label, int &err);

  /** FieldFilenameFormatのXXXXFile出力 */
  bool Write(FILE *fp, const unsigned tab);

  /** ファイル名の生成 */
  string GenerateFileName(const string DirPath,
                          const int nStep,
                          const int nId);

  /** ファイル有無判定 */
  bool FileExist(const string DirPath, const int nStep, const int nId);

};

#endif
