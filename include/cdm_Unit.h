#ifndef _CDM_UNIT_H_
#define _CDM_UNIT_H_

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
 * @file   cdm_Unit.h
 * @brief  cdm_UnitElem & cdm_Unit Class Header
 * @author aics
 */

class cdm_UnitElem {

public:

  std::string Name;         ///<単位の種類名(Length,Velovity,,,)
  std::string Unit;         ///<単位のラベル(m,m/s,Pa,,,,)
  double      reference;    ///<規格化に用いたスケール
  double      difference;   ///<差
  bool        BsetDiff;     ///<differenceの有無（false:なし true:あり）
  //std::string BaseName;     ///<代表値のラベル(L0,P0,,,,)
  //double      BaseValue;    ///<代表値
  //std::string DiffName;     ///<差があるときの名前、空の時DiffValue未定義
  //double      DiffValue;    ///<差の値

  /** コンストラクタ */
  cdm_UnitElem();

  /** コンストラクタ */
  cdm_UnitElem(const std::string _Name,
               const std::string _Unit,
               const double _reference,
               const double _difference,
               const bool _BsetDiff);

  /** デストラクタ */
  ~cdm_UnitElem();

  /**
   * @brief Unit要素の読込み
   * @param [in]  tpCntl      cdm_TextParserクラス
   * @param [in]  label_leaf
   * @return error code
   */
   CDM::E_CDM_ERRORCODE
   Read(cdm_TextParser tpCntl,
        const std::string label_leaf);

  /**
   * @brief DFIファイル:Unit要素を出力する
   * @param [in]  fp      ファイルポインタ
   * @param [in]  tab     インデント
   * @return error code
   */
   CDM::E_CDM_ERRORCODE
   Write(FILE* fp, const unsigned tab);

};


/** index.dfi ファイルの Unit */
class cdm_Unit {

public:

  map<std::string,cdm_UnitElem> UnitList;

  /** コンストラクタ **/
  cdm_Unit();

  /** デストラクタ **/
  ~cdm_Unit();

  /**
   * @brief read Unit(inde.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl);

  /**
   * @brief 該当するUnitElemの取り出し
   * @param [in]  Name 取り出す単位の種類
   * @param [out] unit 取得したcdm_UnitElemクラス
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  GetUnitElem(const std::string Name,
              cdm_UnitElem &unit);

  /**
   * @brief 単位の取り出し("m","cm",,,,,)
   * @param [in]  Name 取り出す単位の種類
   * @param [out] unit 単位文字列
   * @param [out] ref  reference
   * @param [out] diff difference
   * @param [out] bSetDiff difference有無フラグ true:あり、false:なし
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  GetUnit(const std::string Name,
          std::string &unit,
          double &ref,
          double &diff,
          bool &bSetDiff);

  /**
   * @brief DFIファイル:Unit要素を出力する
   * @param [in]  fp      ファイルポインタ
   * @param [in]  tab     インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, const unsigned tab);

  /**
   * @brief Uuitをセットする
   * @param [in] Name       追加する単位系("Length","Velocity",,,,)
   * @param [in] Unit       単位ラベル("M","CM","MM","M/S",,,)
   * @param [in] reference  規格化したスケール値
   * @param [in] difference 差の値
   * @param [in] BsetDiff   differenceの有無
   */
  void
  AddUnit(const std::string Name,
          const std::string Unit,
          const double reference,
          const double difference= 0.0,
          const bool BsetDiff=false);
};

#endif // _CDM_UNIT_H_
