#ifndef _CDM_TEXTPARSER_H_
#define _CDM_TEXTPARSER_H_

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
 * @file   cdm_TextParser.h
 * @brief  TextParser Control class Header
 * @author aics
 */
#include "cdm_Define.h"

#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "string.h"

#include "TextParser.h"

using namespace std;

class cdm_TextParser {

private:
	TextParser* tp;  ///< テキストパーサ

public:
  /** コンストラクタ */
	cdm_TextParser(){};

  /**　デストラクタ */
	~cdm_TextParser(){};


  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（整数型）
   * @param [in] label 取得するベクトルのラベル（絶対パス）
   * @param [out] vec   ベクトル格納配列ポインタ
   * @param [in]  nvec  ベクトルサイズ
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   */
  bool GetVector(const std::string label, int *vec, const int nvec, bool checkPath=true);


  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（実数型）
   * @param [in]  label  取得するベクトルのラベル（絶対パス）
   * @param [out] vec    ベクトル格納配列ポインタ
   * @param [in]  nvec   ベクトルサイズ
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   */
  bool GetVector(const std::string label, double *vec, const int nvec, bool checkPath=true);


  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（文字列型）
   * @param [in]  label  取得するベクトルのラベル（絶対パス）
   * @param [out] vec    ベクトル格納配列ポインタ
   * @param [in]  nvec   ベクトルサイズ
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   */
  bool GetVector(const std::string label, std::string *vec, const int nvec, bool checkPath=true);


  /**
   * @brief TextParser入力ファイルから変数を取得する（整数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   * @param [out] ct    変数格納ポインタ
   */
  bool GetValue(const std::string label, int *ct, bool checkPath=true);


  /**
   * @brief TextParser入力ファイルから変数を取得する（実数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   * @param [out] ct    変数格納ポインタ
   */
  bool GetValue(const std::string label, double *ct, bool checkPath=true);


  /**
   * @brief TextParser入力ファイルから変数を取得する（文字列型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   * @param [out] ct    変数格納ポインタ
   */
  bool GetValue(const std::string label, std::string *ct, bool checkPath=true);


  /**
   * @brief ラベルの有無をチェック
   * @param [in] label チェックするラベル（絶対パス）
   * @param [in]  checkPath true:ラベル存在チェックをフルパスで行う
   */
  bool chkLabel(const std::string label, bool checkPath=true);


  /**
   * @brief ノードの有無をチェック
   * @param [in] label チェックするノード（絶対パス）
   */
  bool chkNode(const std::string label);


  /**
   * @brief ノード以下のnnode番目の文字列を取得する
   * @param [in]  label ノードの絶対パス
   * @param [in]  nnode 取得する文字列が現れる順番
   * @param [out] ct    取得した文字列
   */
  bool GetNodeStr(const std::string label, const int nnode, std::string *ct);


  /**
   * @brief ノード以下のラベルの数を数える
   * @param [in] label ラベルを数えるノードの絶対パス
   * @retval ラベルの数（エラー、もしくはない場合は-1を返す）
   */
  int countLabels(const std::string label);

  /**
   * @brief TextParserLibraryのインスタンス生成
   * @return エラーコード
   */
  void getTPinstance();


  /**
   * @brief TextParserオブジェクトに入力ファイルをセットする
   * @param [in] filename 入力ファイル名
   * @retval エラーコード
   */
  int readTPfile(const std::string filename);


  /** テキストパーサーの内容を破棄 */
  int remove()
  {
    return tp->remove();
  }

#if 1
  TextParser* getTPPtr()
  {
    return tp;
  }
#endif
};

#endif // _CDM_TXETPARSER_H_
