#ifndef _CDM_FILEINFO_H_
#define _CDM_FILEINFO_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_FileInfo.h
 * @brief  cdm_FileInfo Class Header
 * @author aics    
 */


/** index.dfi ファイルの FileInfo */
class cdm_FileInfo {

public:

  //FCONV 20140116.s
  CDM::E_CDM_DFITYPE    DFIType;          ///<dfi種別
  CDM::E_CDM_OUTPUT_FNAME FieldFilenameFormat; ///<ファイル命名基準
  //FCONV 20140116.e
 
  std::string           DirectoryPath;    ///<フィールデータの存在するディレクトリパス
                                          ///< index.dfiからの相対パスまたは絶対パス
  CDM::E_CDM_ONOFF      TimeSliceDirFlag; ///<TimeSlice on or off
  std::string           Prefix;           ///<ファイル接頭文字
  CDM::E_CDM_FORMAT     FileFormat;       ///<ファイルフォーマット "bov","sph",,,
  int                   GuideCell;        ///<仮想セルの数
  CDM::E_CDM_DTYPE      DataType;         ///<配列のデータタイプ "float",,,,
  CDM::E_CDM_ENDIANTYPE Endian;           ///<エンディアンタイプ "big","little"
  CDM::E_CDM_ARRAYSHAPE ArrayShape;       ///<配列形状
  int                   NumVariables;     ///<変数の個数
  vector<std::string>   VariableName;     ///<変数名

  /** コンストラクタ **/
  cdm_FileInfo();

  /** 
   * @brief コンストラクタ 
   * @param [in] _DFIType          dfi種別
   * @param [in] _FieldFilenameFormat ファイル命名基準
   * @param [in] _DirectoryPath    ディレクトリパス
   * @param [in] _TimeSliceDirFlag TimeSlice on or off
   * @param [in] _Prefix           ファイル接頭文字
   * @param [in] _FileFormat       ファイルフォーマット
   * @param [in] _GuideCell        仮想セルの数
   * @param [in] _DataType         配列のデータタイプ
   * @param [in] _Endian           エンディアンタイプ
   * @param [in] _ArrayShape       配列形状
   * @param [in] _NumVariables     変数の個数
   */
  cdm_FileInfo(const CDM::E_CDM_DFITYPE _DFIType,
               const CDM::E_CDM_OUTPUT_FNAME _FieldFilenameFormat,
               const std::string _DirectoryPath, 
               const CDM::E_CDM_ONOFF _TimeSliceDirFlag, 
               const std::string _Prefix, 
               const CDM::E_CDM_FORMAT _FileFormat,
               const int _GuideCell, 
               const CDM::E_CDM_DTYPE _DataType, 
               const CDM::E_CDM_ENDIANTYPE _Endian, 
               const CDM::E_CDM_ARRAYSHAPE _ArrayShape, 
               const int _NumVariables);

  /** デストラクタ **/
  ~cdm_FileInfo();

public:

  /**
   * @brief read FileInfo(inde.dfi)
   * @param [in]   tpCntl  cdm_TextParserクラス 
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Read(cdm_TextParser tpCntl); 

  /**
   * @brief DFIファイル:FileInfo要素を出力する
   * @param [in]  fp      ファイルポインタ
   * @param [in]  tab     インデント
   * @return error code
   */
  CDM::E_CDM_ERRORCODE
  Write(FILE* fp, 
        const unsigned tab); 

  /**
   * @brief 変数名をセットする
   * @param [in] pvari    変数位置 0:u, 1:v, 2:w
   * @param [in] variName 変数名 "u","v","w",,,
   */
  void setVariableName(int pvari, std::string variName); 

  /**
   * @brief 変数名を取得する
   * @param [in] pvari    変数位置 0:u, 1:v, 2:w
   * @return 変数名　変数名が無い場合は空白が返される
   */
  std::string getVariableName(int pvari); 

};

#endif // _cdm_FILEINFO_H_
