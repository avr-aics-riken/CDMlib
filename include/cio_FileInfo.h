#ifndef _CDM_FILEINFO_H_
#define _CDM_FILEINFO_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
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
  int                   Component;        ///<成分数
  vector<std::string>   ComponentVariable;///<成分名

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
   * @param [in] _Component        成分数
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
               const int _Component);

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
   * @brief 成分名をセットする
   * @param [in] pcomp    成分位置 0:u, 1:v, 2:w
   * @param [in] compName 成分名 "u","v","w",,,
   */
  void setComponentVariable(int pcomp, std::string compName); 

  /**
   * @brief 成分名を取得する
   * @param [in] pcomp    成分位置 0:u, 1:v, 2:w
   * @return 成分名　成分名が無い場合は空白が返される
   */
  std::string getComponentVariable(int pcomp); 

};

#endif // _cdm_FILEINFO_H_
