#ifndef _CIO_FILEINFO_H_
#define _CIO_FILEINFO_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_FileInfo.h
 * @brief  cio_FileInfo Class Header
 * @author aics    
 */


/** index.dfi ファイルの FileInfo */
class cio_FileInfo {

public:

  //FCONV 20140116.s
  CIO::E_CIO_DFITYPE    DFIType;          ///<dfi種別
  CIO::E_CIO_OUTPUT_FNAME FieldFilenameFormat; ///<ファイル命名基準
  //FCONV 20140116.e
 
  std::string           DirectoryPath;    ///<フィールデータの存在するディレクトリパス
                                          ///< index.dfiからの相対パスまたは絶対パス
  CIO::E_CIO_ONOFF      TimeSliceDirFlag; ///<TimeSlice on or off
  std::string           Prefix;           ///<ファイル接頭文字
  CIO::E_CIO_FORMAT     FileFormat;       ///<ファイルフォーマット "bov","sph",,,
  int                   GuideCell;        ///<仮想セルの数
  CIO::E_CIO_DTYPE      DataType;         ///<配列のデータタイプ "float",,,,
  CIO::E_CIO_ENDIANTYPE Endian;           ///<エンディアンタイプ "big","little"
  CIO::E_CIO_ARRAYSHAPE ArrayShape;       ///<配列形状
  int                   Component;        ///<成分数
  vector<std::string>   ComponentVariable;///<成分名

  /** コンストラクタ **/
  cio_FileInfo();

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
  cio_FileInfo(const CIO::E_CIO_DFITYPE _DFIType,
               const CIO::E_CIO_OUTPUT_FNAME _FieldFilenameFormat,
               const std::string _DirectoryPath, 
               const CIO::E_CIO_ONOFF _TimeSliceDirFlag, 
               const std::string _Prefix, 
               const CIO::E_CIO_FORMAT _FileFormat,
               const int _GuideCell, 
               const CIO::E_CIO_DTYPE _DataType, 
               const CIO::E_CIO_ENDIANTYPE _Endian, 
               const CIO::E_CIO_ARRAYSHAPE _ArrayShape, 
               const int _Component);

  /** デストラクタ **/
  ~cio_FileInfo();

public:

  /**
   * @brief read FileInfo(inde.dfi)
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  Read(cio_TextParser tpCntl); 

  /**
   * @brief DFIファイル:FileInfo要素を出力する
   * @param [in]  fp      ファイルポインタ
   * @param [in]  tab     インデント
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
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

#endif // _cio_FILEINFO_H_
