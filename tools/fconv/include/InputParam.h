#ifndef _INPUTPARAM_H_
#define _INPUTPARAM_H_

/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   InputParam.h
 * @brief  InputParam Class Header
 * @author aics
 * @date   2013/11/13
 */

#include "cpm_ParaManager.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <errno.h>

#ifndef _WIN32
#include <dirent.h>
#else
#include "sph_win32_util.h"   // for win32
#endif

#include "limits.h" // for UBUNTU

#include "conv_Define.h"
#include "TextParser.h"
#include "cdm_DFI.h"

using namespace std;

class InputParam {

public:

  //dfi情報構造体
  struct dfi_info{
    std::string in_dfi_name;    ///<読込みdfiファイル名
    cdm_DFI* in_dfi;            ///<読込みdfiポインタ
    std::string out_dfi_name;   ///<出力dfiファイル名
    std::string out_proc_name;  ///<出力procファイル名
  };

  vector<dfi_info> m_dfiList;  ///<入出力するDFI名等が格納

  cpm_ParaManager* m_paraMngr;

protected:
  
  bool m_output_dfi_on;      ///<dfi出力指示フラグ
  bool m_bgrid_interp_flag;  ///<節点への補間フラグ
  bool m_out_ftype_crd_on;   ///<座標データの出力形式指示フラグ
  bool m_cropIndexStart_on;  ///<入力領域のスタート指示フラグ
  bool m_cropIndexEnd_on;    ///<入力領域のエンド指示フラグ
  int m_cropIndexStart[3];   ///<入力領域のスタート位置（2014対応予定)
  int m_cropIndexEnd[3];     ///<入力領域のエンド位置（2014対応予定)
  int m_outputDiv[3];        ///<出力分割数 MxNで有効
  int m_thin_count;          ///<間引き数
  int m_outputGuideCell;     ///<出力するガイドセル数
  std::string m_outdir_name; ///<出力先ディレクトリー名

  CDM::E_CDM_DTYPE              m_output_data_type;     ///<出力実数タイプ byte,short,int,float,double
  CDM::E_CDM_FORMAT             m_out_format;           ///<出力ファイルフォーマット sph,bov,avs,plot3d,vtk
  CDM::E_CDM_FILE_TYPE          m_out_file_type;        ///<出力形式 ascii,binary,FortranBinary
  CDM::E_CDM_FILE_TYPE          m_out_file_type_coord;  ///<座標データの出力形式(AVS形式) ascii,binary
  E_CONV_OUTPUT_CONV_TYPE       m_conv_type;            ///<convertタイプ Mx1 MxN MxM
  CDM::E_CDM_ARRAYSHAPE         m_outputArrayShape;     ///<出力配列形状
  CDM::E_CDM_OUTPUT_FNAME       m_outputFilenameFormat; ///<出力ファイル名命名順
  E_CONV_OUTPUT_MULTI_FILE_CAST m_multiFileCasting;     ///<並列処理時のファイル割振り方法
  
public:
  /** コンストラクタ */
  InputParam(cpm_ParaManager* paraMngr);
  
  /**　デストラクタ */
  ~InputParam();
  
public:

  /**
   * @brief CPMのポインタをコピー
   * @param [in] paraMngr  cpm_ParaManagerクラス 
   * @return  エラーコード
   */
  //bool importCPM(cpm_ParaManager* paraMngr);

  /**
   * @brief 入力ファイルの読み込み
   * @param [in] input_file_name 入力TPファイル名
   */
  bool Read(std::string input_file_name); 

  /**
   * @brief 入力パラメータのエラーチェック
   * @return error code
   */
  bool InputParamCheck(); 

  /**
   * @brief 入力パラメータのログ出力
   * @param[in] fp 出力ファイルポインタ
   */
  void PrintParam(FILE* fp); 

  /**
   * @brief dfiListの取り出し
   */
  vector<dfi_info> Get_dfiList() { return m_dfiList; }; 

  /**
   * @brief dfi出力指示フラグの取り出し
   */
  bool Get_Outputdfi_on() { return m_output_dfi_on; };
 
  /**
   * @brief 節点への補間フラグの取り出し
   */
  bool Get_Interp_flag() { return m_bgrid_interp_flag; };

  /**
   * @brief input dfiファイル名リストの取り出し
   */
  //vector<std::string> Get_IndfiNameList() { return m_in_dfi_name; };

  /**
   * @brief output dfiファイル名リストの取り出し
   */
  //vector<std::string> Get_OutdfiNameList() { return m_out_dfi_name; };

  /**
   * @brief output procファイル名リストの取り出し
   */
  //vector<std::string> Get_OutprocNameList() { return m_out_proc_name; };

  /**
   * @brief output procファイル名リストをセットする
   * @param[in] out_proc_name セットするprocファイル名リスト
   */
  /*
  bool Set_OutprocNameList(vector<std::string> out_proc_name)
  { if( out_proc_name.size() < 1 ) return false;
    for(int i=0; i<out_proc_name.size(); i++) {
      m_out_proc_name.push_back(out_proc_name[i]);
    }
    return true;
  }
  */
  /**
   * @brief コンバートタイプの取り出し
   */
  E_CONV_OUTPUT_CONV_TYPE Get_ConvType() { return m_conv_type; };

  /**
   * @brief 出力分割数の取り出し
   */
  int* Get_OutputDivision() { return m_outputDiv; };

  /**
   * @brief 出力ファイルフォーマットの取り出し
   */
  CDM::E_CDM_FORMAT Get_OutputFormat() { return m_out_format; };

  /**
   * 出力ファイルフォーマットの取り出し(文字列）
   */
  std::string Get_OutputFormat_string()
  {
    if( Get_OutputFormat() == CDM::E_CDM_FMT_SPH ) return "sph";
    if( Get_OutputFormat() == CDM::E_CDM_FMT_BOV ) return "bov";
    if( Get_OutputFormat() == CDM::E_CDM_FMT_AVS ) return "avs";
    if( Get_OutputFormat() == CDM::E_CDM_FMT_VTK ) return "vtk";
    if( Get_OutputFormat() == CDM::E_CDM_FMT_PLOT3D ) return "plot3d";
    return "";
  }; 


  /**
   * @brief 出力タイプの取り出し
   */
  CDM::E_CDM_DTYPE
  Get_OutputDataType() { return m_output_data_type; };

  /**
   * @brief 出力タイプの取り出し(文字列)
   */
  std::string Get_OutputDataType_string()
  {
    if( Get_OutputDataType() == CDM::E_CDM_INT8    ) return "Int8";
    if( Get_OutputDataType() == CDM::E_CDM_UINT8   ) return "UInt8";
    if( Get_OutputDataType() == CDM::E_CDM_INT16   ) return "Int16";
    if( Get_OutputDataType() == CDM::E_CDM_UINT16  ) return "UInt16";
    if( Get_OutputDataType() == CDM::E_CDM_INT32   ) return "Int32";
    if( Get_OutputDataType() == CDM::E_CDM_UINT32  ) return "UInt32";
    if( Get_OutputDataType() == CDM::E_CDM_INT64   ) return "Int64";
    if( Get_OutputDataType() == CDM::E_CDM_UINT64  ) return "UInt64";
    if( Get_OutputDataType() == CDM::E_CDM_FLOAT32 ) return "Float32";
    if( Get_OutputDataType() == CDM::E_CDM_FLOAT64 ) return "Float64";
    return "";
  };

  /**
   * @brief 出力形式の取り出し（ascii,binary,FortranBinary)
   */
  CDM::E_CDM_FILE_TYPE
  Get_OutputFileType() { return m_out_file_type; };

  /**
   * @brief 座標データの出力形式の取り出し（ascii,binary)
   * @details AVS形式で利用
   */
  CDM::E_CDM_FILE_TYPE
  Get_OutputFileTypeCoord() { return m_out_file_type_coord; };

  /**
   * @brief 出力先ディレクトリ名の取り出し
   */
  std::string Get_OutputDir() { return m_outdir_name; };

  /**
   * @brief 間引き数の取り出し
   */
  int Get_ThinOut() { return m_thin_count; };

  /**
   * @brief 出力配列形状の取り出し
   */
  CDM::E_CDM_ARRAYSHAPE
  Get_OutputArrayShape() { return m_outputArrayShape; };

  /**
   * @brief 出力配列形状のセット
   */
  void Set_OutputArrayShape( CDM::E_CDM_ARRAYSHAPE outputArrayShape ) 
       { m_outputArrayShape = outputArrayShape; }; 

  /**
   * @brief 出力ファイル名命名順の取り出し
   */
  CDM::E_CDM_OUTPUT_FNAME
  Get_OutputFilenameFormat() { return m_outputFilenameFormat; }; 

  /**
   * @brief 出力ガイドセル数の取り出し
   */
  int Get_OutputGuideCell() { return m_outputGuideCell; };

  /**
   * @brief 出力ガイドセル数の更新
   * @param[in] outgc 出力ガイドセル数
   */
  void Set_OutputGuideCell(int outgc) { m_outputGuideCell = outgc; }; 

  /**
   * @brief 並列処理時のファイル割振り方法の取り出し
   */
  E_CONV_OUTPUT_MULTI_FILE_CAST
  Get_MultiFileCasting() { return m_multiFileCasting; }; 

  /**
   * @brief 入力領域のスタート指示フラグの取り出し
   */
  bool Get_CropIndexStart_on() { return m_cropIndexStart_on; } 

  /**
   * @brief 入力領域のエンド指示フラグの取り出し
   */
  bool Get_CropIndexEnd_on() { return m_cropIndexEnd_on; } 

  /**
   * @brief 入力領域のスタート位置取り出し
   */
  int* Get_CropIndexStart() { return m_cropIndexStart; };  

  /**
   * @brief 入力領域のエンド位置取り出し
   */
  int* Get_CropIndexEnd() { return m_cropIndexEnd; }; 

  /**
   * @brief 入力領域のスタート位置の更新
   * @param[in] sta
   */
  void Set_CropIndexStart(int sta[3])
  {
    for(int i=0; i<3; i++) m_cropIndexStart[i]=sta[i];
  };

  /**
   * @brief 入力領域のエンド位置の更新
   * @param[in] end
   */
  void Set_CropIndexEnd(int end[3])
  {
    for(int i=0; i<3; i++) m_cropIndexEnd[i]=end[i];
  }; 
 
};

#endif // _INPUTPARAM_H_
