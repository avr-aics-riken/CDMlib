/*
 * fconv (File Converter)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   InputParam.C
 * @brief  InputParam Class
 * @author aics
 */

#include "InputParam.h"

// #################################################################
// コンストラクタ
InputParam::InputParam(cpm_ParaManager* paraMngr)
{

  //CPMのポインターをセット
  m_paraMngr=paraMngr;

  m_thin_count=1;
  m_out_file_type=CDM::E_CDM_FILE_TYPE_BINARY;
  m_out_file_type_coord=CDM::E_CDM_FILE_TYPE_BINARY;
  m_outputDiv[0]=-1;m_outputDiv[1]=-1;m_outputDiv[2]=-1;
  m_outputArrayShape=CDM::E_CDM_ARRAYSHAPE_UNKNOWN;
  m_outputFilenameFormat=CDM::E_CDM_FNAME_STEP_RANK;
  m_conv_type=E_CONV_OUTPUT_UNKNOWN;
  m_outputGuideCell=0;
  m_output_data_type = CDM::E_CDM_DTYPE_UNKNOWN;
  //m_in_dfi_name.clear();
  //m_out_dfi_name.clear();
  //m_out_proc_name.clear();
  m_dfiList.clear();
  m_output_dfi_on = false;
  m_bgrid_interp_flag = false;
  m_out_ftype_crd_on = false;

  m_cropIndexStart_on=false;
  m_cropIndexEnd_on  =false;
  m_multiFileCasting = E_CONV_OUTPUT_CAST_UNKNOWN;

}

// #################################################################
// デストラクタ
InputParam::~InputParam()
{

}

// #################################################################
// CPMのポインタをコピー
/*
bool InputParam::importCPM(cpm_ParaManager* paraMngr)
{
  if( !paraMngr ) return false;
  m_paraMngr = paraMngr;
  return true;
}
*/

// #################################################################
// 入力ファイルの読込み
bool InputParam::Read(std::string input_file_name)
{
  FILE* fp=NULL;
  string str;
  string label,label_base,label_leaf;

  // TextParserのインスタンス
  TextParser tpCntl;

  // 入力ファイルをセット
  int err = tpCntl.read(input_file_name);

  // node数の取得
  int nnode=0;
  label_base = "/ConvData";
  if( tpCntl.chkNode(label_base) ) {
    nnode = tpCntl.countLabels(label_base);
  } else Exit(0);

  vector<std::string> in_dfi_name;
  vector<std::string> out_dfi_name;
  vector<std::string> out_proc_name;
  in_dfi_name.clear();
  out_dfi_name.clear();
  out_proc_name.clear();

  int ncnt=0;
  label_base = "/ConvData";
  ncnt++;
  for(int i=0; i<nnode; i++) {
    if( !tpCntl.getNodeStr(label_base, ncnt, str) ) {
      printf("\tParsing error : No Elem name\n");
      Exit(0);
    }
    // 読込みdfiファイル名の読込み
    if( !strcasecmp(str.substr(0,8).c_str(), "InputDFI") ) {
      label = label_base+"/"+str;
      if ( !(tpCntl.getInspectedValue(label, str )) ) {
        printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      //dfiファイル名を格納
      in_dfi_name.push_back(str);
      ncnt++;
      continue;
    } else 
   
    //出力dfiファイル名の読込み
    if( !strcasecmp(str.substr(0,9).c_str(), "OutputDFI") ) {
      label = label_base+"/"+str;
      if ( !(tpCntl.getInspectedValue(label, str )) ) {
        printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      //dfiファイル名を格納
      out_dfi_name.push_back(str);
      ncnt++;
      continue;
    } else

    //出力procファイル名の読込み
    if( !strcasecmp(str.substr(0,13).c_str(), "OutputProcDFI") ) {
      label = label_base+"/"+str;
      if ( !(tpCntl.getInspectedValue(label, str )) ) {
        printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      //dfiファイル名を格納
      out_proc_name.push_back(str);
      ncnt++;
      continue;
    } else

    // コンバートタイプの読込み
    if( !strcasecmp(str.c_str(),"ConvType") ) {
      label = "/ConvData/ConvType";
      if( (tpCntl.getInspectedValue(label,str)) ) {
        if     ( !strcasecmp(str.c_str(), "Mx1") )  m_conv_type = E_CONV_OUTPUT_Mx1;
        else if( !strcasecmp(str.c_str(), "MxN") )  m_conv_type = E_CONV_OUTPUT_MxN;
        else if( !strcasecmp(str.c_str(), "MxM") )  m_conv_type = E_CONV_OUTPUT_MxM;
        else
        {
          Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
      }
      ncnt++;
      continue;
    } else

    // 出力分割数の読込み
    if( !strcasecmp(str.c_str(),"OutputDivision") ) {
      int vec[3];
      label = "/ConvData/OutputDivision";
      if( (tpCntl.getInspectedVector(label, vec, 3)) ) {
        if( vec[0]<1 || vec[1]<1 || vec[2]<1 ) {
           Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
           Exit(0);
        }
        m_outputDiv[0]=vec[0];
        m_outputDiv[1]=vec[1];
        m_outputDiv[2]=vec[2];
      }
      ncnt++;
      continue;
    } else

    // 出力ファイルフォーマットの読込み
    if( !strcasecmp(str.c_str(),"OutputFormat") ) {
      label = "/ConvData/OutputFormat";
      if ( !(tpCntl.getInspectedValue(label, str )) )
      {
         Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
         Exit(0);
      }
      if     ( !strcasecmp(str.c_str(), "sph" ) )    m_out_format = CDM::E_CDM_FMT_SPH;
      else if( !strcasecmp(str.c_str(), "bov" ) )    m_out_format = CDM::E_CDM_FMT_BOV;
      else if( !strcasecmp(str.c_str(), "avs" ) )    m_out_format = CDM::E_CDM_FMT_AVS;
      else if( !strcasecmp(str.c_str(), "plot3d" ) ) m_out_format = CDM::E_CDM_FMT_PLOT3D;
      else if( !strcasecmp(str.c_str(), "vtk" ) )    m_out_format = CDM::E_CDM_FMT_VTK;
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for  '%s'\n", label.c_str());
        Exit(0);
      }
      ncnt++;
      continue;
    } else

     // 出力データタイプの読込み
    if( !strcasecmp(str.c_str(),"OutputDataType") ) {
      label = "/ConvData/OutputDataType";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_output_data_type = CDM::E_CDM_DTYPE_UNKNOWN;
      } else {
        if     ( !strcasecmp(str.c_str(), "UInt8"  ) ) m_output_data_type = CDM::E_CDM_UINT8;
        else if( !strcasecmp(str.c_str(),  "Int8"  ) ) m_output_data_type = CDM::E_CDM_INT8;
        else if( !strcasecmp(str.c_str(), "UInt16" ) ) m_output_data_type = CDM::E_CDM_UINT16;
        else if( !strcasecmp(str.c_str(),  "Int16" ) ) m_output_data_type = CDM::E_CDM_INT16;
        else if( !strcasecmp(str.c_str(), "UInt32" ) ) m_output_data_type = CDM::E_CDM_UINT32;
        else if( !strcasecmp(str.c_str(),  "Int32" ) ) m_output_data_type = CDM::E_CDM_INT32;
        else if( !strcasecmp(str.c_str(), "UInt64" ) ) m_output_data_type = CDM::E_CDM_UINT64;
        else if( !strcasecmp(str.c_str(),  "Int64" ) ) m_output_data_type = CDM::E_CDM_INT64;
        else if( !strcasecmp(str.c_str(), "Float32") ) m_output_data_type = CDM::E_CDM_FLOAT32;
        else if( !strcasecmp(str.c_str(), "Float64") ) m_output_data_type = CDM::E_CDM_FLOAT64;
        else
        {
          printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
      }
      ncnt++;
      continue;
    }

    // 出力形式の読込み
    if( !strcasecmp(str.c_str(),"OutputFileType") ) {
      label = "/ConvData/OutputFileType";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_out_file_type = CDM::E_CDM_FILE_TYPE_BINARY;
      } else {
        if     ( !strcasecmp(str.c_str(), "binary") )        m_out_file_type=CDM::E_CDM_FILE_TYPE_BINARY;
        else if( !strcasecmp(str.c_str(), "ascii") )         m_out_file_type=CDM::E_CDM_FILE_TYPE_ASCII;
        else if( !strcasecmp(str.c_str(), "FortranBinary") ) m_out_file_type=CDM::E_CDM_FILE_TYPE_FBINARY;
        else
        {
           printf("\tInvalid keyword is described for '%s'\n", label.c_str());
           Exit(0);
        }
      }
      ncnt++;
      continue;
    } else

    // 座標データの出力形式の読込み(AVSおよびVTK形式)
    if( !strcasecmp(str.c_str(),"OutputFileTypeCoord") ) {
      label = "/ConvData/OutputFileTypeCoord";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_out_file_type_coord = CDM::E_CDM_FILE_TYPE_BINARY;
      } else {
        if     ( !strcasecmp(str.c_str(), "binary") ) m_out_file_type_coord=CDM::E_CDM_FILE_TYPE_BINARY;
        else if( !strcasecmp(str.c_str(), "ascii") )  m_out_file_type_coord=CDM::E_CDM_FILE_TYPE_ASCII;
        else
        {
           printf("\tInvalid keyword is described for '%s'\n", label.c_str());
           Exit(0);
        }
        m_out_ftype_crd_on = true;
      }
      ncnt++;
      continue;
    } else

    // 節点への補間オプション(AVSおよびVTK形式)
    if( !strcasecmp(str.c_str(),"OutputInterpolation") ) {
      label = "/ConvData/OutputInterpolation";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_bgrid_interp_flag = false;
      } else {
        if     ( !strcasecmp(str.c_str(), "true") ) m_bgrid_interp_flag = true;
        else if( !strcasecmp(str.c_str(), "false") ) m_bgrid_interp_flag = false;
        else
        {
           printf("\tInvalid keyword is described for '%s'\n", label.c_str());
           Exit(0);
        }
      }
      ncnt++;
      continue;
    } else

    // 出力先ディレクトリの読込み
    if( !strcasecmp(str.c_str(),"OutputDir") ) {
      label = "/ConvData/OutputDir";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      } else {
        m_outdir_name = str;
      }
      ncnt++;
      continue;
    } else

    //間引き数の読込み
    if( !strcasecmp(str.c_str(),"ThinningOut") ) {
      int ict;
      label = "/ConvData/ThinningOut";
      if( !(tpCntl.getInspectedValue(label, ict )) ) {
         m_thin_count=1;
      } else {
        if( ict < 0 ) {
           printf("\tInvalid keyword is described for '%s'\n", label.c_str());
           Exit(0);
        }
        m_thin_count = ict;
      }
      ncnt++;
      continue;
    } else

    //出力配列形状の読込み
    if( !strcasecmp(str.c_str(),"OutputArrayShape") ) {
      label = "/ConvData/OutputArrayShape";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_outputArrayShape = CDM::E_CDM_ARRAYSHAPE_UNKNOWN;  
      } else {
        if     ( !strcasecmp(str.c_str(), "ijkn") ) m_outputArrayShape = CDM::E_CDM_IJKN;
        else if( !strcasecmp(str.c_str(), "nijk") ) m_outputArrayShape = CDM::E_CDM_NIJK;
        else {
          printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
      }
      ncnt++;
      continue;
    } else

    //出力ファイル名命名順の読込み
    if( !strcasecmp(str.c_str(),"OutputFilenameFormat") ) {
      label = "/ConvData/OutputFilenameFormat";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_outputFilenameFormat = CDM::E_CDM_FNAME_STEP_RANK;
      } else {
        if     ( !strcasecmp(str.c_str(), "step_rank") ) m_outputFilenameFormat = CDM::E_CDM_FNAME_STEP_RANK;
        else if( !strcasecmp(str.c_str(), "rank_step") ) m_outputFilenameFormat = CDM::E_CDM_FNAME_RANK_STEP;
        else {
          printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
      }
      ncnt++;
      continue;
    } else

    //出力ガイドセル数
    if( !strcasecmp(str.c_str(),"OutputGuideCell") ) {
      int ict;
      label = "/ConvData/OutputGuideCell";
      if( !(tpCntl.getInspectedValue(label, ict )) ) {
        m_outputGuideCell=0;
      } else {
        if( ict < 0 ) {
          printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
        m_outputGuideCell=ict;
      }
      ncnt++;
      continue;
    } else

    //並列処理時のファイル割振り方法
    if( !strcasecmp(str.c_str(),"MultiFileCasting") ) {
      label = "/ConvData/MultiFileCasting";
      if( !(tpCntl.getInspectedValue(label, str )) ) {
        m_multiFileCasting = E_CONV_OUTPUT_STEP;
      } else {
        if     ( !strcasecmp(str.c_str(), "step") ) m_multiFileCasting = E_CONV_OUTPUT_STEP;
        else if( !strcasecmp(str.c_str(), "rank") ) m_multiFileCasting = E_CONV_OUTPUT_RANK;
        else {
          printf("\tInvalid keyword is described for '%s'\n", label.c_str());
          Exit(0);
        }
      }
      ncnt++;
      continue;
    } else 

    //入力領域のスタート位置)
    if( !strcasecmp(str.c_str(),"CropIndexStart") ) { 
      int vec[3];
      label = "/ConvData/CropIndexStart";
      if( (tpCntl.getInspectedVector(label, vec, 3)) ) {
        m_cropIndexStart[0]=vec[0];
        m_cropIndexStart[1]=vec[1];
        m_cropIndexStart[2]=vec[2];
        m_cropIndexStart_on=true;
      } else {
        printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
      ncnt++;
      continue;
    } else

    //入力領域のエンド位置
    if( !strcasecmp(str.c_str(),"CropIndexEnd") ) {
      int vec[3];
      label = "/ConvData/CropIndexEnd";
      if( (tpCntl.getInspectedVector(label, vec, 3)) ) {
        m_cropIndexEnd[0]=vec[0];
        m_cropIndexEnd[1]=vec[1];
        m_cropIndexEnd[2]=vec[2];
        m_cropIndexEnd_on=true;
      } else {
        printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
      ncnt++;
      continue;
    } else {
      printf("\tParsing error : No Elem name : %s\n",str.c_str());
      Exit(0);
    } 
  }

  // TextParserの破棄
  tpCntl.remove();

  //入力dfiファイル名の登録
  if( in_dfi_name.size() > 0 ) {
    dfi_info dInfo;
    for(int i=0; i<in_dfi_name.size(); i++){
      dInfo.in_dfi_name  = in_dfi_name[i];
      dInfo.in_dfi       = NULL;
      dInfo.out_dfi_name = "";
      dInfo.out_proc_name= "";
      m_dfiList.push_back(dInfo);
    }
    in_dfi_name.clear();
  }

  //出力dfiファイル名の登録
  if( out_dfi_name.size() > 0 ) {
    if( out_dfi_name.size() != m_dfiList.size() ) {
      printf("\tmismatch outout dfi file names\n");
      out_dfi_name.clear();
      return false;
    }
    for(int i=0; i<m_dfiList.size(); i++) {
       m_dfiList[i].out_dfi_name = out_dfi_name[i];
    }
    out_dfi_name.clear();
    m_output_dfi_on = true;
  } else m_output_dfi_on = false;


  //出力procファイル名の登録
  if( out_proc_name.size() > 0 ) {
    if( out_proc_name.size() != m_dfiList.size() ) {
      printf("\tmismatch outout proc file names\n");
      out_proc_name.clear();
      return false;
    }
    for(int i=0; i<m_dfiList.size(); i++) {
       m_dfiList[i].out_proc_name = out_proc_name[i];
    }
    out_proc_name.clear();
  }

  //入力パラメータのチェック
  if( !InputParamCheck() ) return false;

  return true;

}

// #################################################################
// 入力パラメータのチェック
bool InputParam::InputParamCheck()
{

  bool ierr=true;

  //読込むDFIファイル名指示のチェック
  if( m_dfiList.size() < 1 ) {
    printf("\tundefined input dfi file name\n");
    ierr=false;
  }

  //コンバートタイプのチェック
  if( m_conv_type == E_CONV_OUTPUT_UNKNOWN ) {
    printf("\tundefine Converter Type\n");
    ierr=false;
  }

  //出力先ディレクトリのチェック
  if( m_outdir_name.empty() ) {
    printf("\tundefine OutputDir\n");
    ierr=false;
  }

  //出力形式のチェック
  if( m_out_file_type == CDM::E_CDM_FILE_TYPE_ASCII ) {
    if( m_out_format != CDM::E_CDM_FMT_PLOT3D &&
        m_out_format != CDM::E_CDM_FMT_VTK ) {
      printf("\tCan't Converter OutputFileType ascii.\n");
      ierr=false;
    }
  }

  //座標データの出力形式のチェック (AVSおよびVTK形式)
  //AVS,VTK以外で座標データの出力形式指示があった場合はメッセージを出力する
  if( m_out_ftype_crd_on ) {
    if( m_out_format != CDM::E_CDM_FMT_AVS &&
        m_out_format != CDM::E_CDM_FMT_VTK ) {
      printf("\tCan't use OutputFileTypeCoord. %s\n",Get_OutputFormat_string().c_str());
    }
  }

  //節点への補間オプションのチェック (AVSおよびVTK形式)
  //AVS,VTK以外で節点への補間オプションがONに指定された場合はOFFにし、メッセージを出力する
  if( m_bgrid_interp_flag ) {
    if( m_out_format != CDM::E_CDM_FMT_AVS &&
        m_out_format != CDM::E_CDM_FMT_VTK ) {
      m_bgrid_interp_flag = false;
      printf("\tCan't use OutputInterpolation. %s\n",Get_OutputFormat_string().c_str());
    }
  }

  //出力配列形状のチェック
  //BOV以外での出力配列形状指示は無効とし、自動的に対応する配列形状で出力
  //なので、指定があった場合はメッセージを出力する
  if( m_outputArrayShape != CDM::E_CDM_ARRAYSHAPE_UNKNOWN ) {
    if( m_out_format != CDM::E_CDM_FMT_BOV ) printf("\tCan't OutputArrayShape.\n");
  }

  //出力配列形状のセット
  if     ( m_out_format == CDM::E_CDM_FMT_SPH ) m_outputArrayShape = CDM::E_CDM_NIJK;
  else if( m_out_format == CDM::E_CDM_FMT_AVS ) m_outputArrayShape = CDM::E_CDM_NIJK;
  else if( m_out_format == CDM::E_CDM_FMT_VTK ) m_outputArrayShape = CDM::E_CDM_NIJK;
  else if( m_out_format == CDM::E_CDM_FMT_PLOT3D ) m_outputArrayShape = CDM::E_CDM_IJKN;
  else if( m_out_format == CDM::E_CDM_FMT_BOV && 
           m_outputArrayShape == CDM::E_CDM_ARRAYSHAPE_UNKNOWN ) 
       m_outputArrayShape = CDM::E_CDM_NIJK;

  //未対応のデータ型への変換チェック
  if( m_output_data_type != CDM::E_CDM_DTYPE_UNKNOWN ) {
    //(sph,plot3d)
    if( m_out_format == CDM::E_CDM_FMT_SPH || m_out_format == CDM::E_CDM_FMT_PLOT3D ) {
      if( m_output_data_type != CDM::E_CDM_FLOAT32 &&
          m_output_data_type != CDM::E_CDM_FLOAT64 ) {
        printf("\tCan't Converter OutputDataType %s.\n",
                Get_OutputDataType_string().c_str());
        ierr=false;
      }
    }

    //(avs)
    if( m_out_format == CDM::E_CDM_FMT_AVS ) {
      if( m_output_data_type != CDM::E_CDM_INT8    &&
          m_output_data_type != CDM::E_CDM_INT16   &&
          m_output_data_type != CDM::E_CDM_INT32   &&
          m_output_data_type != CDM::E_CDM_FLOAT32 &&
          m_output_data_type != CDM::E_CDM_FLOAT64 ) {
        printf("\tCan't Converter OutputDataType %s.\n",
                Get_OutputDataType_string().c_str());
        ierr=false;
      }
    }
  }

  //DFI出力のチェック、出力するDFIファイル名のチェック
  if( m_output_dfi_on ) {
    //DFI出力がSPH，BOV以外で指定された場合はエラー
    //if( m_out_format != CDM::E_CDM_FMT_SPH && m_out_format != CDM::E_CDM_FMT_BOV ) {
    if( m_out_format != CDM::E_CDM_FMT_SPH && m_out_format != CDM::E_CDM_FMT_BOV && m_out_format != CDM::E_CDM_FMT_PLOT3D) {
      printf("\tCan't output dfi OutputFormat. %s\n",Get_OutputFormat_string().c_str());
      ierr=false;
    }
    //出力DFI名のチェック、パスが指定されている場合はエラー
    for(int i=0; i<m_dfiList.size(); i++) {
      std::string inPath = CDM::cdmPath_DirName(m_dfiList[i].out_dfi_name);
      if( inPath != "" && inPath !="./" ) {
        printf("\tIllegal OutputDFI : %s\n",m_dfiList[i].out_dfi_name.c_str());
        ierr=false;
      }
      if( m_dfiList[i].out_proc_name != "" ) {
        inPath = CDM::cdmPath_DirName(m_dfiList[i].out_proc_name);
        if( inPath != "" && inPath !="./" ) {
          printf("\tIllegal OutputProc : %s\n",m_dfiList[i].out_proc_name.c_str());
          ierr=false;
        }
      } else {
      //出力するprocファイル名が省略された場合、出力するdfiファイル名から生成
        std::string proc = CDM::ExtractPathWithoutExt(m_dfiList[i].out_dfi_name);
        std::string fname = proc+"_proc.dfi";
        m_dfiList[i].out_proc_name = fname;
      }
    }
  }

  //出力ガイドセル数のチェック
  if( m_outputGuideCell > 0 ) {
    //sph,bov以外は出力指定不可
    if( m_out_format != CDM::E_CDM_FMT_SPH && m_out_format != CDM::E_CDM_FMT_BOV && m_out_format != CDM::E_CDM_FMT_PLOT3D ) {
    //if( m_out_format != CDM::E_CDM_FMT_SPH && m_out_format != CDM::E_CDM_FMT_BOV ) {
      printf("\tCan't output guide cell : %s\n",Get_OutputFormat_string().c_str());
      ierr=false;
    }
    //間引きありとガイドセル出力を両方指定は不可
    if( m_thin_count > 1 ) {
      printf("\tCan't output guide cell and thinning out\n");
      ierr=false;
    }
  }

  //ファイル割振り方法のチェック
  if( m_multiFileCasting == E_CONV_OUTPUT_RANK && m_conv_type == E_CONV_OUTPUT_Mx1 ) {
    printf("\tCan't multi file casting type rank\n");
  }

  //入力領域指示のチェック
  if( m_conv_type == E_CONV_OUTPUT_MxM ) {
    if( m_cropIndexStart_on || m_cropIndexEnd_on ) {
      printf("\tCan't define CropIndex option exec type MxM\n");
      ierr=false;
    }
  }

  return ierr;

}

// #################################################################
// 入力パラメータのログ出力
void InputParam::PrintParam(FILE* fp)
{
   fprintf(fp,"\n");
   fprintf(fp,"*** fconv file info ***\n");
   fprintf(fp,"\n");

   for(int i=0; i<m_dfiList.size(); i++) {
     fprintf(fp,"\tInputDFI[@]          : \"%s\"\n",m_dfiList[i].in_dfi_name.c_str());
   }

   if( m_output_dfi_on ) {
     for(int i=0; i<m_dfiList.size(); i++) {
       fprintf(fp,"\tOutputDFI[@]         : \"%s\"\n",m_dfiList[i].out_dfi_name.c_str());
     }
     for(int i=0; i<m_dfiList.size(); i++) {
       fprintf(fp,"\tOutputProcDFI[@]     : \"%s\"\n",m_dfiList[i].out_proc_name.c_str());
     }
   }

   if( m_conv_type == E_CONV_OUTPUT_Mx1 ) {
     fprintf(fp,"\tConvType             : \"Mx1\"\n");
   } else if( m_conv_type == E_CONV_OUTPUT_MxM ) {
     fprintf(fp,"\tConvType             : \"MxM\"\n");
   } else if( m_conv_type == E_CONV_OUTPUT_MxN ) {
     fprintf(fp,"\tConvType             : \"MxN\"\n");
   } else {
     Exit(0);
   }

   fprintf(fp,"\tOutputDivision       : (%d,%d,%d)\n",m_outputDiv[0],m_outputDiv[1],m_outputDiv[2]);

   if( m_out_format == CDM::E_CDM_FMT_SPH ) {
     fprintf(fp,"\tOutputFormat         : \"sph\"\n");
   }else if( m_out_format == CDM::E_CDM_FMT_BOV ) {
     fprintf(fp,"\tOutputFormat         : \"bov\"\n");
   }else if( m_out_format == CDM::E_CDM_FMT_AVS ) {
     fprintf(fp,"\tOutputFormat         : \"avs\"\n");
   }else if( m_out_format == CDM::E_CDM_FMT_PLOT3D ) {
     fprintf(fp,"\tOutputFormat         : \"plot3d\"\n");
   }else if( m_out_format == CDM::E_CDM_FMT_VTK ) {
     fprintf(fp,"\tOutputFormat         : \"vtk\"\n");
   }else {
     Exit(0);
   }

   if( m_output_data_type == CDM::E_CDM_INT8 ) {
     fprintf(fp,"\tOutputDataType       : \"Int8\"\n");
   }else if( m_output_data_type == CDM::E_CDM_INT16 ) {
     fprintf(fp,"\tOutputDataType       : \"Int16\"\n");
   }else if( m_output_data_type == CDM::E_CDM_INT32 ) {
     fprintf(fp,"\tOutputDataType       : \"Int32\"\n");
   }else if( m_output_data_type == CDM::E_CDM_INT64 ) {
     fprintf(fp,"\tOutputDataType       : \"Int64\"\n");
   }else if( m_output_data_type == CDM::E_CDM_UINT8 ) {
     fprintf(fp,"\tOutputDataType       : \"UInt8\"\n");
   }else if( m_output_data_type == CDM::E_CDM_UINT16 ) {
     fprintf(fp,"\tOutputDataType       : \"UInt16\"\n");
   }else if( m_output_data_type == CDM::E_CDM_UINT32 ) {
     fprintf(fp,"\tOutputDataType       : \"UInt32\"\n");
   }else if( m_output_data_type == CDM::E_CDM_UINT64 ) {
     fprintf(fp,"\tOutputDataType       : \"UInt64\"\n");
   }else if( m_output_data_type == CDM::E_CDM_FLOAT32 ) {
     fprintf(fp,"\tOutputDataType       : \"Float32\"\n");
   }else if( m_output_data_type == CDM::E_CDM_FLOAT64 ) {
     fprintf(fp,"\tOutputDataType       : \"Float64\"\n");
   } else {
     fprintf(fp,"\tOutputDataType       : \"undefine\"\n");
   }

   if( m_out_file_type == CDM::E_CDM_FILE_TYPE_DEFAULT ) {
     fprintf(fp,"\tOutputFileType       : \"undefine\"\n");
   } else if( m_out_file_type == CDM::E_CDM_FILE_TYPE_ASCII ) {
     fprintf(fp,"\tOutputFileType       : \"ascii\"\n");
   } else if( m_out_file_type == CDM::E_CDM_FILE_TYPE_BINARY ) {
     fprintf(fp,"\tOutputFileType       : \"binary\"\n");
   } else if( m_out_file_type == CDM::E_CDM_FILE_TYPE_FBINARY ) {
     fprintf(fp,"\tOutputFileType       : \"Fortran Binary\"\n");
   }

   if( m_out_ftype_crd_on ) {
     if( m_out_file_type_coord == CDM::E_CDM_FILE_TYPE_DEFAULT ) {
       fprintf(fp,"\tOutputFileTypeCoord  : \"undefine\"\n");
     } else if( m_out_file_type_coord == CDM::E_CDM_FILE_TYPE_ASCII ) {
       fprintf(fp,"\tOutputFileTypeCoord  : \"ascii\"\n");
     } else if( m_out_file_type_coord == CDM::E_CDM_FILE_TYPE_BINARY ) {
       fprintf(fp,"\tOutputFileTypeCoord  : \"binary\"\n");
     }
   }

   if( m_bgrid_interp_flag ) {
     fprintf(fp,"\tOutputInterpolation  : \"true\"\n");
   } else {
     fprintf(fp,"\tOutputInterpolation  : \"false\"\n");
   }

   fprintf(fp,"\tOutputdir            : \"%s\"\n",m_outdir_name.c_str());
   fprintf(fp,"\tThinning             : %d\n",m_thin_count);

   if( m_outputArrayShape == CDM::E_CDM_IJKN ) {
     fprintf(fp,"\tOutputArrayShape     : \"ijkn\"\n");
   }else if( m_outputArrayShape == CDM::E_CDM_NIJK ) {
     fprintf(fp,"\tOutputArrayShape     : \"nijk\"\n");
   }else {
     fprintf(fp,"\tOutputArrayShape     : \"undefine\"\n");
   }

   if( m_outputFilenameFormat == CDM::E_CDM_FNAME_STEP_RANK ) { 
     fprintf(fp,"\tOutputFilenameFormat : \"step_rank\"\n");
   }else if( m_outputFilenameFormat == CDM::E_CDM_FNAME_RANK_STEP ) { 
     fprintf(fp,"\tOutputFilenameFormat : \"rank_step\"\n");
   }else {
     fprintf(fp,"\tOutputFilenameFormat : \"undefine\"\n");
   }

   fprintf(fp,"\tOutputGuideCell      : %d\n",m_outputGuideCell);

   if( m_multiFileCasting == E_CONV_OUTPUT_STEP ) {
     fprintf(fp,"\tMultiFileCasting     : \"step\"\n");
   }else if( m_multiFileCasting == E_CONV_OUTPUT_RANK ) {
     fprintf(fp,"\tMultiFileCasting     : \"rank\"\n");
   }else {
     fprintf(fp,"\tMultiFileCasting     : \"undefine\"\n");
   }

   if( m_cropIndexStart_on ) {
     fprintf(fp,"\tCropIndexStart       : (%d,%d,%d)\n",m_cropIndexStart[0],
                                                        m_cropIndexStart[1],
                                                        m_cropIndexStart[2]);
   }
   if( m_cropIndexEnd_on ) {
     fprintf(fp,"\tCropIndexEnd         : (%d,%d,%d)\n",m_cropIndexEnd[0],
                                                        m_cropIndexEnd[1],
                                                        m_cropIndexEnd[2]);
   }

   fprintf(fp,"\n");
}
