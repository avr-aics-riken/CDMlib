/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_FileInfo.C
 * @brief  cdm_FileInfo Class
 * @author aics    
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K

// #################################################################
// コンストラクタ
cdm_FileInfo::cdm_FileInfo()
{
  DFIType          = CDM::E_CDM_DFITYPE_UNKNOWN;
  FieldFilenameFormat = CDM::E_CDM_FNAME_DEFAULT;
  DirectoryPath    ="";
  TimeSliceDirFlag =CDM::E_CDM_OFF;
  Prefix           ="";
  FileFormat       =CDM::E_CDM_FMT_UNKNOWN;
  GuideCell        =0;
  DataType         =CDM::E_CDM_DTYPE_UNKNOWN;
  Endian           =CDM::E_CDM_ENDIANTYPE_UNKNOWN;
  ArrayShape       =CDM::E_CDM_ARRAYSHAPE_UNKNOWN;
  NumVariables     =0;
}

// #################################################################
// コンストラクタ
cdm_FileInfo::cdm_FileInfo(const CDM::E_CDM_DFITYPE _DFIType,
                           const CDM::E_CDM_OUTPUT_FNAME _FieldFilenameFormat,
                           const std::string _DirectoryPath, 
                           const CDM::E_CDM_ONOFF _TimeSliceDirFlag, 
                           const std::string _Prefix,
                           const CDM::E_CDM_FORMAT _FileFormat,
                           const int _GuideCell, 
                           const CDM::E_CDM_DTYPE _DataType, 
                           const CDM::E_CDM_ENDIANTYPE _Endian,
                           const CDM::E_CDM_ARRAYSHAPE _ArrayShape, 
                           const int _NumVariables)
{
//FCONV 20140116.s
  DFIType          =_DFIType;
  FieldFilenameFormat =_FieldFilenameFormat;
//FCONV 20140116.e
  DirectoryPath    =_DirectoryPath;
  Prefix           =_Prefix;
  TimeSliceDirFlag =_TimeSliceDirFlag;
  FileFormat       =_FileFormat;
  GuideCell        =_GuideCell;
  DataType         =_DataType;
  Endian           =_Endian;
  ArrayShape       =_ArrayShape;
  NumVariables     =_NumVariables;
}

// デストラクタ
cdm_FileInfo::~cdm_FileInfo()
{

}

// #################################################################
// 変数名をセットする
void cdm_FileInfo::setVariableName(int pvari,
                                   std::string variName)
{

  if( VariableName.size()>pvari+1 ) {
    VariableName[pvari]=variName;
  } else {
    for(int i=VariableName.size(); i<pvari+1; i++) {
      VariableName.push_back(variName);
    }
  }
}

// #################################################################
// 変数名を取得する
std::string cdm_FileInfo::getVariableName(int pvari)
{
  std::string VariName="";
  if(VariableName.size()<pvari+1) return VariName;
  return VariableName[pvari];
}

// #################################################################
// FileInfo 読込み関数
CDM::E_CDM_ERRORCODE
cdm_FileInfo::Read(cdm_TextParser tpCntl) 
{

  std::string str;
  std::string label,label_base,label_leaf,label_leaf_leaf;
  int ct;

  int ncnt=0;

//FCONV 20140116.s
  //DFIType
  label = "/FileInfo/DFIType";
  if( !(tpCntl.GetValue(label, &str )) ) {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_DFITYPE;
  } else {
    if( !strcasecmp(str.c_str(),D_CDM_DFITYPE_CARTESIAN ) ) {
      DFIType = CDM::E_CDM_DFITYPE_CARTESIAN;
    } else if( !strcasecmp(str.c_str(),D_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) ) {
      DFIType = CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN;
    } else {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_DFITYPE;
    }
    ncnt++;
  }


  //FieldFilenameFormat
  label = "/FileInfo/FieldFilenameFormat";
  if( !(tpCntl.GetValue(label, &str )) ) {
    FieldFilenameFormat = CDM::E_CDM_FNAME_DEFAULT;
  } else {
    if( !strcasecmp(str.c_str(),"step_rank" ) ) {
      FieldFilenameFormat = CDM::E_CDM_FNAME_STEP_RANK;
    }else if( !strcasecmp(str.c_str(),"rank_step" ) ) {
      FieldFilenameFormat = CDM::E_CDM_FNAME_RANK_STEP;
    }else {
      printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
      return CDM::E_CDM_ERROR_READ_DFI_FIELDFILENAMEFORMAT;
    }
    ncnt++;
  }


  //Directorypath
  label = "/FileInfo/DirectoryPath";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_DIRECTORYPATH;
  }
  DirectoryPath=str;

  ncnt++;

  //TimeSilceDirectory
  label = "/FileInfo/TimeSliceDirectory";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_TIMESLICEDIRECTORY;
  }

  if( !strcasecmp(str.c_str(),"on" ) ) {
    TimeSliceDirFlag=CDM::E_CDM_ON;
  } else if( !strcasecmp(str.c_str(),"off" ) ) {
    TimeSliceDirFlag=CDM::E_CDM_OFF;
  } else {
    printf("\tCDM Parsing error : fail to get '%s'\n",str.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_TIMESLICEDIRECTORY;
  }

  ncnt++;

  //Prefix
  label = "/FileInfo/Prefix";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_PREFIX;
  }
  Prefix=str;

  ncnt++;

  //FileFormat
  label = "/FileInfo/FileFormat";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_FILEFORMAT;
  }
  if( !strcasecmp(str.c_str(), D_CDM_EXT_SPH ) ) {
    //Check DFIType
    if( DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
      printf("\tCDM error : Non_Uniform_Cartesian is not supported in SPH File Format.\n");
      return CDM::E_CDM_ERROR_READ_DFI_FILEFORMAT;
    }
    FileFormat=CDM::E_CDM_FMT_SPH;
    ArrayShape=CDM::E_CDM_NIJK;
  }
  else if( !strcasecmp(str.c_str(), D_CDM_EXT_BOV ) ) {
    //Check DFIType
    if( DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
      printf("\tCDM error : Non_Uniform_Cartesian is not supported in BOV File Format.\n");
      return CDM::E_CDM_ERROR_READ_DFI_FILEFORMAT;
    }
    FileFormat=CDM::E_CDM_FMT_BOV;
    ArrayShape=CDM::E_CDM_IJKN;
  }
  else if( !strcasecmp(str.c_str(), "plot3d" ) ) {
    FileFormat=CDM::E_CDM_FMT_PLOT3D;
    ArrayShape=CDM::E_CDM_IJKN;
  }
  else FileFormat=CDM::E_CDM_FMT_UNKNOWN;

  ncnt++;

  //GuidCell
  label = "/FileInfo/GuideCell";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_GUIDECELL;
  }
  GuideCell=ct;

  ncnt++;

  //DataType
  label = "/FileInfo/DataType";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_DATATYPE;
  }
  DataType=cdm_DFI::ConvDatatypeS2E(str);

  ncnt++;

  //Endian
  label = "/FileInfo/Endian";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_ENDIAN;
  }
  if( !strcasecmp(str.c_str(),"little" ) ) {
    Endian=CDM::E_CDM_LITTLE;
  }else if( !strcasecmp(str.c_str(),"big" ) ) {
    Endian=CDM::E_CDM_BIG;
  }else {
    printf("\tCDM Parsing error : fail to get '%s'\n",str.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_ENDIAN;
  } 

  ncnt++;

  //NumVariables
  label = "/FileInfo/NumVariables";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
    return CDM::E_CDM_ERROR_READ_DFI_NUMVARIABLES;
  }
  NumVariables=ct;

  ncnt++;

  //Variable
  int nvari=0;
  label_leaf_leaf = "/FileInfo/Variable";
  if ( tpCntl.chkNode(label_leaf_leaf) )  //があれば
  {
    nvari = tpCntl.countLabels(label_leaf_leaf);
  }

  ncnt++;

  //フィールドデータの変数の個数と登録された変数名の個数の一致確認
  if ( NumVariables != nvari) {
    printf("\tCDM Parsing error : Number of valiable names\n");
    return CDM::E_CDM_ERROR_UNMATCH_NUM_OF_VARIABLES;
  }

  label_leaf = "/FileInfo";

  if( nvari>0 ) {
    for(int i=0; i<nvari; i++) {
      if(!tpCntl.GetNodeStr(label_leaf,ncnt+i,&str))
      {
        printf("\tCDM Parsing error : No Elem name\n");
        return CDM::E_CDM_ERROR_READ_DFI_NO_MINMAX;
      }
      if( !strcasecmp(str.substr(0,8).c_str(), "variable") ) {
        label_leaf_leaf = label_leaf+"/"+str;

        label = label_leaf_leaf + "/name";
        if ( !(tpCntl.GetValue(label, &str )) ) {
          printf("\tCDM Parsing error : fail to get '%s'\n",label.c_str());
          return CDM::E_CDM_ERROR_READ_DFI_MIN;
        }
        else {
          VariableName.push_back(str);
        }
      }
    }
  }

//Check for SPH format 20141022.s
  if( FileFormat == CDM::E_CDM_FMT_SPH ) {
    if( NumVariables != 1 && NumVariables != 3 ) {
      printf("\tCDM error sph file undefined except for number of valiables 1 or 3.\n");
      return CDM::E_CDM_ERROR_READ_FILEINFO;
    }
  }
//Check for SPH format 20141022.e

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// DFIファイル:FileInfo要素を出力する
CDM::E_CDM_ERRORCODE
cdm_FileInfo::Write(FILE* fp, 
                    const unsigned tab) 
{

  fprintf(fp, "FileInfo {\n");
  fprintf(fp, "\n");

//FCONV 20140116.s
  _CDM_WRITE_TAB(fp, tab+1);
  if(        DFIType == CDM::E_CDM_DFITYPE_CARTESIAN ) {
    fprintf(fp, "DFIType            = \"Cartesian\"\n");
  } else if( DFIType == CDM::E_CDM_DFITYPE_NON_UNIFORM_CARTESIAN ) {
    fprintf(fp, "DFIType            = \"Non_Uniform_Cartesian\"\n");
  }
//FCONV 20140116.s

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "DirectoryPath      = \"%s\"\n", DirectoryPath.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  if(        TimeSliceDirFlag == CDM::E_CDM_OFF ) {
    fprintf(fp, "TimeSliceDirectory = \"off\"\n");
  } else if( TimeSliceDirFlag == CDM::E_CDM_ON ) {
    fprintf(fp, "TimeSliceDirectory = \"on\"\n");
  }

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "Prefix             = \"%s\"\n", Prefix.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  if(        FileFormat == CDM::E_CDM_FMT_SPH ) {
    fprintf(fp, "FileFormat         = \"sph\"\n");
  } else if( FileFormat == CDM::E_CDM_FMT_BOV ) {
    fprintf(fp, "FileFormat         = \"bov\"\n");
  } else if( FileFormat == CDM::E_CDM_FMT_PLOT3D ) {
    fprintf(fp, "FileFormat         = \"plot3d\"\n");
  }  else if( FileFormat == CDM::E_CDM_FMT_VTK ) {
    fprintf(fp, "FileFormat         = \"vtk\"\n");
  }

//FCONV 20140116.s
  _CDM_WRITE_TAB(fp, tab+1);
  if( FieldFilenameFormat == CDM::E_CDM_FNAME_RANK_STEP ) {
    fprintf(fp, "FieldFilenameFormat= \"rank_step\"\n");
  } else {
    fprintf(fp, "FieldFilenameFormat= \"step_rank\"\n");
  }
//FCONV 20140116.e

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "GuideCell          = %d\n", GuideCell);

  _CDM_WRITE_TAB(fp, tab+1);
  std::string Dtype = cdm_DFI::ConvDatatypeE2S((CDM::E_CDM_DTYPE)DataType);
  fprintf(fp, "DataType           = \"%s\"\n",Dtype.c_str());

  _CDM_WRITE_TAB(fp, tab+1);
  if( Endian == CDM::E_CDM_LITTLE) {
    fprintf(fp, "Endian             = \"little\"\n");
  } else {
    fprintf(fp, "Endian             = \"big\"\n");
  }

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "NumVariables       = %d\n",NumVariables);

/*
  if( VariableName.size()>0 ) {
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",VariableName[0].c_str());
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",VariableName[1].c_str());
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",VariableName[2].c_str());
  }
*/
  for(int i=0; i<VariableName.size(); i++) {
    _CDM_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",VariableName[i].c_str());
  }

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}

