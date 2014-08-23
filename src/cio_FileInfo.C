/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_FileInfo.C
 * @brief  cio_FileInfo Class
 * @author aics    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K

// #################################################################
// コンストラクタ
cio_FileInfo::cio_FileInfo()
{
  DFIType          = CIO::E_CIO_DFITYPE_UNKNOWN;
  FieldFilenameFormat = CIO::E_CIO_FNAME_DEFAULT;
  DirectoryPath    ="";
  TimeSliceDirFlag =CIO::E_CIO_OFF;
  Prefix           ="";
  FileFormat       =CIO::E_CIO_FMT_UNKNOWN;
  GuideCell        =0;
  DataType         =CIO::E_CIO_DTYPE_UNKNOWN;
  Endian           =CIO::E_CIO_ENDIANTYPE_UNKNOWN;
  ArrayShape       =CIO::E_CIO_ARRAYSHAPE_UNKNOWN;
  Component        =0;
}

// #################################################################
// コンストラクタ
cio_FileInfo::cio_FileInfo(const CIO::E_CIO_DFITYPE _DFIType,
                           const CIO::E_CIO_OUTPUT_FNAME _FieldFilenameFormat,
                           const std::string _DirectoryPath, 
                           const CIO::E_CIO_ONOFF _TimeSliceDirFlag, 
                           const std::string _Prefix,
                           const CIO::E_CIO_FORMAT _FileFormat,
                           const int _GuideCell, 
                           const CIO::E_CIO_DTYPE _DataType, 
                           const CIO::E_CIO_ENDIANTYPE _Endian,
                           const CIO::E_CIO_ARRAYSHAPE _ArrayShape, 
                           const int _Component)
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
  Component        =_Component;
}

// デストラクタ
cio_FileInfo::~cio_FileInfo()
{

}

// #################################################################
// 成分名をセットする
void cio_FileInfo::setComponentVariable(int pcomp,
                                        std::string compName)
{

  if( ComponentVariable.size()>pcomp+1 ) {
    ComponentVariable[pcomp]=compName;
  } else {
    for(int i=ComponentVariable.size(); i<pcomp+1; i++) {
      ComponentVariable.push_back(compName);
    }
  }
}

// #################################################################
// 成分名を取得する
std::string cio_FileInfo::getComponentVariable(int pcomp)
{
  std::string CompName="";
  if(ComponentVariable.size()<pcomp+1) return CompName;
  return ComponentVariable[pcomp];
}

// #################################################################
// FileInfo 読込み関数
CIO::E_CIO_ERRORCODE
cio_FileInfo::Read(cio_TextParser tpCntl) 
{

  std::string str;
  std::string label,label_base,label_leaf,label_leaf_leaf;
  int ct;

  int ncnt=0;

//FCONV 20140116.s
  //DFIType
  label = "/FileInfo/DFIType";
  if( !(tpCntl.GetValue(label, &str )) ) {
    DFIType = CIO::E_CIO_DFITYPE_CARTESIAN;
  } else {
    if( !strcasecmp(str.c_str(),"Cartesian" ) ) {
      DFIType = CIO::E_CIO_DFITYPE_CARTESIAN;
    } else {
      printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
      return CIO::E_CIO_ERROR_READ_DFI_DFITYPE;
    }
    ncnt++;
  }


  //FieldFilenameFormat
  label = "/FileInfo/FieldFilenameFormat";
  if( !(tpCntl.GetValue(label, &str )) ) {
    FieldFilenameFormat = CIO::E_CIO_FNAME_DEFAULT;
  } else {
    if( !strcasecmp(str.c_str(),"step_rank" ) ) {
      FieldFilenameFormat = CIO::E_CIO_FNAME_STEP_RANK;
    }else if( !strcasecmp(str.c_str(),"rank_step" ) ) {
      FieldFilenameFormat = CIO::E_CIO_FNAME_RANK_STEP;
    }else {
      printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
      return CIO::E_CIO_ERROR_READ_DFI_FIELDFILENAMEFORMAT;
    }
    ncnt++;
  }


  //Directorypath
  label = "/FileInfo/DirectoryPath";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_DIRECTORYPATH;
  }
  DirectoryPath=str;

  ncnt++;

  //TimeSilceDirectory
  label = "/FileInfo/TimeSliceDirectory";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_TIMESLICEDIRECTORY;
  }

  if( !strcasecmp(str.c_str(),"on" ) ) {
    TimeSliceDirFlag=CIO::E_CIO_ON;
  } else if( !strcasecmp(str.c_str(),"off" ) ) {
    TimeSliceDirFlag=CIO::E_CIO_OFF;
  } else {
    printf("\tCIO Parsing error : fail to get '%s'\n",str.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_TIMESLICEDIRECTORY;
  }

  ncnt++;

  //Prefix
  label = "/FileInfo/Prefix";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_PREFIX;
  }
  Prefix=str;

  ncnt++;

  //FileFormat
  label = "/FileInfo/FileFormat";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_FILEFORMAT;
  }
  if( !strcasecmp(str.c_str(),"sph" ) ) {
    FileFormat=CIO::E_CIO_FMT_SPH;
  }
  else if( !strcasecmp(str.c_str(),"bov" ) ) {
    FileFormat=CIO::E_CIO_FMT_BOV;
  }
  else FileFormat=CIO::E_CIO_FMT_UNKNOWN;

  ncnt++;

  //GuidCell
  label = "/FileInfo/GuideCell";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_GUIDECELL;
  }
  GuideCell=ct;

  ncnt++;

  //DataType
  label = "/FileInfo/DataType";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_DATATYPE;
  }
  DataType=cio_DFI::ConvDatatypeS2E(str);

  ncnt++;

  //Endian
  label = "/FileInfo/Endian";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_ENDIAN;
  }
  if( !strcasecmp(str.c_str(),"little" ) ) {
    Endian=CIO::E_CIO_LITTLE;
  }else if( !strcasecmp(str.c_str(),"big" ) ) {
    Endian=CIO::E_CIO_BIG;
  }else {
    printf("\tCIO Parsing error : fail to get '%s'\n",str.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_ENDIAN;
  } 

  ncnt++;

  //ArrayShape  
  label = "/FileInfo/ArrayShape";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_ARRAYSHAPE;
  }
  if( !strcasecmp(str.c_str(),"ijkn" ) ) {
    ArrayShape=CIO::E_CIO_IJKN;
  } else if( !strcasecmp(str.c_str(),"nijk" ) ) {
    ArrayShape=CIO::E_CIO_NIJK;
  }else {
    printf("\tCIO Parsing error : fail to get '%s'\n",str.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_ARRAYSHAPE;
  }

  ncnt++;

  //Componet  
  label = "/FileInfo/Component";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_COMPONENT;
  }
  Component=ct;

  ncnt++;

  //Component Variable
  int ncomp=0;
  label_leaf_leaf = "/FileInfo/Variable";
  if ( tpCntl.chkNode(label_leaf_leaf) )  //があれば
  {
    ncomp = tpCntl.countLabels(label_leaf_leaf);
  }

  ncnt++;

  label_leaf = "/FileInfo";

  if( ncomp>0 ) {
    for(int i=0; i<ncomp; i++) {
      if(!tpCntl.GetNodeStr(label_leaf,ncnt+i,&str))
      {
        printf("\tCIO Parsing error : No Elem name\n");
        return CIO::E_CIO_ERROR_READ_DFI_NO_MINMAX;
      }
      if( !strcasecmp(str.substr(0,8).c_str(), "variable") ) {
        label_leaf_leaf = label_leaf+"/"+str;

        label = label_leaf_leaf + "/name";
        if ( !(tpCntl.GetValue(label, &str )) ) {
          printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
          return CIO::E_CIO_ERROR_READ_DFI_MIN;
        }
        else {
          ComponentVariable.push_back(str);
        }
      }
    }
  }

//FCONV 20140131.s
  if( FileFormat == CIO::E_CIO_FMT_SPH ) {
    if( Component > 1 && ArrayShape == CIO::E_CIO_IJKN ) {
      printf("\tCIO error sph file undefined ijkn component>1.\n");
      return CIO::E_CIO_ERROR_READ_DFI_ARRAYSHAPE;
    }
  }
//FCONV 20140131.e

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// DFIファイル:FileInfo要素を出力する
CIO::E_CIO_ERRORCODE
cio_FileInfo::Write(FILE* fp, 
                    const unsigned tab) 
{

  fprintf(fp, "FileInfo {\n");
  fprintf(fp, "\n");

//FCONV 20140116.s
  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "DFIType            = \"Cartesian\"\n");
//FCONV 20140116.s

  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "DirectoryPath      = \"%s\"\n", DirectoryPath.c_str());

  _CIO_WRITE_TAB(fp, tab+1);
  if(        TimeSliceDirFlag == CIO::E_CIO_OFF ) {
    fprintf(fp, "TimeSliceDirectory = \"off\"\n");
  } else if( TimeSliceDirFlag == CIO::E_CIO_ON ) {
    fprintf(fp, "TimeSliceDirectory = \"on\"\n");
  }

  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "Prefix             = \"%s\"\n", Prefix.c_str());

  _CIO_WRITE_TAB(fp, tab+1);
  if(        FileFormat == CIO::E_CIO_FMT_SPH ) {
    fprintf(fp, "FileFormat         = \"sph\"\n");
  } else if( FileFormat == CIO::E_CIO_FMT_BOV ) {
    fprintf(fp, "FileFormat         = \"bov\"\n");
  }

//FCONV 20140116.s
  _CIO_WRITE_TAB(fp, tab+1);
  if( FieldFilenameFormat == CIO::E_CIO_FNAME_RANK_STEP ) {
    fprintf(fp, "FieldFilenameFormat= \"rank_step\"\n");
  } else {
    fprintf(fp, "FieldFilenameFormat= \"step_rank\"\n");
  }
//FCONV 20140116.e

  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "GuideCell          = %d\n", GuideCell);

  _CIO_WRITE_TAB(fp, tab+1);
  std::string Dtype = cio_DFI::ConvDatatypeE2S((CIO::E_CIO_DTYPE)DataType);
  fprintf(fp, "DataType           = \"%s\"\n",Dtype.c_str());

  _CIO_WRITE_TAB(fp, tab+1);
  if( Endian == CIO::E_CIO_LITTLE) {
    fprintf(fp, "Endian             = \"little\"\n");
  } else {
    fprintf(fp, "Endian             = \"big\"\n");
  }

  _CIO_WRITE_TAB(fp, tab+1);
  if( ArrayShape == CIO::E_CIO_IJKN ) {
    fprintf(fp, "ArrayShape         = \"ijkn\"\n");
  } else {
    fprintf(fp, "ArrayShape         = \"nijk\"\n");
  }

  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "Component          = %d\n",Component);

/*
  if( ComponentVariable.size()>0 ) {
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",ComponentVariable[0].c_str());
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",ComponentVariable[1].c_str());
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",ComponentVariable[2].c_str());
  }
*/
  for(int i=0; i<ComponentVariable.size(); i++) {
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Variable[@]{ name  = \"%s\" }\n",ComponentVariable[i].c_str());
  }

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CIO::E_CIO_SUCCESS;

}

