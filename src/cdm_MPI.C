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
 * @file   cdm_MPI.C
 * @brief  cdm_MPI Class
 * @author aics
 */

#include "cdm_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cdm_MPI::cdm_MPI()
{
  NumberOfRank=0;
  NumberOfGroup=1;
}

// #################################################################
// コンストラクタ
cdm_MPI::cdm_MPI(const int _NumberOfRank, int _NumberOfGroup)
{
  NumberOfRank=_NumberOfRank;
  NumberOfGroup=_NumberOfGroup;
}

// #################################################################
// デストラクタ
cdm_MPI::~cdm_MPI()
{

}

// #################################################################
// DFIファイル：MPI要素の読込み
CDM::E_CDM_ERRORCODE
cdm_MPI::Read(cdm_TextParser tpCntl,
              const cdm_Domain* domain)
{

  std::string str;
  std::string label;
  int ct;

  //NumberOfRank
  label = "/MPI/NumberOfRank";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    ct = domain->GlobalDivision[0]*domain->GlobalDivision[1]*domain->GlobalDivision[2];
  }
  else {
    NumberOfRank = ct;
  }

  //NumberOfGroup
  label = "/MPI/NumberOfGroup";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    ct = 1;
  }
  else {
    NumberOfGroup = ct;
  }

  return CDM::E_CDM_SUCCESS;

}

// #################################################################
// DFIファイル:MPI要素を出力する
CDM::E_CDM_ERRORCODE
cdm_MPI::Write(FILE* fp, const unsigned tab)
{

  fprintf(fp, "MPI {\n");
  fprintf(fp, "\n");

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "NumberOfRank   = %d\n", NumberOfRank);

  _CDM_WRITE_TAB(fp, tab+1);
  fprintf(fp, "NumberOfGroup  = %d\n", 1);

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CDM::E_CDM_SUCCESS;

}
