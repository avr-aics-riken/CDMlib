/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_MPI.C
 * @brief  cio_MPI Class
 * @author aics    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cio_MPI::cio_MPI()
{
  NumberOfRank=0;
  NumberOfGroup=1;
}

// #################################################################
// コンストラクタ
cio_MPI::cio_MPI(const int _NumberOfRank, int _NumberOfGroup)
{
  NumberOfRank=_NumberOfRank;
  NumberOfGroup=_NumberOfGroup;
}

// #################################################################
// デストラクタ
cio_MPI::~cio_MPI()
{

}

// #################################################################
// DFIファイル：MPI要素の読込み
CIO::E_CIO_ERRORCODE
cio_MPI::Read(cio_TextParser tpCntl, 
              const cio_Domain domain) 
{

  std::string str;
  std::string label;
  int ct;

  //NumberOfRank
  label = "/MPI/NumberOfRank";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    ct = domain.GlobalDivision[0]*domain.GlobalDivision[1]*domain.GlobalDivision[2];
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

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// DFIファイル:MPI要素を出力する
CIO::E_CIO_ERRORCODE
cio_MPI::Write(FILE* fp, const unsigned tab)
{

  fprintf(fp, "MPI {\n");
  fprintf(fp, "\n");

  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "NumberOfRank   = %d\n", NumberOfRank);

  _CIO_WRITE_TAB(fp, tab+1);
  fprintf(fp, "NumberOfGroup  = %d\n", 1);  

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CIO::E_CIO_SUCCESS;

}

