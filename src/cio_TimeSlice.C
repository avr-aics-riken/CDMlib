/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_TimeSlice.C 
 * @brief  cio_Slice Class
 * @author aics    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


// #################################################################
// コンストラクタ
cio_Slice::cio_Slice()
{

  step = 0;
  time = 0.0;
  AveragedStep = 0;
  AveragedTime = 0.0;
  VectorMin = 0.0;
  VectorMax = 0.0;

  Min.clear();
  Max.clear();

}


// #################################################################
// デストラクタ
cio_Slice::~cio_Slice()
{

}

// #################################################################
// TimeSliceの読込み
CIO::E_CIO_ERRORCODE
cio_Slice::Read(cio_TextParser tpCntl,
                        std::string label_leaf) 
{
#if 0
  std::string str;
  std::string label,label_leaf_leaf;

  int ct;
  double dt;

  int ncnt=0;

  //Step
  label = label_leaf + "/Step";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_STEP;
  }
  else {
    step=ct;
  }

  ncnt++;

  //Time
  label = label_leaf + "/Time";
  if ( !(tpCntl.GetValue(label, &dt )) ) {
    printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_TIME;
  }
  else {
    time= dt;
  }

  ncnt++;

  //AveragedStep
  label = label_leaf + "/AveragedStep";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    AveragedStep=-1;
  }
  else {
    AveragedStep= ct;
    ncnt++;
  }

  //AveragedTime
  label = label_leaf + "/AveragedTime";
  if ( !(tpCntl.GetValue(label, &dt )) ) {
    AveragedTime=0.0;
  }
  else {
    AveragedTime= dt;
    ncnt++;
  }

  //VectorMinMax/Min
  label = label_leaf + "/VectorMinMax/Min";
  if ( (tpCntl.GetValue(label, &dt )) )
  {
    VectorMin=dt;
    ncnt++;
  }

  //VectorMinMax/Max
  label = label_leaf + "/VectorMinMax/Max";
  if ( (tpCntl.GetValue(label, &dt )) )
  {
    VectorMax=dt;
  }

  //MinMax
  int ncomp=0;
  label_leaf_leaf = label_leaf + "/MinMax";
  if ( tpCntl.chkNode(label_leaf_leaf) )  //があれば
  {
    ncomp = tpCntl.countLabels(label_leaf_leaf);
  }

  ncnt++;

  Min.clear();
  Max.clear();

  for ( int j=0; j<ncomp; j++ ) {

    if(!tpCntl.GetNodeStr(label_leaf,j+ncnt,&str))
    {
      printf("\tCIO Parsing error : No Elem name\n");
      return CIO::E_CIO_ERROR_READ_DFI_NO_MINMAX;
    }
    if( strcasecmp(str.substr(0,6).c_str(), "minmax") ) continue;
    label_leaf_leaf = label_leaf+"/"+str;

    label = label_leaf_leaf + "/Min";
    if ( !(tpCntl.GetValue(label, &dt )) ) {
      printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
      return CIO::E_CIO_ERROR_READ_DFI_MIN;
    }
    else {
      Min.push_back(dt);
    }

    label = label_leaf_leaf + "/Max";
    if ( !(tpCntl.GetValue(label, &dt )) ) {
      printf("\tCIO Parsing error : fail to get '%s'\n",label.c_str());
      return CIO::E_CIO_ERROR_READ_DFI_MAX;
    }
    else {
      Max.push_back(dt);
    }

  }
#else
  std::string str;
  std::string label;
  int ct;
  double dt;

  // カレントを移動
  TextParser *tp = tpCntl.getTPPtr();
  tp->changeNode(label_leaf);

  //Step
  label = "Step";
  if ( !(tpCntl.GetValue(label, &ct, false )) ) {
    printf("\tCIO Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_STEP;
  }
  else {
    step=ct;
  }

  //Time
  label = "Time";
  if ( !(tpCntl.GetValue(label, &dt, false )) ) {
    printf("\tCIO Parsing error : fail to get '%s/%s'\n",label_leaf.c_str(),label.c_str());
    return CIO::E_CIO_ERROR_READ_DFI_TIME;
  }
  else {
    time= dt;
  }

  //AveragedStep
  label = "AveragedStep";
  if ( !(tpCntl.GetValue(label, &ct, false )) ) {
    AveragedStep=-1;
  }
  else {
    AveragedStep= ct;
  }

  //AveragedTime
  label = "AveragedTime";
  if ( !(tpCntl.GetValue(label, &dt, false )) ) {
    AveragedTime=0.0;
  }
  else {
    AveragedTime= dt;
  }

  //VectorMinMax/Min
  label = "VectorMinMax/Min";
  if ( (tpCntl.GetValue(label, &dt, false )) )
  {
    VectorMin=dt;
  }

  //VectorMinMax/Max
  label = "VectorMinMax/Max";
  if ( (tpCntl.GetValue(label, &dt, false )) )
  {
    VectorMax=dt;
  }

  // 子のラベルを取得
  vector<std::string> labels;
  tp->getNodes(labels,1);

  // 子のMinMaxを読み込み
  Min.clear();
  Max.clear();
  for( size_t i=0;i<labels.size();i++ )
  {
    // MinMax要素かどうか確認
    std::string label = labels[i];
    if( strcasecmp(label.substr(0,6).c_str(), "MinMax") ) continue;

    // MinMaxに移動
    std::string leaf = label_leaf + "/" + label;
    tp->changeNode(leaf);

    // Min
    label = "Min";
    if ( !(tpCntl.GetValue(label, &dt, false )) ) {
      printf("\tCIO Parsing error : fail to get '%s/%s'\n",leaf.c_str(),label.c_str());
      return CIO::E_CIO_ERROR_READ_DFI_MIN;
    }
    else {
      Min.push_back(dt);
    }

    label = "Max";
    if ( !(tpCntl.GetValue(label, &dt, false )) ) {
      printf("\tCIO Parsing error : fail to get '%s/%s'\n",leaf.c_str(),label.c_str());
      return CIO::E_CIO_ERROR_READ_DFI_MAX;
    }
    else {
      Max.push_back(dt);
    }
  }
#endif

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// TimeSliceを出力する 
CIO::E_CIO_ERRORCODE
cio_Slice::Write(FILE* fp,
                 const unsigned tab)
{

  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "Step = %u\n",step);

  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "Time = %e\n",time);

  if( !avr_mode ) {
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "AveragedStep = %u\n",AveragedStep);
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "AveragedTime = %e\n",AveragedTime);
  }

  if( Min.size()>1 ) {
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "VectorMinMax {\n");
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Min = %e\n",VectorMin);
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Max = %e\n",VectorMax);
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "}\n");
  }

  for(int j=0; j<Min.size(); j++){
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "MinMax[@] {\n");
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Min = %e\n",Min[j]);
    _CIO_WRITE_TAB(fp, tab+1);
    fprintf(fp, "Max = %e\n",Max[j]);
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "}\n");
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// コンストラクタ
cio_TimeSlice::cio_TimeSlice()
{
  SliceList.clear();
} 

// #################################################################
// デストラクタ
cio_TimeSlice::~cio_TimeSlice()
{

}

// #################################################################
// TimeSliceの読込み
CIO::E_CIO_ERRORCODE
cio_TimeSlice::Read(cio_TextParser tpCntl)
{
#if 0
  std::string str;
  std::string label_base,label_leaf;

  cio_Slice slice;

  int nnode=0;

  CIO::E_CIO_ERRORCODE iret;

  //TimeSlice
  nnode=0;
  label_base = "/TimeSlice";
  if ( tpCntl.chkNode(label_base) )  //があれば
  {
    nnode = tpCntl.countLabels(label_base);
  }

  for (int i=0; i<nnode; i++) {

    int ncnt=0;

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tCIO Parsing error : No Elem name\n");
      return CIO::E_CIO_ERROR_READ_DFI_NO_SLICE;
    }
    if( strcasecmp(str.substr(0,5).c_str(), "Slice") ) continue;
    label_leaf=label_base+"/"+str;

    //Slice要素の読込み
    iret = slice.Read(tpCntl,label_leaf);

    if( iret == CIO::E_CIO_SUCCESS ) {
      SliceList.push_back(slice); 
    } else return iret;

  }
#else
  CIO::E_CIO_ERRORCODE iret;

  // TP
  TextParser *tp = tpCntl.getTPPtr();
  if( !tp )
  {
    return CIO::E_CIO_ERROR_TEXTPARSER;
  }

  // TimeSlice要素の存在チェック
  std::string label_base = "/TimeSlice";
  if( !tpCntl.chkNode(label_base) )
  {
    printf("\tCIO Parsing error : No Elem name [%s]\n", label_base.c_str());
    return CIO::E_CIO_ERROR_READ_TIMESLICE;
  }

  // TimeSliceに移動
  tp->changeNode(label_base);

  // 子のラベルを取得
  vector<std::string> labels;
  tp->getNodes(labels,1);

  // 子のSliceを読み込み
  for( size_t i=0;i<labels.size();i++ )
  {
    // Slice要素かどうか確認
    std::string label = labels[i];
    if( strcasecmp(label.substr(0,5).c_str(), "Slice") ) continue;

    //Slice要素の読込み
    cio_Slice slice;
    std::string leaf = label_base + "/" + label;
    if( (iret = slice.Read(tpCntl,leaf)) == CIO::E_CIO_SUCCESS )
    {
      SliceList.push_back(slice);
    }
    else
    {
      return iret;
    }
  }

  // 元のノードに戻る
  tp->changeNode("/");
#endif

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// TimeSliceを出力する 
CIO::E_CIO_ERRORCODE
cio_TimeSlice::Write(FILE* fp,
                     const unsigned tab)
{

  fprintf(fp, "TimeSlice {\n");
  fprintf(fp, "\n");

  for(int i=0; i<SliceList.size(); i++) {

    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "Slice[@] {\n");

    //Slice要素の出力
    if( SliceList[i].Write(fp,tab+1) != CIO::E_CIO_SUCCESS) return CIO::E_CIO_ERROR;

    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "}\n");

  }

  fprintf(fp, "}\n\n");
  fclose(fp);

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// DFIに出力されているminmaxの合成値を取得
CIO::E_CIO_ERRORCODE 
cio_TimeSlice::getVectorMinMax(const unsigned step,
                               double &vec_min,
                               double &vec_max)
{
  for(int i=0;SliceList.size(); i++) {
    if( (int)step == SliceList[i].step ) {
      vec_min=SliceList[i].VectorMin;
      vec_max=SliceList[i].VectorMax;
      return CIO::E_CIO_SUCCESS;
    }
  }

  return CIO::E_CIO_ERROR;
}

// #################################################################
// DFIに出力されているminmaxとminmaxの合成値を取得
CIO::E_CIO_ERRORCODE 
cio_TimeSlice::getMinMax(const unsigned step,
                         const int compNo,
                         double &min_value,
                         double &max_value)
{

  for(int i=0;SliceList.size(); i++) {
    if( (int)step == SliceList[i].step ) {
      min_value=SliceList[i].Min[compNo];
      max_value=SliceList[i].Max[compNo];
      return CIO::E_CIO_SUCCESS;
    }
  }

  return CIO::E_CIO_ERROR;

}
// #################################################################
// SliceListへの追加
void cio_TimeSlice::AddSlice(int step,
                             double time,
                             double *minmax,
                             int Ncomp,
                             bool avr_mode,
                             int step_avr,
                             double time_avr)
{

  cio_Slice slice;

  slice.step = step;
  slice.time = time;

  //minmaxのセット
  if( minmax ) {
    //成分が１個の場合
    if( Ncomp == 1 ) {
      slice.Min.push_back(minmax[0]);
      slice.Max.push_back(minmax[1]);
    } else {
    //成分が複数個の場合
      for(int i=0; i<Ncomp; i++) {
        slice.Min.push_back(minmax[i*2]);
        slice.Max.push_back(minmax[i*2+1]);
      }
      slice.VectorMin=minmax[6];
      slice.VectorMax=minmax[7];
    }
  }

  //averageのセット
  slice.avr_mode = avr_mode;
  if( !avr_mode ) {
    slice.AveragedStep=step_avr;
    slice.AveragedTime=time_avr;
  } else {
    slice.AveragedStep=0;
    slice.AveragedTime=0.0;
  }

  SliceList.push_back(slice);

}
