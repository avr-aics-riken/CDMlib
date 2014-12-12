#ifndef _CONV_PLOT3D_INLINE_H_
#define _CONV_PLOT3D_INLINE_H_

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
 * @file   conv_plot3d_inline.h
 * @brief  convOutput_PLOT3D クラスのinline関数ヘッダーファイル
 * @author aics
 * @date   2013/11/7
 */

#ifndef CONV_NO_INLINE
 #define CONV_INLINE inline
#else
 #define CONV_INLINE
#endif

/**
 * @brief xyzファイルの出力
 * @param [in] prefix ファイル接頭文字 
 * @param [in] step   ステップ
 * @param [in] rank   ランク
 * @param [in] guide  ガイドセル数
 * @param [in] origin 基点座標
 * @param [in] pitch  ピッチ
 * @param [in] size   セルサイズ
 * @param [in] x      x方向座標ワーク
 * @param [in] y      y方向座標ワーク
 * @param [in] z      z方向座標ワーク
 */
template<class T1, class T2>
CONV_INLINE
void 
convOutput_PLOT3D::OutputPlot3D_xyz(std::string prefix,
                                    int step, 
                                    int rank, 
                                    int guide, 
                                    T1* origin, 
                                    T1* pitch, 
                                    int* size, 
                                    T2* x, 
                                    T2* y, 
                                    T2* z)
{
  //value
  int ngrid=1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int gc_out = 0;//plot3dは常にガイドセルを出力しない

  //間引き数の取得
  int thin_count = m_InputCntl->Get_ThinOut();
  
  int *iblank=NULL;//dummy
  int id,jd,kd;//出力サイズ
  id=size[0]+1;//+2*gc_out
  jd=size[1]+1;//+2*gc_out
  kd=size[2]+1;//+2*gc_out
  
  //間引きのための処理
  int irest=(id-1)%thin_count;
  int jrest=(jd-1)%thin_count;
  int krest=(kd-1)%thin_count;
  id=(id-1)/thin_count;
  jd=(jd-1)/thin_count;
  kd=(kd-1)/thin_count;
  //id=id+1;
  //jd=jd+1;
  //kd=kd+1;
  //格子点への補間は行わず、双対セルとして扱うため、+1は不要
  if(irest!=0) id=id+1;
  if(jrest!=0) jd=jd+1;
  if(krest!=0) kd=kd+1;
  
  // ガイドセル出力があった場合オリジナルポイントを調整しておく
  T2 l_org[3], l_pit[3];
  for (int i=0; i<3; i++)
  {
    l_org[i] = (T2)origin[i] + (T2)pitch[i]*(T2)gd;
    l_pit[i] = (T2)pitch[i];
  }
    
  // 出力ファイル名
  std::string tmp;
  //std::string t_prefix=prefix+"_Grid";
  int fnameformat = m_InputCntl->Get_OutputFilenameFormat();
  tmp = m_InputCntl->Get_OutputDir() +"/"+ 
        cdm_DFI::Generate_FileName(prefix,
                                   rank,
                                   //step,
                                   -1,
                                   "xyz",
                                   (CDM::E_CDM_OUTPUT_FNAME)fnameformat,
                                   false,
                                   CDM::E_CDM_OFF);
    
  //open file
  FILE*fp;
  if( m_InputCntl->Get_OutputFileType() == CDM::E_CDM_FILE_TYPE_ASCII ) {
    if( (fp = fopen(tmp.c_str(), "wa")) == NULL ) {
      printf("\tCan't open file.(%s)\n",tmp.c_str());
      Exit(0);
    }
  } else {
    if( (fp = fopen(tmp.c_str(), "wb")) == NULL ) {
      printf("\tCan't open file.(%s)\n",tmp.c_str());
      Exit(0);
    }
  }
    
  //write block data
  //WriteNgrid(fp,ngrid);//if multi grid
  WriteBlockData(fp,id,jd,kd);
    
  for(int k=0;k<kd;k++){
  for(int j=0;j<jd;j++){
  for(int i=0;i<id;i++){
    size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
    x[ip]=l_org[0]+(T2)thin_count*l_pit[0]*(T2)i;//-pitch[0]*(float)gc_out;
    y[ip]=l_org[1]+(T2)thin_count*l_pit[1]*(T2)j;//-pitch[1]*(float)gc_out;
    z[ip]=l_org[2]+(T2)thin_count*l_pit[2]*(T2)k;//-pitch[2]*(float)gc_out;
  }}}
    
  //x direction modify
  if(irest!=0 && (id-2)>=0 ){
    for(int k=0;k<kd;k++){
    for(int j=0;j<jd;j++){
      size_t ip = _F_IDX_S3D(id, j+1, k+1, id, jd, kd, 0);
      x[ip]=l_org[0]+(T2)thin_count*l_pit[0]*(T2)(id-2)+(T2)irest*l_pit[0];//-pitch[0]*(float)gc_out;
    }}
  }
    
  //y direction modify
  if(jrest!=0 && (jd-2)>=0 ){
    for(int k=0;k<kd;k++){
    for(int i=0;i<id;i++){
      size_t ip = _F_IDX_S3D(i+1, jd, k+1, id, jd, kd, 0);
      y[ip]=l_org[1]+(T2)thin_count*l_pit[1]*(T2)(jd-2)+(T2)jrest*l_pit[1];//-pitch[1]*(float)gc_out;
    }}
  }
    
  //z direction modify
  if(krest!=0 && (kd-2)>=0 ){
    for(int j=0;j<jd;j++){
    for(int i=0;i<id;i++){
      size_t ip = _F_IDX_S3D(i+1, j+1, kd, id, jd, kd, 0);
      z[ip]=l_org[2]+(T2)thin_count*l_pit[2]*(T2)(kd-2)+(T2)krest*l_pit[2];//-pitch[2]*(float)gc_out;
    }}
  }
    
  //z direction modify
  if(krest!=0){
    for(int k=kd-1;k<kd;k++){
    for(int j=0;j<jd;j++){
    for(int i=0;i<id;i++){
      size_t ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      z[ip]=l_org[2]+(T2)krest*l_pit[2]*(T2)k;//-pitch[2]*(float)gc_out;
    }}}
  }
  
  //write
  if(!WriteXYZData(fp, id, jd, kd, ngrid, x, y, z)) printf("\terror WriteXYZData\n");
  
  //close file
  fclose(fp);
    
}

/**
 * @brief gridデータ出力
 * @param [in] fp         出力ファイルポインタ
 * @param [in] id         i方向サイズ
 * @param [in] jd         j方向のサイズ
 * @param [in] kd         k方向のサイズ
 * @param [in] ngrid      1
 * @param [in] x          x座標値
 * @param [in] y          y座標値
 * @param [in] z          z座標値
 */
template<class T>
CONV_INLINE
bool
convOutput_PLOT3D::WriteXYZData(FILE* fp,
                                int id,
                                int jd,
                                int kd,
                                int ngrid,
                                T*  x,
                                T*  y,
                                T*  z)
{

  size_t sz = (size_t)id*(size_t)jd*(size_t)kd;
  unsigned int dmy;

  int s12 =(size_t)id*(size_t)jd;
  int ns12=s12/10;

  switch (m_InputCntl->Get_OutputFileType()) {

    //Fortran Binary 出力
    case CDM::E_CDM_FILE_TYPE_FBINARY:
      dmy = sizeof(T)*sz*3;
      WriteDataMarker(dmy,fp,true);
      fwrite(x, sizeof(T), sz, fp);
      fwrite(y, sizeof(T), sz, fp);
      fwrite(z, sizeof(T), sz, fp);
      WriteDataMarker(dmy,fp,true);
      break;

    //ascii 出力
    case CDM::E_CDM_FILE_TYPE_ASCII:
      WriteXYZ_FORMATTED(fp, id, jd, kd, x);
      WriteXYZ_FORMATTED(fp, id, jd, kd, y);
      WriteXYZ_FORMATTED(fp, id, jd, kd, z);
      break;

    //C Binary 出力
    case CDM::E_CDM_FILE_TYPE_BINARY:
      fwrite(x, sizeof(T), sz, fp);
      fwrite(y, sizeof(T), sz, fp);
      fwrite(z, sizeof(T), sz, fp);
      break;
    default:
      return false;
      break;
  }
  return true;
}

/**
 * @brief Formatted 出力
 * @param [in] fp 出力ファイルポインタ
 * @param [in] id i方向サイズ
 * @param [in] jd j方向サイズ
 * @param [in] kd k方向サイズ
 * @param [in] x   出力座標値配列
 */
template<class T>
CONV_INLINE
void
convOutput_PLOT3D::WriteXYZ_FORMATTED(FILE* fp,
                                      int id,
                                      int jd,
                                      int kd,
                                      T* x)
{

  int s12 =(size_t)id*(size_t)jd;
  int ns12=s12/10;

/*  
  for(int k=0; k<kd; k++) {
    //x-y面の出力
    for(int i=0; i<ns12; i++) {
      for(int ii=0; ii<10; ii++) {
        fprintf(fp,"%15.6E",x[(k*id*jd)+(i*10)+ii]);
      }
      fprintf(fp,"\n");
    }
    //余りの出力
    if( s12%10 > 0 ) {
      for(int i=0; i<s12%10; i++) fprintf(fp,"%15.6E",x[k*id*jd+(ns12*10)+i]);
    }
    fprintf(fp,"\n");
   }
*/
   for(int i=0; i<id*jd*kd; i++) {
     fprintf(fp,"%15.6E\n",x[i]);
   }

}
 
#endif // _CONV_PLOT3D_INLINE_H_ 
