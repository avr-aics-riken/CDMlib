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

#include "fub_DATA.h"

// #################################################################
// コンストラクタ
fub_DATA::fub_DATA(string fubCname, string fubFname)
{

  size[0]=size[1]=size[2]=0;

  xyz=NULL;
  f_xyz=NULL;

  for (int i=0; i<6; i++) nID[i]=NULL;

  fname_Coord = fubCname;
  fname_Field = fubFname;

  endian="";

  fileID=0;

}

// #################################################################
// エンディアンチェック
bool fub_DATA::isMatchEndianSbdmMagick(int ident)
{

/*
  char magick_c[] = "SBDM";
  int  magick_i=0;

  //check match
  magick_i = (magick_c[3]<<24) + (magick_c[2]<<16) + (magick_c[1]<<8) + magick_c[0];
  if( magick_i == ident )
  {
    return true;
  }
*/

  if( ident < 0 || ident > 100 )
  {
    BSWAP32(ident);
    if( ident > 0 && ident < 100 ) return false;
  }

  return true;

}

// #################################################################
// read fub file
bool fub_DATA::readFub(int vc, string &dftype)
{

  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  if( cdumy[0] == 0x01 ) endian="little";
  if( cdumy[0] == 0x00 ) endian="big";

  string fname = fname_Coord;

  FILE *fp;
  // fubファイルオープン
  if( !(fp=fopen(fname.c_str(),"rb")) ) {
    printf("Error open fname : %s\n",fname_Coord.c_str());
    return false;
  }

  unsigned int dmy;

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return false;
  }

  bool echeck = isMatchEndianSbdmMagick(dmy);

  if( !echeck ) {
    if( endian == "little" ) {
      endian = "big";
    } else {
      endian = "little";
    }
  }

  if( !echeck ) BSWAP32(dmy);

  if( (fread(size, sizeof(int), 4 , fp)) == 4 ) {
    if( !echeck ) {
      BSWAP32(size[0]);
      BSWAP32(size[1]);
      BSWAP32(size[2]);
      BSWAP32(size[3]);
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return false;
  }
  if( !echeck ) BSWAP32(dmy);

  int ix,iy,iz,n;
  ix=size[0];
  iy=size[1];
  iz=size[2];
  n =size[3];

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return 0;
  }
  if( !echeck ) BSWAP32(dmy);

  int Leng = ix*iy*iz*n;
  if( (int)dmy/Leng == 4 ) {
    f_xyz = new float[Leng];
    dftype = "Float32";
  } else {
    xyz = new double[Leng];
    dftype = "Float64";
  }

  int ret;
  if( dftype == "Float32" ) {

    if( (ret=fread(f_xyz, sizeof(float),ix*iy*iz*n, fp)) != ix*iy*iz*n ) {
      printf("Error read float data ret : %d\n",ret);
      return 0;
    }
    if( !echeck ) DBSWAPVEC(f_xyz,(ix*iy*iz*n));

  } else if( dftype == "Float64" ) {

    if( (ret=fread(xyz, sizeof(double),ix*iy*iz*n, fp)) != ix*iy*iz*n ) {
      printf("Error read double data ret : %d\n",ret);
      return 0;
    }
    if( !echeck ) DBSWAPVEC(xyz,(ix*iy*iz*n));

  }

  //読み込んだブロックのコーナー点をセット
  setCorner(vc, dftype);

  //読み込んだブロックの最小/最大座標値をセット（実セル、実ノード）
  setMinMax(vc, dftype);

  return true;

}

// #################################################################
// フィールドデータの成分数取得
int fub_DATA::getNval()
{
  int nVal=0;

  string fname = fname_Field;

  FILE *fp;
  //file open
  if( !(fp=fopen(fname.c_str(),"rb")) ) {
    printf("Error open fname : %s\n",fname_Coord.c_str());
    return false;
  }

  unsigned int dmy;

  //fortran record の読み込み
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) {
    fclose(fp);
    return false;
  }

  //endian 判定
  bool echeck = isMatchEndianSbdmMagick(dmy);

  //ensdian 変換
  if( !echeck ) BSWAP32(dmy);

  //data数
  int n = dmy/sizeof(dmy);

  //dataの読み込み
  int ndata[n];
  if( (fread(ndata, sizeof(dmy), n , fp)) == n ) {
    if( !echeck ) {
      for(int i=0; i<n; i++) BSWAP32(ndata[i]);
    }
    nVal = ndata[n-1];
  }

  fclose(fp);
  return nVal;
}

// #################################################################
// コーナー座標のセット
void fub_DATA::setCorner(int vc, string dftype)
{

   int NI,NJ,NK;
   NI=size[0]-2*vc;
   NJ=size[1]-2*vc;
   NK=size[2]-2*vc;

   int iend,jend,kend;
   if( vc>0 ) {
     iend = NI;
     jend = NJ;
     kend = NK;
   } else if( vc == 0 ) {
     iend = NI-1;
     jend = NJ-1;
     kend = NK-1;
   }

   if( dftype == "Float64" ) {
     //xyz(0,0,0)
     corner[_P000].x = xyz[_IDX_S4D(0,0,0,0,NI,NJ,NK,vc)];
     corner[_P000].y = xyz[_IDX_S4D(0,0,0,1,NI,NJ,NK,vc)];
     corner[_P000].z = xyz[_IDX_S4D(0,0,0,2,NI,NJ,NK,vc)];

     //xyz(iend,0,0)
     corner[_P100].x = xyz[_IDX_S4D(iend,0,0,0,NI,NJ,NK,vc)];
     corner[_P100].y = xyz[_IDX_S4D(iend,0,0,1,NI,NJ,NK,vc)];
     corner[_P100].z = xyz[_IDX_S4D(iend,0,0,2,NI,NJ,NK,vc)];

     //xyz(iend,jend,0)
     corner[_P110].x = xyz[_IDX_S4D(iend,jend,0,0,NI,NJ,NK,vc)];
     corner[_P110].y = xyz[_IDX_S4D(iend,jend,0,1,NI,NJ,NK,vc)];
     corner[_P110].z = xyz[_IDX_S4D(iend,jend,0,2,NI,NJ,NK,vc)];

     //xyz(0,jend,0)
     corner[_P010].x = xyz[_IDX_S4D(0,jend,0,0,NI,NJ,NK,vc)];
     corner[_P010].y = xyz[_IDX_S4D(0,jend,0,1,NI,NJ,NK,vc)];
     corner[_P010].z = xyz[_IDX_S4D(0,jend,0,2,NI,NJ,NK,vc)];

     //xyz(0,0,kend)
     corner[_P001].x = xyz[_IDX_S4D(0,0,kend,0,NI,NJ,NK,vc)];
     corner[_P001].y = xyz[_IDX_S4D(0,0,kend,1,NI,NJ,NK,vc)];
     corner[_P001].z = xyz[_IDX_S4D(0,0,kend,2,NI,NJ,NK,vc)];

     //xyz(iend,0,kend)
     corner[_P101].x = xyz[_IDX_S4D(iend,0,kend,0,NI,NJ,NK,vc)];
     corner[_P101].y = xyz[_IDX_S4D(iend,0,kend,1,NI,NJ,NK,vc)];
     corner[_P101].z = xyz[_IDX_S4D(iend,0,kend,2,NI,NJ,NK,vc)];

     //xyz(iend,jend,kend)
     corner[_P111].x = xyz[_IDX_S4D(iend,jend,kend,0,NI,NJ,NK,vc)];
     corner[_P111].y = xyz[_IDX_S4D(iend,jend,kend,1,NI,NJ,NK,vc)];
     corner[_P111].z = xyz[_IDX_S4D(iend,jend,kend,2,NI,NJ,NK,vc)];

     //xyz(0,jend,kend)
     corner[_P011].x = xyz[_IDX_S4D(0,jend,kend,0,NI,NJ,NK,vc)];
     corner[_P011].y = xyz[_IDX_S4D(0,jend,kend,1,NI,NJ,NK,vc)];
     corner[_P011].z = xyz[_IDX_S4D(0,jend,kend,2,NI,NJ,NK,vc)];

   } else if( dftype == "Float32" ) {
     //xyz(0,0,0)
     corner[_P000].x = (double)f_xyz[_IDX_S4D(0,0,0,0,NI,NJ,NK,vc)];
     corner[_P000].y = (double)f_xyz[_IDX_S4D(0,0,0,1,NI,NJ,NK,vc)];
     corner[_P000].z = (double)f_xyz[_IDX_S4D(0,0,0,2,NI,NJ,NK,vc)];

     //xyz(iend,0,0)
     corner[_P100].x = (double)f_xyz[_IDX_S4D(iend,0,0,0,NI,NJ,NK,vc)];
     corner[_P100].y = (double)f_xyz[_IDX_S4D(iend,0,0,1,NI,NJ,NK,vc)];
     corner[_P100].z = (double)f_xyz[_IDX_S4D(iend,0,0,2,NI,NJ,NK,vc)];

     //xyz(iend,jend,0)
     corner[_P110].x = (double)f_xyz[_IDX_S4D(iend,jend,0,0,NI,NJ,NK,vc)];
     corner[_P110].y = (double)f_xyz[_IDX_S4D(iend,jend,0,1,NI,NJ,NK,vc)];
     corner[_P110].z = (double)f_xyz[_IDX_S4D(iend,jend,0,2,NI,NJ,NK,vc)];

     //xyz(0,jend,0)
     corner[_P010].x = (double)f_xyz[_IDX_S4D(0,jend,0,0,NI,NJ,NK,vc)];
     corner[_P010].y = (double)f_xyz[_IDX_S4D(0,jend,0,1,NI,NJ,NK,vc)];
     corner[_P010].z = (double)f_xyz[_IDX_S4D(0,jend,0,2,NI,NJ,NK,vc)];

     //xyz(0,0,kend)
     corner[_P001].x = (double)f_xyz[_IDX_S4D(0,0,kend,0,NI,NJ,NK,vc)];
     corner[_P001].y = (double)f_xyz[_IDX_S4D(0,0,kend,1,NI,NJ,NK,vc)];
     corner[_P001].z = (double)f_xyz[_IDX_S4D(0,0,kend,2,NI,NJ,NK,vc)];

     //xyz(iend,0,kend)
     corner[_P101].x = (double)f_xyz[_IDX_S4D(iend,0,kend,0,NI,NJ,NK,vc)];
     corner[_P101].y = (double)f_xyz[_IDX_S4D(iend,0,kend,1,NI,NJ,NK,vc)];
     corner[_P101].z = (double)f_xyz[_IDX_S4D(iend,0,kend,2,NI,NJ,NK,vc)];

     //xyz(iend,jend,kend)
     corner[_P111].x = (double)f_xyz[_IDX_S4D(iend,jend,kend,0,NI,NJ,NK,vc)];
     corner[_P111].y = (double)f_xyz[_IDX_S4D(iend,jend,kend,1,NI,NJ,NK,vc)];
     corner[_P111].z = (double)f_xyz[_IDX_S4D(iend,jend,kend,2,NI,NJ,NK,vc)];

     //xyz(0,jend,kend)
     corner[_P011].x = (double)f_xyz[_IDX_S4D(0,jend,kend,0,NI,NJ,NK,vc)];
     corner[_P011].y = (double)f_xyz[_IDX_S4D(0,jend,kend,1,NI,NJ,NK,vc)];
     corner[_P011].z = (double)f_xyz[_IDX_S4D(0,jend,kend,2,NI,NJ,NK,vc)];
   }

   return;
}

// #################################################################
// コーナー座標のセット
void fub_DATA::setMinMax(int vc, string dftype)
{

   int NI,NJ,NK;
   NI=size[0]-2*vc;
   NJ=size[1]-2*vc;
   NK=size[2]-2*vc;

   //min,maxの初期値セット
   if( dftype == "Float64" ) {

     min[0] = xyz[_IDX_S4D(0,0,0,0,NI,NJ,NK,vc)];
     min[1] = xyz[_IDX_S4D(0,0,0,1,NI,NJ,NK,vc)];
     min[2] = xyz[_IDX_S4D(0,0,0,2,NI,NJ,NK,vc)];
     max[0] = min[0];
     max[1] = min[1];
     max[2] = min[2];

     for(int k=0; k<NK; k++) {
     for(int j=0; j<NJ; j++) {
     for(int i=0; i<NI; i++) {
       if( min[0]>xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)] ) min[0]=xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)];
       if( max[0]<xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)] ) max[0]=xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)];
       if( min[1]>xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)] ) min[1]=xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)];
       if( max[1]<xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)] ) max[1]=xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)];
       if( min[2]>xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)] ) min[2]=xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)];
       if( max[2]<xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)] ) max[2]=xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)];
     }}}

   } else if( dftype == "Float32") {

     min[0] = (double)f_xyz[_IDX_S4D(0,0,0,0,NI,NJ,NK,vc)];
     min[1] = (double)f_xyz[_IDX_S4D(0,0,0,1,NI,NJ,NK,vc)];
     min[2] = (double)f_xyz[_IDX_S4D(0,0,0,2,NI,NJ,NK,vc)];
     max[0] = min[0];
     max[1] = min[1];
     max[2] = min[2];

     for(int k=0; k<NK; k++) {
     for(int j=0; j<NJ; j++) {
     for(int i=0; i<NI; i++) {
       if( min[0]>(double)f_xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)] ) min[0]=(double)f_xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)];
       if( max[0]<(double)f_xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)] ) max[0]=(double)f_xyz[_IDX_S4D(i,j,k,0,NI,NJ,NK,vc)];
       if( min[1]>(double)f_xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)] ) min[1]=(double)f_xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)];
       if( max[1]<(double)f_xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)] ) max[1]=(double)f_xyz[_IDX_S4D(i,j,k,1,NI,NJ,NK,vc)];
       if( min[2]>(double)f_xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)] ) min[2]=(double)f_xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)];
       if( max[2]<(double)f_xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)] ) max[2]=(double)f_xyz[_IDX_S4D(i,j,k,2,NI,NJ,NK,vc)];
     }}}

   }

   return;

}
// #################################################################
// 隣接ノード検索
void fub_DATA::checkNeighbor(fub_DATA *data2, int &cnt, double tol)
{
  // -X方向
  if( nID[X_MINUS]==NULL && data2->nID[X_PLUS]==NULL) {
    if( checkNeighborMX(data2,tol) ) {
      cnt++;
      return;
    }
  }

  // +X方向
  if( nID[X_PLUS]==NULL && data2->nID[X_MINUS]==NULL) {
    if( checkNeighborPX(data2,tol) ) {
      cnt++;
      return;
    }
  }

  // -Y方向
  if( nID[Y_MINUS]==NULL && data2->nID[Y_PLUS]==NULL) {
    if( checkNeighborMY(data2,tol) ) {
      cnt++;
      return;
    }
  }

  // +Y方向
  if( nID[Y_PLUS]==NULL && data2->nID[Y_MINUS]==NULL) {
    if( checkNeighborPY(data2,tol) ) {
      cnt++;
      return;
    }
  }

  // -Z方向
  if( nID[Z_MINUS]==NULL && data2->nID[Z_PLUS]==NULL) {
    if( checkNeighborMZ(data2,tol) ) {
      cnt++;
      return;
    }
  }

  // +Z方向
  if( nID[Z_PLUS]==NULL && data2->nID[Z_MINUS]==NULL) {
    if( checkNeighborPZ(data2,tol) ) {
      cnt++;
      return;
    }
  }

  return;
}

// #################################################################
// -X方向ノード検索
bool fub_DATA::checkNeighborMX(fub_DATA *data2, double tol)
{
  if( Coord::Same(corner[_P000], data2->corner[_P100], tol) ){
    this->nID[X_MINUS]=data2;
    data2->nID[X_PLUS]=this;
    return true;
  }
  return false;
}

// #################################################################
// +X方向ノード検索
bool fub_DATA::checkNeighborPX(fub_DATA *data2, double tol)
{
  if( Coord::Same(data2->corner[_P000], corner[_P100], tol) ) {
    this->nID[X_PLUS]=data2;
    data2->nID[X_MINUS]=this;
    return true;
  }
  return false;
}

// #################################################################
// -Y方向ノード検索
bool fub_DATA::checkNeighborMY(fub_DATA *data2, double tol)
{
  if( Coord::Same(corner[_P000], data2->corner[_P010], tol) ) {
    this->nID[Y_MINUS]=data2;
    data2->nID[Y_PLUS]=this;
    return true;
  }
  return false;
}

// #################################################################
// +Y方向ノード検索
bool fub_DATA::checkNeighborPY(fub_DATA *data2, double tol)
{
  if( Coord::Same(data2->corner[_P000], corner[_P010], tol) ) {
    this->nID[Y_PLUS]=data2;
    data2->nID[Y_MINUS]=this;
    return true;
  }
  return false;
}

// #################################################################
// -Z方向ノード検索
bool fub_DATA::checkNeighborMZ(fub_DATA *data2, double tol)
{
  if( Coord::Same(corner[_P000], data2->corner[_P001], tol) ) {
    this->nID[Z_MINUS]=data2;
    data2->nID[Z_PLUS]=this;
    return true;
  }
  return false;
}

// #################################################################
// +Z方向ノード検索
bool fub_DATA::checkNeighborPZ(fub_DATA *data2, double tol)
{
  if( Coord::Same(data2->corner[_P000], corner[_P001], tol) ) {
    this->nID[Z_PLUS]=data2;
    data2->nID[Z_MINUS]=this;
    return true;
  }
  return false;
}


// #################################################################
// ノードマップの生成
void fub_DATA::setFubDataMap(int ipos, int jpos, int kpos, int idiv,
         int jdiv, int kdiv, fub_DATA **fubDataMap)
{

  fub_DATA *root = fubDataMap[_IDX_S3D(ipos,jpos,kpos,idiv,jdiv,kdiv,0)];

  // +X
  if( ipos < idiv && root->nID[X_PLUS] ) {
    if( !fubDataMap[_IDX_S3D(ipos+1,jpos,kpos,idiv,jdiv,kdiv,0)] ) {

      fubDataMap[_IDX_S3D(ipos+1,jpos,kpos,idiv,jdiv,kdiv,0)] = root->nID[X_PLUS];

      setFubDataMap(ipos+1,jpos,kpos,idiv,jdiv,kdiv,fubDataMap);

    }
  }

  // +Y
  if( jpos < jdiv && root->nID[Y_PLUS] ) {
    if( !fubDataMap[_IDX_S3D(ipos,jpos+1,kpos,idiv,jdiv,kdiv,0)] ) {

      fubDataMap[_IDX_S3D(ipos,jpos+1,kpos,idiv,jdiv,kdiv,0)] = root->nID[Y_PLUS];

      setFubDataMap(ipos,jpos+1,kpos,idiv,jdiv,kdiv,fubDataMap);

    }
  }

  // +Z
  if( kpos < kdiv && root->nID[Z_PLUS] ) {
    if( !fubDataMap[_IDX_S3D(ipos,jpos,kpos+1,idiv,jdiv,kdiv,0)] ) {

      fubDataMap[_IDX_S3D(ipos,jpos,kpos+1,idiv,jdiv,kdiv,0)] = root->nID[Z_PLUS];

      setFubDataMap(ipos,jpos,kpos+1,idiv,jdiv,kdiv,fubDataMap);

    }
  }

  return;

}

// #################################################################
// proc.dfi ファイル出力
void fub_DATA::WriteProc(fub_DATA *firstNode, int idiv, int jdiv, int kdiv,
                         int xHead[3], int yHead[3], int zHead[3], int vc,
                         fub_DATA **fubDataMap)
{

  //ブロックの最小、最大座標値からリージョンを計算
  double reg[3];
  double min[3],max[3];

  for(int i=0; i<3; i++) {
    min[i]=fubDataMap[0]->min[i];
    max[i]=fubDataMap[0]->max[i];
    reg[i]=0.0e0;
  }
  for(int i=1; i<idiv*jdiv*kdiv; i++) {
  for(int j=0; j<3; j++) {
    if( min[j]>fubDataMap[i]->min[j] ) min[j]=fubDataMap[i]->min[j];
    if( max[j]<fubDataMap[i]->max[j] ) max[j]=fubDataMap[i]->max[j];
  }}

  for(int i=0; i<3; i++) reg[i] = max[i]-min[i];

  FILE* fp;
  if( !(fp=fopen("proc.dfi", "w")) )
  {
    printf("Error open fnaem : proc.dfi\n");
    return;
  }

  //Domain
  fprintf(fp, "Domain {\n");
  fprintf(fp, "\n");
  fprintf(fp, "  GlobalOrigin        = (%e, %e, %e)\n",
          firstNode->corner[_P000].x,firstNode->corner[_P000].y,firstNode->corner[_P000].z);
  fprintf(fp, "  GlobalRegion        = (%e, %e, %e)\n",
          reg[0],reg[1],reg[2]);
  fprintf(fp, "  GlobalVoxel         = (%d, %d, %d)\n",
          xHead[idiv]-1,yHead[jdiv]-1,zHead[kdiv]-1);
  fprintf(fp, "  GlobalDivision      = (%d, %d, %d)\n",
          idiv,jdiv,kdiv);
  fprintf(fp, "  ActiveSubdomainFile = \"\"\n");
  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  //MPI
  fprintf(fp, "MPI {\n");
  fprintf(fp, "\n");
  fprintf(fp, "  NumberOfRank   = %d\n", idiv*jdiv*kdiv);
  fprintf(fp, "  NumberOfGroup  = %d\n", 1);
  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  //Process
  int rank=0;
  fprintf(fp, "Process {\n");
  fprintf(fp, "\n");
  for(int k=0; k<kdiv; k++) {
  for(int j=0; j<jdiv; j++) {
  for(int i=0; i<idiv; i++) {
    fprintf(fp, "  Rank[@] {\n");
    fprintf(fp, "\n");
    //fprintf(fp, "    ID        = %d\n", rank);
    fprintf(fp, "    ID        = %d\n",
      fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->fileID);
    fprintf(fp, "    HostName  = \"%s\"\n", "");
    fprintf(fp, "    VoxelSize = (%d, %d, %d)\n",
      fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->size[0]-2*vc,
      fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->size[1]-2*vc,
      fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->size[2]-2*vc);
    fprintf(fp, "    HeadIndex = (%d, %d, %d)\n",xHead[i],yHead[j],zHead[k]);
    fprintf(fp, "    TailIndex = (%d, %d, %d)\n",
      xHead[i]+fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->size[0]-2*vc-1,
      yHead[j]+fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->size[1]-2*vc-1,
      zHead[k]+fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->size[2]-2*vc-1);
    fprintf(fp, "\n");
    fprintf(fp, "    CellID    = %d\n", 0);
    fprintf(fp, "    BCflagID  = %d\n", 0);
    fprintf(fp, "\n");
/*
    fprintf(fp, "    CoordinateFileName = \"%s\"\n",
      fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->fname_Coord.c_str());
    fprintf(fp, "    FieldDataFileName  = \"%s\"\n",
      fubDataMap[_IDX_S3D(i,j,k,idiv,jdiv,kdiv,0)]->fname_Field.c_str());
    fprintf(fp, "\n");
*/
    fprintf(fp, "  }\n");
    rank++;
  }}}
  fprintf(fp, "\n");
  fprintf(fp, "}\n");

  fclose(fp);
  return;
}
