#ifndef _FUB_DATA_H_
#define _FUB_DATA_H_

/*
 *
 *
 *
 *
 */

#include <stdio.h>
#include <string>

// Text parser
#include "TextParser.h"
//
//#include "cdm_Define.h"
#include "cdm_endianUtil.h"
#include "cdm_TextParser.h"
#include "cdm_FieldFileNameFormat.h"
//
#include "fub_Define.h"
//
#include "Coord.h"

using namespace std;

class fub_DATA {

public:

  int size[4];
  double *xyz;
  float  *f_xyz;

  Coord corner[8];

  fub_DATA *nID[6];

  string fname_Coord;
  string fname_Field;

  double min[3],max[3];

  string endian;

  int fileID;

public:

  /** コンストラクタ */
  fub_DATA(string fubCname, string fubFname);

  /** デストラクタ */
  ~fub_DATA(){
     if( xyz != NULL ) delete [] xyz;
     if( f_xyz != NULL ) delete [] f_xyz;
   };

public:

  /** read fub file */
   bool readFub(int vc, string &dftype);

  /** フィールドデータの成分数取得 */
   int getNval();

  /** エンディアンチェック */
   static bool isMatchEndianSbdmMagick(int ident);

  /** コーナー8座標値のセット*/
   void setCorner(int vc, string dftype);

  /** ブロックの最小/最大座標値をセット（実セル、実ノード）*/
   void setMinMax(int vc, string dftype);

  /** 隣接ノード検索 */
   void checkNeighbor(fub_DATA *data2, int &cnt, double tol);

  /** -X方向のノード検索 */
   bool checkNeighborMX(fub_DATA *data2, double tol);
 
  /** +X方向のノード検索 */
   bool checkNeighborPX(fub_DATA *data2, double tol);
 
  /** -Y方向のノード検索 */
   bool checkNeighborMY(fub_DATA *data2, double tol);
 
  /** +Y方向のノード検索 */
   bool checkNeighborPY(fub_DATA *data2, double tol);
 
  /** -Z方向のノード検索 */
   bool checkNeighborMZ(fub_DATA *data2, double tol);
 
  /** +Z方向のノード検索 */
   bool checkNeighborPZ(fub_DATA *data2, double tol);

  /** ノードマップの生成 */
   void static setFubDataMap(int ipos, int jpos, int kpos, int idiv, int jdiv,
                      int kdiv, fub_DATA **fubDataMap);

  /** proc.dfi 出力 */
   void static WriteProc(fub_DATA *firstNode, int idiv, int jdiv, int kdiv,
                      int xHead[3], int yHead[3], int zHead[3], int vc,
                      fub_DATA **fubDataMap);

};

#endif // _FUB_DATA_H_
