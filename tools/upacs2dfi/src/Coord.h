#ifndef _COORD_H_
#define _COORD_H_
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

class Coord {

public:

  double x,y,z;

public:

  /** コンストラクタ */
  Coord(){
    x=y=z=0.0e0;
  };

  /** デストラクタ */
  ~Coord(){};

public:

  //2点間の距離
  static double Distance(Coord a, Coord b)
  {
    return sqrt((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y)+(b.z-a.z)*(b.z-a.z));
  };

  //same判定
  static bool Same(Coord a, Coord b, double tol)
  {
    if( Distance(a,b) > tol ) return false;
    return true;
  };

};
#endif // _COORD_H_
