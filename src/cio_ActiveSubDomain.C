/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file  cio_ActiveSubDomain.C
 * @brief cio_ActiveSubDomain class 関数
 * @author aics
 */

#include "cio_ActiveSubDomain.h"

/////////////////////////////////////////////////////////////////////
// デフォルトコンストラクタ
cio_ActiveSubDomain::cio_ActiveSubDomain()
{
  clear();
}

/////////////////////////////////////////////////////////////////////
// コンストラクタ
cio_ActiveSubDomain::cio_ActiveSubDomain( int pos[3] )
{
  SetPos(pos);
}

/////////////////////////////////////////////////////////////////////
// デストラクタ
cio_ActiveSubDomain::~cio_ActiveSubDomain()
{
}

/////////////////////////////////////////////////////////////////////
// 情報のクリア
void cio_ActiveSubDomain::clear()
{
  m_pos[0]=0; 
  m_pos[1]=0; 
  m_pos[2]=0; 
}

/////////////////////////////////////////////////////////////////////
// 位置情報のセット
void cio_ActiveSubDomain::SetPos( int pos[3] )
{
  m_pos[0] = pos[0];
  m_pos[1] = pos[1];
  m_pos[2] = pos[2];
}

/////////////////////////////////////////////////////////////////////
// 位置情報の取得
const int* cio_ActiveSubDomain::GetPos() const
{
  return m_pos;
}

/////////////////////////////////////////////////////////////////////
// 比較演算子
bool cio_ActiveSubDomain::operator==(cio_ActiveSubDomain dom)
{
  if( m_pos[0] != dom.m_pos[0] ) return false;
  if( m_pos[1] != dom.m_pos[1] ) return false;
  if( m_pos[2] != dom.m_pos[2] ) return false;
  return true;
}

/////////////////////////////////////////////////////////////////////
// 比較演算子
bool cio_ActiveSubDomain::operator!=(cio_ActiveSubDomain dom)
{
  if( m_pos[0] == dom.m_pos[0] ) return false;
  if( m_pos[1] == dom.m_pos[1] ) return false;
  if( m_pos[2] == dom.m_pos[2] ) return false;
  return true;
}



