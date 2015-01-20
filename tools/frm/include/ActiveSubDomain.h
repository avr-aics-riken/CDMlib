/*
 * frm (File Rank Mapper)
 *
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * Staging
 * @file ActiveSubDomain.h
 * @brif ActiveSubDomain Class Header
 * @author aics
 */

#ifndef _ACTIVESUBDOMAIN_
#define _ACTIVESUBDOMAIN_

/** ActiveSubDomian class */
class ActiveSubDomain 
{

public:

  /** デフォルトコンストラクタ */
  ActiveSubDomain();

  /** コンストラクタ
   * @param[in] pos  領域分割内での位置
   */ 
  ActiveSubDomain( int pos[3]);
  
  /** デストラクタ */
  virtual ~ActiveSubDomain();

  /** 情報のクリア */
  virtual void clear();

  /** 位置のセット
   * @param[in] pos  領域分割内での位置
   */ 
  void SetPos( int pos[3] );

  /** 位置の取得
   * @return 位置情報整数配列のポインタ
   */
  const int* GetPos() const;

  /** 比較演算子
   * @param[in] dom   比較対象の活性サブドメイン情報
   * @retval    true  同じ位置情報を持つ
   * @retval    false 違う位置情報を持つ
   */
  bool operator==(ActiveSubDomain dom);

  /** 比較演算子
   * @param[in] dom   比較対象の活性サブドメイン情報
   * @retval    true  違う位置情報を持つ
   * @retval    false 同じ位置情報を持つ
   */
  bool operator!=(ActiveSubDomain dom);

private:
  int m_pos[3]; ///<領域分割内での位置
   

};

#endif
