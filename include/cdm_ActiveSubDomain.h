/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/*
 * @file ActiveSubDomain.h
 * @brif ActiveSubDomain Class Header
 * @author aics
 */

#ifndef _CDM_ACTIVESUBDOMAIN_
#define _CDM_ACTIVESUBDOMAIN_

/** ActiveSubDomian class */
class cdm_ActiveSubDomain 
{

public:

  /** デフォルトコンストラクタ */
  cdm_ActiveSubDomain();

  /** コンストラクタ
   * @param[in] pos  領域分割内での位置
   */ 
  cdm_ActiveSubDomain( int pos[3]);
  
  /** デストラクタ */
  virtual ~cdm_ActiveSubDomain();

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
  bool operator==(cdm_ActiveSubDomain dom);

  /** 比較演算子
   * @param[in] dom   比較対象の活性サブドメイン情報
   * @retval    true  違う位置情報を持つ
   * @retval    false 同じ位置情報を持つ
   */
  bool operator!=(cdm_ActiveSubDomain dom);

private:
  int m_pos[3]; ///<領域分割内での位置
   

};

#endif
