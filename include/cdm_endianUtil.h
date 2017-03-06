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
 * @file   cdm_endianUtil.h
 * @brief  エンディアンユーティリティマクロ・関数ファイル
 * @author aics
 */
#ifndef _CDM_ENDIAN_UTIL_H_
#define _CDM_ENDIAN_UTIL_H_

#include <stdio.h>
#include <sys/types.h>
#ifndef _WIN32
#include <unistd.h>           // for linux
#else
#include "cdm_win32_util.h"   // for win32
#endif
#include <iostream>
#include <string>

#ifndef _CDM_NO_INLINE_
  #define CDM_INLINE inline
#else
  #define CDM_INLINE
#endif

#ifdef SUPER_UX
  template<class X> CDM_INLINE void BSWAP16(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp;
    tmp = _x_v[0];
    _x_v[0] = _x_v[1];
    _x_v[1] = tmp;
  }
#else // SUPER_UX
  #ifndef BSWAP16
  #define BSWAP_X_16(x) \
    ( (((x) & 0xff00) >>  8) \
    | (((x) & 0x00ff) <<  8) )
  #define BSWAP16(x) { \
    register unsigned short& _x_v = (unsigned short&)(x); \
    _x_v = BSWAP_X_16(_x_v);}
  #endif // BSWAP16
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X> CDM_INLINE void BSWAP32(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp[2];
    tmp[0] = _x_v[0];
    tmp[1] = _x_v[1];
    _x_v[0] = _x_v[3];
    _x_v[1] = _x_v[2];
    _x_v[2] = tmp[1];
    _x_v[3] = tmp[0];
  }
#else // SUPER_UX
  #ifndef BSWAP32
  #define BSWAP_X_32(x) \
    ( (((x) & 0xff000000) >> 24) \
    | (((x) & 0x00ff0000) >>  8) \
    | (((x) & 0x0000ff00) <<  8) \
    | (((x) & 0x000000ff) << 24) )
  #define BSWAP32(x) \
    {register unsigned int& _x_v = (unsigned int&)(x); \
     _x_v = BSWAP_X_32(_x_v);}
  #endif // BSWAP32
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X> CDM_INLINE void BSWAP64(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp[4];
    tmp[0] = _x_v[0];
    tmp[1] = _x_v[1];
    tmp[2] = _x_v[2];
    tmp[3] = _x_v[3];
    _x_v[0] = _x_v[7];
    _x_v[1] = _x_v[6];
    _x_v[2] = _x_v[5];
    _x_v[3] = _x_v[4];
    _x_v[4] = tmp[3];
    _x_v[5] = tmp[2];
    _x_v[6] = tmp[1];
    _x_v[7] = tmp[0];
  }
#else // SUPER_UX
  #ifndef BSWAP64
  #define BSWAP_X_64(x) \
    ( (((x) & 0xff00000000000000ull) >> 56) \
    | (((x) & 0x00ff000000000000ull) >> 40) \
    | (((x) & 0x0000ff0000000000ull) >> 24) \
    | (((x) & 0x000000ff00000000ull) >>  8) \
    | (((x) & 0x00000000ff000000ull) <<  8) \
    | (((x) & 0x0000000000ff0000ull) << 24) \
    | (((x) & 0x000000000000ff00ull) << 40) \
    | (((x) & 0x00000000000000ffull) << 56) )
  #define BSWAP64(x) \
    {register unsigned long long& _x_v = (unsigned long long&)(x); \
     _x_v = BSWAP_X_64(_x_v);}
  #endif // BSWAP64
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X, class Y> CDM_INLINE void SBSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned short _x_v = (unsigned short)a[_i];
      BSWAP16(_x_v);
      a[_i] = _x_v;
    }
  }
#else // SUPER_UX
  #ifndef SBSWAPVEC
  #define SBSWAPVEC(a,n) do{\
    for(register unsigned int _i=0;_i<(n);_i++){BSWAP16(a[_i]);}\
  }while(0)
  #endif // SBSWAPVEC
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X, class Y> CDM_INLINE void BSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned int _x_v = (unsigned int)a[_i];
      BSWAP32(_x_v);
      a[_i] = _x_v;
    }
  }
#else // SUPER_UX
  #ifndef BSWAPVEC

  #define BSWAPVEC(a,n) do{\
    for(register unsigned int _i=0;_i<(n);_i++){BSWAP32(a[_i]);}\
  }while(0)
  #endif // BSWAPVEC
#endif // SUPER_UX

#ifdef SUPER_UX
  template<class X, class Y> CDM_INLINE void DBSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned long long _x_v = (unsigned long long)a[_i];
      BSWAP64(_x_v);
      a[_i] = _x_v;
    }
  }
#else // SUPER_UX
  #ifndef DBSWAPVEC
  #define DBSWAPVEC(a,n) do{\
    for(register unsigned int _i=0;_i<(n);_i++){BSWAP64(a[_i]);}\
  }while(0)
  #endif // DBSWAPVEC
#endif // SUPER_UX

#endif // _CDM_ENDIAN_UTIL_H_
