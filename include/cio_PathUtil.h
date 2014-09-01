/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef _CDM_PATHUTIL_H_
#define _CDM_PATHUTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

#define MAXPATHLEN 512

namespace CDM
{
  inline char cdmPath_getDelimChar()
  {
#ifdef WIN32
    return '\\';
#else
    return '/';
#endif
  }

  inline std::string cdmPath_getDelimString()
  {
    const char dc = CDM::cdmPath_getDelimChar();
    char rs[2] = {dc, '\0'};
    return rs;
  }

  inline bool cdmPath_hasDrive(const std::string& path) {
    if ( path.size() < 2 ) return false;
    char x = path[0];
    if ( ((x >= 'A' && x <= 'Z' ) || (x >= 'a' && x <= 'z')) &&
   path[1] == ':' )
      return true;
    return false;
  }

  inline std::string vfvPath_emitDrive(std::string& path)
  {
    // returns drive (ex. 'C:')
    if ( ! cdmPath_hasDrive(path) ) return std::string();
    std::string driveStr = path.substr(0, 2);
    path = path.substr(2);
    return driveStr;
  }

  // true  : Absolute Path(絶対パス)
  // false : Relative Path(相対パス)
  inline bool cdmPath_isAbsolute(const std::string& path)
  {
    std::string xpath(path);
    vfvPath_emitDrive(xpath);
    char c1, c2;
    c1 = xpath[0];
    c2 = cdmPath_getDelimChar();
    return (c1 == c2);
  }

  inline std::string cdmPath_DirName(const std::string& path,
                                     const char dc = cdmPath_getDelimChar()) {
    char* name = strdup( path.c_str() );
    char* p = name;

    for ( ; ; ++p ) {
      if ( ! *p ) {
        if ( p > name ) {
          char rs[2] = {dc, '\0'};
          return rs;
        } else {
          char rs[3] = {'.', dc, '\0'};
          return rs;
        }
      }
      if ( *p != dc ) break;
    }

    for ( ; *p; ++p );
    while ( *--p == dc ) continue;
    *++p = '\0';

    while ( --p >= name )
      if ( *p == dc ) break;
    ++p;
    if ( p == name )
    {
      char rs[3] = {'.', dc, '\0'};
      return rs;
    }

    while ( --p >= name )
      if ( *p != dc ) break;
    ++p;

    *p = '\0';
    if( p == name ) {
      char rs[2] = {dc, '\0'};
      return rs;
    } else {
      std::string s( name );
      free( name );
      if( !CDM::cdmPath_isAbsolute(s) )
      {
        const char *q = s.c_str();
        if( q[0] != '.' && q[1] != '/' )
        {
          char rs[3] = {'.', dc, '\0'};
          s = std::string(rs) + s;
        }
      }
      return s;
    }
  }

  inline std::string cdmPath_FileName(const std::string& path,
                                      const std::string& addext = std::string(""),
                                      const char dc = cdmPath_getDelimChar()) {
    char* name = strdup( path.c_str() );
    char* p = name;

    for ( ; ; ++p ) {
      if ( ! *p ) {
        if ( p > name ) {
          char rs[2] = {dc, '\0'};
          return rs;
        } else
          return "";
      }
      if ( *p != dc ) break;
    }

    for ( ; *p; ++p ) continue;
    while ( *--p == dc ) continue;
    *++p = '\0';

    while ( --p >= name )
      if ( *p == dc ) break;
    ++p;

    bool add = false;
    if ( addext.length() > 0 ) {
      const int suffixlen = addext.length();
      const int stringlen = strlen( p );
      if ( suffixlen < stringlen ) {
        const int off = stringlen - suffixlen;
        if ( strcasecmp( p + off, addext.c_str()) != 0 )
          add = true;
      }
      else
      {
          add = true;
      }
    }

    std::string s( p );
    if( add ) s += addext;

    free( name );
    return s;
  }

  inline std::string cdmPath_ConnectPath( std::string dirName, std::string fname )
  {
    std::string path = dirName;

    const char *p = dirName.c_str();
    if( p[strlen(p)-1] != CDM::cdmPath_getDelimChar() )
    {
      path += CDM::cdmPath_getDelimString();
    }

    path += fname;

    return path;
  }

  //ファイル名から拡張子を削除
  inline std::string ExtractPathWithoutExt(const std::string &fn)
  {
    std::string::size_type pos;
    if((pos = fn.find_last_of(".")) == std::string::npos){
        return fn;
    }
 
    return fn.substr(0, pos);
  } 

};

#endif /* _CDM_PATHUTIL_H_ */

