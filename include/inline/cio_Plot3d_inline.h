#ifndef _CIO_PLOT3D_INLINE_H_
#define _CIO_PLOT3D_INLINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_Plot3d_inlune.h
 * @brief  cio_DFI_PLOT3D template Header
 * @author aics    
 */

#ifdef CIO_INLINE
 #undef CIO_INLINE
#endif

#ifndef CIO_NO_INLINE
 #define CIO_INLINE inline
#else
 #define CIO_INLINE
#endif

// #################################################################
// xyz 座標値を計算して出力
template<class T>
CIO_INLINE
void
cio_DFI_PLOT3D::write_XYZ(FILE* fp, T* org, T* pit, int sz[3])
{

  int ngrid=1;
  T xyz;

  //ascii
  if( m_output_type == CIO::E_CIO_OUTPUT_TYPE_ASCII ) {
    fprintf(fp,"%5d\n",ngrid);
    fprintf(fp,"%5d%5d%5d\n",sz[0],sz[1],sz[2]);

    //x
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[0]+pit[0]*(T)i;
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //y
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[1]+pit[1]*(T)j;
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //z
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[2]+pit[2]*(T)k;
      fprintf(fp,"%15.6E\n",xyz);
    }}}

  //Fortran Binary
  } else if( m_output_type == CIO::E_CIO_OUTPUT_TYPE_FBINARY ) {
    unsigned int dmy;
    dmy = sizeof(int);
    fwrite(&dmy, sizeof(int), 1, fp);
    fwrite(&ngrid, sizeof(int), 1, fp);
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sizeof(int)*3;
    fwrite(&dmy, sizeof(int), 1, fp);
    fwrite(&sz[0], sizeof(int), 1, fp);
    fwrite(&sz[1], sizeof(int), 1, fp);
    fwrite(&sz[2], sizeof(int), 1, fp);
    fwrite(&dmy, sizeof(int), 1, fp);

    dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*3);
    fwrite(&dmy, sizeof(int), 1, fp);
    //x
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[0]+pit[0]*(T)i;
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //y
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[1]+pit[1]*(T)j;
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //z
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[2]+pit[2]*(T)k;
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}
    fwrite(&dmy, sizeof(int), 1, fp);

  //Binary
  } else {
    fwrite(&ngrid, sizeof(int), 1, fp);
    fwrite(&sz[0], sizeof(int), 1, fp);
    fwrite(&sz[1], sizeof(int), 1, fp);
    fwrite(&sz[2], sizeof(int), 1, fp);

    //x
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[0]+pit[0]*(T)i;
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //y
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[1]+pit[1]*(T)j;
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //z
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = org[2]+pit[2]*(T)k;
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

  }

}

// #################################################################
// fund data 出力
template<class T>
CIO_INLINE
void
cio_DFI_PLOT3D::write_Func(FILE* fp, cio_TypeArray<T>* data, const int sz[3],
                           int ncomp)
{

  //IJNK
  if( data->getArrayShape() == CIO::E_CIO_IJKN ) {
    //ascii
    if( m_output_type == CIO::E_CIO_OUTPUT_TYPE_ASCII ) {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fprintf(fp,"%15.6E\n",data->val(i,j,k,n));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CIO::E_CIO_OUTPUT_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*ncomp);
      fwrite(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->val(i,j,k,n), sizeof(T), 1, fp);
      }}}}
      fwrite(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->val(i,j,k,n), sizeof(T), 1, fp);
      }}}}
    }

  //NIJK
  } else {
    //ascii
    if( m_output_type == CIO::E_CIO_OUTPUT_TYPE_ASCII ) {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fprintf(fp,"%15.6E\n",data->val(n,i,j,k));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CIO::E_CIO_OUTPUT_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*ncomp);
      fwrite(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
      fwrite(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
    }
  }

}

#endif // _CIO_PLOT3D_INLINE_H_
