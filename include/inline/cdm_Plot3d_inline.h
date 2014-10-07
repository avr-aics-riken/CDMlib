#ifndef _CDM_PLOT3D_INLINE_H_
#define _CDM_PLOT3D_INLINE_H_

/*
 * CDMlib - Cartesian Data Management library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cdm_Plot3d_inlune.h
 * @brief  cdm_DFI_PLOT3D template Header
 * @author aics    
 */

#ifdef CDM_INLINE
 #undef CDM_INLINE
#endif

#ifndef CDM_NO_INLINE
 #define CDM_INLINE inline
#else
 #define CDM_INLINE
#endif

// #################################################################
// func data 読込み
template<class T>
CDM_INLINE
CDM::E_CDM_ERRORCODE
cdm_DFI_PLOT3D::read_Func(FILE* fp,
                         cdm_TypeArray<T>* dataS,
                         cdm_TypeArray<T>* dataB,
                         int head[3],
                         bool matchEndian)
{

  const int *szS = dataS->getArraySizeInt();
  const int *szB = dataB->getArraySizeInt();
  int ncompS = dataS->getNcomp();
  int ncompB = dataB->getNcomp();

  //１層ずつ読み込み
  int hzB = head[2];

  //IJNK
  if( dataS->getArrayShape() == CDM::E_CDM_IJKN ) {
    //ascii
    if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
      for(int n=0; n<ncompS; n++) {
      for(int k=0; k<szS[2]; k++) {
         //headインデクスをずらす
        head[2]=hzB+k;
        dataB->setHeadIndex(head);
         //一層読み込み
        for(int j=0; j<szB[1]; j++) {
        for(int i=0; i<szB[0]; i++) {
          fscanf(fp,"%f\n",&dataB->val(i,j,0,0));
        }}
        //コピー
        dataB->copyArrayNcomp(dataS,n);
      }}

    //Fortran Binary
    } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
      unsigned int dmy;
      fread(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<ncompS; n++) {
      for(int k=0; k<szS[2]; k++) {
        //headインデクスをずらす
        head[2]=hzB+k;
        dataB->setHeadIndex(head);
        //一層読み込み
        size_t ndata = dataB->getArrayLength();
        if( dataB->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
        //コピー
        dataB->copyArrayNcomp(dataS,n);
      }}
      fread(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<ncompS; n++) {
      for(int k=0; k<szS[2]; k++) {
        //headインデクスをずらす
        head[2]=hzB+k;
        dataB->setHeadIndex(head);
        //一層読み込み
        size_t ndata = dataB->getArrayLength();
        if( dataB->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
        //コピー
        dataB->copyArrayNcomp(dataS,n);
      }}
    }

  //NIJK
  } /* else {
    //ascii
    if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fscanf(fp,"%f\n",&data->val(n,i,j,k));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*ncomp);
      fread(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fread(&data->val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
      fread(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fread(&data->val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
    }
  }
  */
  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// xyz 座標値を計算して出力
template<class T>
CDM_INLINE
void
cdm_DFI_PLOT3D::write_XYZ(FILE* fp, int sz[3])
{

  int ngrid=1;
  T xyz;

  //ascii
  if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
    fprintf(fp,"%5d\n",ngrid);
    fprintf(fp,"%5d%5d%5d\n",sz[0],sz[1],sz[2]);

    //x
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellX(i));
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //y
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellY(j));
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //z
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellZ(k));
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //iblank
    if (DFI_Domain->iblank != NULL) {
      for(int i=0; i<sz[0]*sz[1]*sz[2]; i++) {
        fprintf(fp,"%2d\n",DFI_Domain->iblank[i]);
      }
    }

  //Fortran Binary
  } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
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
    if (DFI_Domain->iblank != NULL) {
      dmy += sizeof(int)*(sz[0]*sz[1]*sz[2]);
    }
    fwrite(&dmy, sizeof(int), 1, fp);
    //x
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellX(i));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //y
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellY(j));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //z
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellZ(k));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //iblank
    if (DFI_Domain->iblank != NULL) {
      fwrite(DFI_Domain->iblank, sizeof(int), sz[0]*sz[1]*sz[2], fp);
    }
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
      xyz = (T)(DFI_Domain->CellX(i));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //y
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellY(j));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //z
    for(int k=0; k<sz[2]; k++) {
    for(int j=0; j<sz[1]; j++) {
    for(int i=0; i<sz[0]; i++) {
      xyz = (T)(DFI_Domain->CellZ(k));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //iblank
    if (DFI_Domain->iblank != NULL) {
      fwrite(DFI_Domain->iblank, sizeof(int), sz[0]*sz[1]*sz[2], fp);
    }

  }

}

// #################################################################
// fund data 出力
template<class T>
CDM_INLINE
void
cdm_DFI_PLOT3D::write_Func(FILE* fp, cdm_TypeArray<T>* data, const int sz[3],
                           int ncomp)
{

  //IJKN
  if( data->getArrayShape() == CDM::E_CDM_IJKN ) {
    //ascii
    if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fprintf(fp,"%15.6E\n",data->_val(i,j,k,n));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*ncomp);
      fwrite(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(i,j,k,n), sizeof(T), 1, fp);
      }}}}
      fwrite(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(i,j,k,n), sizeof(T), 1, fp);
      }}}}
    }

  //NIJK
  } else {
    //ascii
    if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_ASCII ) {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fprintf(fp,"%15.6E\n",data->_val(n,i,j,k));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CDM::E_CDM_OUTPUT_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*ncomp);
      fwrite(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
      fwrite(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<ncomp; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
    }
  }

}

#endif // _CDM_PLOT3D_INLINE_H_
