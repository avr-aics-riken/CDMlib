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
                          int nz,
                          bool matchEndian)
{

  const int *szB = dataB->getArraySizeInt();
  int nvariS = dataS->getNvari();

  //１層ずつ読み込み
  int hzB = head[2];

  //IJNK
  if( dataS->getArrayShape() == CDM::E_CDM_IJKN ) {
    //ascii
    if( m_input_type == CDM::E_CDM_FILE_TYPE_ASCII ) {
      double temp;
#ifdef CDM_ENABLE_BUFFER_TUNING
      for(int n=0; n<nvariS; n++) {
        for(int k=0; k<szS[2]; k++) {
          for(int j=0; j<szB[1]; j++) {
            for(int i=0; i<szB[0]; i++) {
              fscanf(fp,"%lf\n",&temp);
              dataS->val(i,j,k,n) = (T)temp;
            }
		  }
        }
      }
#else
      for(int n=0; n<nvariS; n++) {
      for(int k=0; k<nz; k++) {
        //headインデクスをずらす
        head[2]=hzB+k;
        dataB->setHeadIndex(head);
        //一層読み込み
        size_t ndata = dataB->getArrayLength();
        if( szB[0]*szB[1] != ndata) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
        for(int j=0; j<szB[1]; j++) {
        for(int i=0; i<szB[0]; i++) {
          fscanf(fp,"%lf\n",&temp);
          dataB->val(i,j,0,0) = (T)temp;
        }}
        //コピー
        dataB->copyArrayNvari(dataS,n);
      }}
#endif
    //Fortran Binary
    } else if( m_input_type == CDM::E_CDM_FILE_TYPE_FBINARY ) {
      unsigned int dmy;
      fread(&dmy, sizeof(int), 1, fp);
#ifdef CDM_ENABLE_BUFFER_TUNING
      size_t ndata = dataS->getArrayLength();
      if( dataS->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
#else
      for(int n=0; n<nvariS; n++) {
      for(int k=0; k<nz; k++) {
        //headインデクスをずらす
        head[2]=hzB+k;
        dataB->setHeadIndex(head);
        //一層読み込み
        size_t ndata = dataB->getArrayLength();
        if( dataB->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
        //コピー
        dataB->copyArrayNvari(dataS,n);
      }}
#endif
      fread(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
#ifdef CDM_ENABLE_BUFFER_TUNING
      size_t ndata = dataS->getArrayLength();
      if( dataS->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
#else
      for(int n=0; n<nvariS; n++) {
      for(int k=0; k<nz; k++) {
        //headインデクスをずらす
        head[2]=hzB+k;
        dataB->setHeadIndex(head);
        //一層読み込み
        size_t ndata = dataB->getArrayLength();
        if( dataB->readBinary(fp,matchEndian) != ndata ) return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
        //コピー
        dataB->copyArrayNvari(dataS,n);
      }}
#endif
    }

  //NIJK (Plot3dの配列形状はIJKN)
  } else {
    return CDM::E_CDM_ERROR_READ_FIELD_DATA_RECORD;
  }

  return CDM::E_CDM_SUCCESS;
}

// #################################################################
// xyz 座標値を計算して出力
template<class T>
CDM_INLINE
void
cdm_DFI_PLOT3D::write_XYZ(FILE* fp, int sz[3], int head[3], const int* iblank)
{

  int ngrid=1;
  T xyz;
  int gc = DFI_Finfo.GuideCell;
  int sz3d_gc = (sz[0]+2*gc)*(sz[1]+2*gc)*(sz[2]+2*gc);

  //ascii
  if( m_output_type == CDM::E_CDM_FILE_TYPE_ASCII ) {
    //fprintf(fp,"%5d\n",ngrid);
    fprintf(fp,"%5d%5d%5d\n",sz[0],sz[1],sz[2]);

    //x
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellX(i+head[0]-1));
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //y
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellY(j+head[1]-1));
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //z
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellZ(k+head[2]-1));
      fprintf(fp,"%15.6E\n",xyz);
    }}}

    //iblank
    if (iblank != NULL) {
      for(int i=0; i<sz3d_gc; i++) {
        fprintf(fp,"%2d\n",iblank[i]);
      }
    }

  //Fortran Binary
  } else if( m_output_type == CDM::E_CDM_FILE_TYPE_FBINARY ) {
    unsigned int dmy;
    //dmy = sizeof(int);
    //fwrite(&dmy, sizeof(int), 1, fp);
    //fwrite(&ngrid, sizeof(int), 1, fp);
    //fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sizeof(int)*3;
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sz[0]+2*gc;
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sz[1]+2*gc;
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sz[2]+2*gc;
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sizeof(int)*3;
    fwrite(&dmy, sizeof(int), 1, fp);

    dmy = sizeof(T)*(sz3d_gc*3);
    if (iblank != NULL) {
      dmy += sizeof(int)*sz3d_gc;
    }
    fwrite(&dmy, sizeof(int), 1, fp);
    //x
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellX(i+head[0]-1));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //y
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellY(j+head[1]-1));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //z
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellZ(k+head[2]-1));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //iblank
    if (iblank != NULL) {
      fwrite(iblank, sizeof(int), sz3d_gc, fp);
    }
    fwrite(&dmy, sizeof(int), 1, fp);

  //Binary
  } else {
    unsigned int dmy;
    //fwrite(&ngrid, sizeof(int), 1, fp);
    dmy = sz[0]+2*gc;
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sz[1]+2*gc;
    fwrite(&dmy, sizeof(int), 1, fp);
    dmy = sz[2]+2*gc;
    fwrite(&dmy, sizeof(int), 1, fp);

    //x
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellX(i+head[0]-1));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //y
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellY(j+head[1]-1));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //z
    for(int k=-gc; k<sz[2]+gc; k++) {
    for(int j=-gc; j<sz[1]+gc; j++) {
    for(int i=-gc; i<sz[0]+gc; i++) {
      xyz = (T)(DFI_Domain->CellZ(k+head[2]-1));
      fwrite(&xyz, sizeof(T), 1, fp);
    }}}

    //iblank
    if (iblank != NULL) {
      fwrite(iblank, sizeof(int), sz3d_gc, fp);
    }

  }

}

// #################################################################
// fund data 出力
template<class T>
CDM_INLINE
void
cdm_DFI_PLOT3D::write_Func(FILE* fp, cdm_TypeArray<T>* data, const int sz[3],
                           int nvari)
{

  //IJKN
  if( data->getArrayShape() == CDM::E_CDM_IJKN ) {
    //ascii
    if( m_output_type == CDM::E_CDM_FILE_TYPE_ASCII ) {
      for(int n=0; n<nvari; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fprintf(fp,"%15.6E\n",data->_val(i,j,k,n));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CDM::E_CDM_FILE_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*nvari);
      fwrite(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<nvari; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(i,j,k,n), sizeof(T), 1, fp);
      }}}}
      fwrite(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<nvari; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(i,j,k,n), sizeof(T), 1, fp);
      }}}}
    }

  //NIJK
  } else {
    //ascii
    if( m_output_type == CDM::E_CDM_FILE_TYPE_ASCII ) {
      for(int n=0; n<nvari; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fprintf(fp,"%15.6E\n",data->_val(n,i,j,k));
      }}}}

    //Fortran Binary
    } else if( m_output_type == CDM::E_CDM_FILE_TYPE_FBINARY ) {
      unsigned int dmy;
      dmy = sizeof(T)*(sz[0]*sz[1]*sz[2]*nvari);
      fwrite(&dmy, sizeof(int), 1, fp);
      for(int n=0; n<nvari; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
      fwrite(&dmy, sizeof(int), 1, fp);

    //binary
    } else {
      for(int n=0; n<nvari; n++) {
      for(int k=0; k<sz[2]; k++) {
      for(int j=0; j<sz[1]; j++) {
      for(int i=0; i<sz[0]; i++) {
        fwrite(&data->_val(n,i,j,k), sizeof(T), 1, fp);
      }}}}
    }
  }

}

#endif // _CDM_PLOT3D_INLINE_H_
