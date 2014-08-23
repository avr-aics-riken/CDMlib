!########################################################################
!
! CDMlib - Cartesian Data Management library
!
! Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!########################################################################

!! *********************************************************************
!! szS srcの実ボクセル数
!! gcS srcの仮想セル数
!! szD dstの実ボクセル数(=szS*2)
!! gcD dstの仮想セル数(=gcS*2)
!! nc  成分数
!! src 粗データの配列
!! dst 密データの配列

  do k=1-gcS,szS(3)+gcS
    kk=(k-1)*2+1
  do j=1-gcS,szS(2)+gcS
    jj=(j-1)*2+1
  do i=1-gcS,szS(1)+gcS
    ii=(i-1)*2+1
  do n=1,nc

    q = src(n,i,j,k)

    dst(n,ii  ,jj  ,kk  ) = q
    dst(n,ii+1,jj  ,kk  ) = q
    dst(n,ii  ,jj+1,kk  ) = q
    dst(n,ii+1,jj+1,kk  ) = q
    dst(n,ii  ,jj  ,kk+1) = q
    dst(n,ii+1,jj  ,kk+1) = q
    dst(n,ii  ,jj+1,kk+1) = q
    dst(n,ii+1,jj+1,kk+1) = q

  enddo
  enddo
  enddo
  enddo

