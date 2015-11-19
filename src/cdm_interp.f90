!########################################################################
!
! CDMlib - Cartesian Data Management library
! 
! Copyright (c) 2013-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!########################################################################
!
!> @file    cdm_interp.f90
!! @brief   粗データsrcを密データdstに補間するルーチン群
!! @author  kero

!  *********************************************************************
!! @brief 粗データsrcを密データdstに補間する(i,j,k,n) real版
!! @param szS srcの実ボクセル数
!! @param gcS srcの仮想セル数
!! @param szD dstの実ボクセル数(=szS*2)
!! @param gcD dstの仮想セル数(=gcS*2)
!! @param nc  変数の個数
!! @param src 粗データの配列
!! @param dst 密データの配列
!! @todo
!!    - 粗ボクセル(i,j,k,n)位置に対応する密ボクセルは
!!    - dst(i  ,j  ,k  ,n)
!!    - dst(i+1,j  ,k  ,n)
!!    - dst(i  ,j+1,k  ,n)
!!    - dst(i+1,j+1,k  ,n)
!!    - dst(i  ,j  ,k+1,n)
!!    - dst(i+1,j  ,k+1,n)
!!    - dst(i  ,j+1,k+1,n)
!!    - dst(i+1,j+1,k+1,n)
subroutine cdm_interp_ijkn_r4(szS,gcS,szD,gcD,nc,src,dst)
  implicit none
  ! arguments
  integer :: szS(3),gcS,szD(3),gcD,nc
  real*4,dimension(1-gcS:szS(1)+gcS,1-gcS:szS(2)+gcS,1-gcS:szS(3)+gcS,nc) :: src
  real*4,dimension(1-gcD:szD(1)+gcD,1-gcD:szD(2)+gcD,1-gcD:szD(3)+gcD,nc) :: dst
  ! local variables
  integer :: i,j,k,n
  integer :: ii,jj,kk
  real*4  :: q

  include 'cdm_interp_ijkn.h'

  return
end subroutine cdm_interp_ijkn_r4

!  *********************************************************************
!! @brief 粗データsrcを密データdstに補間する(i,j,k,n) double版
!! @param szS srcの実ボクセル数
!! @param gcS srcの仮想セル数
!! @param szD dstの実ボクセル数(=szS*2)
!! @param gcD dstの仮想セル数(=gcS*2)
!! @param nc  変数の個数
!! @param src 粗データの配列
!! @param dst 密データの配列
!! @todo
!!    - 粗ボクセル(i,j,k,n)位置に対応する密ボクセルは
!!    - dst(i  ,j  ,k  ,n)
!!    - dst(i+1,j  ,k  ,n)
!!    - dst(i  ,j+1,k  ,n)
!!    - dst(i+1,j+1,k  ,n)
!!    - dst(i  ,j  ,k+1,n)
!!    - dst(i+1,j  ,k+1,n)
!!    - dst(i  ,j+1,k+1,n)
!!    - dst(i+1,j+1,k+1,n)
subroutine cdm_interp_ijkn_r8(szS,gcS,szD,gcD,nc,src,dst)
  implicit none
  ! arguments
  integer :: szS(3),gcS,szD(3),gcD,nc
  real*8,dimension(1-gcS:szS(1)+gcS,1-gcS:szS(2)+gcS,1-gcS:szS(3)+gcS,nc) :: src
  real*8,dimension(1-gcD:szD(1)+gcD,1-gcD:szD(2)+gcD,1-gcD:szD(3)+gcD,nc) :: dst
  ! local variables
  integer :: i,j,k,n
  integer :: ii,jj,kk
  real*8  :: q

  include 'cdm_interp_ijkn.h'

  return
end subroutine cdm_interp_ijkn_r8

!  *********************************************************************
!! @brief 粗データsrcを密データdstに補間する(n,i,j,k) real版
!! @param szS srcの実ボクセル数
!! @param gcS srcの仮想セル数
!! @param szD dstの実ボクセル数(=szS*2)
!! @param gcD dstの仮想セル数(=gcS*2)
!! @param nc  変数の個数
!! @param src 粗データの配列
!! @param dst 密データの配列
!! @todo
!!    - 粗ボクセル(n,i,j,k)位置に対応する密ボクセルは
!!    - dst(n,i  ,j  ,k  )
!!    - dst(n,i+1,j  ,k  )
!!    - dst(n,i  ,j+1,k  )
!!    - dst(n,i+1,j+1,k  )
!!    - dst(n,i  ,j  ,k+1)
!!    - dst(n,i+1,j  ,k+1)
!!    - dst(n,i  ,j+1,k+1)
!!    - dst(n,i+1,j+1,k+1)
subroutine cdm_interp_nijk_r4(szS,gcS,szD,gcD,nc,src,dst)
  implicit none
  ! arguments
  integer :: szS(3),gcS,szD(3),gcD,nc
  real*4,dimension(nc,1-gcS:szS(1)+gcS,1-gcS:szS(2)+gcS,1-gcS:szS(3)+gcS) :: src
  real*4,dimension(nc,1-gcD:szD(1)+gcD,1-gcD:szD(2)+gcD,1-gcD:szD(3)+gcD) :: dst
  ! local variables
  integer :: i,j,k,n
  integer :: ii,jj,kk
  real*4  :: q

  include 'cdm_interp_nijk.h'

  return
end subroutine cdm_interp_nijk_r4

!  *********************************************************************
!! @brief 粗データsrcを密データdstに補間する(n,i,j,k) double版
!! @param szS srcの実ボクセル数
!! @param gcS srcの仮想セル数
!! @param szD dstの実ボクセル数(=szS*2)
!! @param gcD dstの仮想セル数(=gcS*2)
!! @param nc  変数の個数
!! @param src 粗データの配列
!! @param dst 密データの配列
!! @todo
!!    - 粗ボクセル(n,i,j,k)位置に対応する密ボクセルは
!!    - dst(n,i  ,j  ,k  )
!!    - dst(n,i+1,j  ,k  )
!!    - dst(n,i  ,j+1,k  )
!!    - dst(n,i+1,j+1,k  )
!!    - dst(n,i  ,j  ,k+1)
!!    - dst(n,i+1,j  ,k+1)
!!    - dst(n,i  ,j+1,k+1)
!!    - dst(n,i+1,j+1,k+1)
subroutine cdm_interp_nijk_r8(szS,gcS,szD,gcD,nc,src,dst)
  implicit none
  ! arguments
  integer :: szS(3),gcS,szD(3),gcD,nc
  real*8,dimension(nc,1-gcS:szS(1)+gcS,1-gcS:szS(2)+gcS,1-gcS:szS(3)+gcS) :: src
  real*8,dimension(nc,1-gcD:szD(1)+gcD,1-gcD:szD(2)+gcD,1-gcD:szD(3)+gcD) :: dst
  ! local variables
  integer :: i,j,k,n
  integer :: ii,jj,kk
  real*8  :: q

  include 'cdm_interp_nijk.h'

  return
end subroutine cdm_interp_nijk_r8

