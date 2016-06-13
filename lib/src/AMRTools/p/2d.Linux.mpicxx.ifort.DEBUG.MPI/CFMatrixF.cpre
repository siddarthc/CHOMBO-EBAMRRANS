      subroutine GETCOARSEFINEINTERPMATRIX(
     & cfmat
     & ,icfmatlo0,icfmatlo1
     & ,icfmathi0,icfmathi1
     & ,ncfmatcomp
     & ,nref
     & ,nnbrs
     & ,npwrs
     & ,degree
     & ,nfine
     & ,nbrsv
     & ,inbrsvhi0
     & ,ifboxlo0,ifboxlo1
     & ,ifboxhi0,ifboxhi1
     & ,idegboxlo0,idegboxlo1
     & ,idegboxhi0,idegboxhi1
     & )
      implicit none
      integer ncfmatcomp
      integer icfmatlo0,icfmatlo1
      integer icfmathi0,icfmathi1
      REAL*8 cfmat(
     & icfmatlo0:icfmathi0,
     & icfmatlo1:icfmathi1,
     & 0:ncfmatcomp-1)
      integer nref
      integer nnbrs
      integer npwrs
      integer degree
      integer nfine
      integer inbrsvhi0
      integer nbrsv(
     & 0:inbrsvhi0)
      integer ifboxlo0,ifboxlo1
      integer ifboxhi0,ifboxhi1
      integer idegboxlo0,idegboxlo1
      integer idegboxhi0,idegboxhi1
      integer i0,i1
      integer ipwr, sump, ifine, inbr, idim, ic
      integer pwrs(npwrs, 2)
      integer nbrs(nnbrs, 2)
      integer fine(nfine, 2)
      REAL*8 mat(nfine, nnbrs)
      ipwr = 0
      do i1 = idegboxlo1,idegboxhi1
      do i0 = idegboxlo0,idegboxhi0
         sump = i0 +i1
         if ((sump .gt. 0) .and. (sump .le. degree))then
            ipwr = ipwr + 1
            pwrs(ipwr, 1) = i0
            pwrs(ipwr, 2) = i1
         endif
      enddo
      enddo
      ic = 0
      do inbr = 1, nnbrs
         do idim = 1, 2
            nbrs(inbr, idim) = nbrsv(ic)
            ic = ic + 1
         end do
      end do
      ifine = 0
      do i1 = ifboxlo1,ifboxhi1
      do i0 = ifboxlo0,ifboxhi0
         ifine = ifine + 1
         fine(ifine, 1) = i0
         fine(ifine, 2) = i1
      enddo
      enddo
      call coarsefineleastsquares(nref, nnbrs, npwrs, nfine,
     & nbrs, pwrs, fine, mat)
      ifine = 0
      do i1 = icfmatlo1,icfmathi1
      do i0 = icfmatlo0,icfmathi0
         ifine = ifine + 1
         do inbr = 1, nnbrs
            cfmat(i0,i1, inbr-1) = mat(ifine, inbr)
         end do
      enddo
      enddo
      return
      end
      subroutine APPLYCOARSEFINEINTERP(
     & fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,stmat
     & ,istmatlo0,istmatlo1
     & ,istmathi0,istmathi1
     & ,nstmatcomp
     & ,fbase
     & ,cbase
     & ,nbrsv
     & ,inbrsvhi0
     & ,ifboxlo0,ifboxlo1
     & ,ifboxhi0,ifboxhi1
     & )
      implicit none
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer nstmatcomp
      integer istmatlo0,istmatlo1
      integer istmathi0,istmathi1
      REAL*8 stmat(
     & istmatlo0:istmathi0,
     & istmatlo1:istmathi1,
     & 0:nstmatcomp-1)
      integer fbase(0:1)
      integer cbase(0:1)
      integer inbrsvhi0
      integer nbrsv(
     & 0:inbrsvhi0)
      integer ifboxlo0,ifboxlo1
      integer ifboxhi0,ifboxhi1
      integer ibase0,ibase1
      integer icoar0,icoar1
      integer ifine0,ifine1
      integer ncomp, stsize, n, stpt, stbase
      REAL*8 val
      ncomp = nfinecomp
      stsize = (inbrsvhi0 + 1) / 2
      do n = 0, ncomp-1
      do ibase1 = ifboxlo1,ifboxhi1
      do ibase0 = ifboxlo0,ifboxhi0
            ifine0 = ibase0 + fbase(0)
            ifine1 = ibase1 + fbase(1)
            val = 0
            stbase = 0
            do stpt = 0, stsize-1
               icoar0 = cbase(0) + nbrsv(stbase )
               icoar1 = cbase(1) + nbrsv(stbase + 1)
               val = val +
     & stmat(ibase0,ibase1, stpt) *
     & coarse(icoar0,icoar1, n)
               stbase = stbase + 2
            end do
            fine(ifine0,ifine1, n) = val
      enddo
      enddo
      end do
      return
      end
