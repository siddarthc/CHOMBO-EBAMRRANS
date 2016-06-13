#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

        subroutine REGPROLONG(
     &           phithislevel
     &           ,iphithislevello0,iphithislevello1
     &           ,iphithislevelhi0,iphithislevelhi1
     &           ,correctcoarse
     &           ,icorrectcoarselo0,icorrectcoarselo1
     &           ,icorrectcoarsehi0,icorrectcoarsehi1
     &           ,icoarboxlo0,icoarboxlo1
     &           ,icoarboxhi0,icoarboxhi1
     &           ,irefboxlo0,irefboxlo1
     &           ,irefboxhi0,irefboxhi1
     &           ,reftocoar
     &           )

      implicit none
      integer iphithislevello0,iphithislevello1
      integer iphithislevelhi0,iphithislevelhi1
      REAL_T phithislevel(
     &           iphithislevello0:iphithislevelhi0,
     &           iphithislevello1:iphithislevelhi1)
      integer icorrectcoarselo0,icorrectcoarselo1
      integer icorrectcoarsehi0,icorrectcoarsehi1
      REAL_T correctcoarse(
     &           icorrectcoarselo0:icorrectcoarsehi0,
     &           icorrectcoarselo1:icorrectcoarsehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer irefboxlo0,irefboxlo1
      integer irefboxhi0,irefboxhi1
      integer reftocoar
        integer iic,jjc
        integer iie,jje
        integer iif,jjf
        real_t fineval, coarval
        
      do jjc = icoarboxlo1,icoarboxhi1
      do iic = icoarboxlo0,icoarboxhi0

        
      do jje = irefboxlo1,irefboxhi1
      do iie = irefboxlo0,irefboxhi0

        
        iif  =  reftocoar*iic  + iie
        jjf  =  reftocoar*jjc  + jje
        fineval = phithislevel(iif,jjf)
        coarval = correctcoarse(iic,jjc)
        phithislevel(iif,jjf) = fineval + coarval
        
      enddo
      enddo
        
      enddo
      enddo
        return
        end
        subroutine PROLONGADDSLOPE(
     &           phithislevel
     &           ,iphithislevello0,iphithislevello1
     &           ,iphithislevelhi0,iphithislevelhi1
     &           ,correctcoarse
     &           ,icorrectcoarselo0,icorrectcoarselo1
     &           ,icorrectcoarsehi0,icorrectcoarsehi1
     &           ,icoarboxlo0,icoarboxlo1
     &           ,icoarboxhi0,icoarboxhi1
     &           ,irefboxlo0,irefboxlo1
     &           ,irefboxhi0,irefboxhi1
     &           ,idir
     &           ,dxf
     &           ,dxc
     &           ,reftocoar
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer iphithislevello0,iphithislevello1
      integer iphithislevelhi0,iphithislevelhi1
      REAL_T phithislevel(
     &           iphithislevello0:iphithislevelhi0,
     &           iphithislevello1:iphithislevelhi1)
      integer icorrectcoarselo0,icorrectcoarselo1
      integer icorrectcoarsehi0,icorrectcoarsehi1
      REAL_T correctcoarse(
     &           icorrectcoarselo0:icorrectcoarsehi0,
     &           icorrectcoarselo1:icorrectcoarsehi1)
      integer icoarboxlo0,icoarboxlo1
      integer icoarboxhi0,icoarboxhi1
      integer irefboxlo0,irefboxlo1
      integer irefboxhi0,irefboxhi1
      integer idir
      REAL_T dxf
      REAL_T dxc
      integer reftocoar
        integer iic,jjc
        integer iie,jje
        integer iif,jjf
        integer ioff,joff
        real_t coarslope, hival, loval, phifold, phifnew
        real_t dxcoar, dxfine, coarloc, fineloc, midval
        integer iindexc, iindexf
        dxcoar = dxc
        dxfine = dxf
        ioff = chf_id(0,idir)
                  joff = chf_id(1,idir)
        
      do jjc = icoarboxlo1,icoarboxhi1
      do iic = icoarboxlo0,icoarboxhi0

        iindexc = ioff*iic + joff*jjc
        coarloc = dxcoar*(iindexc + half)
        hival   = correctcoarse(iic + ioff,jjc + joff)
        loval   = correctcoarse(iic - ioff,jjc - joff)
        midval  = correctcoarse(iic,jjc)
        coarslope = (hival -loval)/(two*dxcoar)
        
      do jje = irefboxlo1,irefboxhi1
      do iie = irefboxlo0,irefboxhi0

        
        iif  =  reftocoar*iic  + iie
        jjf  =  reftocoar*jjc  + jje
        iindexf = ioff*iif + joff*jjf
        fineloc = dxfine*(iindexf + half)
        phifold = phithislevel(iif,jjf)
        phifnew = phifold + coarslope*(fineloc - coarloc)
        phithislevel(iif,jjf) = phifnew
        
      enddo
      enddo
        
      enddo
      enddo
        return
        end
