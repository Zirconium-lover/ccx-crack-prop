!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2024 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine eqspacednodes(co,istartfront,iendfront,nnfront,
     &     ifrontprop,nk,nfront,ifronteq,charlen,
     &     istartfronteq,iendfronteq,nfronteq,acrackglob,ier,
     &     iendcrackfro,iincglob,iinc,dnglob,ncyctot)
!     
!     Determine equally spaced nodes along the propagated crack front.
!     
!     acrackglob, iincglob and dnglob are the crack length, the
!     increment number and the total number of cycles at the END
!     of the present increment and are therefore attached to the
!     propagated front; they get values assigned in the present routine.
!     
      implicit none
!     
!     FRONT_MULT must match the allocation strategy in crackpropagation.c
!     
      integer, parameter :: FRONT_MULT   = 5
!     
!     Optional per-segment cap to avoid very large insertion on a single
!     edge (set large if you want to rely only on nklim).
!     
      integer, parameter :: NODESNUM_MAX = 100
!     
      integer i,k,m,n1,n2,mm,j
      integer istartfront(*),iendfront(*),iendcrackfro(*)
      integer nnfront,ifrontprop(*),nodesnum,ier,icrack,nk,nfront
      integer ifronteq(*),istartfronteq(*),iendfronteq(*),nfronteq
      integer iincglob(*),iinc,ncyctot,nklim,nk0
!     
      real*8 co(3,*),dist,charlen(*),x1,x2
      real*8 acrackglob(*),dnglob(*)
!     
      ier = 0
      nk0 = nk
      nklim = nk0 + FRONT_MULT*nfront
!     
!     Loop over all fronts
!     
      icrack  = 1
      nfronteq = 0
!     
      do i=1,nnfront
        istartfronteq(i)=nfronteq+1
!     
!       Determine crack index for this front
!     
        if(iendcrackfro(icrack).lt.istartfront(i)) then
          icrack=icrack+1
        endif
!     
!       Check characteristic length
!     
        if(charlen(icrack).le.0.d0) then
          write(*,*) '*ERROR in eqspacednodes: charlen<=0'
          ier=1
          return
        endif
        if(charlen(icrack).lt.1.d-12) charlen(icrack)=1.d-12
!     
!       First node of the equivalent front: keep position of the first
!       propagated node (for consistent numbering)
!     
        nk = nk + 1
        if (nk.gt.nklim) then
          write(*,*) '*ERROR in eqspacednodes: nk>nklim'
          ier = 1
          return
        endif
        n1 = ifrontprop(istartfront(i))
        do k=1,3
          co(k,nk)=co(k,n1)
        enddo
        acrackglob(nk)=acrackglob(n1)
        iincglob(nk)=iinc+1
        dnglob(nk)=1.d0*ncyctot
        ifronteq(istartfronteq(i))=nk
!     
!       Insert nodes along segments (n1 -> n2), excluding n1 and including n2
!     
        m = 0
        do mm=istartfront(i)+1,iendfront(i)
          n2 = ifrontprop(mm)
          dist=dsqrt((co(1,n2)-co(1,n1))**2+
     &               (co(2,n2)-co(2,n1))**2+
     &               (co(3,n2)-co(3,n1))**2)
!         If two adjacent nodes almost coincide: skip segment
          if(dist.lt.1.d-12) then
            n1 = n2
            cycle
          endif
!         Number of inserted nodes on this segment
          nodesnum=nint(dist/charlen(icrack))
          if(nodesnum.lt.1) nodesnum=1
          if(nodesnum.gt.NODESNUM_MAX) nodesnum=NODESNUM_MAX
!         Soft memory cap: do not exceed nklim
          if(nk+nodesnum.gt.nklim) then
            nodesnum = nklim - nk
            if(nodesnum.lt.1) then
              write(*,*) '*ERROR in eqspacednodes: nk>nklim'
              ier=1
              return
            endif
          endif
!         Create nodes (including a copy of n2 at j=nodesnum)
          do j=1,nodesnum
            x2 = dble(j)/dble(nodesnum)
            x1 = 1.d0 - x2
            nk = nk + 1
            if (nk.gt.nklim) then
              write(*,*) '*ERROR in eqspacednodes: nk>nklim'
              ier=1
              return
            endif
            co(1,nk)=x1*co(1,n1)+x2*co(1,n2)
            co(2,nk)=x1*co(2,n1)+x2*co(2,n2)
            co(3,nk)=x1*co(3,n1)+x2*co(3,n2)
            acrackglob(nk)=x1*acrackglob(n1)+x2*acrackglob(n2)
            iincglob(nk)=iinc+1
            dnglob(nk)=1.d0*ncyctot
            ifronteq(istartfronteq(i)+m+j)=nk
          enddo
          m  = m + nodesnum
          n1 = n2
        enddo
!     
!       End of this front
!     
        nfronteq = nfronteq + m + 1
        iendfronteq(i) = nfronteq
      enddo
!     
      return
      end
