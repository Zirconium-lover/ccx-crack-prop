!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2024 Guido Dhondt
!     
      subroutine stressintensity(nfront,ifrontrel,stress,xt,xn,xa,
     &     dk1,dk2,dk3,xkeq,phi,psi,acrack,shape,nstep)
!     
!     calculate the stress intensity factors along the crack fronts
!     
      implicit none
!
! --- Debug switches ---
! Print one short line per m (step in the stress history)
      logical, parameter :: debug_inc    = .true.
! Print full tables (VERY verbose)
      logical, parameter :: debug_tables = .false.
! Print diagnostics only when "bad" states are detected
      logical, parameter :: debug_bad    = .true.
! Limit number of "bad" prints (avoid huge logs)
      integer, parameter :: bad_print_max = 50
! ----------------------
!
      integer nfront,ifrontrel(*),i,noderel,nstep,m
      integer bad_print_count
!
      real*8 s(3,3),stress(6,nstep,*),xt(3,*),xn(3,*),xa(3,*),
     &     dk1(nstep,*),dk2(nstep,*),dk3(nstep,*),xkeq(nstep,*),
     &     phi(nstep,*),psi(nstep,*),pi,c2,c3,c4,term,
     &     acrack(*),ratio,constant,t(3),shape(3,*),c1
!
! --- Diagnostics locals ---
      logical bad
      real*8 xn1,xn2,xn3,xa1,xa2,xa3,xt1,xt2,xt3
      real*8 nx,na,nt,smax,ac,sh1,sh2,sh3
      real*8 dk1raw,dk2raw,dk3raw
!
      pi=4.d0*datan(1.d0)
      c2=70.d0*pi/180.d0
      c3=78.d0*pi/180.d0
      c4=33.d0*pi/180.d0
!
      bad_print_count=0
!
! --- compute K1,K2,K3 from stresses ---
      do i=1,nfront
!
        noderel=ifrontrel(i)
!
        do m=1,nstep
!
          s(1,1)=stress(1,m,noderel)
          s(1,2)=stress(4,m,noderel)
          s(1,3)=stress(6,m,noderel)
          s(2,1)=s(1,2)
          s(2,2)=stress(2,m,noderel)
          s(2,3)=stress(5,m,noderel)
          s(3,1)=s(1,3)
          s(3,2)=s(2,3)
          s(3,3)=stress(3,m,noderel)
!
!         traction vector on the crack plane (t = s * xn)
!
          t(1)=s(1,1)*xn(1,i)+s(1,2)*xn(2,i)+s(1,3)*xn(3,i)
          t(2)=s(2,1)*xn(1,i)+s(2,2)*xn(2,i)+s(2,3)*xn(3,i)
          t(3)=s(3,1)*xn(1,i)+s(3,2)*xn(2,i)+s(3,3)*xn(3,i)
!
!         stress intensity factors in local crack front basis
!         (raw, before shape*sqrt(pi*a))
!
          dk1raw=t(1)*xn(1,i)+t(2)*xn(2,i)+t(3)*xn(3,i)
          dk2raw=t(1)*xa(1,i)+t(2)*xa(2,i)+t(3)*xa(3,i)
          dk3raw=t(1)*xt(1,i)+t(2)*xt(2,i)+t(3)*xt(3,i)
!
!         subsurface circular crack scaling
!
          ac=acrack(i)
          sh1=shape(1,i)
          sh2=shape(2,i)
          sh3=shape(3,i)
!
!         protect sqrt for non-positive acrack (nonphysical -> diagnostic)
!
          if(ac.le.0.d0) then
            constant=0.d0
          else
            constant=dsqrt(pi*ac)
          endif
!
          dk1(m,i)=dk1raw*sh1*constant
          dk2(m,i)=dk2raw*sh2*constant
          dk3(m,i)=dk3raw*sh3*constant
!
!         Diagnostics: catch abnormal scaling / broken local basis
!
          if(debug_bad) then
!
            xn1=xn(1,i)
            xn2=xn(2,i)
            xn3=xn(3,i)
            xa1=xa(1,i)
            xa2=xa(2,i)
            xa3=xa(3,i)
            xt1=xt(1,i)
            xt2=xt(2,i)
            xt3=xt(3,i)
!
            nx=dsqrt(xn1*xn1+xn2*xn2+xn3*xn3)
            na=dsqrt(xa1*xa1+xa2*xa2+xa3*xa3)
            nt=dsqrt(xt1*xt1+xt2*xt2+xt3*xt3)
!
            smax=dmax1(dabs(stress(1,m,noderel)),
     &                 dabs(stress(2,m,noderel)))
            smax=dmax1(smax,dabs(stress(3,m,noderel)))
            smax=dmax1(smax,dabs(stress(4,m,noderel)))
            smax=dmax1(smax,dabs(stress(5,m,noderel)))
            smax=dmax1(smax,dabs(stress(6,m,noderel)))
!
!           "Bad" triggers (tune if needed):
!           - norms of basis vectors far from 1
!           - acrack negative or huge
!           - shape huge
!           - stresses huge
!
            bad=.false.
            if(nx.gt.1.d2) bad=.true.
            if(na.gt.1.d2) bad=.true.
            if(nt.gt.1.d2) bad=.true.
            if(ac.le.0.d0) bad=.true.
            if(ac.gt.1.d6) bad=.true.
            if(dabs(sh1).gt.1.d3) bad=.true.
            if(dabs(sh2).gt.1.d3) bad=.true.
            if(dabs(sh3).gt.1.d3) bad=.true.
            if(smax.gt.1.d8) bad=.true.
!
            if(bad.and.(bad_print_count.lt.bad_print_max)) then
              bad_print_count=bad_print_count+1
              write(*,*) '*** BAD_STRESSINTENSITY ***'
              write(*,2001) i,m,noderel
              write(*,2002) nx,na,nt
              write(*,2003) ac,sh1,sh2,sh3
              write(*,2004) smax
              write(*,2005) dk1raw,dk2raw,dk3raw
              write(*,2006) dk1(m,i),dk2(m,i),dk3(m,i)
              call flush(6)
            endif
!
          endif
!
        enddo
      enddo
!
! --- compute Keq, phi, psi (Richard) ---
      do i=1,nfront
        do m=1,nstep
!
          if(dk1(m,i).ge.0.d0) then
            xkeq(m,i)=(dk1(m,i)+dsqrt(dk1(m,i)*dk1(m,i)+
     &           5.3361*dk2(m,i)*dk2(m,i)+
     &           4.d0*dk3(m,i)*dk3(m,i)))/2.d0
!
            if(xkeq(m,i).gt.1.d-20) then
              term=dk1(m,i)+dabs(dk2(m,i))+dabs(dk3(m,i))
              if(term.le.1.d-30) then
                phi(m,i)=0.d0
                psi(m,i)=0.d0
              else
                if(dabs(dk2(m,i)).gt.1.d-30) then
                  ratio=dabs(dk2(m,i))/term
                  phi(m,i)=-c2*ratio*(2.d0-ratio)*
     &                     dk2(m,i)/dabs(dk2(m,i))
                else
                  phi(m,i)=0.d0
                endif
!
                if(dabs(dk3(m,i)).gt.1.d-30) then
                  ratio=dabs(dk3(m,i))/term
                  psi(m,i)=-ratio*(c3-c4*ratio)*
     &                     dk3(m,i)/dabs(dk3(m,i))
                else
                  psi(m,i)=0.d0
                endif
              endif
            else
              phi(m,i)=0.d0
              psi(m,i)=0.d0
            endif
!
          else
            xkeq(m,i)=-(-dk1(m,i)+dsqrt(dk1(m,i)*dk1(m,i)+
     &           5.3361*dk2(m,i)*dk2(m,i)+
     &           4.d0*dk3(m,i)*dk3(m,i)))/2.d0
            phi(m,i)=0.d0
            psi(m,i)=0.d0
          endif
!
        enddo
      enddo
!
! --- Minimal debug output: one line per m (does not spam tables) ---
      if(debug_inc) then
        do m=1,nstep
          write(*,3000) m,nstep
          call flush(6)
        enddo
      endif
!
! --- Optional verbose tables ---
      if(debug_tables) then
        write(*,*) 'stressintensity k1 k2 k3'
        write(*,*)
        c1=1.d0/dsqrt(1000.d0)
        do m=1,nstep
          do i=1,nfront
            write(*,100) m,i,c1*dk1(m,i),c1*dk2(m,i),c1*dk3(m,i)
          enddo
          write(*,*)
        enddo
!
        write(*,*) 'stressintensity keq phi psi'
        write(*,*)
        do m=1,nstep
          do i=1,nfront
            write(*,100) m,i,c1*xkeq(m,i),phi(m,i)*180.d0/pi,
     &           psi(m,i)*180.d0/pi
          enddo
          write(*,*)
        enddo
      endif
!
 100  format(2i10,3(1x,e11.4))
 2001 format(1x,'i=',i8,' m=',i8,' noderel=',i10)
 2002 format(1x,'|xn|=',e12.4,' |xa|=',e12.4,' |xt|=',e12.4)
 2003 format(1x,'acrack=',e12.4,' sh=',3(1x,e12.4))
 2004 format(1x,'smax=',e12.4)
 2005 format(1x,'dk_raw=',3(1x,e12.4))
 2006 format(1x,'dk_fin=',3(1x,e12.4))
 3000 format(1x,'stressintensity: m=',i6,' / nstep=',i6)
!
      return
      end