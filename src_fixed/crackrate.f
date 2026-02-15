!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2024 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
      subroutine crackrate(nfront,ifrontrel,xkeq,phi,ifront,
     &     dadn,ncyc,icritic,datarget,crcon,temp,ncrtem,
     &     crconloc,ncrconst,xk1,xk2,xk3,nstep,acrack,
     &     wk1,wk2,wk3,xkeqmin,xkeqmax,dkeq,domstep,domphi,
     &     param,nparam,law,ier,r)
!
      implicit none
!
      character*132 param(*)
!
      integer i,nfront,ifrontrel(*),noderel,icritic,ncyc,ncrconst,
     &     ncrtem,nstep,m,ifront(*),nparam,law,ier
!
      real*8 datarget,damax,xkeq(nstep,*),phi(nstep,*),
     &     acrack(*),dadn(*),crcon(0:ncrconst,*),t1l,
     &     crconloc(*),dadnref,dkref,xm,epsilon,dkth,delta,dkc,
     &     xk1(nstep,*),xk2(nstep,*),xk3(nstep,*),temp(nstep,*),
     &     wk1(*),wk2(*),wk3(*),w,r(*),fr,
     &     xkeqmin(*),xkeqmax(*),dkeq(*),domstep(*),domphi(*)
!
      integer kcreached
      real*8 dkeqeff, tiny, jumpfactor, da_i
!
      tiny       = 1.d-6
      jumpfactor = 5.d0
!
      ier     = 0
      icritic = 0
      kcreached = 0
!
      if(datarget.le.0.d0) then
        write(*,*) '*WARNING in crackrate: datarget<=0, reset'
        datarget = 1.d-12
      endif
!
      damax=0.d0
      do i=1,nfront
        xkeqmin(i)=1.d30
        xkeqmax(i)=-1.d30
        wk1(i)=0.d0
        wk2(i)=0.d0
        wk3(i)=0.d0
!
        do m=1,nstep
          if(xkeq(m,i).gt.xkeqmax(i)) then
            xkeqmax(i)=xkeq(m,i)
            domphi(i)=phi(m,i)
            domstep(i)=1.d0*m
          endif
          xkeqmin(i)=min(xkeqmin(i),xkeq(m,i))
          if(dabs(xk1(m,i)).gt.dabs(wk1(i))) wk1(i)=xk1(m,i)
          if(dabs(xk2(m,i)).gt.dabs(wk2(i))) wk2(i)=xk2(m,i)
          if(dabs(xk3(m,i)).gt.dabs(wk3(i))) wk3(i)=xk3(m,i)
        enddo
!
        if(nstep.eq.1) then
          xkeqmin(i)=0.d0
        endif
        dkeq(i)=xkeqmax(i)-xkeqmin(i)
!
        if(dabs(xkeqmax(i)).lt.1.d-10) then
          if(xkeqmax(i).lt.0.d0) then
            r(i)=1.d10
          else
            r(i)=-1.d10
          endif
        else
          r(i)=xkeqmin(i)/xkeqmax(i)
        endif
!
        noderel=ifrontrel(i)
        t1l=temp(1,noderel)
        call materialdata_crack(crcon,ncrconst,ncrtem,t1l,crconloc)
!
        dadnref=crconloc(1)
        dkref  =crconloc(2)
        xm     =crconloc(3)
        epsilon=crconloc(4)
        dkth   =crconloc(5)
        delta  =crconloc(6)
        dkc    =crconloc(7)
        w      =crconloc(8)
!
        if(dabs(1.d0-r(i)).lt.1.d-10) then
          fr=0.d0
        else
          fr=1.d0/(1.d0-r(i))**((1.d0-w)*xm)
        endif
!
!       Kc reached: do NOT stop, but mark and use effective dk for formula
!
        dkeqeff = dkeq(i)
        if(dkeq(i).ge.dkc) then
          kcreached = 1
          write(*,*) '*WARNING in crackrate: Kc is reached'
          write(*,*) '         original K-range: ',dkeq(i)
          dkeqeff = dkc*(1.d0-tiny)
          write(*,*) '         dk used in law : ',dkeqeff
          write(*,*)
        endif
!
        if(dkeqeff.le.dkth) then
          dadn(i)=0.d0
        else
          dadn(i)=dadnref*(dkeqeff/dkref)**xm*
     &         (1.d0-dexp(epsilon*(1.d0-dkeqeff/dkth)))/
     &         (1.d0-dexp(delta*(xkeqmax(i)/dkc-1.d0)))*fr
        endif
!
        damax=max(dadn(i),damax)
      enddo
!
!     cycles selection: keep jumps when Kc is reached
!
      if(damax.le.0.d0) then
        icritic=-1
        ncyc=1
        return
      endif
!
      if(kcreached.eq.1) then
        ncyc = 1
      else
        ncyc = nint(datarget/damax)
        if(ncyc.le.0) ncyc=1
      endif
!
!     soft cap of per-increment jump: da_i <= jumpfactor*datarget
!
      do i=1,nfront
        da_i = dadn(i)*dble(ncyc)
        if(da_i.gt.jumpfactor*datarget) then
          dadn(i) = (jumpfactor*datarget)/dble(ncyc)
        endif
      enddo
!
      icritic = 0
      return
      end
