!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2024 Guido Dhondt
!     
      subroutine crackprop(ifrontrel,ibounnod,domphi,da,co,costruc,nk,
     &     xa,xn,nnfront,istartfront,iendfront,doubleglob,integerglob,
     &     isubsurffront,dadn,ncyc,ifrontprop,nstep,acrack,acrackglob,
     &     datarget,iincglob,iinc,dnglob,ncyctot,ier)
!
!     Changes included:
!       1) Relaxed "outside" criterion for end-nodes: disttol_loc=1.d-12
!       2) hloc uses averaged spacing (reduces "accordion" front)
!       3) forbid backward advance: da>=0 and if cos(domphi)<0 => domphi=0
!       4) (optional speed-up) Cgeom increased (tune if needed)
!
      implicit none
!
      integer i,j,k,m,node,ifrontrel(*),nk,nnfront,istartfront(*),
     &     iendfront(*),integerglob(*),nktet,netet,ne,nkon,nfaces,
     &     nfield,nselect,imastset,iselect(6),nterms,nelem,
     &     ialset(1),iendset(1),istartset(1),konl(20),loopa,
     &     noderel,ibounnod(*),isubsurffront(*),ifrontprop(*),
     &     ncyc,nstep,nodep,nodeq,iincglob(*),iinc,
     &     ncyctot,ier
!
      real*8 da(*),domphi(*),co(3,*),costruc(3,*),xa(3,*),xn(3,*),
     &     doubleglob(*),coords(3),ratio(20),dist,pi,theta,ctheta,
     &     stheta,acrack(*),c(3,3),r0(3),r(3),dadn(*),dnglob(*),
     &     value(1),acrackglob(*),datarget,p(3),q(3),dd,al,
     &     hprev,hnext,hloc,damaxgeom,Cgeom,epsout,den,epsden,
     &     domphilim,disttol_loc,bestdist,bestcoords(3)
!
      nktet=integerglob(1)
      netet=integerglob(2)
      ne=integerglob(3)
      nkon=integerglob(4)
      nfaces=integerglob(5)
      nfield=13
!
      nselect=0
      imastset=0
      loopa=8
!
      pi=4.d0*datan(1.d0)
      theta=5.d0*pi/180.d0
      domphilim=0.5d0*pi-1.d-8
!
!     geometric cap factor (fraction of local front spacing)
!     (tune: 0.30 slow, 0.60-0.80 faster)
!
      Cgeom=0.70d0
!
      ier=0
!
      do k=1,nnfront
        do i=istartfront(k),iendfront(k)
!
          noderel=ifrontrel(i)
          node=ibounnod(noderel)
!
!         base increment (force non-negative)
!
          da(i)=dadn(i)*dble(ncyc)
          da(i)=dabs(da(i))
!
!         clamp domphi, and forbid backward along xa
!
          if(domphi(i).gt.domphilim) domphi(i)=domphilim
          if(domphi(i).lt.-domphilim) domphi(i)=-domphilim
          if(dcos(domphi(i)).lt.0.d0) domphi(i)=0.d0
!
!         minimum propagation (classic CCX style: 1% of datarget)
!
          if(da(i).lt.0.01d0*datarget) then
            da(i)=0.01d0*datarget
            domphi(i)=0.d0
          endif
!
!         local front spacing hloc (AVERAGED to reduce oscillations)
!
          hprev=1.d30
          hnext=1.d30
          if(i.gt.istartfront(k)) then
            nodep=ibounnod(ifrontrel(i-1))
            hprev=dsqrt((co(1,node)-co(1,nodep))**2+
     &                  (co(2,node)-co(2,nodep))**2+
     &                  (co(3,node)-co(3,nodep))**2)
          endif
          if(i.lt.iendfront(k)) then
            nodeq=ibounnod(ifrontrel(i+1))
            hnext=dsqrt((co(1,node)-co(1,nodeq))**2+
     &                  (co(2,node)-co(2,nodeq))**2+
     &                  (co(3,node)-co(3,nodeq))**2)
          endif
!
          if((hprev.lt.1.d20).and.(hnext.lt.1.d20)) then
            hloc=0.5d0*(hprev+hnext)
          elseif(hprev.lt.1.d20) then
            hloc=hprev
          elseif(hnext.lt.1.d20) then
            hloc=hnext
          else
            hloc=datarget
          endif
          hloc=dmax1(hloc,1.d-12)
!
!         geometric cap da <= Cgeom * hloc
!         for endpoints keep a floor to avoid collapse
!
          damaxgeom=Cgeom*hloc
          if((i.eq.istartfront(k)).or.(i.eq.iendfront(k))) then
            damaxgeom=dmax1(damaxgeom,0.05d0*datarget)
          endif
          if(da(i).gt.damaxgeom) then
            da(i)=damaxgeom
          endif
!
!         push outside: adaptive epsout
!
          epsout=dmax1(1.2d-6,1.d-3*datarget,1.d-2*hloc)
          da(i)=dsqrt(epsout**2+da(i)**2)
!
          nk=nk+1
          ifrontprop(i)=nk
!
          acrackglob(nk)=acrack(i)+da(i)
          iincglob(nk)=iinc+1
          dnglob(nk)=1.d0*ncyctot
!
          do j=1,3
            co(j,nk)=co(j,node)+(xa(j,i)*dcos(domphi(i))+
     &           xn(j,i)*dsin(domphi(i)))*da(i)
          enddo
!
!         subsurface front: no surface rotation needed
!
          if(isubsurffront(k).eq.1) cycle
!
!         endpoint handling: rotate propagated end nodes outside
!
          if(i.eq.istartfront(k)) then
!
            nodep=node
            nodeq=ibounnod(ifrontrel(i+1))
            do j=1,3
              p(j)=co(j,nodep)
              q(j)=co(j,nodeq)
              r(j)=costruc(j,noderel)
            enddo
            dd=(r(1)-q(1))**2+(r(2)-q(2))**2+(r(3)-q(3))**2
            den=(r(1)-q(1))*(p(1)-q(1))+
     &          (r(2)-q(2))*(p(2)-q(2))+
     &          (r(3)-q(3))*(p(3)-q(3))
            epsden=1.d-20*dmax1(dd,1.d0)
            if(dabs(den).gt.epsden) then
              al=dd/den
            else
              al=1.d0
            endif
            if((al.lt.0.d0).or.(al.gt.1.d0)) al=1.d0
            do j=1,3
              costruc(j,noderel)=q(j)+al*(p(j)-q(j))
            enddo
!
            do j=1,3
              r(j)=(xa(j,i)*dcos(domphi(i))+
     &             xn(j,i)*dsin(domphi(i)))*da(i)
              co(j,nk)=costruc(j,noderel)+r(j)
              r0(j)=r(j)
            enddo
!
            ctheta=dcos(theta)
            stheta=-dsin(theta)
!
            bestdist=-1.d30
            do j=1,3
              bestcoords(j)=co(j,nk)
            enddo
!
!           relaxed criterion for end-node
!
            disttol_loc=1.d-12
!
            ier=1
            do m=1,72
              do j=1,3
                coords(j)=co(j,nk)
              enddo
              call basis(doubleglob(1),doubleglob(netet+1),
     &             doubleglob(2*netet+1),doubleglob(3*netet+1),
     &             doubleglob(4*netet+1),doubleglob(5*netet+1),
     &             integerglob(6),integerglob(netet+6),
     &             integerglob(2*netet+6),doubleglob(6*netet+1),
     &             integerglob(3*netet+6),nktet,netet,
     &             doubleglob(4*nfaces+6*netet+1),nfield,
     &             doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &             integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &             integerglob(2*ne+7*netet+6),
     &             integerglob(nkon+2*ne+7*netet+6),coords(1),
     &             coords(2),coords(3),value,ratio,iselect,
     &             nselect,istartset,iendset,ialset,imastset,
     &             integerglob(nkon+2*ne+8*netet+6),nterms,konl,
     &             nelem,loopa,dist)
!
              if(dist.gt.bestdist) then
                bestdist=dist
                do j=1,3
                  bestcoords(j)=co(j,nk)
                enddo
              endif
              if(dist.ge.disttol_loc) then
                ier=0
                exit
              endif
!
              c(1,1)=ctheta+(1-ctheta)*(xn(1,i)**2)
              c(1,2)=-stheta*xn(3,i)+(1-ctheta)*xn(1,i)*xn(2,i)
              c(1,3)=stheta*xn(2,i)+(1-ctheta)*xn(1,i)*xn(3,i)
              c(2,1)=stheta*xn(3,i)+(1-ctheta)*xn(2,i)*xn(1,i)
              c(2,2)=ctheta+(1-ctheta)*(xn(2,i)**2)
              c(2,3)=-stheta*xn(1,i)+(1-ctheta)*xn(2,i)*xn(3,i)
              c(3,1)=-stheta*xn(2,i)+(1-ctheta)*xn(3,i)*xn(1,i)
              c(3,2)=stheta*xn(1,i)+(1-ctheta)*xn(3,i)*xn(2,i)
              c(3,3)=ctheta+(1-ctheta)*(xn(3,i)**2)
              do j=1,3
                r0(j)=r(j)
              enddo
              do j=1,3
                r(j)=c(j,1)*r0(1)+c(j,2)*r0(2)+c(j,3)*r0(3)
              enddo
              do j=1,3
                co(j,nk)=costruc(j,noderel)+r(j)
              enddo
            enddo
!
            if(ier.eq.1) then
              do j=1,3
                co(j,nk)=bestcoords(j)
              enddo
              write(*,'(a)') '*WARN crackprop: start end inside'
              write(*,'(a)') ' use best candidate'
              write(*,9000) iinc+1,k,i,bestdist
              ier=0
            endif
!
          elseif(i.eq.iendfront(k)) then
!
            nodep=ibounnod(ifrontrel(i-1))
            nodeq=node
            do j=1,3
              p(j)=co(j,nodep)
              q(j)=co(j,nodeq)
              r(j)=costruc(j,noderel)
            enddo
            dd=(r(1)-q(1))**2+(r(2)-q(2))**2+(r(3)-q(3))**2
            den=(r(1)-q(1))*(p(1)-q(1))+
     &          (r(2)-q(2))*(p(2)-q(2))+
     &          (r(3)-q(3))*(p(3)-q(3))
            epsden=1.d-20*dmax1(dd,1.d0)
            if(dabs(den).gt.epsden) then
              al=dd/den
            else
              al=0.d0
            endif
            if((al.lt.0.d0).or.(al.gt.1.d0)) al=0.d0
            do j=1,3
              costruc(j,noderel)=q(j)+al*(p(j)-q(j))
            enddo
!
            do j=1,3
              r(j)=(xa(j,i)*dcos(domphi(i))+
     &             xn(j,i)*dsin(domphi(i)))*da(i)
              co(j,nk)=costruc(j,noderel)+r(j)
              r0(j)=r(j)
            enddo
!
            ctheta=dcos(theta)
            stheta=dsin(theta)
!
            bestdist=-1.d30
            do j=1,3
              bestcoords(j)=co(j,nk)
            enddo
!
!           relaxed criterion for end-node
!
            disttol_loc=1.d-12
!
            ier=1
            do m=1,72
              do j=1,3
                coords(j)=co(j,nk)
              enddo
              call basis(doubleglob(1),doubleglob(netet+1),
     &             doubleglob(2*netet+1),doubleglob(3*netet+1),
     &             doubleglob(4*netet+1),doubleglob(5*netet+1),
     &             integerglob(6),integerglob(netet+6),
     &             integerglob(2*netet+6),doubleglob(6*netet+1),
     &             integerglob(3*netet+6),nktet,netet,
     &             doubleglob(4*nfaces+6*netet+1),nfield,
     &             doubleglob(nstep*13*nktet+4*nfaces+6*netet+1),
     &             integerglob(7*netet+6),integerglob(ne+7*netet+6),
     &             integerglob(2*ne+7*netet+6),
     &             integerglob(nkon+2*ne+7*netet+6),coords(1),
     &             coords(2),coords(3),value,ratio,iselect,
     &             nselect,istartset,iendset,ialset,imastset,
     &             integerglob(nkon+2*ne+8*netet+6),nterms,konl,
     &             nelem,loopa,dist)
!
              if(dist.gt.bestdist) then
                bestdist=dist
                do j=1,3
                  bestcoords(j)=co(j,nk)
                enddo
              endif
              if(dist.ge.disttol_loc) then
                ier=0
                exit
              endif
!
              c(1,1)=ctheta+(1-ctheta)*(xn(1,i)**2)
              c(1,2)=-stheta*xn(3,i)+(1-ctheta)*xn(1,i)*xn(2,i)
              c(1,3)=stheta*xn(2,i)+(1-ctheta)*xn(1,i)*xn(3,i)
              c(2,1)=stheta*xn(3,i)+(1-ctheta)*xn(2,i)*xn(1,i)
              c(2,2)=ctheta+(1-ctheta)*(xn(2,i)**2)
              c(2,3)=-stheta*xn(1,i)+(1-ctheta)*xn(2,i)*xn(3,i)
              c(3,1)=-stheta*xn(2,i)+(1-ctheta)*xn(3,i)*xn(1,i)
              c(3,2)=stheta*xn(1,i)+(1-ctheta)*xn(3,i)*xn(2,i)
              c(3,3)=ctheta+(1-ctheta)*(xn(3,i)**2)
              do j=1,3
                r0(j)=r(j)
              enddo
              do j=1,3
                r(j)=c(j,1)*r0(1)+c(j,2)*r0(2)+c(j,3)*r0(3)
              enddo
              do j=1,3
                co(j,nk)=costruc(j,noderel)+r(j)
              enddo
            enddo
!
            if(ier.eq.1) then
              do j=1,3
                co(j,nk)=bestcoords(j)
              enddo
              write(*,'(a)') '*WARN crackprop: end end inside'
              write(*,'(a)') ' use best candidate'
              write(*,9000) iinc+1,k,i,bestdist
              ier=0
            endif
!
          endif
!
        enddo
      enddo
!
      return
!
 9000 format(1x,'inc=',i6,' front=',i6,' i=',i6,' dist=',e11.4)
!
      end