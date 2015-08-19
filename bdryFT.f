!      Computes Fourier coeffs from outer plasma boundary r,z
!      Writes rmnc,zmnc, rmns,zmns etc in format to go into vmec input file
!      reads from stdin, writes to stdout, User should pipe from & to files at call

!      Note: Assumes toroidal symmetry (so good only for tokamaks)
!            i.e., LTHREED=F, LASYM=F (symmetric top-down) or T (asymmetric top-down)
!            So only relevant coeffs are from column 1 or 2 below, i.e.,
!            rmncc,rmnsc, zmncc,zmnsc   which we will call (for brevity)
!            rmnc, rmns,  zmnc, zmns
!            (for toroidal symm, no coeffs for s=sin in v direction)

!       From readin.f in vmec/input_output:
!       rbcc,rbss,rbcs,rbsc
!                boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
!       zbcc,zbss,zbcs,zbsc
!                boundary Fourier coefficient arrays for Z
!
!       XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC
!
!       STACKING ORDER DEPENDS ON LASYM AND LTHREED. EACH COMPONENT XCC, XSS, XSC, XCS
!       HAS SIZE = mns. (PHIFAC, MSE TAKE UP 1 INDEX EACH AT END OF ARRAY)
!
!         LTHREED=F,      LTHREED=F,      LTHREED=T,      LTHREED=T 
!         LASYM=F         LASYM=T         LASYM=F         LASYM=T
!
!          rmnc           rmnc             rmnc            rmnc           
!          zmns           rmns             rmnss           rmnss
!          lmnsc          zmns             zmns            rmns
!                         zmnc             zmncs           rmncs
!                         lmnsc            lmnsc           zmns
!                         lmncc            lmncs           zmncs
!                                                          zmnc
!                                                          zmnss
!                                                          lmnsc
!                                                          lmncs
!                                                          lmncc
!                                                          lmnss
!
      
!      lasym=1: top-down asymmetric, lasym=0: symmetric

       implicit none 
       integer, parameter :: mxp = 10000, mfc = 128
       double precision r(0:mxp),z(0:mxp),thpol(0:mxp)   ! input r,z
       double precision rr(0:mxp),zz(0:mxp),rs(0:mxp),zs(0:mxp)         ! output r,z
       double precision rcntr,zcntr,dth,thi,thm,cthi,cthm,sthi,sthm
       double precision rmnc(0:mfc),rmns(0:mfc)
       double precision zmnc(0:mfc),zmns(0:mfc)
       double precision pi,twopi,rerrmx,rerrav,zerrmx,zerrav
       integer n,i,m,k,nf,lasym
       character*1 cdum

       pi=dacos(-1.d0)
       twopi=2.d0*pi

!      Reading sequence compatible with input for xcurve (descur)
!      so vmcBD file made by corsica toVmec can be read by both

!      read n = number of bdry points
       read(*,*) n
!      read in r(i),z(i) of plasma boundary
       do i=1,n
         read(*,*) r(i),z(i)
       enddo

!      read plasma axis R,Z (not needed by descur) & symmetry type
       read(*,*) rcntr, zcntr, lasym  
       write(*,*) "# rcntr, zcntr, lasym ", rcntr, zcntr
       write(*,*) "# lasym ", lasym

       do i=1,n
         r(i)=r(i)-rcntr
         z(i)=z(i)-zcntr
       enddo

!      Calculate poloidal angle at each r,z (0 < th < 2*pi)
       do i=1,n
         thpol(i)=datan2(z(i),r(i))
         if (thpol(i).lt.0) thpol(i)=thpol(i)+twopi
       enddo

!      Close curve ends
       r(0)=r(n)
       z(0)=z(n)      
       thpol(0)=thpol(n)

!      Read number of fourier modes needed 
       read(*,*) nf     ! Ignores comments after these numbers
       write(*,*) "# nf ", nf

!      Calculate nf coeffs (no symmetry assumed)
       rmnc=0
       rmns=0
       zmnc=0
       zmns=0

       do k=0,nf      
         do i=1,n
           m=i-1
           dth=thpol(i)-thpol(m)
           if (dth.gt. pi) dth=dth-twopi
           if (dth.lt.-pi) dth=dth+twopi
           dth=dabs(dth/2.d0)

           thi=k*thpol(i)
           cthi=dcos(thi)
           sthi=dsin(thi)

           thm=k*thpol(m)
           cthm=dcos(thm)
           sthm=dsin(thm)           

           rmnc(k)=rmnc(k)+dth*(r(m)*cthm +r(i)*cthi)
           rmns(k)=rmns(k)+dth*(r(m)*sthm +r(i)*sthi)

           zmnc(k)=zmnc(k)+dth*(z(m)*cthm +z(i)*cthi)
           zmns(k)=zmns(k)+dth*(z(m)*sthm +z(i)*sthi)
         enddo
       enddo

       rmnc(0)=rmnc(0)*0.5d0
       rmns(0)=rmns(0)*0.5d0
       zmnc(0)=zmnc(0)*0.5d0
       zmns(0)=zmns(0)*0.5d0

       rmnc=rmnc/pi
       rmns=rmns/pi
       zmnc=zmnc/pi
       zmns=zmns/pi

!      Calculate errors by recalculating rr,zz from fourier coeffs & thpol(i)
       rerrmx=0.d0
       zerrmx=0.d0
       rerrav=0.d0
       zerrav=0.d0
       rr=0.d0
       zz=0.d0
       do i=1,n
         do k=0,nf
           thi=k*thpol(i)
           cthi=dcos(thi)
           sthi=dsin(thi)
           rr(i)=rr(i)+rmnc(k)*cthi+rmns(k)*sthi
           rs(i)=rs(i)+rmnc(k)*cthi   ! Symmetrized shape
           zz(i)=zz(i)+zmnc(k)*cthi+zmns(k)*sthi
           zs(i)=zs(i)+zmns(k)*sthi   ! Symmetrized shape
         enddo
         rerrmx=dmax1(rerrmx,dabs(rr(i)-r(i)))
         zerrmx=dmax1(zerrmx,dabs(zz(i)-z(i)))
         rerrav=rerrav+dabs(rr(i)-r(i))
         zerrav=zerrav+dabs(zz(i)-z(i))
       enddo
       rerrav=rerrav/n
       zerrav=zerrav/n
      
!      Write answers to stdout (can be piped wherever)
 
       write(*,*) "# ", rerrmx, zerrmx, " Max errors r,z"
       write(*,*) "# ", rerrav, zerrav, " Ave errors r,z"
 
       rmnc(0)=rmnc(0)+rcntr
       zmns(0)=zmns(0)+zcntr

       write(*,*) "# lasym=1 fourier coeffs for vmec"
       do k=0,nf
         WRITE(*,'(" RBC(0,",i3,")=",e20.12,
     1            ", RBS(0,",i3,")=",e20.12,
     2            ", ZBC(0,",i3,")=",e20.12,
     3            ", ZBS(0,",i3,")=",e20.12 )')
     4     k,rmnc(k), k,rmns(k), k,zmnc(k), k,zmns(k)
        enddo

       write(*,'(a,2e20.12)') 
     1  "# In Out r, z, rcntr, zcntr=",rcntr,zcntr
       do i=1,n
           write(*,'(5(g14.6,1x))')
     1     r(i)+rcntr,z(i)+zcntr,rr(i)+rcntr,zz(i)+zcntr,thpol(i)
       enddo

       write(*,*) "# lasym=0 fourier coeffs for vmec"
       do k=0,nf
         WRITE(*,'(" RBC(0,",i3,")=",e20.12,
     3            ", ZBS(0,",i3,")=",e20.12 )')
     4    k,rmnc(k), k,zmns(k)
        enddo

       write(*,*) "# Symmetrized & original r,z"
       do i=1,n
           write(*,'(5(g14.6,1x))')
     1     rs(i)+rcntr,zs(i)+zcntr,rr(i)+rcntr,zz(i)+zcntr,thpol(i)
       enddo

!      Warn if not enough theta points n > 10*nf
       if( n .lt. 10*nf ) then
         write(*,*) "# Num of pts ",n," < 10*nf ",nf
         write(*,*) "# Fourier transform probably bad"
       endif

       stop
       end


