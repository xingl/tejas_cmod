!  Cubic spline version of vmpol (which came from vmgs)
!-----------------------
!  Read vmec output file wout_ext and write a gs2 input file gs2_ext.nc
!  Both use NETCDF data structures (why, I do not know)
!  P. Valanju April 2005

!To compile this code, type ncmpf vmgs
!           This links it to vmec and netcdf libraries. (see ~/bin/ncmpf)

!To run it, type vmgs
!  It will ask you for wout extension and number of theta points
!  Enter for example:  test  64

!  Reads from vmec wout_ext.nc or wout_ext.txt file (whichever exists)
!  Writes to gs2_ext.nc or gs2_ext.txt file (whichever is specified)
!  nc => NETCDF, txt => ascii text file

! ************** Caution *******************
! Note: because originally this was meant for gs2, and CORSICA needs cgs,
!       so this routine writes presf[in cgs (dynes/cm^2)] = 10*presf[in Pascal]
!       I do not change it here because many old scripts will break, but
!       when writing new ones (say for GENE), please keep this in mind.
!       Also, since gs2 needed cgs, there may be cm instead of m in this one.
! ******************************************

! Auxiliary gnuplottable files written:
!   fort.11  : R[m], Z[m], |B|[Tesla] vs (psi,nr)  Plot flux surfaces
!   fort.12  : i, a(i), rCenter(i), Shafr_shift(i), m_th_out, m_th_in
!   fort.13  : i, psi_tor(vmec)[Wb], psi_pol(gs2)[Wb], iota
!
!     Strategy: 
!     Either: Change all ncdgt (get data) to ncdpt (put data) in eqin
!             This way we can use existing gs2 without recompiling it 
!     Or    : Dump NETCDF, just use Write/Read mirror routines.
!             This will need adding the Read routine to gs2i
!             and recompiling it
!     I would rather Use second strategy because NETCDF is too complex
!               and over-powerful (do not need it)
!     This is consistent with the rest of gs2:
!             onle geq.f90 and peq.f90 use NETCDF calls

!-----------------------
      PROGRAM vmec2000_to_gs2
      USE read_wout_mod
      IMPLICIT NONE
      CHARACTER(len=128) :: wout_ext
      INTEGER :: nt, ierr

!     Read wout_ext and number of theta points from user
      print*, 'Enter wout_extension and number of theta points > '      
      read(5,*) wout_ext, nt

      call read_wout_file(wout_ext, ierr)
      print*, 'IERR = ', ierr

!     Write rmn, zmn for next vmec input run

!     Write gs2 input file
      call eqout_vmec_gs2_nc(wout_ext, nt)
     
      END PROGRAM vmec2000_to_gs2
!-----------------------
! Writes a gs2 file in toq (geq) or transp (peq) format
! Read answers from VMEC in double precision (wout gives real(rprec)
! Write them for gs2 in real*4 single precision (eqin reads real*4)

!     This version writes to a netcdf file that can be read in with
!     existing eqin routine in gs2

      SUBROUTINE eqout_vmec_gs2_nc( gs2_ext, nt )

!     USE statements MUST go before others (fortrancrap!)   
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod
      USE vparams, ONLY: mu0
      USE vmec_utils
      USE stel_constants, ONLY: zero

      IMPLICIT NONE

!     netcdf junk
      include 'netcdf.inc'
      integer :: ncid, id, ifail, istat
      integer :: rdim, zdim, psidim, sdim, tdim, nrtdim
      integer :: start(2), cnt(2), dims(2)
      real*4, DIMENSION(:),ALLOCATABLE :: workk 
      real*4 work1

      CHARACTER*(256) :: aegis_file

!     VMEC stuff
      character(Len=*) :: gs2_ext
      character*1 cdum
      REAL(rprec),DIMENSION(:,:),ALLOCATABLE :: R_psi, Z_psi, B_psi,
     1 vep, vep_arr, B_psi_half
      REAL(rprec),DIMENSION(:),ALLOCATABLE :: psi,rho_d,
     1 R_cntr,q,shat,pprim, fp, beta, phi_h,
     2 pFlxMtk,pMtk,denMtk,TeMtk,TiMtk,
     3 pFlxN,pPlas,dens,Telec,Tion
      REAL(rprec) :: Pi,FourPiInv, rDum, a_min_fac,surfArea,
     1 zdist,zd_in,zd_out, rcntr,shft,temp, tempp, R_diff, q95,
     2 th_step,theta,xm_theta,cos_theta,sin_theta,
     3 vol_tot,beta_vav
      INTEGER :: i, j, k, m, n, m_in, m_out, nr, nt, nrt
      
      INTEGER, PARAMETER :: nvep=10
      REAL(rprec) eta_vep(nvep), al_vep(nvep), Br, Bphi, Bz,
     1 betaN_troy, betaN_geo, beta_geo, press_vav, Bcntr
      DATA eta_vep/1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4., 1.e6/
      al_vep(1:nvep)=1.0/(1.0+eta_vep(1:nvep))
           
!.............................
!     Notes: P. Valanju April 2005 pvalanju@mail.utexas.edu
!     This subroutine writes a generic NetCDF equilibrium file
!     containing the axisymmetric magnetic field geometry in flux 
!     coordinates that can be read by eqin routine in gs2/geo/geq.f90

!     Note: Need to calculate R(nr,nt), Z_psi(nr,nt) on many theta=nt
!           for each flux surface=nr 
!           from VMEC Fourier coefficients at each psi=nr
!..........................
!      Double precision constants
!       Pi    = atan2(0.0d0,-1.0d0)
       Pi = 2 * dacos(0.0d0)
       FourPiInv =1.0d0/(4*Pi)
	
!..................................
!      First do calculations to turn VMEC output info to gs2 input info

!      1. Convert vmec fourier coeffs to flux surface shapes in R, Z, and |B|
!      rmnc(mnmax,ns) to  R_psi(nr,nt),  (nr,nt) = (rho,theta)
!      zmns(mnmax,ns) to  Z_psi(nr,nt),
!      bmnc(mnmax_nyq,ns) to B_psi(nr,nt)
!       PRINT*, ' Calculating R,Z,B(psi,theta) from Fourier coeffs'

       nr = ns
       nrt = nr*nt
       ALLOCATE ( R_psi(nr,nt), Z_psi(nr,nt), B_psi(nr,nt), 
     1  B_psi_half(nr,nt), fp(nr),  beta(nr),
     2  rho_d(nr), R_cntr(nr), q(nr), shat(nr), pprim(nr),
     3  phi_h(nr), vep(nr,nvep),vep_arr(nr,nvep) )
       ALLOCATE ( workk(nrt) )

!...................................
!      1. Calculate FS r,z. Note rmn, zmn are on full mesh, but
!         Note: bmnc,s are on half mesh
       write(6,*) 'Rmn, Bmn(1,1) = ', rmnc(1,1), bmnc(1,1)
       R_psi(1:nr,1:nt)=0.0d0
       Z_psi(1:nr,1:nt)=0.0d0
       B_psi(1:nr,1:nt)=0.0d0
       B_psi_half(1:nr,1:nt)=0.0d0
       IF (lasym) then  ! non-axisymmetric
         print*,"lasym is TRUE ",lasym
         th_step = 2.0d0*Pi/(nt-1)
       else
         th_step = Pi/(nt-1)
       end if
       do i = 1, nr    !over all radial grid points
        do j = 1, nt
	 theta = th_step * (j-1)
         do m = 1, mnmax
          xm_theta = xm(m) * theta
          cos_theta = dcos( xm_theta )
          sin_theta = dsin( xm_theta )
          R_psi(i,j) = R_psi(i,j) + rmnc(m,i) * cos_theta
          Z_psi(i,j) = Z_psi(i,j) + zmns(m,i) * sin_theta
          B_psi_half(i,j) = B_psi_half(i,j) + bmnc(m,i) * cos_theta
	  IF (lasym) then
            R_psi(i,j) = R_psi(i,j)+rmns(m,i) * sin_theta
            Z_psi(i,j) = Z_psi(i,j)+zmnc(m,i) * cos_theta
	    B_psi_half(i,j) = B_psi_half(i,j) + bmns(m,i) * sin_theta
	  end if
         end do
        end do
       end do

!..................................
!      Change ALL half-mesh arrays to full mesh NOW
!      Because they will be used later

!      To do this with cubic splines needs phi=toroidal flux as x-array on half and full grids
!      wout gives phi on full grid, so first get phi_h on half-grid
!      Note: This cannot be a spline interp because phi is the x-axis for interpolation
!            This will work for non-uniform grids also
       phi_h(1) = 0.0d0   ! Dummy - never used
       Do i = 2, nr
         phi_h(i) = 0.5d0*( phi(i-1) + phi(i) )
       end do      

!      Change B_psi from half-mesh to full mesh (bmnc,s are on half grid)
       temp=0.0d0
       do j = 1, nt
         call 
     1   half_to_full_grid(nr,phi_h,B_psi_half(1,j), phi,B_psi(1,j))
	 temp=temp+B_psi(1,j)
       end do
       B_psi(1,1:nt)=temp/nt

!..................................
!      2. Calculate q from iota on full mesh
       q(1:nr) = 1.0d0/iotaf(1:nr)

!..................................
!      3. rho_d(1:nr) == half diameters of flux surfaces at elevation cm
!                        of Z_mag on the radial grid (minor radius)
!         We allow for a vertical shift of magnetic axis
!         These are all on full mesh.
       rho_d(1)  = 0.0d0
       R_cntr(1) = R_psi(1,1)
       do i = 2, nr
        zd_in  = 1.0d20
	zd_out = 1.0d20
        do j = 1, nt
	 zdist  = DABS( Z_psi(i,j)-Z_psi(1,1) )
	 if( R_psi(i,j) .gt. R_psi(1,1) ) then
	  if ( zd_out .gt. zdist ) then
	   zd_out = zdist
	   m_out  = j
	  end if
	 else
	  if ( zd_in .gt. zdist ) then
	   zd_in = zdist
	   m_in   = j
	  end if
	 end if
	end do
	rho_d(i)  = 0.5d0*(R_psi(i,m_out)-R_psi(i,m_in))
        R_cntr(i) = 0.5d0*(R_psi(i,m_out)+R_psi(i,m_in))
!        write(99,*) i, m_in, m_out, rho_d(i)
       end do
!      Note: remember the last m_in and m_out calculated at nr
!            for later use in calculating rho.
       
!..................................
!      4. Calculate beta on full mesh (presf is on full mesh):
!      This beta matches the one used in gs2
       beta = 2 *mu0 * presf *((R_cntr(nr)/rbtor)**2)
!       beta = 2 *mu0 * presf /(b0**2)
!      Do NOT use vol aver beta for gs2: gives error in shat       
!       call half_to_full_grid( nr, phi_h,beta_vol, phi,beta)

!..................................
!      5. Since psi for gs2 is the poloidal flux,
!      while vmec spits out only toroidal flux, so convert
!      psi(n) = Sum[ iota(i) * (phi(i+1)-phi(i)) ]/(2 Pi)
!      Note: iotaf and phi are both on full mesh
!       PRINT*, ' Converting phi_toroidal (vmec) to psi_poloidal (gs2)'
       ALLOCATE ( psi(1:nr) )
       psi(1) = 0.0d0
       do i = 2, nr      !Trapezoid: can be unequally spaced points
        j=i-1
        psi(i)=psi(j) +FourPiInv*(iotaf(i)+iotaf(j))*(phi(i)-phi(j))
       end do
      
!..................................
!     6. Calculate shat, pprim, vep, beta all on full mesh
      do i = 2, nr-1
       j=i+1
       m=i-1
       shat(i)=rho_d(i)*iotaf(i)*(q(j)-q(m))/(rho_d(j)-rho_d(m))
       pprim(i)=(beta(j)-beta(m))/(rho_d(j)-rho_d(m))
       do k = 1, nvep
        vep_arr(i,k)=pprim(i)/(beta(i)**al_vep(k))
       end do
       pprim(i)=pprim(i)/beta(i)
      end do
      do k = 1, nvep
       vep_arr(1,k)=vep_arr(2,k)
       vep_arr(nr,k)=vep_arr(nr-1,k)
      end do

      do i = 2, nr-1
       j=i+1
       m=i-1
       temp=(rho_d(j)-rho_d(m))
       do k = 1, nvep
        vep(i,k)=(beta(i)**(al_vep(k)-1))*
     1     ( vep_arr(j,k)-vep_arr(m,k) )/temp
       end do
      end do
      shat(1)=0
      shat(nr)=shat(nr-1)
      pprim(1)=0
      pprim(nr)=pprim(nr-1)
      do k = 1, nvep
       vep(1,k)=vep(2,k)
       vep(nr,k)=vep(nr-1,k)
      end do

!...........................................................
!     Write to files for gnu plotting VMEC 2-D & 3-D arrays:

      open(unit=31,file='ExB_'//TRIM(gs2_ext))
      write(31,'(a,10(1x,g12.4))') '# alp',(al_vep(k),k=1,nvep)
      write(31,'(a,10(1x,g12.4))') '# eta',(eta_vep(k),k=1,nvep)
      do i = 1, nr
       write(31,'(11(1x,g12.4))')rho_d(i)/rho_d(nr),(vep(i,k),k=1,nvep)
!       write(331,*) i, rho_d(i), beta(i), vep(i,3), vep(i,7)
      end do
      close(unit=31)

      vol_tot  = SUM( vp(2:nr) )

!      beta_vav = SUM( beta(2:nr) * vp(2:nr) )
!      beta_vav = ABS( beta_vav / vol_tot)
!      Get beta_vav from threed1
       call system(
     1 'grep -i "^ *beta total *=" threed1.'
     2//TRIM(gs2_ext)//" | awk '{print $NF}' > tmp")
       open(29, file='tmp')
       read(29,*) beta_vav
       close(29)

      Bcntr=b0*R_cntr(1)/R_cntr(nr)
      betaN_troy=dabs(beta_vav/(ITor*1.0D-6/(rho_d(nr)*Bcntr)))

!      Get press_vav from threed1
       call system(
     1 'grep -i "^ *pressure *=" threed1.'
     2//TRIM(gs2_ext)//" | awk '{print $NF}' > tmp")
       open(29, file='tmp')
       read(29,*) press_vav
       close(29)
      beta_geo=2*mu0 *press_vav /(Bcntr**2)
      betaN_geo=dabs(beta_geo/(ITor*1.0D-6/(rho_d(nr)*Bcntr)))

      vol_tot  = vol_tot *4*Pi*Pi/(nr-1)

!     Find q95 (PMV 28 Jan 2015)
      q95=3.0
      temp=1.0E20
      do i = 1, nr
       tempp = abs( 0.95 - abs(psi(i)/psi(nr)) )
       if( tempp .lt. temp ) then
         temp = tempp
         q95 = q(i)
       end if
      end do

!     Write r,z coeffs of outermost surface for use in input file      
      open(unit=22,file='rzmn_'//TRIM(gs2_ext))
      write(22,*) 'RAXIS(0)=', Rmajor,', ZAXIS(0)=0.0'
      do i = 1, mnmax 
        m = xm(i)
	n = xn(i)
        write(22,'(a,i4,a,i4,a,g22.14,a,i4,a,i4,a,g22.14)')
     1   'RBC(',n,',',m,')=',rmnc(i,nr),
     2   ', ZBS(',n,',',m,')=',zmns(i,nr)  
	IF (lasym) 
     1   write(22,'(a,i4,a,i4,a,g22.14,a,i4,a,i4,a,g22.14)')
     2   'RBS(',n,',',m,')=',rmns(i,nr),
     3   ', ZBC(',n,',',m,')=',zmnc(i,nr)  
      end do
      close(unit=22)
      
      open(unit=12,file='gp2_'//TRIM(gs2_ext))
      write(12,'(a,3(g12.4,1x))')
     1 '# B0, Bcntr = ', b0, Bcntr
      write(12,'(a,4(g12.4,1x),i)')
     1 '# Itor [A], R0, Rctr, a, nr = ', 
     2 ITor,100*R_psi(1,1),100*R_cntr(nr),100*rho_d(nr),nr
      write(12,'(a,3(g12.4,1x))')
     1 '# Phi_pol, Psi_tor max, q95 = ', phi(nr), psi(nr), q95
      write(12,'(a,5(g12.4,1x))') 
     1'# beta_Ax, beta_avg, betaN_Troy, betaN_geo, Vol = ',
     2 beta(1), beta_vav, betaN_troy, betaN_geo, vol_tot
      write(12,'(10a)')    '#T 1.r/a   ','2.Phi_Tor  ','3.Psi_pol  ',
     1 '4.q-prof   ','5.s-hat    ','6.beta  ','7.Pres[PA] ',
     1 '8.Pprime   ','9.<j*B>    ','10.Volume  '

!     Added 15 Oct 2013, for jbs loop
      open(unit=41,file='gPOL_'//TRIM(gs2_ext))
      write(41,'(6a20)') '# r/a     ','PhiTor    ',
     1  'PsiPol    ','<j.B>    ','<jTor>  ','vp'
      do i = 1, nr
       write(41,'(6(g20.12,1x))') rho_d(i)/rho_d(nr),
     1  phi(i)/phi(nr),psi(i)/psi(nr),
     2  dabs(jdotb(i)),dabs(jcurv(i)),vp(i)
      enddo
      close(unit=41)

!     Added 5 Dec 2013, for jbs_mtk_geom_in
      open(unit=42,file='jbs_mtk_geom_in_'//TRIM(gs2_ext))
      a_min_fac=sqrt(dabs((phi(nr)-phi(1))/(Pi*b0)))/rho_d(nr)
!     Get outermost area from threed1.new
       call system(
     1 'grep "Normal Surface Area" threed1.'
     2//TRIM(gs2_ext)//" | awk '{print $5}' > tmp")
       open(29, file='tmp')
       read(29,*) surfArea
       close(29)
!      write(42,*) '#   nr-1,  poloidal psi(nr)-psi(1),   a_min_fac,  surfArea'
      write(42,'(a,i5,3g20.12)') '# ',nr-1,
     1 psi(nr)-psi(1), a_min_fac, surfArea
!      write(42,*) b0, phi(nr)-phi(1), rho_d(nr)
!      write(42,*) '#    R(m),      a(m),       q,      r*BT'
      do i = 1, nr
        write(42,'(4(g20.12,1x))')
     1  R_cntr(i), rho_d(i), q(i), rbtor0
      enddo
      close(unit=42)

!     Added 26 Jun 2014, to send fixed bdry eqlb to corsica
      open(unit=43,file='toCorsica_'//TRIM(gs2_ext))
!      write plasma boundary R,Z (cm) at nt theta angles
       write(43,*) 2*nt-1
        do j = 1, nt
	 write(43,*) 100*R_psi(nr,j),100*Z_psi(nr,j)
	end do
        do j = nt-1, 1, -1
	 write(43,*) 100*R_psi(nr,j),-100*Z_psi(nr,j)
	end do
!      Write plasma profiles. Corsica cgs press= 10*vmec MKS presf
       write(43,*) nr, ITor/1.0D6
        do i = 1, nr
         write(43,*) psi(i)/psi(nr),
     1    10.0*presf(i), dabs(jdotb(i))
	end do
      close(unit=43)

!     Added 15 Oct 2013, for jbs loop

      open(unit=44,file='prof_mtk_out',status='OLD',iostat=i)
      if(i.eq.0) then
!       Read 1.pFlxN, 2.p, 3.n, 4.T_i, 5.T_e from prof_mtk_out, and interpolate to nr grid
        read(44,*) cdum, j
        read(44,*) cdum
!        j=j+1
        ALLOCATE ( 
     1    pFlxMtk(j),pMtk(j),denMtk(j),TeMtk(j),TiMtk(j),
     2    pFlxN(nr),pPlas(nr),dens(nr),Telec(nr),Tion(nr) )
        do i=1,j
          read(44,*) pFlxMtk(i), pMtk(i),
     1        denMtk(i), TeMtk(i), TiMtk(i)
        enddo
        pFlxN(1:nr)=psi(1:nr)/psi(nr)  ! Normalized polFlx
        call spline_arr(j,pFlxMtk,denMtk, nr,pFlxN,dens)
        call spline_arr(j,pFlxMtk,pMtk, nr,pFlxN,pPlas)
        call spline_arr(j,pFlxMtk,TeMtk, nr,pFlxN,Telec)
        call spline_arr(j,pFlxMtk,TiMtk, nr,pFlxN,Tion)
!       Write Tion, Telec, dens etc info into file for gene
        open(unit=45,file='toGene_'//TRIM(gs2_ext))
        write(45,'(a,6a18)') '#','PhiTor','PsiPol','p(mks)',
     1    'dens(10^19)','Te(kev)','Ti(kev)'
        do i = 1, nr
          write(45,'(6(g20.12,1x))')
     1     phi(i)/phi(nr),pFlxN(i),
     2     pPlas(i),dens(i),Telec(i),Tion(i)
          enddo
        close(unit=45)
      else
        write(*,*) "No prof_mtk_out file found, no toGene written"
      endif

!     Added 30 Apr 2010, for Bmod contour gnuplot
      open(unit=23,file='Bmod_'//TRIM(gs2_ext))
        write(23,'(a)') '# Bmod contour gnuplot'
        do i = 2, nr
	  write(23,*) '# ir = ', i
          do j = 1, nt
	    write(23,'(3(g12.4,2x),i)') 
     1         sqrt(phi(i)/phi(nr)),th_step*(j-1),B_psi(i,j)
	  end do
	  write(23,*) ' '
        end do
      close(unit=23)

      open(unit=21,file='gp3_'//TRIM(gs2_ext))
      write(21,*)  '#T R(i,j),   ',' z(i,j)    ','  Bmod(i,j)','     j'
      do i = 1, nr
       write(12,'(10(g12.4,1x))') rho_d(i)/rho_d(nr),
     1  phi(i)/phi(nr), psi(i)/psi(nr),
     1  q(i), shat(i), beta(i), 10*presf(i),
     2  pprim(i), jdotb(i), vp(i)*4*Pi*Pi/(nr-1)
!     2  pprim(i), jcurv(i), vp(i)*4*Pi*Pi/(nr-1)
!       write(13,'(10(g12.4,1x))') rho_d(i),
!     1  q(i), q(i),
!     1  q(i), q(i), q(i), q(i),
!     2  R_cntr(i)-R_cntr(nr), q(i), R_cntr(i)
	write(21,*) '# ir = ', i
        do j = 1, nt
	 write(21,'(3(g12.4,2x),i)') R_psi(i,j),Z_psi(i,j),B_psi(i,j),j-1
!	 write(212,*) R_psi(i,j),Z_psi(i,j)
	end do
	write(21,*) ' '
       end do
      close(unit=12)
!      close(unit=13)
      close(unit=21)
      
!     Write radial profiles at z=0 (added Apr 15, 2011, PMV)
      open(unit=14,file='rp_'//TRIM(gs2_ext))
      open(unit=15,file='B_Ratio_test')
      write(14,*)  '#T r/a,  <j*B>,  beta, q,  Major R, Mod B,  Bz'
      do i = nr, 2, -1

!      PMV Aug 2011 added to calc 3 components of B, print Bz=Bpol
       call GetBcyl_WOUT(R_psi(i,nt), zero, zero, Br, Bphi, Bz)       
       write(15,*) i,"Bratio=", sqrt(Br**2+Bphi**2+Bz**2)/B_psi(i,nt)

       write(14,'(7(g12.4,2x))') rho_d(i)/rho_d(nr), jdotb(i), 
     1  beta(i), q(i), R_psi(i,nt), B_psi(i,nt), Bz
      end do

      rDum=0.5d0*(R_psi(1,1)+R_psi(1,nt))
      call GetBcyl_WOUT(rDum, zero, zero, Br, Bphi, Bz)       
      write(15,*) i,"Bratio=", sqrt(Br**2+Bphi**2+Bz**2)/B_psi(i,nt)
      write(14,'(7(g12.4,2x))') 
     1 rho_d(1)/rho_d(nr), jdotb(1),beta(1), q(1), 
     2 0.5d0*(R_psi(1,1)+R_psi(1,nt)),0.5d0*(B_psi(1,1)+B_psi(1,nt)),
     3 Bz

      do i = 2, nr

       call GetBcyl_WOUT(R_psi(i,1), zero, zero, Br, Bphi, Bz)       
       write(15,*) i,"Bratio=", sqrt(Br**2+Bphi**2+Bz**2)/B_psi(i,1)

       write(14,'(7(g12.4,2x))') rho_d(i)/rho_d(nr), jdotb(i), 
     1  beta(i), q(i), R_psi(i,1), B_psi(i,1), Bz
      end do

      close(unit=14)
      close(unit=15)

!...........................
!     create netCDF file: will be opened in 'define mode', in which
!     the variables in the file are defined, but not read or written
      ncid = nccre ( 'gs2_'//TRIM(gs2_ext)//'.nc', NCCLOB, ifail)
!      print*,'Created file gs2_'//TRIM(gs2_ext)//'.nc ', ncid
      
!     Also open a parallel output text file 
      open(unit=7,file='gs2_'//TRIM(gs2_ext)//'.txt')

!........................
!     create dimensions (array sizes)
      sdim= 0
      psidim = ncddef (ncid, 'psi', nr, ifail)
!      print*,'psidim ',psidim
      tdim = ncddef (ncid, 'theta', nt, ifail)
!      print*,'tdim ',tdim
!      print*,'Created dimensions (array sizes)'

!........................
!     create scalar fields
      id = ncvdef (ncid, 'psi_0', NCFLOAT, 0, sdim, ifail)
      id = ncvdef (ncid, 'psi_a', NCFLOAT, 0, sdim, ifail)
      id = ncvdef (ncid, 'B_T', NCFLOAT, 0, sdim, ifail)
      id = ncvdef (ncid, 'I_0', NCFLOAT, 0, sdim, ifail)

      id = ncvdef (ncid, 'Rmag', NCFLOAT, 0, sdim, ifail)
      id = ncvdef (ncid, 'Zmag', NCFLOAT, 0, sdim, ifail)
      id = ncvdef (ncid, 'a', NCFLOAT, 0, sdim, ifail)
!      print*,'Created scalar fields'

!........................
!     create 1d fields
      id = ncvdef (ncid, 'rho', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'psi', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'JoverR', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'psibar', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'fp', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'q', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'beta', NCFLOAT, 1, psidim, ifail)
      id = ncvdef (ncid, 'pressure', NCFLOAT, 1, psidim, ifail)
!      print*,'Created 1-d fields'

!........................
!     create 2d fields
      dims(1) = psidim
      dims(2) = tdim
      id = ncvdef (ncid, 'R_psi', NCFLOAT, 2, dims, ifail)
      id = ncvdef (ncid, 'Z_psi', NCFLOAT, 2, dims, ifail)
      id = ncvdef (ncid, 'B_psi', NCFLOAT, 2, dims, ifail)
!      print*,'Created 2-d fields'

!........................
!     Add self-documenting attributes to variables
      id = ncvid (ncid, 'psi_0', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1       20, 'psi on magnetic axis', ifail)

      id = ncvid (ncid, 'psi_a', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      22, 'psi on plasma boundary', ifail)

      id = ncvid (ncid, 'B_T', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      42, 'vacuum toroidal magnetic field at R_center', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      5, 'Gauss', ifail)

      id = ncvid (ncid, 'I_0', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      20, 'total plasma current', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      1, 'A', ifail)

      id = ncvid (ncid, 'Rmag', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      27, 'R position of magnetic axis', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      2, 'cm', ifail)

      id = ncvid (ncid, 'Zmag', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      27, 'Z position of magnetic axis', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      2, 'cm', ifail)

      id = ncvid (ncid, 'a', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      12, 'minor radius', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      2, 'cm', ifail)

      id = ncvid (ncid, 'rho', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1 59,'half diameter of flux surface at elevation of magnetic axis',
     1      ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      1, 'a', ifail)

      id = ncvid (ncid, 'JoverR', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1 11,'< J*R_0/R >',ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,1, '?', ifail)

      id = ncvid (ncid, 'psi', ifail)
      call ncaptc (ncid, id, 'long_name', NCCHAR,
     1      13, 'poloidal flux', ifail)

      id = ncvid (ncid, 'psibar', ifail)

      id = ncvid (ncid, 'fp', ifail)

      id = ncvid (ncid, 'q', ifail)

      id = ncvid (ncid, 'beta', ifail)

      id = ncvid (ncid, 'pressure', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      11, 'pressure(0)', ifail)

      id = ncvid (ncid, 'R_psi', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      2, 'cm', ifail)

      id = ncvid (ncid, 'Z_psi', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      2, 'cm', ifail)

      id = ncvid (ncid, 'B_psi', ifail)
      call ncaptc (ncid, id, 'units', NCCHAR,
     1      2, '??', ifail)

!     leave 'define mode'
      call ncendf (ncid, ifail)

!........................................................
!     Write data to netCDF file that has now been set up (formatted)
!     Convert ALL VMEC in MKS, Tesla to cm Gauss, flux -> 1.0e8 *flux
!     Simultaneously write them into the text file unit = 7

!........................
!     Write scalars:

      id = ncvid (ncid, 'psi_0', ifail)
      work1 = -1.0e8*psi(1)
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'psi_0 = ', work1

      id = ncvid (ncid, 'psi_a', ifail)
      work1 = -1.0e8*psi(nr)
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'psi_a = ', work1

      id = ncvid (ncid, 'B_T', ifail)
!      work1 = Abs(1.0e4*b0)
      work1 = Abs(1.0e4*rbtor/R_cntr(nr))
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'B_T = ', work1, b0, R_cntr(1), R_cntr(nr)

      id = ncvid (ncid, 'I_0', ifail)
      work1 = Itor
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'I_0 = ', work1

      id = ncvid (ncid, 'Rmag', ifail)
      work1 = 100*R_psi(1,1)
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'Rmag = ', work1

      id = ncvid (ncid, 'Zmag', ifail)
      work1 = 100*Z_psi(1,1)
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'Zmag = ', work1

      id = ncvid (ncid, 'a', ifail)
      work1 = 100*rho_d(nr)
      call ncvpt1 (ncid, id, 1, work1, ifail)
      write(7,*)'a = ', work1

      write(6,*)'Wrote scalars'

!........................
!     Write vectors(1:nr), all on full mesh:

! 1.  Normalized rho (goes from 0 to 1)
!     Note: use the last m_in and m_out calculated above at nr
      id = ncvid (ncid, 'rho', ifail)
      R_diff = R_psi(nr,m_out)-R_psi(nr,m_in)
      do i=1,nr
         workk(i) = (R_psi(i,m_out) - R_psi(i,m_in))/R_diff
      enddo
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'rho = ', workk(1:nr)
      
!  2. J over R
      id = ncvid (ncid, 'JoverR', ifail)
      workk(1:nr)=0.     !Set to zero for now
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'JoverR = ', workk(1:nr)

!  3. Poloidal flux psi in gs2 units
      id = ncvid (ncid, 'psi', ifail)
      workk(1:nr) = -1.0e8*psi(1:nr)
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'psi = ', workk(1:nr)

!  4. Normalized psi going from 0 to 1
      id = ncvid (ncid, 'psibar', ifail)
      do i=1,nr
         workk(i) = (psi(i)-psi(1))/(psi(nr)-psi(1))
      enddo
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'psibar = ', workk(1:nr)

!  5. fprime = bvco on full mesh
      id = ncvid (ncid, 'fp', ifail)
      workk(1:nr) = -1.0e6*fp(1:nr)
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'fp = ', workk(1:nr)

!  6. q profile
      id = ncvid (ncid, 'q', ifail)
      workk(1:nr) = q(1:nr)
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'q = ', workk(1:nr)

!  7. beta (NOT beta_vol)
      id = ncvid (ncid, 'beta', ifail)
      workk(1:nr) = beta(1:nr)
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'beta = ', workk(1:nr)

!  8. Pressure:  eiktest needs pressure normalized to 1 at center
!      workk(1:nr) = 10.0*presf(1:nr) would be Pressure in cgs units
      id = ncvid (ncid, 'pressure', ifail)
      workk(1:nr) = presf(1:nr)/presf(1)    !Normalized pressure
      call ncvpt (ncid, id, 1, nr, workk, ifail)
      write(7,*)'pressure = ', workk(1:nr)

      write(6,*)'Wrote vectors'

!      write(112,*)    '#T 1.r/a    ','2.eqpsi  ','3.Psi_bar   ',
!     1 '4.fprime    ','5.qsf       ','6.beta   ','7.Presssur  '
!      do i = 1, nr
!       write(112,'(7(g11.5,1x))') rho_d(i),-1.0e8*psi(i),
!     1  psi(i)/psi(nr),-1.0e6*fp(i),q(i),beta(i),presf(i)/presf(1)
!      end do

!........................
!     Write 2-D matrices:

      start(1) = 1
      start(2) = 1
      cnt(1) = nr
      cnt(2) = nt
      nrtdim=nr*nt
      id = ncvid (ncid, 'R_psi', ifail)
      do j=1,nt
         do i=1,nr
            workk(1+i-1+nr*(j-1)) = 100.0*R_psi(i,j)
         enddo
      enddo
      call ncvpt (ncid, id, start, cnt, workk, ifail)
      write(7,*)'R_psi = ', workk(1:nrtdim)

      id = ncvid (ncid, 'Z_psi', ifail)
      do j=1,nt
         do i=1,nr
            workk(1+i-1+nr*(j-1)) = 100.0*Z_psi(i,j)
         enddo
      enddo
      call ncvpt (ncid, id, start, cnt, workk, ifail)
      write(7,*)'Z_psi = ', workk(1:nrtdim)

      id = ncvid (ncid, 'B_psi', ifail)
      do j=1,nt
         do i=1,nr
            workk(1+i-1+nr*(j-1)) = ABS(1.0e4*B_psi(i,j))
         enddo
      enddo
      call ncvpt (ncid, id, start, cnt, workk, ifail)
!      print*,'Wrote 2-D matrices'
      write(7,*)'B_psi = ', workk(1:nrtdim)

      call ncclos (ncid, ifail)
      
      close(unit=7)

! ----------------------
! Write aegis input file
      IF (lasym) then
        write(*,*) 'No aegis file written: LASYM=T: NOT AXISYMMETRIC!'
      else
        aegis_file = "aegis_gs_" // TRIM(input_extension) // ".txt"
        OPEN (unit=51,FILE=aegis_file,FORM='FORMATTED',iostat=istat)
        IF (istat .ne. 0) STOP 'Error writing aegis output file'

        WRITE (51,*) 2*nt-1, nr   
      
        WRITE (51,*) psi(1:nr)         !pol flux, full mesh (included 2*pi factor)
        WRITE (51,*) phi(1:nr)         !R*BT, full mesh
        WRITE (51,*) presf(1:nr)/mu0   !pressure, full mesh (MKS units)
        WRITE (51,*) q(1:nr)           !q, full mesh      
           
        do i = 1, nr
          WRITE (51,*) ( R_psi(i,j),j=1,nt)
          WRITE (51,*) ( R_psi(i,j),j=nt-1,1,-1)
          WRITE (51,*) ( Z_psi(i,j),j=1,nt)
          WRITE (51,*) (-Z_psi(i,j),j=nt-1,1,-1)
        end do 

        CLOSE (unit=51)
      end if

      END SUBROUTINE eqout_vmec_gs2_nc

!-------------------------------------------------------------------------
!      Spline-subroutines to go between half and full mesh
!      Note: These work for any x-arr (tor or pol flux, non-uniform, etc)
!            xh(1),yh(1) is alwyas dummy in vmec

       SUBROUTINE half_to_full_grid( nr, xh,yh, xf,yf )
       USE stel_kinds, ONLY: rprec
       integer :: i, nr
       REAL(rprec) :: xh(nr),yh(nr), xf(nr),yf(nr)

       call spline_arr(nr-1,xh(2:nr),yh(2:nr), nr,xf,yf)
       
       END SUBROUTINE        
!-----------------------

       SUBROUTINE full_to_half_grid( nr, xf,yf, xh,yh )
       USE stel_kinds, ONLY: rprec
       integer :: i, nr
       REAL(rprec) :: xh(nr),yh(nr), xf(nr),yf(nr)

       call spline_arr(nr,xf,yf, nr-1,xh(2:nr),yh(2:nr))
       yh(1) = 0.0d0   ! Dummy - never used
       
       END SUBROUTINE        
! --------------------------------------------------------------------------
!     From one array get another by spline interpolation
!     PMV Nov 2013
      subroutine spline_arr(n,x,y, m,xx,yy)

      USE stel_kinds
      implicit none
     
      integer, intent(in) :: n,m
      real(rprec), dimension(n), intent(in) :: x,y
      real(rprec), dimension(m), intent(in) :: xx,yy
!     local -------------------
      real(rprec) :: yp1,ypp1, ypn,yppn
      real(rprec), dimension(n)  :: sc
      integer :: i

!     fix boundary derivatives by quadratic fit
!     left part
      call quad_fit(x(1),y(1),x(2),y(2),x(3),y(3),yp1,ypp1)
!     right part
      call quad_fit(x(n),y(n),x(n-1),y(n-1),x(n-2),y(n-2),ypn,yppn)

!     Calculate spline coeffs sc for (x,y) array only one time
      call spline_nr(x,y,n,yp1,ypn,sc)

!     Interpolate to xx array point by point
      do i = 1, m
        call splint_ep(x,y,sc,n, xx(i),yy(i))
      end do

      return
      end subroutine spline_arr

! ---------------------------------------
      subroutine splint_ep(xx,yy,sc,n, x,y)

!     Given spline coeffs, find value at x
!     If x outside range, extrapolate using the last cubic spline

      USE stel_kinds
      implicit none
     
      integer, intent(in) :: n
      real(rprec), intent(in) :: x
      real(rprec), intent(out) :: y
      real(rprec), dimension(n), intent(in) :: xx, yy, sc
!     local -------------------
      real(rprec) :: a, b, h
      integer :: k, khi, klo

      if(x < xx(1) ) then
!     x before left point
        klo=1
        khi=2
      else if(x > xx(n) ) then
!     x after right point
        klo=n-1
        khi=n
      else
!     x within range, use splint_nr(xx,yy,sc,n,x,y)
        call splint_nr(xx,yy,sc,n,x,y)
        return
      endif

!     Extrapolate last cubic spline to points beyond range      
      h=xx(khi)-xx(klo)
      if(h ==0.d+0)
     >    stop "spline_cubic: bad xx input! xx(i) have to be distinct!"
      a=(xx(khi)-x)/h
      b=(x-xx(klo))/h
      y=a*yy(klo)+b*yy(khi) +
     >  ((a**3-a)*sc(klo)+(b**3-b)*sc(khi))*(h**2)/6.d+0

      return
      end subroutine splint_ep

! ---------------------------------------
      subroutine quad_fit(x1,y1,x2,y2,x3,y3,b,c)

      USE stel_kinds
      implicit none
      real(rprec), intent(in) :: x1,y1,x2,y2,x3,y3
      real(rprec), intent(out) :: b,c
      real(rprec) :: d
     
!     Fit a quadratic to 3 points (xi,yi)
!     y=y1+b(x-x1)+c(x-x1)**2  solve for 
!     slope y'=b, y''=2c at x=x1
!
!     Solves for b,c
!     y2-y1=b(x2-x1)+c(x2-x1)^2
!     y3-y1=b(x3-x1)+c(x3-x1)^2
!     Det d=(x2-x1)(x3-x1)^2 -(x3-x1)(x2-x1)^2
!           =(x2-x1)(x3-x1)(x3-x1-x2+x1)
!           =(x1-x2)(x2-x3)(x3-x1)
     
      d=(x1-x2)*(x2-x3)*(x3-x1)
      b=( (y2-y1)*(x3-x1)**2-(y3-y1)*(x2-x1)**2)/d
      c=(-(y2-y1)*(x3-x1)   +(y3-y1)*(x2-x1)   )/d

      return
      end subroutine quad_fit

! -----------------------------------------
! Following routines are exxctly from LIBSTELL, brought here for compiling with ncmpf
! ---------------------------------------
      subroutine spline_cubic(x,y,xx,yy,n,iflag)
      USE stel_kinds
      implicit none
!----------------------------------------------------------------------
! iflag = 0   normal return
! iflag =-1   x-request outside of bounds
! iflag =-2   xx arrays with two equal entries or not in correct order
!----------------------------------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: n
      integer, intent(inout)  :: iflag
      real(rprec), intent(in)  :: x
      real(rprec), dimension(n), intent(in)  :: xx,yy
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec), intent(out) :: y
      real(rprec), dimension(n)  :: y2, dxx
      real(rprec) :: yp1, ypn
      real(rprec) :: c
!-----------------------------------------------

      iflag = 0  !initialization
      if(x < xx(1) .or. x > xx(n)) then
        iflag = -1
        y=0.d+0
        return
      endif
      dxx(1:n-1)=xx(2:n)-xx(1:n-1)
      if(minval(dxx(1:n-1)) <= 0.d+0) then
        iflag=-2
        return
      endif

! fix boundary derivatives by quadratic fit
! left part
      c=((yy(3)-yy(1))/(xx(3)-xx(1))-(yy(2)-yy(1))/(xx(2)-xx(1)))
     >  /(xx(3)-xx(2))
      yp1=(yy(2)-yy(1))/(xx(2)-xx(1))-c*(xx(2)-xx(1))
! right part
      c=((yy(n-2)-yy(n))/(xx(n-2)-xx(n)) 
     >  -(yy(n-1)-yy(n))/(xx(n-1)-xx(n)))
     >  /(xx(n-2)-xx(n-1))
      ypn=(yy(n-1)-yy(n))/(xx(n-1)-xx(n))-c*(xx(n-1)-xx(n))

      call spline_nr(xx,yy,n,yp1,ypn,y2)
      call splint_nr(xx,yy,y2,n,x,y)

      return

      end subroutine spline_cubic

      subroutine spline_nr(x,y,n,yp1,ypn,y2)

      USE stel_kinds
      implicit none
     
! taken from numerical recipes f77 and recoded.
! Given the arrays x(1:n) and y(1:n) containing the tabulated function
! with x(1) < x(2) <...< x(n) and given values yp1 and ypn for the first
! derivative of the interpolating function at points q and n, respectively,
! this routine returns an array y2(1:n) of length n which contains the
! second derivatives of the interpolating function at the tabulated points x(i).
! If yp1 and/or ypn are equatl to 1ed+30 or larger, the routine is signaled
! to set the correspoinding boundary condition for a natural spline with zero
! derivative on that boundary.
! nmax is the largest anticipated value of n.
      integer, intent(in) :: n
      integer, parameter :: nmax=500
      real(rprec), intent(in) :: yp1, ypn
      real(rprec), dimension(n), intent(in) :: x, y
      real(rprec), dimension(n), intent(out) :: y2
      integer :: i,k
      real(rprec) :: p, qn, sig, un
      real(rprec), dimension(n) :: u

      if(yp1 > .99d+30) then
        y2(1)=0.d+0
        u(1) =0.d+0
      else
        y2(1)=-0.5d+0
        u(1)=(3.d+0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d+0
        y2(i)=(sig-1.d+0)/p
        u(i)=(6.d+0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     >       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn > .99d+30)then
        qn=0.d+0
        un=0.d+0
      else
        qn=0.5d+0
        un=(3.d+0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d+0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end subroutine spline_nr

      subroutine splint_nr(xx,yy,y2a,n,x,y)
! Given the arrays xx(1:n) and yy(1:n) of length n, which tabulate a function
! with the xx(i)'s in order), and given the array y2a(1:n), which is the
! output from spline above, and given a value of x, this routine returns
! a cubic-spline interpolated value y.
      USE stel_kinds
      implicit none
     
      integer, intent(in) :: n
      real(rprec), intent(in) :: x
      real(rprec), intent(out) :: y
      real(rprec), dimension(n), intent(in) :: xx, yy, y2a
!- local -------------------
      real(rprec) :: a, b, h
      integer :: k, khi, klo

      klo=1
      khi=n
      do
        if (khi-klo <= 1) exit  !inverted num.rec. condition for endless do-loop exit
        k=(khi+klo)/2
        if(xx(k) > x) then
          khi=k
        else
          klo=k
        endif
      enddo

      h=xx(khi)-xx(klo)
      if(h ==0.d+0)
     >    stop "spline_cubic: bad xx input! xx(i) have to be distinct!"
      a=(xx(khi)-x)/h
      b=(x-xx(klo))/h
      y=a*yy(klo)+b*yy(khi) +
     >  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d+0
      return
      end subroutine splint_nr


