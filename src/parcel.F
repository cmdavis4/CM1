  MODULE parcel_module

  implicit none

  private
  public :: parcel_driver,parcel_interp,parcel_write,setup_parcel_vars,getparcelzs
  public :: parcel_massprofile,parcel_resolve_box
  public :: parcel_alloc_local,parcel_append_local,parcel_global_count,parcel_fixup_2d
  public :: parcel_emit,parcel_gather,parcel_gather_positions,parcel_write_par
  public :: np_loc

  !-------------------------------------------------------------------------
  !  Distributed-parcel state (used only when nint(var17).eq.1).
  !  In distributed mode each rank stores ONLY the parcels physically inside
  !  its subdomain, in rows 1..np_loc of pdl, with global parcel id pidl.
  !  In legacy mode (var17=0) these are unused; the stock argument-passed
  !  pdata(nparcels,npvals) is used instead.
  !-------------------------------------------------------------------------
  integer :: np_loc = 0                  ! # parcels owned by this rank
  integer :: np_max = 0                  ! allocated local capacity
  real,    allocatable, target :: pdl(:,:)   ! local parcel data (np_max,npvals)
  integer, allocatable :: pidl(:)        ! global parcel id  (np_max)

  !-------------------------------------------------------------------------
  !  Scalar-point output registry entry (see parcel_interp).  Each entry
  !  pairs a scalar-point source field with the parcel-output column it
  !  feeds.  This lets the scalar interpolation be a single loop instead of
  !  one hand-written, shape-aware tri_interp call per variable -- adding a
  !  new scalar output is then one line, with no grid-shape/weight bookkeeping.
  !-------------------------------------------------------------------------
  type p3d
    real, pointer :: f(:,:,:)
  end type p3d

  CONTAINS

      subroutine parcel_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,zs,    &
                               sigma,sigmaf,znt,rho,ua,va,wa,pdata)
      use input
      use constants
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine updates the parcel locations
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: znt
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(nparcels,npvals), target :: pdata

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d
      integer :: nrkp
      real :: dt2,uu1,vv1,ww1
      real :: z0,rznt,var
      real :: sigdot,sig1,zsp,sig3d

      ! distributed parcels (var17=1): operate on the local arrays (pdl) via
      ! a pointer, so the entire interpolation/RK2 body below is reused
      ! unchanged and legacy/distributed produce identical trajectories.
      real, pointer :: pw(:,:)
      integer :: npp,ngc

      logical, parameter :: debug = .false.

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      ! no slip:
      if( imove.eq.1 )then
        IF(j.le.je)THEN
        do i=ib,ie+1
          ua(i,j,0) = 0.0 - umove
        enddo
        ENDIF
        do i=ib,ie
          va(i,j,0) = 0.0 - vmove
        enddo
      else
        IF(j.le.je)THEN
        do i=ib,ie+1
          ua(i,j,0) = 0.0
        enddo
        ENDIF
        do i=ib,ie
          va(i,j,0) = 0.0
        enddo
      endif
    ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = 0.0
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      IF(j.le.je)THEN
      do i=ib,ie
        wa(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, update it's location:

    dt2 = dt/2.0

    ! select working array: legacy = global pdata (all nparcels);
    ! distributed = this rank's local pdl (np_loc owned parcels).
    IF( nint(var17).eq.0 )THEN
      pw  => pdata
      npp =  nparcels
    ELSE
      pw  => pdl
      npp =  np_loc
    ENDIF

    nploop:  &
    DO np=1,npp

      x3d = pw(np,prx)
      y3d = pw(np,pry)
      if( .not. terrain_flag )then
        z3d = pw(np,prz)
      else
        sig3d = pw(np,prsig)
      endif

      iflag = -100
      jflag = -100
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit1:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit1

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pw(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pw(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                     (jflag.ge.1.and.jflag.le.nj) )THEN

      rkloop:  DO nrkp = 1,2

      IF( nrkp.eq.1 )THEN
        i=iflag
        j=jflag
      ELSE
        iflag = -100
        jflag = -100
        IF(nx.eq.1)THEN
          iflag = 1
        ELSE
          ! cm1r19:
          i = ni+2
          do while( iflag.lt.0 .and. i.gt.0 )
            i = i-1
            if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
              iflag = i
            endif
          enddo
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          jflag = 1
        ELSE
          do j=0,nj+1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) ) jflag=j
          enddo
          ! cm1r19:
          j = nj+2
          do while( jflag.lt.0 .and. j.gt.0 )
            j = j-1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
              jflag = j
            endif
          enddo
        ENDIF
        i=iflag
        j=jflag
      ENDIF

        IF(debug)THEN
        if( i.lt.0 .or. i.gt.(ni+1) .or. j.lt.0 .or. j.gt.(nj+1) )then
          print *,'  myid,i,j = ',myid,i,j
          print *,'  x,x1     = ',x3d,pw(np,prx)
          print *,'  y,y1     = ',y3d,pw(np,pry)
          do i=0,ni+1
            print *,i,abs(xh(i)-x3d),0.5*dx*ruh(i)
          enddo
          do j=0,nj+1
            print *,j,abs(yh(j)-y3d),0.5*dy*rvh(j)
          enddo
          call stopcm1
        endif
        ENDIF

        kflag = 1
        if( .not. terrain_flag )then
          do while( z3d.ge.zf(iflag,jflag,kflag+1) )
            kflag = kflag+1
          enddo
        else
          do while( sig3d.ge.sigmaf(kflag+1) )
            kflag = kflag+1
          enddo
        endif

        IF(debug)THEN
        if( kflag.le.0 .or. kflag.ge.(nk+1) )then
          print *,myid,nrkp
          print *,iflag,jflag,kflag
          print *,pw(np,prx),pw(np,pry),pw(np,prz)
          print *,x3d,y3d,z3d
          print *,uval,vval,wval
          print *,zf(iflag,jflag,kflag),z3d,zf(iflag,jflag,kflag+1)
          print *,'  16667 '
          call stopcm1
        endif
        ENDIF

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        if( .not. terrain_flag )then
          if( z3d.lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
          rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        else
          if( sig3d.lt.sigma(k) )then
            k=k-1
          endif
          rz = ( sig3d-sigma(k) )/( sigma(k+1)-sigma(k) )
        endif

        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )

        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),x3d,xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( .not. terrain_flag )then
          if( z3d.lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( sig3d.lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),y3d,yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        rx = rxs
        ry = rys
        if( .not. terrain_flag )then
          rz = ( z3d-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )
        else
          rz = ( sig3d-sigmaf(k) )/( sigmaf(k+1)-sigmaf(k) )
        endif

        ! saveit:
        rzw = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),z3d,zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa,wval)

        if( terrain_flag )then
          call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa,sigdot)
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,zs,zsp)
          z3d = zsp + sig3d*((zt-zsp)*rzt)
        endif

!----------------------------------------------------------------------
!  uv for parcels below lowest model level:

      IF( bbc.eq.3 )THEN
        ! semi-slip lower boundary condition:
        if( z3d.lt.zh(1,1,1) )then
          ! re-calculate velocities if parcel is below lowest model level:
          !------
          ! u at lowest model level:
          i=iflag
          j=jflag
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,ua(ib,jb,1),uval)
          !------
          ! v at lowest model level:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,va(ib,jb,1),vval)
          !------
          ! z0:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt,z0)
          !------
          ! get u,v from (neutral) log-layer equation:
          rznt = 1.0/z0
          var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
          if( imove.eq.1 )then
            uval = (uval+umove)*var - umove
            vval = (vval+vmove)*var - vmove
          else
            uval = uval*var
            vval = vval*var
          endif
        endif
      ENDIF

!-----------------------------------------------------
!  Update parcel positions:
!-----------------------------------------------------

      ! RK2 scheme:
      IF(nrkp.eq.1)THEN
        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pw(np,prx)+dt*uval
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pw(np,pry)+dt*vval
        ENDIF
        if( terrain_flag )then
          sig3d = pw(np,prsig) + dt*sigdot
          sig1 = sigdot
        else
          z3d = pw(np,prz)+dt*wval
          ww1=wval
        endif
        uu1=uval
        vv1=vval
      ELSE
        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pw(np,prx)+dt2*(uu1+uval)
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pw(np,pry)+dt2*(vv1+vval)
        ENDIF
        if( terrain_flag )then
          sig3d = pw(np,prsig) + dt2*(sig1+sigdot)
          IF( sig3d.lt.0.0 )THEN
            print *,'  parcel is below surface:  np,x3d,y3d,sig3d = ',np,x3d,y3d,sig3d
            sig3d=1.0e-6
          ENDIF
          sig3d=min(sig3d,maxz)
        else
          z3d = pw(np,prz)+dt2*(ww1+wval)
          IF( z3d.lt.0.0 )THEN
            print *,'  parcel is below surface:  np,x3d,y3d,z3d = ',np,x3d,y3d,z3d
            z3d=1.0e-6
          ENDIF
          z3d=min(z3d,maxz)
        endif
      ENDIF


      ENDDO  rkloop

!-----------------------------------------------------
!  Account for boundary conditions (if necessary)
!-----------------------------------------------------

        ! New for cm1r17:  if parcel exits domain,
        ! just assume periodic lateral boundary conditions
        ! (no matter what actual settings are for wbc,ebc,sbc,nbc)
        ! NOTE: in distributed mode (var17=1) the wrap is DEFERRED to
        ! parcel_migrate, so an edge parcel keeps its un-wrapped position
        ! and routes to the correct (periodic) neighbor rank.
      IF( nint(var17).eq.0 )THEN
        if(x3d.lt.minx)then
          x3d=x3d+(maxx-minx)
        endif
        if(x3d.gt.maxx)then
          x3d=x3d-(maxx-minx)
        endif

        if( (y3d.gt.maxy).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d-(maxy-miny)
        endif
        if( (y3d.lt.miny).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d+(maxy-miny)
        endif
      ENDIF

        pw(np,prx)=x3d
        pw(np,pry)=y3d
        if( .not. terrain_flag )then
          pw(np,prz)=z3d
        else
          pw(np,prsig)=sig3d
        endif

#ifdef MPI

      ELSE

        ! set to really small number (so we can use the allreduce command below)
        pw(np,prx) = -1.0e30
        pw(np,pry) = -1.0e30
        if( .not. terrain_flag )then
          pw(np,prz) = -1.0e30
        else
          pw(np,prsig) = -1.0e30
        endif
#endif

      ENDIF  myparcel

    ENDDO  nploop

!----------------------------------------------------------------------
!  communicate data  (for MPI runs)

#ifdef MPI
      IF( nint(var17).eq.0 )THEN
        ! legacy: every rank computed its owned parcels (others sentineled);
        ! merge so all ranks hold all parcels.
        if( .not. terrain_flag )then
          call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,1),3*nparcels,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        else
          call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,1),2*nparcels,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,prsig),nparcels,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        endif
      ELSE
        ! distributed: hand boundary-crossing parcels to neighbor ranks.
        call parcel_migrate(xf,yf)
        ! optional per-step conservation check (var15=1): the global parcel
        ! count must stay == nparcels (no lost or duplicated parcels).
        IF( nint(var15).eq.1 )THEN
          ngc = parcel_global_count()
          if( myid.eq.0 ) print *,'  [parcel] count, np_loc(0), nparcels =',ngc,np_loc,nparcels
          if( ngc.ne.nparcels .and. myid.eq.0 ) print *,'  [parcel] *** COUNT MISMATCH: lost/duplicated parcels! ***'
        ENDIF
      ENDIF
#endif

!----------------------------------------------------------------------
!  get height ASL:

      if( terrain_flag )then
            call getparcelzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,pdata)
            DO np=1,nparcels
              ! get z from sigma:
              ! (see Section 3 of "The governing equations for CM1", 
              !  http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf)
              pw(np,prz) = pw(np,przs) + pw(np,prsig)*((zt-pw(np,przs))*rzt)
            ENDDO
      endif

!----------------------------------------------------------------------

      end subroutine parcel_driver


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_interp(dt,mtime,xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                               zh,mh,rmh,zf,mf,znt,ust,c1,c2,          &
                               zs,sigma,sigmaf,rds,gz,                 &
                               pi0,th0,thv0,qv0,qc0,qi0,rth0,          &
                               dum1,dum2,dum3,dum4,zv  ,qt  ,prs,rho,  &
                               dum7,dum8,buoy,vpg  ,                   &
                               u3d,v3d,w3d,pp3d,th   ,t     ,th3d,q3d, &
                               kmh,kmv,khh,khv,tke3d,pt3d,pdata,       &
                               tdiag,qdiag,                            &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,reqs_p, &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)
      use input
      use constants
      use cm1libs , only : rslf,rsif
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine interpolates model information to the parcel locations
!  (diagnostic only ... not used for model integration)
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      double precision, intent(in) :: mtime
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: znt,ust
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(kb:ke) :: rds
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,thv0,qv0,qc0,qi0,rth0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke), target :: dum1,dum2,dum3,dum4,zv,qt,prs,rho
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum7,dum8
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: w3d,buoy,vpg
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pp3d,th3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke), target :: th,t
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq), target :: q3d
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tke3d
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt), target :: pt3d
      real, intent(inout), dimension(nparcels,npvals), target :: pdata
      real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
      real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag), target :: qdiag
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(cmp,jmp,kmt) :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt) :: tks1,tks2,tkn1,tkn2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tem,tem1
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d,z0,rznt,var

      ! distributed parcels (var17=1): operate on local arrays (pdl) via pointer.
      real, pointer :: pwk(:,:)
      integer :: npp

      ! scalar-point output registry (built once, consumed by the parcel loop)
      type(p3d), dimension(maxvars) :: scl
      integer,   dimension(maxvars) :: scl_slot
      integer :: nscl

      logical, parameter :: debug = .false.

!----------------------------------------------------------------------
!  Get derived variables:

    IF(imoist.eq.1)THEN
      ! with moisture:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
    do k=1,nk

      do j=1,nj
      do i=1,ni
        qt(i,j,k)=q3d(i,j,k,nqv)
      enddo
      enddo
      do n=nql1,nql2
        do j=1,nj
        do i=1,ni
          qt(i,j,k)=qt(i,j,k)+q3d(i,j,k,n)
        enddo
        enddo
      enddo
      IF(iice.eq.1)THEN
        do n=nqs1,nqs2
        do j=1,nj
        do i=1,ni
          qt(i,j,k)=qt(i,j,k)+q3d(i,j,k,n)
        enddo
        enddo
        enddo
      ENDIF
      IF( prth.ge.1 .or. prt.ge.1 .or. prqsl.ge.1 .or. prqsi.ge.1 .or.  prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          th(i,j,k) = (th0(i,j,k)+th3d(i,j,k))
          t(i,j,k) = th(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
      ENDIF
      IF( prb.ge.1 .or. prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k)             &
                           +repsm1*(q3d(i,j,k,nqv)-qv0(i,j,k))  &
                           -(qt(i,j,k)-q3d(i,j,k,nqv)-qc0(i,j,k)-qi0(i,j,k))   )
        enddo
        enddo
      ENDIF
      IF( prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum8(i,j,k) = th(i,j,k)*(1.0+reps*q3d(i,j,k,nqv))/(1.0+qt(i,j,k))
        enddo
        enddo
      ENDIF

    enddo

    ELSE
      ! dry:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    do k=1,nk

      IF( prth.ge.1 .or. prt.ge.1 .or. prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          th(i,j,k)= (th0(i,j,k)+th3d(i,j,k))
          t(i,j,k) = th(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
        enddo
        enddo
      ENDIF
      IF( prb.ge.1 .or. prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum7(i,j,k) = g*( th3d(i,j,k)*rth0(i,j,k) )
        enddo
        enddo
      ENDIF
      IF( prvpg.ge.1 )THEN
        do j=1,nj
        do i=1,ni
          dum8(i,j,k) = th(i,j,k)
        enddo
        enddo
      ENDIF

    enddo

    ENDIF


    IF( prb.ge.1 .or. prvpg.ge.1 )THEN
      do k=2,nk
      do j=1,nj
      do i=1,ni
        buoy(i,j,k) = (c1(1,1,k)*dum7(i,j,k-1)+c2(1,1,k)*dum7(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call    prepcornert(buoy,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,1)
      do j=0,nj+1
      do i=0,ni+1
        buoy(i,j,1) = buoy(i,j,2)+(buoy(i,j,3)-buoy(i,j,2))  &
                                 *(  zf(i,j,1)-  zf(i,j,2))  &
                                 /(  zf(i,j,3)-  zf(i,j,2))
        buoy(i,j,nk+1) = buoy(i,j,nk)+(buoy(i,j,nk  )-buoy(i,j,nk-1))  &
                                     *(  zf(i,j,nk+1)-  zf(i,j,nk  ))  &
                                     /(  zf(i,j,nk  )-  zf(i,j,nk-1))
      enddo
      enddo
    ENDIF
    IF( prvpg.ge.1 )THEN
    if( .not. terrain_flag )then
      do k=2,nk
      tem1 = rdz*cp*mf(1,1,k)
      do j=1,nj
      do i=1,ni
        vpg(i,j,k) = -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))  &
                          *(c2(1,1,k)*dum8(i,j,k)+c1(1,1,k)*dum8(i,j,k-1))
      enddo
      enddo
      enddo
    else
      do k=2,nk
      tem1 = rds(k)*cp
      do j=1,nj
      do i=1,ni
        vpg(i,j,k) = -tem1*(pp3d(i,j,k)-pp3d(i,j,k-1))*gz(i,j)  &
                          *(c2(1,1,k)*dum8(i,j,k)+c1(1,1,k)*dum8(i,j,k-1))
      enddo
      enddo
      enddo
    endif
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call    prepcornert(vpg ,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,  &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,1)
      ! cmr18:  at top/bottom boundaries, vpg + buoy = 0
      do j=0,nj+1
      do i=0,ni+1
        vpg(i,j,1) = -buoy(i,j,1)
        vpg(i,j,nk+1) = -buoy(i,j,nk+1)
      enddo
      enddo
    ENDIF

    if(timestats.ge.1) time_parcels=time_parcels+mytime()

!----------------------------------------------------------------------
!  get corner info for MPI runs
!  (may not parallelize correctly if this is not done)

#ifdef MPI
      call getcorneru(u3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcu2(u3d)
      call getcornerv(v3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcv2(v3d)
      call getcornerw(w3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcw2(w3d)
#endif

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v,s the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        u3d(i,j,0) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        v3d(i,j,0) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      ! no slip:
      if( imove.eq.1 )then
        IF(j.le.je)THEN
        do i=ib,ie+1
          u3d(i,j,0) = 0.0 - umove
        enddo
        ENDIF
        do i=ib,ie
          v3d(i,j,0) = 0.0 - vmove
        enddo
      else
        IF(j.le.je)THEN
        do i=ib,ie+1
          u3d(i,j,0) = 0.0
        enddo
        ENDIF
        do i=ib,ie
          v3d(i,j,0) = 0.0
        enddo
      endif
    ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
      do i=ib,ie+1
        u3d(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        v3d(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
      do i=ib,ie+1
        u3d(i,j,nk+1) = 0.0
      enddo
      ENDIF
      do i=ib,ie
        v3d(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      IF(j.le.je)THEN
      do i=ib,ie
        w3d(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      if( prth.ge.1 )then
        call prepcorners(th ,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prt.ge.1 )then
        call prepcorners(t  ,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prprs.ge.1 )then
        call prepcorners(prs,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if( prrho.ge.1 )then
        call prepcorners(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif
      if(prpt1.ge.1)then
        do n=1,npt
          call prepcorners(pt3d(ib,jb,kb,n),nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        enddo
      endif
      if( prqv.ge.1 )then
        call prepcorners(q3d(ib,jb,kb,nqv),nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
      endif
      if( prq1.ge.1 .or. prnc1.ge.1 )then
        do n = 1,numq
          call prepcorners(q3d(ib,jb,kb,n),nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,0)
        enddo
      endif
      if( prkm.ge.1 )then
        call prepcornert(kmh,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
        call prepcornert(kmv,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
      endif
      if( prkh.ge.1 )then
        call prepcornert(khh,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
        call prepcornert(khv,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
      endif
      if( prtke.ge.1 )then
        call prepcornert(tke3d,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,reqs_p,0)
      endif
      if( prdbz.ge.1 .and. qd_dbz.ge.1 )then
        call prepcorners(qdiag(ibdq,jbdq,kbdq,qd_dbz),  &
                             nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
      endif

!----------------------------------------------------------------------

    IF( prqsl.ge.1 )THEN
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum1(i,j,k) = rslf( prs(i,j,k) , t(i,j,k) )
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcorners(dum1,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    ENDIF
    IF( prqsi.ge.1 )THEN
      do k=1,nk
      do j=1,nj
      do i=1,ni
        dum2(i,j,k) = rsif( prs(i,j,k) , t(i,j,k) )
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_parcels=time_parcels+mytime()
      call prepcorners(dum2,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
    ENDIF

!----------------------------------------------------------------------
!  Get zvort at appropriate C-grid location:
!  (assuming no terrain)
!  cm1r18:  below lowest model level:
!           Use extrapolated velocities for bbc=1,2
!           Use log-layer equations for bbc=3 (see below)

    IF( przv.ge.1)THEN

      do k=0,nk+1
      do j=1,nj+1
      do i=1,ni+1
        zv(i,j,k) = (v3d(i,j,k)-v3d(i-1,j,k))*rdx*uf(i)   &
                   -(u3d(i,j,k)-u3d(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, get interpolated info:

    ! legacy: global pdata (all nparcels);  distributed: local pdl (np_loc).
    IF( nint(var17).eq.0 )THEN
      pwk => pdata
      npp =  nparcels
    ELSE
      pwk => pdl
      npp =  np_loc
    ENDIF

!----------------------------------------------------------------------
!  Build the scalar-point output registry (once, before the parcel loop).
!  Each line pairs a scalar-point source field with its output column
!  (pr* slot).  The interpolation happens once, in the loop below, so to
!  add a new scalar output you add ONE line here -- no grid-shape, weight,
!  or tri_interp bookkeeping.  Fields on other stagger points (u,v,w and
!  zvort) keep their own dedicated blocks below.

    nscl = 0

    IF( imoist.eq.1 )THEN
      if( prdbz.ge.1 .and. qd_dbz.ge.1 )then
        nscl=nscl+1 ; scl(nscl)%f => qdiag(:,:,:,qd_dbz) ; scl_slot(nscl) = prdbz
      endif
      if( prqv.ge.1 )then
        nscl=nscl+1 ; scl(nscl)%f => q3d(:,:,:,nqv) ; scl_slot(nscl) = prqv
      endif
      if( prq1.ge.1 )then
        do n=nql1,nql1+(prq2-prq1)
          nscl=nscl+1 ; scl(nscl)%f => q3d(:,:,:,n) ; scl_slot(nscl) = prq1+(n-nql1)
        enddo
      endif
      if( prnc1.ge.1 )then
        do n=nnc1,nnc1+(prnc2-prnc1)
          nscl=nscl+1 ; scl(nscl)%f => q3d(:,:,:,n) ; scl_slot(nscl) = prnc1+(n-nnc1)
        enddo
      endif
      if( prqsl.ge.1 )then
        nscl=nscl+1 ; scl(nscl)%f => dum1 ; scl_slot(nscl) = prqsl
      endif
      if( prqsi.ge.1 )then
        nscl=nscl+1 ; scl(nscl)%f => dum2 ; scl_slot(nscl) = prqsi
      endif
    ENDIF

    if( prth.ge.1 )then
      nscl=nscl+1 ; scl(nscl)%f => th  ; scl_slot(nscl) = prth
    endif
    if( prt.ge.1 )then
      nscl=nscl+1 ; scl(nscl)%f => t   ; scl_slot(nscl) = prt
    endif
    if( prprs.ge.1 )then
      nscl=nscl+1 ; scl(nscl)%f => prs ; scl_slot(nscl) = prprs
    endif
    if( prrho.ge.1 )then
      nscl=nscl+1 ; scl(nscl)%f => rho ; scl_slot(nscl) = prrho
    endif
    if( prpt1.ge.1 )then
      do n=1,npt
        nscl=nscl+1 ; scl(nscl)%f => pt3d(:,:,:,n) ; scl_slot(nscl) = prpt1+n-1
      enddo
    endif

    nploop2:  &
    DO np=1,npp

      pwk(np,prtime) = mtime

      x3d = pwk(np,prx)
      y3d = pwk(np,pry)
      z3d = pwk(np,prz)

      iflag = -100
      jflag = -100
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit2:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit2

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pwk(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pwk(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myprcl:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                   (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        kflag = 1
        if( .not. terrain_flag )then
          do while( pwk(np,prz).ge.zf(iflag,jflag,kflag+1) )
            kflag = kflag+1
          enddo
        else
          do while( pwk(np,prsig).ge.sigmaf(kflag+1) )
            kflag = kflag+1
          enddo
        endif

        x3d = pwk(np,prx)
        y3d = pwk(np,pry)
        z3d = pwk(np,prz)

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pwk(np,pry).lt.yh(j) )then
          j=j-1
        endif
        if( .not. terrain_flag )then
          if( pwk(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
          rz = ( pwk(np,prz)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
        else
          if( pwk(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
          rz = ( pwk(np,prsig)-sigma(k) )/( sigma(k+1)-sigma(k) )
        endif

        rx = ( pwk(np,prx)-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( pwk(np,pry)-yh(j) )/( yh(j+1)-yh(j) )

        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),pwk(np,prx),xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pwk(np,pry),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pwk(np,prz),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,u3d,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pwk(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( .not. terrain_flag )then
          if( pwk(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( pwk(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = ( pwk(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( pwk(np,pry)-yf(j) )/( yf(j+1)-yf(j) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pwk(np,prx),xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),pwk(np,pry),yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pwk(np,prz),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,v3d,vval)

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( pwk(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pwk(np,pry).lt.yh(j) )then
          j=j-1
        endif

!!!        rx = ( pwk(np,prx)-xh(i) )/( xh(i+1)-xh(i) )
!!!        ry = ( pwk(np,pry)-yh(j) )/( yh(j+1)-yh(j) )
        rx = rxs
        ry = rys
        if( .not. terrain_flag )then
          rz = ( pwk(np,prz)-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )
        else
          rz = ( pwk(np,prsig)-sigmaf(k) )/( sigmaf(k+1)-sigmaf(k) )
        endif

        ! saveit:
        rzw = rz

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pwk(np,prx),xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pwk(np,pry),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),pwk(np,prz),zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,w3d ,wval)
      if(prkm.ge.1)then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kmh,pwk(np,prkm  ))
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kmv,pwk(np,prkm+1))
      endif
      if(prkh.ge.1)then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,khh,pwk(np,prkh  ))
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,khv,pwk(np,prkh+1))
      endif
      if( prtke.ge.1 )then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,tke3d,pwk(np,prtke))
      endif
      if( prb.ge.1 )then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,buoy,pwk(np,prb))
      endif
      if( prvpg.ge.1 )then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,vpg,pwk(np,prvpg))
      endif

!----------------------------------------------------------------------
!  Data on scalar points

        i=iflag
        j=jflag
        k=kflag

        if( pwk(np,prx).lt.xh(i) )then
          i=i-1
        endif
        if( pwk(np,pry).lt.yh(j) )then
          j=j-1
        endif
        if( .not. terrain_flag )then
          if( pwk(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( pwk(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = rxs
        ry = rys
        rz = rzs

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15558: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        ! all scalar-point outputs, driven by the registry built above:
        do n=1,nscl
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,scl(n)%f,pwk(np,scl_slot(n)))
        enddo

!----------------------------------------------------------------------
!  Data on zvort points

      IF( przv.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( .not. terrain_flag )then
          if( pwk(np,prz).lt.zh(iflag,jflag,k) )then
            k=k-1
          endif
        else
          if( pwk(np,prsig).lt.sigma(k) )then
            k=k-1
          endif
        endif

        rx = rxu
        ry = ryv
        rz = rzs

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zv,pwk(np,przv))

      ENDIF

!----------------------------------------------------------------------
!  surface variables  and  uv for parcels below lowest model level:

      IF( prznt.ge.1 .or. prust.ge.1 .or. bbc.eq.3 )THEN
        i=iflag
        j=jflag
        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt,z0)
        if( prznt.ge.1 ) pwk(np,prznt) = z0
        if( prust.ge.1 )  &
        call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,ust,pwk(np,prust))
      ENDIF

      IF( bbc.eq.3 )THEN
        ! semi-slip lower boundary condition:
        if( z3d.lt.zh(1,1,1) )then
          ! re-calculate velocities if parcel is below lowest model level:
          !------
          ! u at lowest model level:
          i=iflag
          j=jflag
          if( y3d.lt.yh(j) )then
            j=j-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,u3d(ib,jb,1),uval)
          !------
          ! v at lowest model level:
          i=iflag
          j=jflag
          if( x3d.lt.xh(i) )then
            i=i-1
          endif
          call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,v3d(ib,jb,1),vval)
          !------
          ! get u,v from (neutral) log-layer equation:
          rznt = 1.0/z0
          var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
          uval = (uval+umove)*var
          vval = (vval+vmove)*var
          !------
          IF( przv.ge.1 )THEN
            do j=jflag-1,jflag+1
            do i=iflag  ,iflag+1
              z0 = 0.5*(znt(i-1,j)+znt(i,j))
              rznt = 1.0/z0
              dum3(i,j,1) = (u3d(i,j,1)+umove)*alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
            enddo
            enddo
            do j=jflag  ,jflag+1
            do i=iflag-1,iflag+1
              z0 = 0.5*(znt(i,j-1)+znt(i,j))
              rznt = 1.0/z0
              dum4(i,j,1) = (v3d(i,j,1)+vmove)*alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
            enddo
            enddo
            do j=jflag,jflag+1
            do i=iflag,iflag+1
              dum7(i,j,1) = (dum4(i,j,1)-dum4(i-1,j,1))*rdx*uf(i)   &
                           -(dum3(i,j,1)-dum3(i,j-1,1))*rdy*vf(j)
            enddo
            enddo
            i=iflag
            j=jflag
            call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 1, 0, 0,dum7(ib,jb,1),pwk(np,przv))
          ENDIF
        endif
      ENDIF


!----------------------------------------------------------------------

        pwk(np,pru)=uval
        pwk(np,prv)=vval
        pwk(np,prw)=wval

#ifdef MPI
      ELSE

        ! set to really small number (so we can use the allreduce command below)
        do n=4,npvals
          pwk(np,n) = -1.0e30
        enddo
#endif

      ENDIF  myprcl

    ENDDO  nploop2

!----------------------------------------------------------------------
!  communicate data
#ifdef MPI
  ! legacy only: gather each rank's owned-parcel diagnostics onto rank 0
  ! (via MAX over the -1e30 sentinels) so rank 0 can write them.  In
  ! distributed mode every rank already holds its own parcels' diagnostics
  ! in pdl; the gather for output is handled by the Phase-5 write path.
  IF( nint(var17).eq.0 )THEN
    do np=4,npvals
      if( myid.eq.0 )then
        call MPI_REDUCE(MPI_IN_PLACE,pdata(1,np),nparcels,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_REDUCE( pdata(1,np),pdata(1,np),nparcels,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      endif
    enddo
  ENDIF
#endif
!----------------------------------------------------------------------

      end subroutine parcel_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_write(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc)
      use input
#ifdef NETCDF
      use writeout_nc_module, only : writepdata_nc
#endif
      implicit none

      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,3) :: ploc

      integer :: i,n,np

!----------------------------------------------------------------------
!  write out data

    IF(myid.eq.0)THEN

      IF(output_format.eq.1)THEN
        ! GrADS format:

        ! write GrADS descriptor file:
        call write_prclctl(name_prcl,desc_prcl,unit_prcl,prec)

        do i=1,maxstring
          string(i:i) = ' '
        enddo

        string = 'cm1out_pdata.dat'

!!!        write(string(15:20),101) prec
!!!101     format(i6.6)

        if(dowr) write(outfile,*) string
        open(unit=61,file=string,form='unformatted',access='direct',   &
             recl=4*npvals*nparcels,status='unknown')

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  pdata prec = ',prec

        write(61,rec=prec) ((pdata(np,n),np=1,nparcels),n=1,npvals)

        close(unit=61)

#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN

        call     writepdata_nc(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,pdata,ploc(1,1))

#endif
      ENDIF
      if(dowr) write(outfile,*)

    ENDIF   ! endif for myid=0

      return
      end subroutine parcel_write


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tri_interp(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,pdata)
      use input
      implicit none

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      real :: pdata

      pdata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

      end subroutine tri_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine get2d(i,j,x3d,y3d,xh,xf,yh,yf,xs,ys,is,js,s,sval)
    use input
    implicit none

    integer, intent(in) :: i,j
    real, intent(in) :: x3d,y3d
    real, intent(in), dimension(ib:ie) :: xh
    real, intent(in), dimension(ib:ie+1) :: xf
    real, intent(in), dimension(jb:je) :: yh
    real, intent(in), dimension(jb:je+1) :: yf

    ! 0 = scalar point
    ! 1 = velocity point
    integer, intent(in) :: xs,ys
    integer, intent(in) :: is,js

    real, intent(in), dimension(ib:ie+is,jb:je+js) :: s
    real, intent(out) :: sval

    real :: wg1,wg2,wg3,wg4
    real :: x13,x23,x33,x43
    real :: w1,w2,w3,w7,rx,ry,rz

    logical, parameter :: debug = .false.

!-----------------------------------------------------------------------
      ! tri-linear interp:

      IF(xs.eq.1)THEN
        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
      ELSE
        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
      ENDIF

      IF(ys.eq.1)THEN
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
      ELSE
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
      ENDIF

        w1=(1.0-rx)*(1.0-ry)
        w2=rx*(1.0-ry)
        w3=(1.0-rx)*ry
        w7=rx*ry

      IF( debug )THEN
        if( rx.lt.-0.000001 .or. rx.gt.1.000001 .or.        &
            ry.lt.-0.000001 .or. ry.gt.1.000001 .or.        &
            (w1+w2+w3+w7).lt.0.999999 .or.  &
            (w1+w2+w3+w7).gt.1.000001       &
          )then
          print *,'  x3d,y3d     = ',x3d,y3d
          print *,'  i,j         = ',i,j
          print *,'  rx,ry       = ',rx,ry
          print *,'  w1,w2,w3,w7 = ',w1,w2,w3,w7
          print *,'  w1+w2+w3+w7 = ',w1+w2+w3+w7
          print *,' 22346 '
          call stopcm1
        endif
      ENDIF

      sval =s(i  ,j  )*w1    &
           +s(i+1,j  )*w2    &
           +s(i  ,j+1)*w3    &
           +s(i+1,j+1)*w7

!-----------------------------------------------------------------------

    end subroutine get2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getparcelzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,pdata)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(inout), dimension(nparcels,npvals) :: pdata

      integer :: i,j,iflag,jflag,np
      real :: x3d,y3d

    zsnploop:  &
    DO np=1,nparcels

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)

      iflag = -100
      jflag = -100

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  zshaveit1:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  zshaveit1

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      zsmyparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                       (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,zs,pdata(np,przs))

      ELSE  zsmyparcel

        pdata(np,przs) = -1.0e30

      ENDIF  zsmyparcel

    ENDDO  zsnploop

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,pdata(1,przs) ,nparcels,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

      end subroutine getparcelzs


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setup_parcel_vars(name_prcl,desc_prcl,unit_prcl,qname)
      use input
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      character(len=3), intent(in), dimension(maxq) :: qname

      integer :: n,n2
      character(len=8) :: text1
      character(len=30) :: text2

      prcl_out = 0

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'x'
      desc_prcl(prcl_out) = 'x position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'y'
      desc_prcl(prcl_out) = 'y position'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'z'
      desc_prcl(prcl_out) = 'z position (above sea level)'
      unit_prcl(prcl_out) = 'm'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'u'
      desc_prcl(prcl_out) = 'u velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'v'
      desc_prcl(prcl_out) = 'v velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'w'
      desc_prcl(prcl_out) = 'w velocity'
      unit_prcl(prcl_out) = 'm/s'

      prcl_out = prcl_out+1
      name_prcl(prcl_out) = 'mtime'
      desc_prcl(prcl_out) = 'model time (seconds since beginning of simulation)'
      unit_prcl(prcl_out) = 's'

      if( prth.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'th'
        desc_prcl(prcl_out) = 'potential temperature'
        unit_prcl(prcl_out) = 'K'
      endif

      if( prt.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 't'
        desc_prcl(prcl_out) = 'temperature'
        unit_prcl(prcl_out) = 'K'
      endif

      if( prprs.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'prs'
        desc_prcl(prcl_out) = 'pressure'
        unit_prcl(prcl_out) = 'Pa'
      endif

      if(prpt1.ge.1)then
        do n=1,npt
          text1='pt      '
          if(n.le.9)then
            write(text1(3:3),155) n
155         format(i1.1)
          elseif(n.le.99)then
            write(text1(3:4),154) n
154         format(i2.2)
          else
            write(text1(3:5),153) n
153         format(i3.3)
          endif

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = 'passive tracer mixing ratio'
          unit_prcl(prcl_out) = 'kg/kg'
        enddo
      endif

      if( prqv.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qv'
        desc_prcl(prcl_out) = 'water vapor mixing ratio'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if(prq1.ge.1)then
        n2 = nql2
        if( iice.eq.1 ) n2 = nqs2
        do n=nql1,n2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = text2
          unit_prcl(prcl_out) = 'kg/kg'
        enddo
      endif

      if(prnc1.ge.1)then
        do n=nnc1,nnc2
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)

          prcl_out = prcl_out+1
          name_prcl(prcl_out) = text1
          desc_prcl(prcl_out) = text2
          unit_prcl(prcl_out) = '1/kg'
        enddo
      endif

      if( prkm.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'kmh'
        desc_prcl(prcl_out) = 'horiz eddy viscosity for momentum'
        unit_prcl(prcl_out) = 'm^2/s'

        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'kmv'
        desc_prcl(prcl_out) = 'vert eddy viscosity for momentum'
        unit_prcl(prcl_out) = 'm^2/s'
      endif

      if( prkh.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'khh'
        desc_prcl(prcl_out) = 'horiz eddy diffusivity for scalars'
        unit_prcl(prcl_out) = 'm^2/s'

        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'khv'
        desc_prcl(prcl_out) = 'vert eddy diffusivity for scalars'
        unit_prcl(prcl_out) = 'm^2/s'
      endif

      if( prtke.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'tke'
        desc_prcl(prcl_out) = 'subgrid tke'
        unit_prcl(prcl_out) = 'm^2/s^2'
      endif

      if( prdbz.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'dbz'
        desc_prcl(prcl_out) = 'reflectivity'
        unit_prcl(prcl_out) = 'dBZ'
      endif

      if( prb.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'b'
        desc_prcl(prcl_out) = 'buoyancy'
        unit_prcl(prcl_out) = 'm/s/s'
      endif

      if( prvpg.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'vpg'
        desc_prcl(prcl_out) = 'vertical perturbation pressure gradient'
        unit_prcl(prcl_out) = 'm/s/s'
      endif

      if( przv.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'zvort'
        desc_prcl(prcl_out) = 'vertical vorticity'
        unit_prcl(prcl_out) = '1/s'
      endif

      if( prrho.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'rho'
        desc_prcl(prcl_out) = 'dry-air density'
        unit_prcl(prcl_out) = 'kg/m^3'
      endif

      if( prqsl.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qsl'
        desc_prcl(prcl_out) = 'saturation mixing ratio wrt liquid'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prqsi.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'qsi'
        desc_prcl(prcl_out) = 'saturation mixing ratio wrt ice'
        unit_prcl(prcl_out) = 'kg/kg'
      endif

      if( prznt.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'znt'
        desc_prcl(prcl_out) = 'surface roughness length'
        unit_prcl(prcl_out) = 'm'
      endif

      if( prust.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'ust'
        desc_prcl(prcl_out) = 'surface friction velocity'
        unit_prcl(prcl_out) = 'm/s'
      endif

      if( przs.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'zs'
        desc_prcl(prcl_out) = 'terrain height'
        unit_prcl(prcl_out) = 'm'
      endif

      if( prsig.ge.1 )then
        prcl_out = prcl_out+1
        name_prcl(prcl_out) = 'sigma'
        desc_prcl(prcl_out) = 'sigma (nondimensional height)'
        unit_prcl(prcl_out) = 'nondimensional'
      endif

!-----------------------------------------------------------------------

      end subroutine setup_parcel_vars


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_prclctl(name_prcl,desc_prcl,unit_prcl,prec)
      use input
      use constants , only : grads_undef
      implicit none

      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl
      integer, intent(in) :: prec

      integer :: i,n,nn
      character(len=16) :: a16

      !---------------------------------------------------------------
      ! This subroutine writes the GrADS descriptor file for parcels
      !---------------------------------------------------------------

    idcheck:  &
    IF( myid.eq.0 )THEN

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string = 'cm1out_parcel.ctl'
      if(dowr) write(outfile,*) string
      open(unit=50,file=string,status='unknown')

      string = 'cm1out_pdata.dat'

      write(50,401) string
      write(50,402) trim(cm1version)
      write(50,403) grads_undef
      write(50,404) nparcels
      write(50,405)
      write(50,406)
      write(50,407) prec

      write(50,408) prcl_out

      DO n = 1 , prcl_out
        a16 = '                '
        nn = len(trim(unit_prcl(n)))
        write(a16(2:15),214) unit_prcl(n)
        write(a16(1:1),201 )       '('
        write(a16(nn+2:nn+2),201 ) ')'
        write(50,409) name_prcl(n),desc_prcl(n),a16
      ENDDO

      write(50,410)

      close(unit=50)

    ENDIF  idcheck

201   format(a1)
214   format(a14)

401   format('dset ^',a)
402   format('title CM1 parcel data output, using version ',a,'; time is generic, see variable mtime for actual times')
403   format('undef ',f10.1)
404   format('xdef ',i10,' linear 1 1')
405   format('ydef          1 linear 0 1')
406   format('zdef          1 linear 0 1')
407   format('tdef ',i10,' linear 00:00Z01JAN0001 1YR')
408   format('vars ',i6)
409   format(a12,' 1 99 ',a40,1x,a16)
410   format('endvars')

      end subroutine write_prclctl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_massprofile(zf,rho0,qv0,cummass,mcol)
      use input
      implicit none

!-----------------------------------------------------------------------
!  Builds the cumulative DRY-air mass profile (per unit horizontal area)
!  at the vertical faces zf, for mass-weighted parcel initialization.
!
!    cummass(k) = integral from the surface up to face zf(k) of
!                 rhod0 dz', where rhod0 = rho0/(1+qv0) is the dry-air
!                 (partial) density.  cummass(1)=0, cummass(nk+1)=mcol.
!    mcol       = total column dry-air mass per unit area (kg/m^2).
!
!  Assumes a horizontally homogeneous base state (flat terrain), so the
!  i=1,j=1 column represents the whole domain.
!-----------------------------------------------------------------------

      real, intent(in),  dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in),  dimension(ib:ie,jb:je,kb:ke)   :: rho0,qv0
      real, intent(out), dimension(1:nk+1)              :: cummass
      real, intent(out)                                 :: mcol

      integer :: k
      real :: rhod0

      cummass(1) = 0.0
      do k=1,nk
        rhod0 = rho0(1,1,k)/(1.0+qv0(1,1,k))
        cummass(k+1) = cummass(k) + rhod0*( zf(1,1,k+1)-zf(1,1,k) )
      enddo
      mcol = cummass(nk+1)

      end subroutine parcel_massprofile


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_resolve_box(zf,cummass,                         &
                                    bx0,bx1,by0,by1,bz0,bz1,mlo,mhi)
      use input
      implicit none

!-----------------------------------------------------------------------
!  Resolve the optional parcel seeding bounding box (var9..var14) into
!  concrete, in-domain bounds plus the cumulative dry-air mass at the z
!  bounds.  SINGLE source of truth, called by BOTH the count derivation
!  (cm1.F, var18=2) and the placement (init3d.F) so the two cannot drift.
!
!    var9 ,var10 = x min,max (m)
!    var11,var12 = y min,max (m)
!    var13,var14 = z min,max (m ASL)
!
!  Per axis, a non-positive-width request (max<=min, e.g. both left at
!  the default 0) falls back to the full domain extent for that axis, so
!  the box is opt-in and per-axis.  Bounds are then CLAMPED to the domain
!  ([minx,maxx], [miny,maxy], and the model column for z); clamping is
!  reported.  A box that lands entirely outside the domain (zero or
!  negative width after clamping) is a fatal error -- otherwise it would
!  silently inflate nparcels (var18=2) and/or place parcels outside the
!  domain (dropped in distributed mode -> conservation mismatch).
!-----------------------------------------------------------------------

      real, intent(in),  dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in),  dimension(1:nk+1)              :: cummass
      real, intent(out)  :: bx0,bx1,by0,by1,bz0,bz1,mlo,mhi

      integer :: k
      real :: frac
      logical :: clamped

      ! per-axis defaults: a non-positive-width request -> full domain.
      bx0 = var9  ; bx1 = var10
      by0 = var11 ; by1 = var12
      bz0 = var13 ; bz1 = var14
      if( bx1.le.bx0 )then ; bx0 = minx       ; bx1 = maxx         ; endif
      if( by1.le.by0 )then ; by0 = miny       ; by1 = maxy         ; endif
      if( bz1.le.bz0 )then ; bz0 = zf(1,1,1)  ; bz1 = zf(1,1,nk+1) ; endif

      ! clamp to the domain; report if a user-supplied bound was moved.
      clamped = .false.
      if( bx0.lt.minx       )then ; bx0 = minx       ; clamped = .true. ; endif
      if( bx1.gt.maxx       )then ; bx1 = maxx       ; clamped = .true. ; endif
      if( by0.lt.miny       )then ; by0 = miny       ; clamped = .true. ; endif
      if( by1.gt.maxy       )then ; by1 = maxy       ; clamped = .true. ; endif
      if( bz0.lt.zf(1,1,1)    )then ; bz0 = zf(1,1,1)    ; clamped = .true. ; endif
      if( bz1.gt.zf(1,1,nk+1) )then ; bz1 = zf(1,1,nk+1) ; clamped = .true. ; endif
      if( clamped .and. dowr )then
        write(outfile,*) '  WARNING: parcel seeding box clamped to the domain:'
        write(outfile,*) '           x =',bx0,bx1
        write(outfile,*) '           y =',by0,by1
        write(outfile,*) '           z =',bz0,bz1
      endif

      ! a degenerate / entirely-outside box is fatal.
      if( bx1.le.bx0 .or. by1.le.by0 .or. bz1.le.bz0 )then
        if(dowr)then
          write(outfile,*)
          write(outfile,*) '  ERROR: parcel seeding box (var9..var14) is degenerate or'
          write(outfile,*) '         lies entirely outside the domain (zero width after'
          write(outfile,*) '         clamping).  resolved bounds:'
          write(outfile,*) '         x =',bx0,bx1
          write(outfile,*) '         y =',by0,by1
          write(outfile,*) '         z =',bz0,bz1
          write(outfile,*)
        endif
        call stopcm1
      endif

      ! cumulative dry-air mass at each z bound (same face-linear interp as
      ! the placement inversion), so a vertical draw confined to [bz0,bz1]
      ! stays mass-weighted: drawing uniformly in [mlo,mhi] reproduces the
      ! column mass pdf restricted to the box.
      mlo = cummass(1)
      mhi = cummass(1)
      do k=1,nk
        if( bz0.ge.zf(1,1,k) .and. bz0.le.zf(1,1,k+1) )then
          frac = ( bz0 - zf(1,1,k) ) / max( zf(1,1,k+1)-zf(1,1,k) , 1.0e-30 )
          mlo  = cummass(k) + frac*( cummass(k+1)-cummass(k) )
          exit
        endif
      enddo
      do k=1,nk
        if( bz1.ge.zf(1,1,k) .and. bz1.le.zf(1,1,k+1) )then
          frac = ( bz1 - zf(1,1,k) ) / max( zf(1,1,k+1)-zf(1,1,k) , 1.0e-30 )
          mhi  = cummass(k) + frac*( cummass(k+1)-cummass(k) )
          exit
        endif
      enddo

      end subroutine parcel_resolve_box


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_local_cell(x3d,y3d,xf,yf,iflag,jflag)
      use input
      implicit none

!-----------------------------------------------------------------------
!  Ownership test (distributed parcels):  returns the LOCAL cell indices
!  (iflag,jflag) of the cell containing parcel (x3d,y3d).  This rank owns
!  the parcel iff, on return, 1<=iflag<=ni .and. 1<=jflag<=nj.
!  Otherwise iflag/jflag are <1 (not on this subdomain).
!
!  This reproduces verbatim the ownership + boundary tie-break logic used
!  in parcel_driver (so distributed and legacy agree on who owns what).
!-----------------------------------------------------------------------

      real, intent(in) :: x3d,y3d
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf
      integer, intent(out) :: iflag,jflag

      integer :: i,j

      iflag = -100
      jflag = -100

      IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.   &
          y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

        IF(nx.eq.1)THEN
          iflag = 1
        ELSE
          i = ni+1
          do while( iflag.lt.0 .and. i.gt.1 )
            i = i-1
            if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) ) iflag = i
          enddo
        ENDIF

        IF(axisymm.eq.1.or.ny.eq.1)THEN
          jflag = 1
        ELSE
          j = nj+1
          do while( jflag.lt.0 .and. j.gt.1 )
            j = j-1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) ) jflag = j
          enddo
        ENDIF

      ENDIF

#ifdef MPI
      ! boundary tie-break: a parcel exactly on the shared east/north face
      ! belongs to the neighbor, not to us (unless we are the last rank):
      IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
          (jflag.ge.1.and.jflag.le.nj) )THEN
        IF( iflag.eq.ni .and. x3d.eq.xf(iflag+1) .and. nodex.gt.1 .and. myi.ne.nodex ) iflag = -1
        IF( jflag.eq.nj .and. y3d.eq.yf(jflag+1) .and. nodey.gt.1 .and. myj.ne.nodey ) jflag = -1
      ENDIF
#endif

      end subroutine parcel_local_cell


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_alloc_local(npv)
      use input
      implicit none
!-----------------------------------------------------------------------
!  Allocate this rank's local parcel arrays (distributed mode).  Capacity
!  starts at headroom*(global/numprocs); it grows on demand if a rank
!  accumulates more than expected (clustering).  Resets np_loc=0.
!-----------------------------------------------------------------------
      integer, intent(in) :: npv
      integer :: navg

      navg = (nparcels + numprocs - 1)/numprocs    ! ceil(nparcels/numprocs)
      np_max = max( 1000, 4*navg )
      if( allocated(pdl)  ) deallocate(pdl)
      if( allocated(pidl) ) deallocate(pidl)
      allocate( pdl(np_max,npv) )
      allocate( pidl(np_max) )
      pdl  = 0.0
      pidl = 0
      np_loc = 0

      end subroutine parcel_alloc_local


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_grow_local(need)
      use input
      implicit none
!-----------------------------------------------------------------------
!  Ensure local capacity for at least 'need' parcels (grow ~1.5x).
!-----------------------------------------------------------------------
      integer, intent(in) :: need
      integer :: newmax
      real,    allocatable :: t2(:,:)
      integer, allocatable :: t1(:)

      if( need.le.np_max ) return
      newmax = max( need, (3*np_max)/2 )
      allocate( t2(np_max,npvals) ) ;  t2(:,:) = pdl(:,:)
      allocate( t1(np_max)        ) ;  t1(:)   = pidl(:)
      deallocate(pdl)  ; deallocate(pidl)
      allocate( pdl(newmax,npvals) ) ;  pdl = 0.0 ;  pdl(1:np_max,:) = t2(:,:)
      allocate( pidl(newmax)       ) ;  pidl = 0  ;  pidl(1:np_max)  = t1(:)
      deallocate(t2) ;  deallocate(t1)
      np_max = newmax

      end subroutine parcel_grow_local


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_append_local(nid,x3d,y3d,z3d,xf,yf,stored)
      use input
      implicit none
!-----------------------------------------------------------------------
!  If parcel (x3d,y3d,z3d) with global id nid is owned by this rank,
!  append it to the local arrays (positions set; other columns left 0).
!  Returns stored=.true. if it was kept here.
!-----------------------------------------------------------------------
      integer, intent(in) :: nid
      real, intent(in) :: x3d,y3d,z3d
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf
      logical, intent(out) :: stored

      integer :: iflag,jflag

      stored = .false.
      call parcel_local_cell(x3d,y3d,xf,yf,iflag,jflag)
      if( iflag.ge.1 .and. iflag.le.ni .and. jflag.ge.1 .and. jflag.le.nj )then
        if( np_loc.ge.np_max ) call parcel_grow_local(np_loc+1)
        np_loc = np_loc + 1
        pdl(np_loc,prx) = x3d
        pdl(np_loc,pry) = y3d
        pdl(np_loc,prz) = z3d
        pidl(np_loc)    = nid
        stored = .true.
      endif

      end subroutine parcel_append_local


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_emit(nid,x3d,y3d,z3d,pdata,xf,yf)
      use input
      implicit none
!-----------------------------------------------------------------------
!  Placement-agnostic parcel sink.  Emit the parcel with global id 'nid'
!  at position (x3d,y3d,z3d).  Every placement scheme (grid, mass-weighted,
!  file, ...) routes new parcels through here, so all of them behave
!  identically in legacy and distributed modes:
!    var17=0 (legacy)      -> written to the global pdata(nid,*).
!    var17=1 (distributed) -> filtered into this rank's local arrays
!                             (kept only if (x3d,y3d) is owned here).
!  In distributed mode pdata may be a size-1 placeholder; the legacy
!  branch that indexes it is never taken, so that is safe.
!-----------------------------------------------------------------------
      integer, intent(in) :: nid
      real, intent(in) :: x3d,y3d,z3d
      real, intent(inout), dimension(nparcels,npvals) :: pdata
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf

      logical :: stored

      IF( nint(var17).eq.0 )THEN
        pdata(nid,prx) = x3d
        pdata(nid,pry) = y3d
        pdata(nid,prz) = z3d
      ELSE
        call parcel_append_local(nid,x3d,y3d,z3d,xf,yf,stored)
      ENDIF

      end subroutine parcel_emit


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_wrap_inplace(coord,pwmin,pwmax)
      use input
      implicit none
!-----------------------------------------------------------------------
!  Periodic wrap of the 'coord' column of all local parcels into
!  [pwmin,pwmax).  Used when there is a single rank in this direction
!  (no migration needed, but parcels still wrap periodically).
!-----------------------------------------------------------------------
      integer, intent(in) :: coord
      real,    intent(in) :: pwmin,pwmax
      integer :: n
      real :: span
      span = pwmax-pwmin
      do n=1,np_loc
        if( pdl(n,coord).ge.pwmax ) pdl(n,coord) = pdl(n,coord)-span
        if( pdl(n,coord).lt.pwmin ) pdl(n,coord) = pdl(n,coord)+span
      enddo
      end subroutine parcel_wrap_inplace


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_migrate_axis(coord,mneg,mpos,pwmin,pwmax,lobnd,hibnd,islast)
#ifdef MPI
      use input
      use mpi
      implicit none
!-----------------------------------------------------------------------
!  One migration pass along a single direction (x or y).
!    coord        = column of pdl holding this direction's position (prx/pry)
!    mneg,mpos    = MPI ranks of the negative/positive neighbors
!    pwmin,pwmax    = global domain bounds in this direction (periodic span)
!    lobnd,hibnd  = this rank's local subdomain face bounds (xf(1),xf(ni+1))
!    islast       = .true. if this is the last rank in this direction
!  Parcels with position < lobnd go to mneg; > hibnd (or == hibnd unless
!  this is the last rank) go to mpos.  Departing parcels are packed and
!  exchanged; received parcels are appended and periodically wrapped into
!  [pwmin,pwmax).  Positions are un-wrapped on entry (wrap deferred from the
!  RK step), so an edge parcel exits to the periodic neighbor correctly.
!-----------------------------------------------------------------------
      integer, intent(in) :: coord,mneg,mpos
      real,    intent(in) :: pwmin,pwmax,lobnd,hibnd
      logical, intent(in) :: islast

      real,    allocatable :: snd(:,:),spd(:,:),rnd(:,:),rpd(:,:)
      integer, allocatable :: sni(:),spi(:),rni(:),rpi(:)
      integer :: n,nsn,nsp,nrn,nrp,nkeep,ireq,need
      integer :: req(8),sta(MPI_STATUS_SIZE,8),st1(MPI_STATUS_SIZE),ie2
      real :: c,span
      logical :: goneg,gopos

      span = pwmax-pwmin

      ! --- count departing parcels ---
      nsn=0 ; nsp=0
      do n=1,np_loc
        c = pdl(n,coord)
        goneg = ( c.lt.lobnd )
        gopos = ( c.gt.hibnd ) .or. ( c.eq.hibnd .and. (.not.islast) )
        if( goneg )then
          nsn=nsn+1
        elseif( gopos )then
          nsp=nsp+1
        endif
      enddo

      allocate( snd(npvals,max(1,nsn)) , sni(max(1,nsn)) )
      allocate( spd(npvals,max(1,nsp)) , spi(max(1,nsp)) )

      ! --- pack departing; compact keepers in place ---
      nsn=0 ; nsp=0 ; nkeep=0
      do n=1,np_loc
        c = pdl(n,coord)
        goneg = ( c.lt.lobnd )
        gopos = ( c.gt.hibnd ) .or. ( c.eq.hibnd .and. (.not.islast) )
        if( goneg )then
          nsn=nsn+1 ; snd(:,nsn)=pdl(n,:) ; sni(nsn)=pidl(n)
        elseif( gopos )then
          nsp=nsp+1 ; spd(:,nsp)=pdl(n,:) ; spi(nsp)=pidl(n)
        else
          nkeep=nkeep+1
          if( nkeep.ne.n )then
            pdl(nkeep,:)=pdl(n,:) ; pidl(nkeep)=pidl(n)
          endif
        endif
      enddo
      np_loc = nkeep

      ! --- exchange counts ---
      ! my pos-going -> mpos ; I receive mneg's pos-going (it entered my neg side)
      call MPI_SENDRECV(nsp,1,MPI_INTEGER,mpos,101, nrn,1,MPI_INTEGER,mneg,101, MPI_COMM_WORLD,st1,ie2)
      ! my neg-going -> mneg ; I receive mpos's neg-going (it entered my pos side)
      call MPI_SENDRECV(nsn,1,MPI_INTEGER,mneg,102, nrp,1,MPI_INTEGER,mpos,102, MPI_COMM_WORLD,st1,ie2)

      allocate( rnd(npvals,max(1,nrn)) , rni(max(1,nrn)) )
      allocate( rpd(npvals,max(1,nrp)) , rpi(max(1,nrp)) )

      ! --- exchange payloads (data + ids), non-blocking ---
      ireq=0
      ireq=ireq+1 ; call MPI_IRECV(rnd,npvals*max(1,nrn),MPI_REAL,   mneg,201,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_IRECV(rni,       max(1,nrn),MPI_INTEGER,mneg,202,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_IRECV(rpd,npvals*max(1,nrp),MPI_REAL,   mpos,203,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_IRECV(rpi,       max(1,nrp),MPI_INTEGER,mpos,204,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_ISEND(spd,npvals*max(1,nsp),MPI_REAL,   mpos,201,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_ISEND(spi,       max(1,nsp),MPI_INTEGER,mpos,202,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_ISEND(snd,npvals*max(1,nsn),MPI_REAL,   mneg,203,MPI_COMM_WORLD,req(ireq),ie2)
      ireq=ireq+1 ; call MPI_ISEND(sni,       max(1,nsn),MPI_INTEGER,mneg,204,MPI_COMM_WORLD,req(ireq),ie2)
      call MPI_WAITALL(ireq,req,sta,ie2)

      ! --- append received, with periodic wrap into [pwmin,pwmax) ---
      need = np_loc + nrn + nrp
      if( need.gt.np_max ) call parcel_grow_local(need)
      do n=1,nrn
        np_loc = np_loc+1
        pdl(np_loc,:) = rnd(:,n) ; pidl(np_loc) = rni(n)
        if( pdl(np_loc,coord).ge.pwmax ) pdl(np_loc,coord)=pdl(np_loc,coord)-span
        if( pdl(np_loc,coord).lt.pwmin ) pdl(np_loc,coord)=pdl(np_loc,coord)+span
      enddo
      do n=1,nrp
        np_loc = np_loc+1
        pdl(np_loc,:) = rpd(:,n) ; pidl(np_loc) = rpi(n)
        if( pdl(np_loc,coord).ge.pwmax ) pdl(np_loc,coord)=pdl(np_loc,coord)-span
        if( pdl(np_loc,coord).lt.pwmin ) pdl(np_loc,coord)=pdl(np_loc,coord)+span
      enddo

      deallocate( snd,sni,spd,spi,rnd,rni,rpd,rpi )
#else
      integer, intent(in) :: coord,mneg,mpos
      real,    intent(in) :: pwmin,pwmax,lobnd,hibnd
      logical, intent(in) :: islast
#endif
      end subroutine parcel_migrate_axis


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_migrate(xf,yf)
      use input
      implicit none
!-----------------------------------------------------------------------
!  Distributed-parcel migration: hand parcels that crossed a subdomain
!  boundary to the owning neighbor rank.  Two passes (x then y) which,
!  together, route corner/diagonal exits correctly.  Replaces the stock
!  per-step global MPI_ALLREDUCE.  Called once per parcel update step.
!-----------------------------------------------------------------------
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je+1) :: yf

      ! x-direction:
      IF( nx.gt.1 )THEN
        IF( nodex.gt.1 )THEN
          call parcel_migrate_axis(prx,mywest,myeast,minx,maxx,xf(1),xf(ni+1),myi.eq.nodex)
        ELSE
          call parcel_wrap_inplace(prx,minx,maxx)
        ENDIF
      ENDIF

      ! y-direction (skip for axisymmetric / 2D x-z):
      IF( ny.gt.1 .and. axisymm.ne.1 )THEN
        IF( nodey.gt.1 )THEN
          call parcel_migrate_axis(pry,mysouth,mynorth,miny,maxy,yf(1),yf(nj+1),myj.eq.nodey)
        ELSE
          call parcel_wrap_inplace(pry,miny,maxy)
        ENDIF
      ENDIF

      end subroutine parcel_migrate


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function parcel_global_count()
      use input
#ifdef MPI
      use mpi
#endif
      implicit none
!-----------------------------------------------------------------------
!  Global total of np_loc across all ranks (verification: must equal
!  nparcels every step -- no lost or duplicated parcels).
!-----------------------------------------------------------------------
      integer :: gtot

      gtot = np_loc
#ifdef MPI
      call MPI_ALLREDUCE(np_loc,gtot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
      parcel_global_count = gtot

      end function parcel_global_count


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_gather(pdata)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none
!-----------------------------------------------------------------------
!  Serial output gather (distributed mode):  assemble the global, id-ordered
!  pdata(nparcels,npvals) on rank 0 from every rank's local pdl/pidl, so the
!  existing rank-0 writer (parcel_write / writepdata_nc / GrADS) can write it
!  unchanged.  Each parcel id in 1..nparcels is owned by exactly one rank, so
!  pdata(id,:) is filled exactly once.  Only rank 0's pdata is meaningful on
!  return.  (The scalable alternative is the var16=1 parallel-write path.)
!-----------------------------------------------------------------------
      real, intent(inout), dimension(nparcels,npvals) :: pdata

#ifdef MPI
      integer :: r,nr,k,gid
      integer :: st(MPI_STATUS_SIZE)
      integer, allocatable :: rids(:)
      real,    allocatable :: rdat(:,:)

      IF( myid.eq.0 )THEN
        ! place rank 0's own local parcels:
        do k=1,np_loc
          pdata(pidl(k),1:npvals) = pdl(k,1:npvals)
        enddo
        ! receive and place every other rank's parcels:
        do r=1,numprocs-1
          call MPI_RECV(nr,1,MPI_INTEGER,r,301,MPI_COMM_WORLD,st,ierr)
          if( nr.gt.0 )then
            allocate( rids(nr) , rdat(npvals,nr) )
            call MPI_RECV(rids,       nr,MPI_INTEGER,r,302,MPI_COMM_WORLD,st,ierr)
            call MPI_RECV(rdat,npvals*nr,MPI_REAL,   r,303,MPI_COMM_WORLD,st,ierr)
            do k=1,nr
              gid = rids(k)
              pdata(gid,1:npvals) = rdat(1:npvals,k)
            enddo
            deallocate( rids , rdat )
          endif
        enddo
      ELSE
        call MPI_SEND(np_loc,1,MPI_INTEGER,0,301,MPI_COMM_WORLD,ierr)
        if( np_loc.gt.0 )then
          allocate( rids(np_loc) , rdat(npvals,np_loc) )
          do k=1,np_loc
            rids(k)        = pidl(k)
            rdat(1:npvals,k) = pdl(k,1:npvals)
          enddo
          call MPI_SEND(rids,       np_loc,MPI_INTEGER,0,302,MPI_COMM_WORLD,ierr)
          call MPI_SEND(rdat,npvals*np_loc,MPI_REAL,   0,303,MPI_COMM_WORLD,ierr)
          deallocate( rids , rdat )
        endif
      ENDIF
#else
      ! single-process build: pdl already holds everything in id order? No --
      ! copy local parcels into pdata by id.
      integer :: k
      do k=1,np_loc
        pdata(pidl(k),1:npvals) = pdl(k,1:npvals)
      enddo
#endif

      end subroutine parcel_gather


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_gather_positions(ploc)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none
!-----------------------------------------------------------------------
!  Gather just the parcel POSITIONS (prx,pry,prz) into the global,
!  id-ordered ploc(nparcels,3) on rank 0 -- for the restart file, which
!  only stores positions.  Unlike parcel_gather (all npvals into pdata),
!  this needs only rank-0's ploc to stay full-size, so pdata can be a
!  size-1 placeholder even on rank 0 in parallel-output mode.  Collective.
!-----------------------------------------------------------------------
      real, intent(inout), dimension(nparcels,3) :: ploc   ! valid on rank 0

#ifdef MPI
      integer :: r,nr,k,gid
      integer :: st(MPI_STATUS_SIZE)
      integer, allocatable :: rids(:)
      real,    allocatable :: rdat(:,:)

      IF( myid.eq.0 )THEN
        do k=1,np_loc
          ploc(pidl(k),1) = pdl(k,prx)
          ploc(pidl(k),2) = pdl(k,pry)
          ploc(pidl(k),3) = pdl(k,prz)
        enddo
        do r=1,numprocs-1
          call MPI_RECV(nr,1,MPI_INTEGER,r,401,MPI_COMM_WORLD,st,ierr)
          if( nr.gt.0 )then
            allocate( rids(nr) , rdat(3,nr) )
            call MPI_RECV(rids,  nr,MPI_INTEGER,r,402,MPI_COMM_WORLD,st,ierr)
            call MPI_RECV(rdat,3*nr,MPI_REAL,   r,403,MPI_COMM_WORLD,st,ierr)
            do k=1,nr
              gid = rids(k)
              ploc(gid,1) = rdat(1,k)
              ploc(gid,2) = rdat(2,k)
              ploc(gid,3) = rdat(3,k)
            enddo
            deallocate( rids , rdat )
          endif
        enddo
      ELSE
        call MPI_SEND(np_loc,1,MPI_INTEGER,0,401,MPI_COMM_WORLD,ierr)
        if( np_loc.gt.0 )then
          allocate( rids(np_loc) , rdat(3,np_loc) )
          do k=1,np_loc
            rids(k)   = pidl(k)
            rdat(1,k) = pdl(k,prx)
            rdat(2,k) = pdl(k,pry)
            rdat(3,k) = pdl(k,prz)
          enddo
          call MPI_SEND(rids,  np_loc,MPI_INTEGER,0,402,MPI_COMM_WORLD,ierr)
          call MPI_SEND(rdat,3*np_loc,MPI_REAL,   0,403,MPI_COMM_WORLD,ierr)
          deallocate( rids , rdat )
        endif
      ENDIF
#endif

      end subroutine parcel_gather_positions


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_write_par(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl)
      use input
      use writeout_nc_module, only : writepdata_nc_par
#ifdef MPI
      use mpi
#endif
      implicit none
!-----------------------------------------------------------------------
!  Parallel-output orchestrator (var16=1).  Redistribute this rank's local
!  parcels (pdl/pidl, spatially scattered ids) into contiguous global-id
!  BLOCKS via Alltoallv -- rank r ends up holding ids [r*bsz+1 .. r*bsz+nblk]
!  placed in id order in obuf -- then every rank writes its block to the one
!  shared file in parallel.  Never touches the global pdata array, so no rank
!  (not even rank 0) holds all parcels: this is the path that scales to 1e8.
!-----------------------------------------------------------------------
      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=40), intent(in), dimension(maxvars) :: name_prcl,desc_prcl,unit_prcl

#ifdef MPI
      integer :: bsz,bstart,nblk,k,dr,j,nsend,nrecv
      integer, allocatable :: scnt(:),rcnt(:),sdsp(:),rdsp(:),pos(:)
      integer, allocatable :: scv(:),rcv(:),sdv(:),rdv(:)
      integer, allocatable :: sid(:),rid(:)
      real,    allocatable :: sdat(:,:),rdat(:,:),obuf(:,:)

      bsz    = (nparcels+numprocs-1)/numprocs           ! ceil(nparcels/numprocs)
      bstart = myid*bsz
      nblk   = max( 0 , min( bsz , nparcels-bstart ) )

      allocate( scnt(0:numprocs-1),rcnt(0:numprocs-1) )
      allocate( sdsp(0:numprocs-1),rdsp(0:numprocs-1),pos(0:numprocs-1) )
      allocate( scv(0:numprocs-1),rcv(0:numprocs-1),sdv(0:numprocs-1),rdv(0:numprocs-1) )

      ! count local parcels destined for each id-block owner:
      scnt = 0
      do k=1,np_loc
        dr = (pidl(k)-1)/bsz
        if( dr.gt.numprocs-1 ) dr = numprocs-1
        scnt(dr) = scnt(dr)+1
      enddo
      call MPI_ALLTOALL(scnt,1,MPI_INTEGER,rcnt,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      sdsp(0)=0 ; rdsp(0)=0
      do k=1,numprocs-1
        sdsp(k) = sdsp(k-1)+scnt(k-1)
        rdsp(k) = rdsp(k-1)+rcnt(k-1)
      enddo
      nsend = sdsp(numprocs-1)+scnt(numprocs-1)
      nrecv = rdsp(numprocs-1)+rcnt(numprocs-1)         ! == nblk

      ! pack send buffers in destination order:
      allocate( sid(max(1,nsend)) , sdat(npvals,max(1,nsend)) )
      pos = sdsp
      do k=1,np_loc
        dr = (pidl(k)-1)/bsz
        if( dr.gt.numprocs-1 ) dr = numprocs-1
        pos(dr) = pos(dr)+1
        sid(pos(dr))    = pidl(k)
        sdat(:,pos(dr)) = pdl(k,:)
      enddo

      ! exchange ids, then the npvals-wide payload (counts/displs scaled):
      allocate( rid(max(1,nrecv)) , rdat(npvals,max(1,nrecv)) )
      call MPI_ALLTOALLV(sid,scnt,sdsp,MPI_INTEGER, rid,rcnt,rdsp,MPI_INTEGER, MPI_COMM_WORLD,ierr)
      scv = scnt*npvals ; rcv = rcnt*npvals ; sdv = sdsp*npvals ; rdv = rdsp*npvals
      call MPI_ALLTOALLV(sdat,scv,sdv,MPI_REAL, rdat,rcv,rdv,MPI_REAL, MPI_COMM_WORLD,ierr)

      ! scatter received parcels into id order within this rank's block:
      allocate( obuf(bsz,npvals) ) ; obuf = 0.0
      do k=1,nrecv
        j = rid(k)-bstart                               ! 1-based index in block
        obuf(j,:) = rdat(:,k)
      enddo

      deallocate( scnt,rcnt,sdsp,rdsp,pos,scv,rcv,sdv,rdv,sid,rid,sdat,rdat )

      call writepdata_nc_par(prec,rtime,qname,name_prcl,desc_prcl,unit_prcl,obuf,bstart,nblk,bsz)

      deallocate( obuf )
#endif
      end subroutine parcel_write_par


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine parcel_fixup_2d()
      use input
      implicit none
!-----------------------------------------------------------------------
!  2D-domain fixups for distributed local parcels (mirror of the legacy
!  pdata fixups): collapse the singleton horizontal dimension to 0.
!-----------------------------------------------------------------------
      integer :: n

      IF(axisymm.eq.1.or.ny.eq.1)THEN
        do n=1,np_loc
          pdl(n,pry) = 0.0
        enddo
      ENDIF
      IF(nx.eq.1)THEN
        do n=1,np_loc
          pdl(n,prx) = 0.0
        enddo
      ENDIF

      end subroutine parcel_fixup_2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE parcel_module
