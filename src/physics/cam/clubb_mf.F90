module clubb_mf

! =============================================================================== !
! Mass-flux module for use with CLUBB                                             !
! Together (CLUBB+MF) they comprise a eddy-diffusivity mass-flux approach (EDMF)  !
! =============================================================================== !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: integrate_mf, &
            clubb_mf_readnl, &
            do_clubb_mf, &
            do_clubb_mf_diag, &
            clubb_mf_nup

  real(r8) :: clubb_mf_L0   = 0._r8
  real(r8) :: clubb_mf_ent0 = 0._r8
  integer, protected :: clubb_mf_nup  = 0
  logical, protected :: do_clubb_mf = .false.
  logical, protected :: do_clubb_mf_diag = .false.

  contains

  subroutine clubb_mf_readnl(nlfile)

  ! =============================================================================== !
  ! MF namelists                                                                    !
  ! =============================================================================== !

    use namelist_utils,  only: find_group_name
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, mpi_integer, mpi_logical

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    character(len=*), parameter :: sub = 'clubb_mf_readnl'

    integer :: iunit, read_status, ierr


    namelist /clubb_mf_nl/ clubb_mf_L0, clubb_mf_ent0, clubb_mf_nup, do_clubb_mf, do_clubb_mf_diag

    if (masterproc) then
      open( newunit=iunit, file=trim(nlfile), status='old' )
      call find_group_name(iunit, 'clubb_mf_nl', status=read_status)
      if (read_status == 0) then
         read(iunit, clubb_mf_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun('clubb_mf_readnl: ERROR reading namelist')
         end if
      end if
      close(iunit)
    end if

    call mpi_bcast(clubb_mf_L0,   1, mpi_real8,   mstrid, mpicom, ierr) 
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_L0")
    call mpi_bcast(clubb_mf_ent0, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ent0")
    call mpi_bcast(clubb_mf_nup,  1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_nup")
    call mpi_bcast(do_clubb_mf,      1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf")
    call mpi_bcast(do_clubb_mf_diag, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_diag")

    if ((.not. do_clubb_mf) .and. do_clubb_mf_diag ) then
       call endrun('clubb_mf_readnl: Error - cannot turn on do_clubb_mf_diag without also turning on do_clubb_mf')
    end if
    

  end subroutine clubb_mf_readnl

  subroutine integrate_mf( nz,      dzm,     zm,      p_zm,      iexner_zm,         & ! input
                           rho_zt,  dzt,     zt,      p_zt,      iexner_zt,         & ! input
                           u,       v,       thl,     qt,        thv,               & ! input
                                             th,      qv,        qc,                & ! input
                                             thl_zm,  qt_zm,                        & ! input
                                             th_zm,   qv_zm,     qc_zm,             & ! input
                                             wthl,    wqt,       pblh,              & ! input
                           upa,                                                     & ! output - plume diagnostics
                           upw,                                                     & ! output - plume diagnostics
                           upqt,                                                    & ! output - plume diagnostics
                           upthl,                                                   & ! output - plume diagnostics
                           upthv,                                                   & ! output - plume diagnostics
                           upth,                                                    & ! output - plume diagnostics
                           upqc,                                                    & ! output - plume diagnostics
                           upbuoy,                                                  & ! output - plume diagnostics
                           upent,                                                   & ! output - plume diagnostics
                           ae,      aw,                                             & ! output - diagnosed fluxes BEFORE mean field update
                           awthl,   awqt,                                           & ! output - diagnosed fluxes BEFORE mean field update
                           awql,    awqi,                                           & ! output - diagnosed fluxes BEFORE mean field update
                           awth,    awqv,                                           & ! output - diagnosed fluxes BEFORE mean field update
                           awu,     awv,                                            & ! output - diagnosed fluxes BEFORE mean field update
                           thflx,   qvflx,                                          & ! output - diagnosed fluxes BEFORE mean field update
                                    qcflx,                                          & ! output - diagnosed fluxes BEFORE mean field update
                           thlflx,  qtflx )                                           ! output - variables needed for solver

  ! ================================================================================= !
  ! Mass-flux algorithm                                                               !
  !                                                                                   !
  ! Provides rtm and thl fluxes due to mass flux ensemble,                            !
  ! which are fed into the mixed explicit/implicit clubb solver as explicit terms     !
  !                                                                                   !
  ! Mass flux variables are computed on edges (i.e. momentum grid):                   !
  ! upa,upw,upqt,...                                                                  !
  ! dry_a,moist_a,dry_w,moist_w, ...                                                  !
  !                                                                                   ! 
  ! In CLUBB (unlike CAM) nlevs of momentum grid = nlevs of thermodynamic grid,       !
  ! due to a subsurface thermodynamic layer. To avoid confusion, below the variables  !  
  ! are grouped by the grid they are on.                                              !
  !                                                                                   !
  ! *note that state on the lowest thermo level is equal to state on the lowest       !
  ! momentum level due to state_zt(1) = state_zt(2), and lowest momentum level        !
  ! is a weighted combination of the lowest two thermodynamic levels.                 !
  !                                                                                   !
  ! ---------------------------------Authors----------------------------------------  !
  ! Marcin Kurowski, JPL                                                              !
  ! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM       !
  ! Additional modifications by Adam Herrington, NCAR                                 !
  ! ================================================================================= !

     use physconst,          only: rair, cpair, gravit, latvap, latice, zvir

     integer,  intent(in)                :: nz
     real(r8), dimension(nz), intent(in) :: u,      v,            & ! thermodynamic grid
                                            thl,    thv,          & ! thermodynamic grid
                                            th,     qv,           & ! thermodynamic grid
                                            qt,     qc,           & ! thermodynamic grid
                                            rho_zt,               & ! thermodynamic grid
                                            dzt,    zt,           & ! thermodynamic grid
                                            p_zt,   iexner_zt,    & ! thermodynamic grid
                                            thl_zm,               & ! momentum grid
                                            th_zm,  qv_zm,        &
                                            qt_zm,  qc_zm,        & ! momentum grid
                                            dzm,    zm,           & ! momentum grid
                                            p_zm,   iexner_zm       ! momentum grid

     real(r8), intent(in)                :: wthl,wqt
     real(r8), intent(inout)             :: pblh

     real(r8),dimension(nz,clubb_mf_nup), intent(out) :: upa,     & ! momentum grid
                                                         upw,     & ! momentum grid
                                                         upqt,    & ! momentum grid
                                                         upthl,   & ! momentum grid
                                                         upthv,   & ! momentum grid
                                                         upth,    & ! momentum grid
                                                         upqc,    & ! momentum grid
                                                         upbuoy,  & ! momentum grid
                                                         upent

     real(r8),dimension(nz), intent(out) :: ae,      aw,          & ! momentum grid
                                            awthl,   awqt,        & ! momentum grid
                                            awql,    awqi,        & ! momentum grid
                                            awth,    awqv,        & ! momentum grid
                                            awu,     awv,         & ! momentum grid
                                            thflx,   qvflx,       & ! momentum grid 
                                            qcflx,                & ! momentum grid
                                            thlflx,  qtflx          ! momentum grid

     ! =============================================================================== !
     ! INTERNAL VARIABLES
     !
     ! sums over all plumes
     real(r8), dimension(nz)              :: moist_th,   dry_th,      & ! momentum grid
                                             awqc,                    & ! momentum grid                     
                                             awthl_conv, awqt_conv,   & ! momentum grid
                                             thl_env,    qt_env         ! thermodynamic grid
     !
     ! updraft properties
     real(r8), dimension(nz,clubb_mf_nup) :: upqv,     upqs,           & ! momentum grid
                                             upql,     upqi,           & ! momentum grid
                                             upu,      upv,            & ! momentum grid 
                                             uplmix                      ! momentum grid
     !
     ! microphyiscs terms
     real(r8), dimension(nz,clubb_mf_nup) :: supqt,    supthl,         & ! thermodynamic grid 
                                             uprr                        ! thermodynamic grid
     !
     ! entrainment profiles
     real(r8), dimension(nz,clubb_mf_nup) :: ent,      entf              ! thermodynamic grid
     integer,  dimension(nz,clubb_mf_nup) :: enti                        ! thermodynamic grid
     ! 
     ! other variables
     integer                              :: k,kstart,i,ih
     real(r8), dimension(clubb_mf_nup)    :: zcb
     real(r8)                             :: zcb_unset,                &
                                             wthv,                     &
                                             wstar,  qstar,   thvstar, & 
                                             sigmaw, sigmaqt, sigmathv,&
                                                     wmin,    wmax,    & 
                                             wlv,    wtv,     wp,      & 
                                             B,                        & ! thermodynamic grid
                                             entexp, entexpu, entw,    & ! thermodynamic grid
                                             lmixt,                    & ! thermodynamic grid
                                             qtovqs, sevap,            & ! thermodynamic grid
                                             betathl,betaqt,           & ! thermodynamic grid        
                                             thln,   thvn,    thn,     & ! momentum grid
                                             qtn,    qsn,              & ! momentum grid
                                             qcn,    qln,     qin,     & ! momentum grid
                                             un,     vn,      wn2,     & ! momentum grid
                                             lmixn,   srfarea,         & ! momentum grid
                                             srfwqtu, srfwthvu,        &
                                             facqtu,  facthvu
     
     ! parameters defining initial conditions for updrafts
     real(r8),parameter                   :: pwmin = 1.5_r8,           &
                                             pwmax = 3._r8

     !
     ! alpha, z-scores after Suselj etal 2019
     real(r8),parameter                   :: alphw   = 0.572_r8,       &
                                             alphqt  = 2.890_r8,       &     
                                             alphthv = 2.890_r8
     !
     ! w' covariance after Suselj etal 2019
     real(r8),parameter                   :: cwqt  = 0.32_r8,          &
                                             cwthv = 0.58_r8
     !
     ! virtual mass coefficients for w-eqn after Suselj etal 2019
     real(r8),parameter                   :: wa = 1.0_r8,              &
                                             wb = 1.5_r8
     !
     ! min values to avoid singularities
     real(r8),parameter                   :: wstarmin = 1.e-3_r8,      &
                                             pblhmin  = 100._r8
     !
     ! to condensate or not to condensate
     logical                              :: do_condensation = .true.
     !
     ! to precip or not to precip
     logical                              :: do_precip = .false.
     !
     ! evaporation efficiency after Suselj etal 2019
     real(r8),parameter                   :: ke = 2.5e-4_r8
     !
     ! fraction of rain detrained into downdrafts
     real(r8),parameter                   :: fdd = 0._r8
     !
     ! fixed entrainment rate (debug only)
     real(r8),parameter                   :: fixent = 0.004_r8
     !
     ! to scale surface fluxes
     logical                              :: scalesrf = .false. 
     !
     ! to dry flux
     logical                              :: dryflux = .false.
     !
     ! to debug flag
     logical                              :: debug  = .false.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! INITIALIZE OUTPUT VARIABLES
     ! outputs - variables needed for solver
     aw        = 0._r8
     awth      = 0._r8
     awthl     = 0._r8
     awqt      = 0._r8
     awqv      = 0._r8
     awqc      = 0._r8
     awql      = 0._r8
     awqi      = 0._r8
     awu       = 0._r8
     awv       = 0._r8
     thlflx    = 0._r8
     qtflx     = 0._r8
     thflx     = 0._r8
     qvflx     = 0._r8
     qcflx     = 0._r8

     ent       = 0._r8
     entf      = 0._r8
     enti      = 0

     ! this is the environmental area - by default 1.
     ae = 1._r8

     ! START MAIN COMPUTATION
     upw   = 0._r8
     upth  = 0._r8
     upthl = 0._r8
     upthv = 0._r8
     upqt  = 0._r8
     upa   = 0._r8
     upu   = 0._r8
     upv   = 0._r8
     upqc  = 0._r8
     upth  = 0._r8
     upql  = 0._r8
     upqi  = 0._r8
     upqv  = 0._r8
     upqs  = 0._r8
     upbuoy= 0._r8
     uplmix= 0._r8
     uprr  = 0._r8
     upent = 0._r8
     supqt = 0._r8
     supthl= 0._r8

     ! unique identifier
     zcb_unset = 9999999._r8
     zcb       = zcb_unset

     pblh = max(pblh,pblhmin)
     wthv = wthl+zvir*thv(1)*wqt

     ! if surface buoyancy is positive then do mass-flux
     if ( wthv > 0._r8 ) then

       if (debug) then
         ! overide stochastic entrainment with fixent
         ent(:,:) = fixent
       else

         ! get entrainment coefficient, dz/L0
         do i=1,clubb_mf_nup
           do k=1,nz
             entf(k,i) = dzm(k) / clubb_mf_L0
           enddo
         enddo

         ! get poisson, P(dz/L0)
         call poisson( nz, clubb_mf_nup, entf, enti, u(2:5))

         ! get entrainment, ent=ent0/dz*P(dz/L0)
         do i=1,clubb_mf_nup
           do k=1,nz
             ent(k,i) = real( enti(k,i))*clubb_mf_ent0/dzm(k)
           enddo
         enddo

       end if

       ! get surface conditions
       wstar   = max( wstarmin, (gravit/thv(1)*wthv*pblh)**(1._r8/3._r8) )
       qstar   = wqt / wstar
       thvstar = wthv / wstar

       sigmaw   = alphw * wstar
       sigmaqt  = alphqt * abs(qstar)
       sigmathv = alphthv * abs(thvstar)

       wmin = sigmaw * pwmin
       wmax = sigmaw * pwmax

       do i=1,clubb_mf_nup
         wlv = wmin + (wmax-wmin) / (real(clubb_mf_nup,r8)) * (real(i-1, r8))
         wtv = wmin + (wmax-wmin) / (real(clubb_mf_nup,r8)) * real(i,r8)

         upw(1,i) = 0.5_r8 * (wlv+wtv)
         upa(1,i) = 0.5_r8 * erf( wtv/(sqrt(2._r8)*sigmaw) ) &
                    - 0.5_r8 * erf( wlv/(sqrt(2._r8)*sigmaw) )

         upu(1,i) = u(1)
         upv(1,i) = v(1)

         upqt(1,i)  = cwqt * upw(1,i) * sigmaqt/sigmaw
         upthv(1,i) = cwthv * upw(1,i) * sigmathv/sigmaw
       enddo

       facqtu=1._r8
       facthvu=1._r8
       if (scalesrf) then 
         ! scale surface fluxes
         srfwqtu = 0._r8
         srfwthvu = 0._r8
         srfarea = 0._r8
         do i=1,clubb_mf_nup
             srfwqtu=srfwqtu+upqt(1,i)*upw(1,i)*upa(1,i)
             srfwthvu=srfwthvu+upthv(1,i)*upw(1,i)*upa(1,i)
             srfarea = srfarea+upa(1,i)
         end do
         facqtu=srfarea*wqt/srfwqtu
         facthvu=srfarea*wthv/srfwthvu
       end if

       do i=1,clubb_mf_nup

         betaqt = (qt(4)-qt(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         if (betaqt*0.5_r8*(dzt(2)+dzt(1)) .le. qt(2)) then
           upqt(1,i)= qt(2)-betaqt*0.5_r8*(dzt(2)+dzt(1))+facqtu*upqt(1,i)
         else
           upqt(1,i)= qt(1)+facqtu*upqt(1,i)
         end if

         betathl = (thv(4)-thv(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         upthv(1,i)= thv(2)-betathl*0.5_r8*(dzt(2)+dzt(1))+facthvu*upthv(1,i)

         upthl(1,i) = upthv(1,i) / (1._r8+zvir*upqt(1,i))
         upth(1,i)  = upthl(1,i)

         ! get cloud, lowest thermodynamic level (assume it has properties of lowest momentum level) 
         if (do_condensation) then
           call condensation_mf(upqt(1,i), upthl(1,i), p_zm(1), iexner_zm(1), &
                                thvn, qcn, thn, qln, qin, qsn, lmixn)
           upthv(1,i) = thvn
           upqc(1,i)  = qcn
           upql(1,i)  = qln
           upqi(1,i)  = qin
           upqs(1,i)  = qsn
           upth(1,i)  = thn
           if (qcn > 0._r8) zcb(i) = zm(1)
         else
           ! assume no cldliq
           upqc(1,i)  = 0._r8
         end if
       end do

       ! get updraft properties 
       do i=1,clubb_mf_nup
         do k=1,nz-1

           ! get microphysics, autoconversion
           if (do_precip .and. upqc(k,i) > 0._r8) then
             call precip_mf(upqs(k,i),upqt(k,i),upw(k,i),dzm(k),zt(k+1)-zcb(i),supqt(k,i))

             supthl(k,i) = -1._r8*lmixn*supqt(k,i)*iexner_zm(k)/cpair
           else
             supqt(k,i)  = 0._r8
             supthl(k,i) = 0._r8
           end if

           ! integrate updraft
           entexp  = exp(-ent(k,i)*dzm(k))
           entexpu = exp(-ent(k,i)*dzm(k)/3._r8)

           qtn  = qt_zm(k) *(1._r8-entexp ) + upqt (k,i)*entexp + supqt(k,i)
           thln = thl_zm(k)*(1._r8-entexp ) + upthl(k,i)*entexp + supthl(k,i)           
           !un   = u(k+1)  *(1._r8-entexpu) + upu  (k,i)*entexpu
           !vn   = v(k+1)  *(1._r8-entexpu) + upv  (k,i)*entexpu

           ! get cloud, thermodynamics levels
           if (do_condensation) then
             call condensation_mf(qtn, thln, p_zt(k+1), iexner_zt(k+1), &
                                  thvn, qcn, thn, qln, qin, qsn, lmixn)
             if (zcb(i).eq.zcb_unset .and. qcn > 0._r8) zcb(i) = zt(k+1)
           else
             thvn = thln*(1._r8+zvir*qtn)
           end if

           ! get buoyancy
           B=gravit*(thvn/thv(k+1)-1._r8)
           if (debug) then
             if ( masterproc ) then
               write(iulog,*) "B(k,i), k, i ", B, k, i
             end if
           end if

           ! get wn^2
           wp = wb*0.5_r8*(ent(k,i)+ent(k+1,i))
           if (wp==0._r8) then
             wn2 = upw(k,i)**2._r8+2._r8*wa*B*dzt(k+1)
           else
             entw = exp(-2._r8*wp*dzt(k+1))
             wn2 = entw*upw(k,i)**2._r8+wa*B/wp*(1._r8-entw)
           end if

           if (wn2>0._r8) then
             upw(k+1,i)   = sqrt(wn2)
             upthv(k+1,i) = thvn
             upthl(k+1,i) = thln
             upqt(k+1,i)  = qtn
             upqc(k+1,i)  = qcn
             upqs(k+1,i)  = qsn
             upu(k+1,i)   = un
             upv(k+1,i)   = vn
             upa(k+1,i)   = upa(k,i)
             upql(k+1,i)  = qln
             upqi(k+1,i)  = qin
             if (qtn - qcn.gt.0._r8) upqv(k+1,i) = qtn - qcn
             uplmix(k+1,i)= lmixn
             upth(k+1,i)  = thn
             upbuoy(k+1,i)= B
             upent(k,i) = ent(k,i)
           else
             exit
           end if
         enddo
       enddo

       ! downward sweep for rain evaporation, snow melting 
       if (do_precip) then
         do i=1,clubb_mf_nup
           do k=nz-1,1,-1
             ! get rain evaporation
             if ((upqs(k,i) + upqs(k-1,i)).le.0._r8) then
               qtovqs = 0._r8
             else
               qtovqs = (upqt(k,i) + upqt(k-1,i))/(upqs(k,i) + upqs(k-1,i))
             end if
             qtovqs = min(1._r8,qtovqs)
             sevap = ke*(1._r8 - qtovqs)*sqrt(max(uprr(k+1,i),0._r8))

             ! get rain rate
             uprr(k,i) = uprr(k+1,i) &
                         - rho_zt(k)*dzt(k)*( supqt(k,i)*(1._r8-fdd) + sevap )

             if (debug) then
               if ( masterproc ) then
                 write(iulog,*) "uprr(k,i), k, i ", uprr(k,i), k, i
               end if
             end if

             ! update source terms
             lmixt = 0.5_r8*(uplmix(k,i)+uplmix(k-1,i))
             supqt(k,i) = supqt(k,i) + sevap
             supthl(k,i) = supthl(k,i) - lmixt*sevap*iexner_zt(k)/cpair
           end do
         end do
       end if

       do k=1,nz
         do i=1,clubb_mf_nup
           ae  (k) = ae  (k) - upa(k,i)
           aw  (k) = aw  (k) + upa(k,i)*upw(k,i)
           awu (k) = awu (k) + upa(k,i)*upw(k,i)*upu(k,i)
           awv (k) = awv (k) + upa(k,i)*upw(k,i)*upv(k,i)
           awthl(k)= awthl(k)+ upa(k,i)*upw(k,i)*upthl(k,i) 
           awth(k) = awth(k) + upa(k,i)*upw(k,i)*upth(k,i)
           awqt(k) = awqt(k) + upa(k,i)*upw(k,i)*upqt(k,i)
           awqv(k) = awqv(k) + upa(k,i)*upw(k,i)*upqv(k,i)
           awql(k) = awql(k) + upa(k,i)*upw(k,i)*upql(k,i)
           awqi(k) = awqi(k) + upa(k,i)*upw(k,i)*upqi(k,i)
           awqc(k) = awqc(k) + upa(k,i)*upw(k,i)*upqc(k,i)
         enddo
       enddo

       if (dryflux) then
         awthl_conv = awth
         awqt_conv = awqv
         thl_env = th
         qt_env = qv
       else
         awthl_conv = awthl       
         awqt_conv = awqt
         thl_env = thl
         qt_env = qt
       end if

       ! use upwinding to compute fluxes

       ! get thl & qt fluxes
       if (scalesrf) then
         kstart = 1
         betathl = (thl_env(4)-thl_env(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         betaqt = (qt_env(4)-qt_env(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         thl_env(1) = thl_env(2)-betathl*0.5_r8*(dzt(2)+dzt(1))
         qt_env(1) = qt_env(2)-betaqt*0.5_r8*(dzt(2)+dzt(1))
         if (qt_env(1).lt.0._r8) qt_env(1) = 0._r8
       else
         kstart = 2
       end if

       do k=kstart,nz
         thlflx(k)= awthl_conv(k) - aw(k)*thl_env(k)
         qtflx(k)= awqt_conv(k) - aw(k)*qt_env(k)
       enddo

       ! get th & qv fluxes
       thl_env = th
       qt_env = qv
       if (scalesrf) then
         betathl = (th(4)-th(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         betaqt = (qv(4)-qv(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         thl_env(1) = thl_env(2)-betathl*0.5_r8*(dzt(2)+dzt(1))
         qt_env(1) = qt_env(2)-betaqt*0.5_r8*(dzt(2)+dzt(1))
         if (qt_env(1).lt.0._r8) qt_env(1) = 0._r8
       end if

       do k=kstart,nz
         thflx(k)= awth(k) - aw(k)*thl_env(k)
         qvflx(k)= awqv(k) - aw(k)*qt_env(k)
       enddo

       ! get qc fluxes
       qt_env = qc
       if (scalesrf) then
         betaqt = (qc(4)-qc(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         qt_env(1) = qt_env(2)-betaqt*0.5_r8*(dzt(2)+dzt(1))
         if (qt_env(1).lt.0._r8) qt_env(1) = 0._r8
       end if

       do k=kstart,nz
         qcflx(k)= awqc(k) - aw(k)*qt_env(k)
       enddo

     end if  ! ( wthv > 0.0 )

  end subroutine integrate_mf

  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi, qs, lmix )
  ! =============================================================================== !
  ! zero or one condensation for edmf: calculates thv and qc                        !
  ! =============================================================================== !
     use physconst,          only: cpair, zvir, h2otrip
     use wv_saturation,      only : qsat

     real(r8),intent(in) :: qt,thl,p,iex
     real(r8),intent(out):: thv,qc,th,ql,qi,qs,lmix

     !local variables
     integer  :: niter,i
     real(r8) :: diff,t,qstmp,qcold,es,wf

     ! max number of iterations
     niter=50
     ! minimum difference
     diff=2.e-5_r8

     qc=0._r8
     t=thl/iex

     !by definition:
     ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
     ! Thl = Th - L/cp*ql/Exner              (2)
     !so:
     ! Th  = Thl + L/cp*ql/Exner             (3)
     ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
     !     = Thl*Exner + L/cp*ql
     do i=1,niter
       wf = get_watf(t)
       !wf = 1._r8
       t = thl/iex+get_alhl(wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal (check!)
       call qsat(t,p,es,qstmp)
       qcold = qc
       qc = max(0.5_r8*qc+0.5_r8*(qt-qstmp),0._r8)
       if (abs(qc-qcold)<diff) exit
     enddo

     wf = get_watf(t)
     !wf = 1._r8 
     t = thl/iex+get_alhl(wf)/cpair*qc

     call qsat(t,p,es,qs)
     qc = max(qt-qs,0._r8)
     thv = (thl+get_alhl(wf)/cpair*iex*qc)*(1._r8+zvir*(qt-qc)-qc)
     lmix = get_alhl(wf)
     th = t*iex
     qi = qc*(1._r8-wf)
     ql = qc*wf

     contains

     function get_watf(t)
       real(r8)            :: t,get_watf,tc
       real(r8), parameter :: &
                              tmax=-10._r8, &
                              tmin=-40._r8

       tc=t-h2otrip

       if (tc>tmax) then
         get_watf=1._r8
       else if (tc<tmin) then
         get_watf=0._r8
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(wf)
     !latent heat of the mixture based on water fraction
       use physconst,        only : latvap , latice
       real(r8) :: get_alhl,wf

       get_alhl = wf*latvap+(1._r8-wf)*(latvap+latice)

     end function get_alhl

  end subroutine condensation_mf

  subroutine precip_mf(qs,qt,w,dz,dzcld,Supqt)
  !**********************************************************************
  ! Precipitation microphysics
  ! By Adam Herrington, after Kay Suselj
  !**********************************************************************

       real(r8),intent(in)  :: qs,qt,w,dz,dzcld
       real(r8),intent(out) :: Supqt
       ! 
       ! local vars
       real(r8)            :: tauwgt, tau,       & ! time-scale vars
                              qstar                ! excess cloud liquid                   

       real(r8),parameter  :: tau0  = 15._r8,    & ! base time-scale
                              zmin  = 300._r8,   & ! small cloud thick
                              zmax  = 3000._r8,  & ! large cloud thick
                              qcmin = 0.00125_r8   ! supersat threshold 

       qstar = qs+qcmin
       
       if (qt > qstar) then
         ! get precip efficiency
         tauwgt = (dzcld-zmin)/(zmax-zmin)
         tauwgt = min(max(tauwgt,0._r8),1._r8)
         tau    = tauwgt/tau0
 
         ! get source for updraft
         Supqt = (qstar-qt)*(1._r8 - exp(-1._r8*tau*dz/w))
       else
         Supqt = 0._r8
       end if

  end subroutine precip_mf

  subroutine poisson(nz,nup,lambda,poi,state)
  !**********************************************************************
  ! Set a unique (but reproduceble) seed for the kiss RNG
  ! Call Poisson deviate
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       integer,                     intent(in)  :: nz,nup
       real(r8), dimension(4),      intent(in)  :: state
       real(r8), dimension(nz,nup), intent(in)  :: lambda
       integer,  dimension(nz,nup), intent(out) :: poi
       integer,  dimension(1,4)                 :: tmpseed
       integer                                  :: i,j
       type(ShrKissRandGen)                     :: kiss_gen

       ! Compute seed
       tmpseed(1,1) = int((state(1) - int(state(1))) * 1000000000._r8)
       tmpseed(1,2) = int((state(2) - int(state(2))) * 1000000000._r8)
       tmpseed(1,3) = int((state(3) - int(state(3))) * 1000000000._r8)
       tmpseed(1,4) = int((state(4) - int(state(4))) * 1000000000._r8)

       ! Set seed
       kiss_gen = ShrKissRandGen(tmpseed)

       do i=1,nz
         do j=1,nup
           call knuth(kiss_gen,lambda(i,j),poi(i,j))
         enddo
       enddo

  end subroutine poisson

  subroutine knuth(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth 
  ! The Art of Computer Programming, v2, 137-138
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(r8),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       ! Local variables
       real(r8), dimension(1,1) :: tmpuni
       real(r8)                 :: puni, explam
       integer                  :: k

       k = 0
       explam = exp(-1._r8*lambda)
       puni = 1._r8
       do while (puni > explam)
         k = k + 1
         call kiss_gen%random(tmpuni)
         puni = puni*tmpuni(1,1)
       end do
       kout = k - 1

  end subroutine knuth

end module clubb_mf
