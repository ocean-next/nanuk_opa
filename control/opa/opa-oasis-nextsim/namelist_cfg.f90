!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OPA  :  1 - run manager      (namrun)
!! namelists    2 - Domain           (namcfg, namzgr, namzgr_sco, namdom, namtsd)
!!              3 - Surface boundary (namsbc, namsbc_ana, namsbc_flx, namsbc_clio, namsbc_core, namsbc_sas
!!                                    namsbc_cpl, namtra_qsr, namsbc_rnf,
!!                                    namsbc_apr, namsbc_ssr, namsbc_alb)
!!              4 - lateral boundary (namlbc, namcla, namobc, namagrif, nambdy, nambdy_tide)
!!              5 - bottom  boundary (nambfr, nambbc, nambbl)
!!              6 - Tracer           (nameos, namtra_adv, namtra_ldf, namtra_dmp)
!!              7 - dynamics         (namdyn_adv, namdyn_vor, namdyn_hpg, namdyn_spg, namdyn_ldf)
!!              8 - Verical physics  (namzdf, namzdf_ric, namzdf_tke, namzdf_kpp, namzdf_ddm, namzdf_tmx, namzdf_tmx_new)
!!              9 - diagnostics      (namnc4, namtrd, namspr, namflo, namhsb, namsto)
!!             10 - miscellaneous    (namsol, nammpp, namctl)
!!             11 - Obs & Assim      (namobs, nam_asminc)
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!======================================================================
!!                   ***  Run management namelists  ***
!!======================================================================
!!   namrun       parameters of the run
!!======================================================================
!
!-----------------------------------------------------------------------
&namrun        !   parameters of the run
!-----------------------------------------------------------------------
   nn_no       =       0   !  job number (no more used...)
   cn_exp      =  '<CONFCASE>'
   nn_it000    =    <IT000> !  first time step
   nn_itend    =    <ITEND> !  last  time step (std 5475)
   nn_date0    =   <NDATE0> !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
   nn_leapy    =       1    !  Leap year calendar (1) or not (0)
   ln_rstart   = .<RSTRT>.  !  start from rest (F) or from a restart file (T)
   nn_euler    =       1    !  = 0 : start with forward time step if ln_rstart=T  #lulu
   nn_rstctl   =  <IRCTL>   !  restart control => activated only if ln_rstart = T
                            !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
                            !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
                            !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
   cn_ocerst_in    = '<CN_OCERST_IN>'
   cn_ocerst_indir = '<CN_OCERST_INDIR>'   !  directory from which to read input ocean restarts
   cn_ocerst_out = 'restart_oce'
   cn_ocerst_outdir = '<CN_OCERST_OUTDIR>'      !  directory in which to write output ocean restarts
   nn_istate   =   0    !  output the initial state (1) or not (0)
   nn_stock    =   35040   !1Y at dt=900  frequency of creation of a restart file (modulo referenced to 1)
   nn_write    =   35040   !1Y at dt=900  frequency of write in the output file   (modulo referenced to nn_it000)
   ln_dimgnnn  = .false.   !  DIMG file format: 1 file for all processors (F) or by processor (T)
   ln_mskland  = .false.   !  mask land points in NetCDF outputs (costly: + ~15%)
   ln_cfmeta   = .false.   !  output additional data to netCDF files required for compliance with the CF metadata standard
   ln_clobber  = .true.    !  clobber (overwrite) an existing file
   nn_chunksz  =       0   !  chunksize (bytes) for NetCDF file (works only with iom_nf90 routines)
/
!
!!======================================================================
!!                      ***  Domain namelists  ***
!!======================================================================
!!   namcfg       parameters of the configuration
!!   namzgr       vertical coordinate
!!   namzgr_sco   s-coordinate or hybrid z-s-coordinate
!!   namdom       space and time domain (bathymetry, mesh, timestep)
!!   namtsd       data: temperature & salinity
!!======================================================================
!
!-----------------------------------------------------------------------
&namcfg     !   parameters of the configuration
!-----------------------------------------------------------------------
   cp_cfg      =  "nanuk"              !  name of the configuration
   cp_cfz      =  "no zoom"            !  name of the zoom of configuration
   jp_cfg      =       025             !  resolution of the configuration
   jpidta      =     528               !  1st lateral dimension ( >= jpi )
   jpjdta      =     603               !  2nd    "         "    ( >= jpj )
   jpkdta      =      75               !  number of levels      ( >= jpk )
   jpiglo      =     528               !  1st dimension of global domain --> i =jpidta
   jpjglo      =     603               !  2nd    -                  -    --> j =jpjdta
   jpizoom     =       1               !  left bottom (i,j) indices of the zoom
   jpjzoom     =       1               !  in data domain indices
   jperio      =       0               !  lateral cond. type (between 0 and 6)
                                       !  = 0 closed                 ;   = 1 cyclic East-West
                                       !  = 2 equatorial symmetric   ;   = 3 North fold T-point pivot
                                       !  = 4 cyclic East-West AND North fold T-point pivot
                                       !  = 5 North fold F-point pivot
                                       !  = 6 cyclic East-West AND North fold F-point pivot
   ln_use_jattr = .false.              !  use (T) the file attribute: open_ocean_jstart, if present
                                       !  in netcdf input files, as the start j-row for reading
/
!-----------------------------------------------------------------------
&namzgr        !   vertical coordinate
!-----------------------------------------------------------------------
   ln_zco      = .false.   !  z-coordinate - full    steps   (T/F)      ("key_zco" may also be defined)
   ln_zps      = .true.    !  z-coordinate - partial steps   (T/F)
   ln_sco      = .false.   !  s- or hybrid z-s-coordinate    (T/F)
   ln_isfcav   = .false.   !  ice shelf cavity               (T/F)
/
!-----------------------------------------------------------------------
&namzgr_sco    !   s-coordinate or hybrid z-s-coordinate
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namdom        !   space and time domain (bathymetry, mesh, timestep)
!-----------------------------------------------------------------------
   nn_bathy    =    1      !  compute (=0) or read (=1) the bathymetry file
   rn_bathy    =    0.     !  value of the bathymetry. if (=0) bottom flat at jpkm1
   nn_closea   =    0      !  remove (=0) or keep (=1) closed seas and lakes (ORCA)
   nn_msh      =    0      !  create (=1) a mesh file or not (=0)
   rn_hmin     =   -3.     !  min depth of the ocean (>0) or min number of ocean level (<0)
   rn_e3zps_min=   25.     !  partial step thickness is set larger than the minimum of
   rn_e3zps_rat=    0.2    !  rn_e3zps_min and rn_e3zps_rat*e3t, with 0<rn_e3zps_rat<1
                           !
   rn_rdt      =  900.     !  time step for the dynamics (and tracer if nn_acc=0)
   rn_atfp     =    0.1    !  asselin time filter parameter
   nn_acc      =    0      !  acceleration of convergence : =1      used, rdt < rdttra(k)
                                 !                          =0, not used, rdt = rdttra
   rn_rdtmin   = 28800.          !  minimum time step on tracers (used if nn_acc=1)
   rn_rdtmax   = 28800.          !  maximum time step on tracers (used if nn_acc=1)
   rn_rdth     =  800.           !  depth variation of tracer time step  (used if nn_acc=1)
   ln_crs      = .false.      !  Logical switch for coarsening module
   jphgr_msh   =       0               !  type of horizontal mesh
                                       !  = 0 curvilinear coordinate on the sphere read in coordinate.nc
                                       !  = 1 geographical mesh on the sphere with regular grid-spacing
                                       !  = 2 f-plane with regular grid-spacing
                                       !  = 3 beta-plane with regular grid-spacing
                                       !  = 4 Mercator grid with T/U point at the equator
   ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
   ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
   ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
   ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
   ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
   ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
   ppsur       =   -3958.951371276829   !  ORCA r4, r2 and r05 coefficients
   ppa0        =    103.9530096000000   ! (default coefficients)
   ppa1        =    2.415951269000000   !
   ppkth       =    15.35101370000000   !
   ppacr       =    7.000000000000000             !
   ppdzmin     =  999999.              !  Minimum vertical spacing
   pphmax      =  999999.              !  Maximum depth
   ldbletanh   =  .TRUE.              !  Use/do not use double tanf function for vertical coordinates
   ppa2        =  100.7609285000000              !  Double tanh function parameters
   ppkth2      =   48.02989372000000             !
   ppacr2      =   13.00000000000             !
/
!-----------------------------------------------------------------------
&namsplit      !   time splitting parameters                            ("key_dynspg_ts")
!-----------------------------------------------------------------------
   ln_bt_fw      =    .FALSE.          !  Forward integration of barotropic equations
   ln_bt_av      =    .TRUE.           !  Time filtering of barotropic variables
   ln_bt_nn_auto =    .TRUE.           !  Set nn_baro automatically to be just below
                                       !  a user defined maximum courant number (rn_bt_cmax)
   nn_baro       =    48               !  Number of iterations of barotropic mode
                                       !  during rn_rdt seconds. Only used if ln_bt_nn_auto=F
   rn_bt_cmax    =    0.7              !  Maximum courant number allowed if ln_bt_nn_auto=T
   nn_bt_flt     =    1                !  Time filter choice
                                       !  = 0 None
                                       !  = 1 Boxcar over   nn_baro barotropic steps
                                       !  = 2 Boxcar over 2*nn_baro     "        "
/
!-----------------------------------------------------------------------
&namcrs        !   Grid coarsening for dynamics output and/or
               !   passive tracer coarsened online simulations
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtsd    !   data : Temperature  & Salinity
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!          !  file name                   ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!          !                              !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_tem_ini  = 'woa09_temperature01_monthly_1deg_t_an_CMA_drowned_Ex_L75' ,     -1    ,  't_an'   ,    .false.   , .true. , 'yearly'   , 'reshape_WOA09_REG1toCREG025_bilin.nc'  ,   ''    ,    ''
   sn_sal_ini  = 'woa09_salinity01_monthly_1deg_s_an_CMA_drowned_Ex_L75'    ,     -1    ,  's_an'   ,    .false.   , .true. , 'yearly'   , 'reshape_WOA09_REG1toCREG025_bilin.nc'  ,   ''    ,    ''
   !
   ! data used for damping ( tradmp)
   sn_tem_dmp  = 'woa09_temperature_monthly_1deg_t_an_CMA_drowned_Ex_L75' ,    -12     ,  't_an'   ,    .true.   , .true. , 'yearly'   , 'reshape_WOA09_REG1toCREG025_bilin.nc'  ,   ''    ,    ''
   sn_sal_dmp  = 'woa09_salinity_monthly_1deg_s_an_CMA_drowned_Ex_L75'    ,    -12     ,  's_an'   ,    .true.   , .true. , 'yearly'   , 'reshape_WOA09_REG1toCREG025_bilin.nc'  ,   ''    ,    ''
   !
   cn_dir        = './'     !  root directory for the location of the runoff files
   ln_tsd_init   = .<TSD_INIT>.   ! Initialisation of ocean T & S with T &S input data (T) or not (F)
   ln_tsd_tradmp = .true.   !  damping of ocean T & S toward T &S input data (T) or not (F)
/
!!======================================================================
!!            ***  Surface Boundary Condition namelists  ***
!!======================================================================
!!   namsbc          surface boundary condition
!!   namsbc_ana      analytical         formulation
!!   namsbc_flx      flux               formulation
!!   namsbc_clio     CLIO bulk formulae formulation
!!   namsbc_core     CORE bulk formulae formulation
!!   namsbc_mfs      MFS  bulk formulae formulation
!!   namsbc_cpl      CouPLed            formulation                     ("key_oasis3")
!!   namsbc_sas      StAndalone Surface module
!!   namtra_qsr      penetrative solar radiation
!!   namsbc_rnf      river runoffs
!!   namsbc_isf      ice shelf melting/freezing
!!   namsbc_apr      Atmospheric Pressure
!!   namsbc_ssr      sea surface restoring term (for T and/or S)
!!   namsbc_alb      albedo parameters
!!======================================================================
!
!-----------------------------------------------------------------------
&namsbc        !   Surface Boundary Condition (surface module)
!-----------------------------------------------------------------------
   nn_fsbc     = 1         !  frequency of surface boundary condition computation
                           !     (also = the frequency of sea-ice model call)
   ln_ana      = .false.   !  analytical formulation                    (T => fill namsbc_ana )
   ln_flx      = .false.   !  flux formulation                          (T => fill namsbc_flx )
   ln_blk_clio = .false.   !  CLIO bulk formulation                     (T => fill namsbc_clio)
   ln_blk_core = .false.   !  CORE bulk formulation                     (T => fill namsbc_core)
   ln_blk_mfs  = .false.   !  MFS bulk formulation                      (T => fill namsbc_mfs )
   ln_cpl      = .false.   !  atmosphere coupled   formulation          ( requires key_oasis3 )
   ln_mixcpl   = .false.   !  forced-coupled mixed formulation          ( requires key_oasis3 )
   nn_components = 1       !  configuration of the opa-sas OASIS coupling
                           !  =0 no opa-sas OASIS coupling: default single executable configuration
                           !  =1 opa-sas OASIS coupling: multi executable configuration, OPA component
                           !  =2 opa-sas OASIS coupling: multi executable configuration, SAS component
   ln_apr_dyn  = .false.   !  Patm gradient added in ocean & ice Eqs.   (T => fill namsbc_apr )
   nn_ice      = 0         !  =0 no ice boundary condition   ,
                           !  =1 use observed ice-cover      ,
                           !  =2 ice-model used                         ("key_lim3" or "key_lim2")
   nn_ice_embd = 1         !  =0 levitating ice (no mass exchange, concentration/dilution effect)
                           !  =1 levitating ice with mass and salt exchange but no presure effect
                           !  =2 embedded sea-ice (full salt and mass exchanges and pressure)
   ln_dm2dc    = .false.   !LOLO!  daily mean to diurnal cycle on short wave
   ln_rnf      = .true.    !  runoffs                                   (T   => fill namsbc_rnf)
   nn_isf      = 0         !  ice shelf melting/freezing                (/=0 => fill namsbc_isf)
                           !  0 =no isf                  1 = presence of ISF
                           !  2 = bg03 parametrisation   3 = rnf file for isf
                           !  4 = ISF fwf specified
                           !  option 1 and 4 need ln_isfcav = .true. (domzgr)
   ln_ssr      = .true.    !  Sea Surface Restoring on T and/or S       (T => fill namsbc_ssr)
   nn_fwb      = 0         !  FreshWater Budget: =0 unchecked
                           !     =1 global mean of e-p-r set to zero at each time step
                           !     =2 annual global mean of e-p-r set to zero
   ln_wave = .false.       !  Activate coupling with wave (either Stokes Drift or Drag coefficient, or both)  (T => fill namsbc_wave)
   ln_cdgw = .false.       !  Neutral drag coefficient read from wave model (T => fill namsbc_wave)
   ln_sdw  = .false.       !  Computation of 3D stokes drift                (T => fill namsbc_wave)
   nn_lsm  = 0             !  =0 land/sea mask for input fields is not applied (keep empty land/sea mask filename field) ,
                           !  =1:n number of iterations of land/sea mask application for input fields (fill land/sea mask filename field)
   nn_limflx = -1          !  LIM3 Multi-category heat flux formulation (use -1 if LIM3 is not used)
                           !  =-1  Use per-category fluxes, bypass redistributor, forced mode only, not yet implemented coupled
                           !  = 0  Average per-category fluxes (forced and coupled mode)
                           !  = 1  Average and redistribute per-category fluxes, forced mode only, not yet implemented coupled
                           !  = 2  Redistribute a single flux over categories (coupled mode only)
/
!-----------------------------------------------------------------------
&namsbc_ana    !   analytical surface boundary condition
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_flx    !   surface boundary condition : flux formulation
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_clio   !   namsbc_clio  CLIO bulk formulae
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_core   !   namsbc_core  CORE bulk formulae
!-----------------------------------------------------------------------
!! WE ARE OPA! WE receive surface fluxes we don't compute them...
/
!-----------------------------------------------------------------------
&namsbc_cpl    !   coupled ocean/atmosphere model                       ("key_oasis3")
!-----------------------------------------------------------------------
!! OPA - oasis - SAS mode doesn't need this!!!
   sn_rcv_rnf    =       'none'                ,    'no'    ,     ''      ,         ''          ,   ''
/
!-----------------------------------------------------------------------
&namsbc_sas    !   analytical surface boundary condition
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtra_qsr    !   penetrative solar radiation
!-----------------------------------------------------------------------
!              !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_chl      ='chlaseawifs_c1m-99-05_smooth_CREG_R025_vh20161117',        -1         , 'CHLA'    ,   .true.     , .true. , 'yearly'  , ''       , ''       , ''

   cn_dir      = './'      !  root directory for the location of the runoff files
   ln_traqsr   = .true.    !  Light penetration (T) or not (F)
   ln_qsr_rgb  = .false.   !  RGB (Red-Green-Blue) light penetration
   ln_qsr_2bd  = .true.    !  2 bands              light penetration
   ln_qsr_bio  = .false.   !  bio-model light penetration
   nn_chldta   =      1    !  RGB : 2D Chl data (=1), 3D Chl data (=2) or cst value (=0)
   rn_abs      =   0.56    !  RGB & 2 bands: fraction of light (rn_si1)
   rn_si0      =   0.35    !  RGB & 2 bands: shortess depth of extinction
   rn_si1      =   23.0    !  2 bands: longest depth of extinction
   ln_qsr_ice  = .true.    !  light penetration for ice-model LIM3
/
!-----------------------------------------------------------------------
&namsbc_rnf    !   runoffs namelist surface boundary condition
!-----------------------------------------------------------------------
!              !  file name                 ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                            !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_rnf      = 'CREG025_runoff_monthly_combined_Dai_Trenberth_Bamber'     ,        -1         , 'runoff',     .true.     , .false. , 'yearly'  , ''       , ''       , ''
   sn_cnf      = 'CREG025_runoff_monthly_combined_Dai_Trenberth_Bamber'     ,         0         , 'socoefr' ,   .false.    , .false. , 'yearly'  , ''       , ''       , ''
   sn_s_rnf    = 'runoffs'                  ,        24         , 'rosaline',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
   sn_t_rnf    = 'runoffs'                  ,        24         , 'rotemper',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
   sn_dep_rnf  = 'runoffs'                  ,         0         , 'rodepth' ,   .false.    , .true. , 'yearly'  , ''       , ''       , ''

   cn_dir       = './RNF/'  !  root directory for the location of the runoff files
   ln_rnf_mouth = .true.    !  specific treatment at rivers mouths
   rn_hrnf      =  50.e0    !  depth over which enhanced vertical mixing is used
   rn_avt_rnf   =   1.e-3   !  value of the additional vertical mixing coef. [m2/s]
   rn_rfact     =   1.e0    !  multiplicative factor for runoff
   ln_rnf_depth = .false.   !  read in depth information for runoff
   ln_rnf_tem   = .false.   !  read in temperature information for runoff
   ln_rnf_sal   = .false.   !  read in salinity information for runoff
   ln_rnf_depth_ini = .false.  ! compute depth at initialisation from runoff file
   rn_rnf_max   = 5.735e-4  !  max value of the runoff climatologie over global domain ( ln_rnf_depth_ini = .true )
   rn_dep_max   = 150.      !  depth over which runoffs is spread ( ln_rnf_depth_ini = .true )
   nn_rnf_depth_file = 0    !  create (=1) a runoff depth file or not (=0)
/
!-----------------------------------------------------------------------
&namsbc_isf    !  Top boundary layer (ISF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_apr    !   Atmospheric pressure used as ocean forcing or in bulk
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_ssr    !   surface boundary condition : sea surface restoring
!-----------------------------------------------------------------------
!              !  file name  ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !             !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_sst      = 'woa09_sst01-12_monthly_1deg_t_an_CMA_drowned_Ex_L75' ,   -12   ,  't_an'   ,  .true.   , .true. , 'yearly'   , 'SST_reshape_WOA09_REG1toCREG025_bilin.nc'       , ''       , ''
   sn_sss      = 'woa09_sss01-12_monthly_1deg_s_an_CMA_drowned_Ex_L75' ,   -12   ,  's_an'   ,  .true.    , .true. , 'yearly'  , 'SSS_reshape_WOA09_REG1toCREG025_bilin.nc'       , ''       , ''
   sn_coast    = 'dist_coast_CREG025_ct20190401'           ,         0    ,  'Tcoast'    ,  .false.   , .true. , 'yearly'  ,  ''      , ''

   cn_dir      = './'      !  root directory for the location of the runoff files
   nn_sstr     =     0     !  add a retroaction term in the surface heat       flux (=1) or not (=0)
   nn_sssr     =     2     !  add a damping     term in the surface freshwater flux (=2)
                           !  or to SSS only (=1) or no damping term (=0)
   rn_dqdt     =   -40.    !  magnitude of the retroaction on temperature   [W/m2/K]
   rn_deds     =  -166.67  !  magnitude of the damping on salinity   [mm/day]
   ln_sssr_bnd =   .true.  !  flag to bound erp term (associated with nn_sssr=2)
   rn_sssr_bnd =   4.e0    !  ABS(Max/Min) value of the damping erp term [mm/day]
   ln_sssr_flt = .true.    ! use filtering of SSS model for sss restoring
   nn_shap_iter   =  100   ! number of iteration of the shapiro filter 
   ln_sssr_msk = .true.    ! use a mask near the coast
   rn_dist    =  150.      ! distance to the coast
   ln_sssr_ice = .false.   ! Apply SSS restoring under sea-ice
/
!-----------------------------------------------------------------------
&namsbc_alb    !   albedo parameters
!-----------------------------------------------------------------------
   nn_ice_alb  =    0   !  parameterization of ice/snow albedo
                        !     0: Shine & Henderson-Sellers (JGR 1985)
                        !     1: "home made" based on Brandt et al. (J. Climate 2005)
                        !                         and Grenfell & Perovich (JGR 2004)
   rn_albice   =  0.58  !  albedo of bare puddled ice (values from 0.49 to 0.58)
                        !     0.53 (default) => if nn_ice_alb=0
                        !     0.50 (default) => if nn_ice_alb=1
/
!-----------------------------------------------------------------------
&namberg       !   iceberg parameters
!-----------------------------------------------------------------------
      ln_icebergs              = .false.
/

!!======================================================================
!!               ***  Lateral boundary condition  ***
!!======================================================================
!!   namlbc        lateral momentum boundary condition
!!   namcla        cross land advection
!!   namagrif      agrif nested grid ( read by child model only )       ("key_agrif")
!!   nambdy        Unstructured open boundaries                         ("key_bdy")
!!   namtide       Tidal forcing at open boundaries                     ("key_bdy_tides")
!!======================================================================
!
!-----------------------------------------------------------------------
&namlbc        !   lateral momentum boundary condition
!-----------------------------------------------------------------------
   rn_shlat    =    2.     !  shlat = 0  !  0 < shlat < 2  !  shlat = 2  !  2 < shlat
                           !  free slip  !   partial slip  !   no slip   ! strong slip
   ln_vorlat   = .false.   !  consistency of vorticity boundary condition with analytical eqs.
/
!-----------------------------------------------------------------------
&namcla        !   cross land advection
!-----------------------------------------------------------------------
   nn_cla      =    0      !  advection between 2 ocean pts separates by land
/
!-----------------------------------------------------------------------
&namagrif      !  AGRIF zoom                                            ("key_agrif")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nam_tide      !   tide parameters (#ifdef key_tide)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nambdy        !  unstructured open boundaries                          ("key_bdy")
!-----------------------------------------------------------------------
    nb_bdy         = 2                    !  number of open boundary sets
    ln_coords_file = .false.,.false.      !  =T : read bdy coordinates from file
    cn_coords_file = '',''                !  bdy coordinates files
    ln_mask_file   = .false.              !  =T : read mask from file
    cn_mask_file   = ''                   !  name of mask file (if ln_mask_file=.TRUE.)
    cn_dyn2d       = 'flather','frs'      !
    nn_dyn2d_dta   = 1,1                  !  = 0, bdy data are equal to the initial state
                                          !  = 1, bdy data are read in 'bdydata   .nc' files
                                          !  = 2, use tidal harmonic forcing data from files
                                          !  = 3, use external data AND tidal harmonic forcing
    cn_dyn3d      =  'frs','frs'          !
    nn_dyn3d_dta  =  1,1                  !  = 0, bdy data are equal to the initial state
                                          !  = 1, bdy data are read in 'bdydata   .nc' files
    cn_tra        =  'frs','frs'          !
    nn_tra_dta    =  1,1                  !  = 0, bdy data are equal to the initial state
                                          !  = 1, bdy data are read in 'bdydata   .nc' files
    !! This is OPA not LIM !!!
    cn_ice_lim      =  'none','none'       !
    nn_ice_lim_dta  =  0,0                !  = 0, bdy data are equal to the initial state
                                          !  = 1, bdy data are read in 'bdydata   .nc' files
    rn_ice_tem      = 270.,270.           !  lim3 only: arbitrary temperature of incoming sea ice
    rn_ice_sal      =  10., 10.           !  lim3 only:      --   salinity           --
    rn_ice_age      =  30., 30.           !  lim3 only:      --   age                --

    ln_tra_dmp    =.false.,.false.        !  open boudaries conditions for tracers
    ln_dyn3d_dmp  =.false.,.false.        !  open boundary condition for baroclinic velocities
    rn_time_dmp   =  1.,1.                ! Damping time scale in days
    rn_time_dmp_out =  10.,10.            ! Outflow damping time scale
    nn_rimwidth   = 10,1                  !  width of the relaxation zone
    ln_vol        = .true.                !  total volume correction (see nn_volctl parameter)
    nn_volctl     = 1                     !  = 0, the total water flux across open boundaries is zero
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition     ("key_bdy")
!-----------------------------------------------------------------------
    ctypebdy ='S'                   ! Open boundary type (W,E,S or N)
    nbdyind  = 2                    ! indice of velocity row or column
                                    ! if ==-1, set obc at the domain boundary
                                    !        , discard start and end indices
    nbdybeg  = 54                   ! indice of segment start
    nbdyend  = 324                  ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition     ("key_bdy")
!-----------------------------------------------------------------------
    ctypebdy ='N'                   ! Open boundary type (W,E,S or N)
    nbdyind  = 601                  ! indice of velocity row or column
                                    ! if ==-1, set obc at the domain boundary
                                    !        , discard start and end indices
    nbdybeg  = 202                  ! indice of segment start
    nbdyend  = 223                  ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_dta      !  open boundaries - external data           ("key_bdy")
!-----------------------------------------------------------------------
!              !   file name    ! frequency (hours) !  variable  ! time interpol. !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                !  (if <0  months)  !    name    !    (logical)   !  (T/F)  ! 'monthly' ! filename ! pairing  ! filename      !
   bn_ssh =    'ORCA025.L75-GJM189-CREG025.L75-BDY_south_SSH' ,  -1   , 'sossheig' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_u3d  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_south_U'   ,  -1   , 'vozocrtx' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_v3d  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_south_V'   ,  -1   , 'vomecrty' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_tem  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_south_T'   ,  -1   , 'votemper' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_sal  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_south_S'   ,  -1   , 'vosaline' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   !
   cn_dir  =    './BDY/'
   ln_full_vel = .true.
/
!-----------------------------------------------------------------------
&nambdy_dta      !  open boundaries - external data           ("key_bdy")
!-----------------------------------------------------------------------
!              !   file name    ! frequency (hours) !  variable  ! time interpol. !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                !  (if <0  months)  !    name    !    (logical)   !  (T/F)  ! 'monthly' ! filename ! pairing  ! filename      !
   bn_ssh  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_north_SSH'    ,  -1   , 'sossheig' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_u3d  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_north_U'      ,  -1   , 'vozocrtx' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_v3d  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_north_V'      ,  -1   , 'vomecrty' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_tem  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_north_T'      ,  -1   , 'votemper' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   bn_sal  =   'ORCA025.L75-GJM189-CREG025.L75-BDY_north_S'      ,  -1   , 'vosaline' ,     .true.     , .false. ,  'yearly'  ,    ''    ,   ''     ,     ''
   !
   ! NO ICE! WE ARE OPA only!!! 
   cn_dir  =    './BDY/'
   ln_full_vel = .true.
/

!-----------------------------------------------------------------------
&nambdy_tide     ! tidal forcing at open boundaries
!-----------------------------------------------------------------------
/
!!======================================================================
!!                 ***  Bottom boundary condition  ***
!!======================================================================
!!   nambfr        bottom friction
!!   nambbc        bottom temperature boundary condition
!!   nambbl        bottom boundary layer scheme                         ("key_trabbl")
!!======================================================================
!
!-----------------------------------------------------------------------
&nambfr        !   bottom friction
!-----------------------------------------------------------------------
   nn_bfr      =    2      !  type of bottom friction :   = 0 : free slip,  = 1 : linear friction
                           !                              = 2 : nonlinear friction
   rn_bfri1    =    4.e-4  !  bottom drag coefficient (linear case)
   rn_bfri2    =    2.e-3  !  bottom drag coefficient (non linear case). Minimum coeft if ln_loglayer=T
   rn_bfri2_max =   1.e-1  !  max. bottom drag coefficient (non linear case and ln_loglayer=T)
   rn_bfeb2    =    2.5e-3 !  bottom turbulent kinetic energy background  (m2/s2)
   rn_bfrz0    =    3.e-3  !  bottom roughness [m] if ln_loglayer=T
   ln_bfr2d    = .false.   !  horizontal variation of the bottom friction coef (read a 2D mask file )
   rn_bfrien   =    50.    !  local multiplying factor of bfr (ln_bfr2d=T)
   rn_tfri1    =    4.e-4  !  top drag coefficient (linear case)
   rn_tfri2    =    2.5e-3 !  top drag coefficient (non linear case). Minimum coeft if ln_loglayer=T
   rn_tfri2_max =   1.e-1  !  max. top drag coefficient (non linear case and ln_loglayer=T)
   rn_tfeb2    =    0.0    !  top turbulent kinetic energy background  (m2/s2)
   rn_tfrz0    =    3.e-3  !  top roughness [m] if ln_loglayer=T
   ln_tfr2d    = .false.   !  horizontal variation of the top friction coef (read a 2D mask file )
   rn_tfrien   =    50.    !  local multiplying factor of tfr (ln_tfr2d=T)

   ln_bfrimp   = .true.    !  implicit bottom friction (requires ln_zdfexp = .false. if true)
   ln_loglayer = .true.   !  logarithmic formulation (non linear case)
/
!-----------------------------------------------------------------------
&nambbc        !   bottom temperature boundary condition
!-----------------------------------------------------------------------
!              !                              !  (if <0  months)  !  
!              !  file name      ! frequency (hours) ! variable   ! time interp.   !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !                 !  (if <0  months)  !   name     !   (logical)    !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
   sn_qgh      ='geothermal_heating.nc',  -12.  , 'heatflow'      ,   .false.      , .true.  , 'yearly'  , ''       , ''       , ''
   !
   cn_dir      = './'      !  root directory for the location of the runoff files
   ln_trabbc   = .false.    !  Apply a geothermal heating at the ocean bottom
   nn_geoflx   =    0      !  geothermal heat flux: = 0 no flux
                           !     = 1 constant flux
                           !     = 2 variable flux (read in geothermal_heating.nc in mW/m2)
   rn_geoflx_cst = 86.4e-3 !  Constant value of geothermal heat flux [W/m2]

/
!-----------------------------------------------------------------------
&nambbl        !   bottom boundary layer scheme
!-----------------------------------------------------------------------
   nn_bbl_ldf  =  1      !  diffusive bbl (=1)   or not (=0)
   nn_bbl_adv  =  1      !  advective bbl (=1/2) or not (=0)
   rn_ahtbbl   =  1000.  !  lateral mixing coefficient in the bbl  [m2/s]
   rn_gambbl   =  10.    !  advective bbl coefficient                 [s]
/

!!======================================================================
!!                        Tracer (T & S ) namelists
!!======================================================================
!!   nameos        equation of state
!!   namtra_adv    advection scheme
!!   namtra_adv_mle   mixed layer eddy param. (Fox-Kemper param.)
!!   namtra_ldf    lateral diffusion scheme
!!   namtra_dmp    T & S newtonian damping
!!======================================================================
!
!-----------------------------------------------------------------------
&nameos        !   ocean physical parameters
!-----------------------------------------------------------------------
   nn_eos      =  0     !  type of equation of state and Brunt-Vaisala frequency
                                 !  =-1, TEOS-10
                                 !  = 0, EOS-80
                                 !  = 1, S-EOS   (simplified eos)
   ln_useCT    = .false.  ! use of Conservative Temp. ==> surface CT converted in Pot. Temp. in sbcssm
/
!-----------------------------------------------------------------------
&namtra_adv    !   advection scheme for tracer
!-----------------------------------------------------------------------
   ln_traadv_cen2   =  .false.   !  2nd order centered scheme
   ln_traadv_tvd    =  .true.    !  TVD scheme
   ln_traadv_muscl  =  .false.   !  MUSCL scheme
   ln_traadv_muscl2 =  .false.   !  MUSCL2 scheme + cen2 at boundaries
   ln_traadv_ubs    =  .false.   !  UBS scheme
   ln_traadv_qck    =  .false.   !  QUICKEST scheme
   ln_traadv_msc_ups=  .false.   !  use upstream scheme within muscl
   ln_traadv_tvd_zts=  .false.  !  TVD scheme with sub-timestepping of vertical tracer advection
/
!-----------------------------------------------------------------------
&namtra_adv_mle !   mixed layer eddy parametrisation (Fox-Kemper param)
!-----------------------------------------------------------------------
   ln_mle    = .false.     ! (T) use the Mixed Layer Eddy (MLE) parameterisation
/
!----------------------------------------------------------------------------------
&namtra_ldf    !   lateral diffusion scheme for tracers
!----------------------------------------------------------------------------------
   !                       !  Operator type:
   ln_traldf_lap    =  .true.   !  laplacian operator
   ln_traldf_bilap  =  .false.  !  bilaplacian operator
   !                       !  Direction of action:
   ln_traldf_level  =  .false.  !  iso-level
   ln_traldf_hor    =  .false.  !  horizontal (geopotential)   (needs "key_ldfslp" when ln_sco=T)
   ln_traldf_iso    =  .true.   !  iso-neutral                 (needs "key_ldfslp")
   !		       	   !  Griffies parameters              (all need "key_ldfslp")
   ln_traldf_grif   =  .false.  !  use griffies triads
   ln_traldf_gdia   =  .false.  !  output griffies eddy velocities
   ln_triad_iso     =  .false.  !  pure lateral mixing in ML
   ln_botmix_grif   =  .false.  !  lateral mixing on bottom
   !                       !  Coefficients
   ! Eddy-induced (GM) advection always used with Griffies; otherwise needs "key_traldf_eiv"
   ! Value rn_aeiv_0 is ignored unless = 0 with Held-Larichev spatially varying aeiv
   !                                  (key_traldf_c2d & key_traldf_eiv & key_orca_r2, _r1 or _r05)
   rn_aeiv_0        =     0.    !  eddy induced velocity coefficient [m2/s]
   rn_aht_0         =   300.    !  horizontal eddy diffusivity for tracers [m2/s]
   rn_ahtb_0        =     0.    !  background eddy diffusivity for ldf_iso [m2/s]
   !                                           (normally=0; not used with Griffies)
   rn_slpmax        =     0.01  !  slope limit
   rn_chsmag        =     1.    !  multiplicative factor in Smagorinsky diffusivity
   rn_smsh          =     1.    !  Smagorinsky diffusivity: = 0 - use only sheer
   rn_aht_m         =  2000.    !  upper limit or stability criteria for lateral eddy diffusivity (m2/s)
/
!-----------------------------------------------------------------------
&namtra_dmp    !   tracer: T & S newtonian damping
!-----------------------------------------------------------------------
   ln_tradmp   =  .true.   !  add a damping termn (T) or not (F)
   nn_hdmp     =   -2      !  horizontal shape =-1, damping in Med and Red Seas only
                           !                   =XX, damping poleward of XX degrees (XX>0)
                           !                      + F(distance-to-coast) + Red and Med Seas
                           !                   =-2, DRAKKAR customization
   nn_zdmp     =    0      !  vertical   shape =0    damping throughout the water column
                           !                   =1 no damping in the mixing layer (kz  criteria)
                           !                   =2 no damping in the mixed  layer (rho crieria)
   nn_file     =    0      !  create a damping.coeff NetCDF file (=1) or not (=0)
   ln_dmpmask  = .false.   !  Read dmp_mask.nc file  when T (between 0 and 1 )
   rn_timsk    =  730.     !  Time scale used for dmp_mask
   cn_resto    = 'resto.nc' ! Name of file containing restoration coefficient field (use dmp_tools to create this)
/

!!======================================================================
!!                      ***  Dynamics namelists  ***
!!======================================================================
!!   namdyn_adv    formulation of the momentum advection
!!   namdyn_vor    advection scheme
!!   namdyn_hpg    hydrostatic pressure gradient
!!   namdyn_spg    surface pressure gradient                            (CPP key only)
!!   namdyn_ldf    lateral diffusion scheme
!!======================================================================
!
!-----------------------------------------------------------------------
&namdyn_adv    !   formulation of the momentum advection
!-----------------------------------------------------------------------
   ln_dynadv_vec = .true.  !  vector form (T) or flux form (F)
   nn_dynkeg     = 1       ! scheme for grad(KE): =0   C2  ;  =1   Hollingsworth correction
   ln_dynadv_cen2= .false. !  flux form - 2nd order centered scheme
   ln_dynadv_ubs = .false. !  flux form - 3rd order UBS      scheme
   ln_dynzad_zts = .false. !  Use (T) sub timestepping for vertical momentum advection
/
!-----------------------------------------------------------------------
&nam_vvl    !   vertical coordinate options
!-----------------------------------------------------------------------
 /
!-----------------------------------------------------------------------
&namdyn_vor    !   option of physics/algorithm (not control by CPP keys)
!-----------------------------------------------------------------------
   ln_dynvor_ene = .false. !  enstrophy conserving scheme
   ln_dynvor_ens = .false. !  energy conserving scheme
   ln_dynvor_mix = .false. !  mixed scheme
   ln_dynvor_een = .true.  !  energy & enstrophy scheme
   ln_dynvor_een_old = .false.  !  energy & enstrophy scheme - original formulation
/
!-----------------------------------------------------------------------
&namdyn_hpg    !   Hydrostatic pressure gradient option
!-----------------------------------------------------------------------
   ln_hpg_zco  = .false.   !  z-coordinate - full steps
   ln_hpg_zps  = .true.    !LOLO! NO VVL!!!  z-coordinate - partial steps (interpolation)
   ln_hpg_sco  = .false.   !  s-coordinate (standard jacobian formulation)
   ln_hpg_isf  = .false.   !  s-coordinate (sco ) adapted to isf
   ln_hpg_djc  = .false.   !  s-coordinate (Density Jacobian with Cubic polynomial)
   ln_hpg_prj  = .false.   !  s-coordinate (Pressure Jacobian scheme)
   ln_dynhpg_imp = .false. !  time stepping: semi-implicit time scheme  (T)
                                 !           centered      time scheme  (F)
/
!-----------------------------------------------------------------------
!namdyn_spg    !   surface pressure gradient   (CPP key only)
!-----------------------------------------------------------------------
!                          !  explicit free surface                     ("key_dynspg_exp")
!                          !  filtered free surface                     ("key_dynspg_flt")
!                          !  split-explicit free surface               ("key_dynspg_ts")

!-----------------------------------------------------------------------
&namdyn_ldf    !   lateral diffusion on momentum
!-----------------------------------------------------------------------
   !                       !  Type of the operator :
   ln_dynldf_lap    =  .false.  !  laplacian operator
   ln_dynldf_bilap  =  .true.   !  bilaplacian operator
   !                       !  Direction of action  :
   ln_dynldf_level  =  .false.  !  iso-level
   ln_dynldf_hor    =  .true.   !  horizontal (geopotential)            (require "key_ldfslp" in s-coord.)
   ln_dynldf_iso    =  .false.  !  iso-neutral                          (require "key_ldfslp")
   !                       !  Coefficient
   rn_ahm_0_lap     =     0.    !  horizontal laplacian eddy viscosity   [m2/s]
   rn_ahmb_0        =     0.    !  background eddy viscosity for ldf_iso [m2/s]
   rn_ahm_0_blp     =   -3.0e11 !  horizontal bilaplacian eddy viscosity [m4/s]
   rn_cmsmag_1      =     3.    !  constant in laplacian Smagorinsky viscosity
   rn_cmsmag_2      =     3     !  constant in bilaplacian Smagorinsky viscosity
   rn_cmsh          =     1.    !  1 or 0 , if 0 -use only shear for Smagorinsky viscosity
   rn_ahm_m_blp     =    -1.e12 !  upper limit for bilap  abs(ahm) < min( dx^4/128rdt, rn_ahm_m_blp)
   rn_ahm_m_lap     = 40000.    !  upper limit for lap  ahm < min(dx^2/16rdt, rn_ahm_m_lap)
/

!!======================================================================
!!             Tracers & Dynamics vertical physics namelists
!!======================================================================
!!    namzdf            vertical physics
!!    namzdf_ric        richardson number dependent vertical mixing     ("key_zdfric")
!!    namzdf_tke        TKE dependent vertical mixing                   ("key_zdftke")
!!    namzdf_kpp        KPP dependent vertical mixing                   ("key_zdfkpp")
!!    namzdf_ddm        double diffusive mixing parameterization        ("key_zdfddm")
!!    namzdf_tmx        tidal mixing parameterization                   ("key_zdftmx")
!!    namzdf_tmx_new    new tidal mixing parameterization               ("key_zdftmx_new")
!!======================================================================
!
!-----------------------------------------------------------------------
&namzdf        !   vertical physics
!-----------------------------------------------------------------------
   rn_avm0     =  1.4e-06  !  vertical eddy viscosity   [m2/s]          (background Kz if not "key_zdfcst")
   rn_avt0     =  1.0e-10  !  vertical eddy diffusivity [m2/s]          (background Kz if not "key_zdfcst")
   nn_avb      =    0      !  profile for background avt & avm (=1) or not (=0)
   nn_havtb    =    0      !  horizontal shape for avtb (=1) or not (=0)
   ln_zdfevd   = .true.    !  enhanced vertical diffusion (evd) (T) or not (F)
   nn_evdm     =    1      !  evd apply on tracer (=0) or on tracer and momentum (=1)
   rn_avevd    =  10.      !  evd mixing coefficient [m2/s]
   ln_zdfnpc   = .false.   !  Non-Penetrative Convective algorithm (T) or not (F)
   nn_npc      =    1            !  frequency of application of npc
   nn_npcp     =  365            !  npc control print frequency
   ln_zdfexp   = .false.   !  time-stepping: split-explicit (T) or implicit (F) time stepping
   nn_zdfexp   =    3            !  number of sub-timestep for ln_zdfexp=T
/
!-----------------------------------------------------------------------
&namzdf_ric    !   richardson number dependent vertical diffusion       ("key_zdfric" )
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_tke    !   turbulent eddy kinetic dependent vertical diffusion  ("key_zdftke")
!-----------------------------------------------------------------------
   rn_ediff    =   0.1     !  coef. for vertical eddy coef. (avt=rn_ediff*mxl*sqrt(e) )
   rn_ediss    =   0.7     !  coef. of the Kolmogoroff dissipation
   rn_ebb      =  67.83    !  coef. of the surface input of tke (=67.83 suggested when ln_mxl0=T)
   rn_ebbice   =   3.75    !  coef. of the surface input of tke under ice
   nn_havti    =    1      !  horizontal shape for avtb (=1) or not (=0) under ice
   rn_emin     =   1.e-10  !  minimum value of tke [m2/s2]
   rn_emin0    =   1.e-4   !  surface minimum value of tke [m2/s2]
   rn_bshear   =   1.e-20  ! background shear (>0) currently a numerical threshold (do not change it)
   nn_mxl      =   3       !  mixing length: = 0 bounded by the distance to surface and bottom
                           !                 = 1 bounded by the local vertical scale factor
                           !                 = 2 first vertical derivative of mixing length bounded by 1
                           !                 = 3 as =2 with distinct disspipative an mixing length scale
   nn_pdl      =   1       !  Prandtl number function of richarson number (=1, avt=pdl(Ri)*avm) or not (=0, avt=avm)
   ln_mxl0     = .true.    !  surface mixing length scale = F(wind stress) (T) or not (F)
   rn_mxl0     =   0.01    !  surface  buoyancy lenght scale minimum value
   ln_lc       = .true.    !  Langmuir cell parameterisation (Axell 2002)
   rn_lc       =   0.15    !  coef. associated to Langmuir cells
   nn_etau     =   1       !  penetration of tke below the mixed layer (ML) due to internal & intertial waves
                           !        = 0 no penetration
                           !        = 1 add a tke source below the ML
                           !        = 2 add a tke source just at the base of the ML
                           !        = 3 as = 1 applied on HF part of the stress    ("key_oasis3")
   rn_efr      =   0.05    !  fraction of surface tke value which penetrates below the ML (nn_etau=1 or 2)
   nn_htau     =   1       !  type of exponential decrease of tke penetration below the ML
                           !        = 0  constant 10 m length scale
                           !        = 1  0.5m at the equator to 30m poleward of 40 degrees
                           !        = 2  0.5m at the equator to 10m poleward of 40 degrees
/
!------------------------------------------------------------------------
&namzdf_kpp    !   K-Profile Parameterization dependent vertical mixing  ("key_zdfkpp", and optionally:
!------------------------------------------------------------------------ "key_kppcustom" or "key_kpplktb")
/
!-----------------------------------------------------------------------
&namzdf_gls                !   GLS vertical diffusion                   ("key_zdfgls")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_ddm    !   double diffusive mixing parameterization             ("key_zdfddm")
!-----------------------------------------------------------------------
   rn_avts     = 1.e-4     !  maximum avs (vertical mixing on salinity)
   rn_hsbfr    = 1.6       !  heat/salt buoyancy flux ratio
/
!-----------------------------------------------------------------------
&namzdf_tmx    !   tidal mixing parameterization                        ("key_zdftmx")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_tmx_new    !   new tidal mixing parameterization                ("key_zdftmx_new")
!-----------------------------------------------------------------------
!              ! file name ! frequency (hours) !   variable   ! time interp.   !  clim   ! 'yearly'/ ! weights  ! rotation ! land/sea mask !
!              !           !  (if <0  months)  !     name     !   (logical)    !  (T/F ) ! 'monthly' ! filename ! pairing  ! filename      !
      sn_tmxn_mpb =  'mixing_power_bot' , 0       , 'field'   ,  .false.        , .true.  , 'yearly'  , ''       , ''       , ''
      sn_tmxn_mpp =  'mixing_power_pyc' , 0       , 'field'   ,  .false.        , .true.  , 'yearly'  , ''       , ''       , ''
      sn_tmxn_mpc =  'mixing_power_cri' , 0       , 'field'   ,  .false.        , .true.  , 'yearly'  , ''       , ''       , ''
      sn_tmxn_dsb =  'decay_scale_bot'  , 0       , 'field'   ,  .false.        , .true.  , 'yearly'  , ''       , ''       , ''
      sn_tmxn_dsc =  'decay_scale_cri'  , 0       , 'field'   ,  .false.        , .true.  , 'yearly'  , ''       , ''       , ''
   nn_zpyc     = 2         !  pycnocline-intensified dissipation scales as N (=1) or N^2 (=2)
   ln_mevar    = .true.    !  variable (T) or constant (F) mixing efficiency
   ln_tsdiff   = .true.    !  account for differential T/S mixing (T) or not (F)
/
!!======================================================================
!!                  ***  Miscellaneous namelists  ***
!!======================================================================
!!   namsol            elliptic solver / island / free surface
!!   nammpp            Massively Parallel Processing                    ("key_mpp_mpi)
!!   namctl            Control prints & Benchmark
!!   namc1d            1D configuration options                         ("key_c1d")
!!   namc1d_uvd        data: U & V currents                             ("key_c1d")
!!   namc1d_dyndmp     U & V newtonian damping                          ("key_c1d")
!!   namsto            Stochastic parametrization of EOS
!!======================================================================
!
!-----------------------------------------------------------------------
&namsol        !   elliptic solver / island / free surface
!-----------------------------------------------------------------------
   nn_solv     =      1    !  elliptic solver: =1 preconditioned conjugate gradient (pcg)
                           !                   =2 successive-over-relaxation (sor)
   nn_sol_arp  =      0    !  absolute/relative (0/1) precision convergence test
   rn_eps      =  1.e-6    !  absolute precision of the solver
   nn_nmin     =    300    !  minimum of iterations for the SOR solver
   nn_nmax     =    800    !  maximum of iterations for the SOR solver
   nn_nmod     =     10    !  frequency of test for the SOR solver
   rn_resmax   =  1.e-10   !  absolute precision for the SOR solver
   rn_sor      =  1.92     !  optimal coefficient for SOR solver (to be adjusted with the domain)
/
!-----------------------------------------------------------------------
&nammpp        !   Massively Parallel Processing                        ("key_mpp_mpi)
!-----------------------------------------------------------------------
   cn_mpi_send =  'I'      !  mpi send/recieve type   ='S', 'B', or 'I' for standard send,
                           !  buffer blocking send or immediate non-blocking sends, resp.
   nn_buffer   =   0       !  size in bytes of exported buffer ('B' case), 0 no exportation
   ln_nnogather=  .false.  !  activate code to avoid mpi_allgather use at the northfold
   jpni        =   <JPNI>  !  jpni   number of processors following i (set automatically if < 1)
   jpnj        =   <JPNJ>  !  jpnj   number of processors following j (set automatically if < 1)
   jpnij       =   <JPNIJ> !  jpnij  number of local domains (set automatically if < 1)
/
!-----------------------------------------------------------------------
&namctl        !   Control prints & Benchmark
!-----------------------------------------------------------------------
   ln_ctl      = .false.   !  trends control print (expensive!)
   nn_print    =    0      !  level of print (0 no extra print)
   nn_ictls    =    0      !  start i indice of control sum (use to compare mono versus
   nn_ictle    =    0      !  end   i indice of control sum        multi processor runs
   nn_jctls    =    0      !  start j indice of control               over a subdomain)
   nn_jctle    =    0      !  end   j indice of control
   nn_isplt    =    1      !  number of processors in i-direction
   nn_jsplt    =    1      !  number of processors in j-direction
   nn_bench    =    0      !  Bench mode (1/0): CAUTION use zero except for bench
                           !     (no physical validity of the results)
   nn_timing   =    0      !  timing by routine activated (=1) creates timing.output file, or not (=0)
/
!-----------------------------------------------------------------------
&namc1d_uvd    !   data: U & V currents                                 ("key_c1d")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namc1d_dyndmp !   U & V newtonian damping                              ("key_c1d")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsto       ! Stochastic parametrization of EOS
!-----------------------------------------------------------------------
/

!!======================================================================
!!                  ***  Diagnostics namelists  ***
!!======================================================================
!!   namnc4       netcdf4 chunking and compression settings             ("key_netcdf4")
!!   namtrd       dynamics and/or tracer trends
!!   namptr       Poleward Transport Diagnostics
!!   namflo       float parameters                                      ("key_float")
!!   namhsb       Heat and salt budgets
!!======================================================================
!
!-----------------------------------------------------------------------
&namnc4        !   netcdf4 chunking and compression settings            ("key_netcdf4")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtrd        !   diagnostics on dynamics and/or tracer trends
!              !       and/or mixed-layer trends and/or barotropic vorticity
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namflo       !   float parameters                                      ("key_float")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namptr       !   Poleward Transport Diagnostic
!-----------------------------------------------------------------------
   ln_diaptr  = .false.    !  Poleward heat and salt transport (T) or not (F)
   ln_subbas  = .false.     !  Atlantic/Pacific/Indian basins computation (T) or not
/
!-----------------------------------------------------------------------
&namhsb       !  Heat and salt budgets
!-----------------------------------------------------------------------
   ln_diahsb  = .false.    !  check the heat and salt budgets (T) or not (F)
/
!-----------------------------------------------------------------------
&nam_diaharm   !   Harmonic analysis of tidal constituents ('key_diaharm')
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namdct        ! transports through sections
!-----------------------------------------------------------------------
/

!!======================================================================
!!            ***  Observation & Assimilation namelists ***
!!======================================================================
!!   namobs       observation and model comparison                      ('key_diaobs')
!!   nam_asminc   assimilation increments                               ('key_asminc')
!!======================================================================
!
!-----------------------------------------------------------------------
&namobs       !  observation usage switch                               ('key_diaobs')
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nam_asminc   !   assimilation increments                               ('key_asminc')
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_wave   ! External fields from wave model
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namdyn_nept  !   Neptune effect (simplified: lateral and vertical diffusions removed)
!-----------------------------------------------------------------------
/

