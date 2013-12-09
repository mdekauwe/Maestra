MODULE maestcom
    IMPLICIT NONE
    
    CHARACTER(LEN=3), PARAMETER :: format_ascii = 'asc' 
    CHARACTER(LEN=3), PARAMETER :: format_binary = 'bin'
    
    ! Maximum dimensions of arrays
    INTEGER, PARAMETER :: MAXT = 10000        ! Maximum no. of trees in plot
    INTEGER, PARAMETER :: MAXLAY = 15         ! Maximum no. of layers for radiation
    INTEGER, PARAMETER :: MAXSOILLAY = 100    ! Maximum no. of layers of soil (RAD)
    INTEGER, PARAMETER :: MAXSP = 10          ! Maximum no. of species (RAD).
    INTEGER, PARAMETER :: MAXP = 10000        ! Maximum no. of gridpoints
    INTEGER, PARAMETER :: MAXC = 3            ! Maximum no. of leaf area distributions
    INTEGER, PARAMETER :: MAXANG = 9          ! Maximum no. of zenith & leaf angles
    INTEGER, PARAMETER :: MAXD = 13           ! For resp prog
    INTEGER, PARAMETER :: MAXDATE = 400      ! Maximum no. of dates for tree or physiol parameters
    !INTEGER, PARAMETER :: maxdate = 5        ! Maximum no. of dates for physiol parameters
    !INTEGER, PARAMETER :: MAXMET = 18         ! Maximum columns in met data file
    ! changed by mgdk, just to try and compile, not 18 but MHET has 20?!!
    INTEGER, PARAMETER :: MAXMET = 20         ! Maximum columns in met data file
    INTEGER, PARAMETER :: MAXHISTO = 200      ! Maximum bins in PAR histogram
    REAL, PARAMETER    :: TOL = 0.02          ! Tolerance for leaf temp iteration
    INTEGER, PARAMETER :: MAXDAY = 900        ! For sumtrees program
    INTEGER, PARAMETER :: MAXHRS = 96         ! Maximum number of time periods in a day (ie 15 mins)
    INTEGER, PARAMETER :: MAXSTP = 10000      ! For ODEINT (Utils.for)
    INTEGER, PARAMETER :: NMAX = 4            ! "
    INTEGER, PARAMETER :: KMAXX = 200         ! "
    

    ! Values of physical constants
    REAL, PARAMETER :: TINY = 1.0E-30       ! "
    REAL, PARAMETER :: DEFWIND = 2.5       ! Default wind speed (m s-1)
    REAL, PARAMETER :: UMOLPERJ = 4.57     ! Conversion from J to umol quanta
    REAL, PARAMETER :: FPAR = 0.5          ! Fraction of global radiation that is PAR
    REAL, PARAMETER :: ABSZERO = -273.15   ! Absolute zero in degrees Celsius
    REAL, PARAMETER :: FREEZE = 273.15     ! Zero degrees Celsius in Kelvin
    REAL, PARAMETER :: TAU = 0.76          ! Transmissivity of atmosphere
    REAL, PARAMETER :: PI = 3.1415927      ! Pi
    REAL, PARAMETER :: TWOPI = 2.0 * PI    ! Two times Pi
    REAL, PARAMETER :: PID2 = PI / 2.0     ! Pi divided by two
    REAL, PARAMETER :: PID180 = PI / 180.0 ! Pi divided by 180 degrees
    REAL, PARAMETER :: AIRMA = 29.e-3      ! mol mass air (kg/mol)
    REAL, PARAMETER :: PATM = 1.0125E5     ! atmospheric pressure - standard condns (Pa)
    REAL, PARAMETER :: CPAIR = 1010.0      ! heat capacity of air (J kg-1 K-1)
    REAL, PARAMETER :: CPH2O = 4.186E06    ! heat capacity of water (J kg-1 K-1)
    REAL, PARAMETER :: CPQUARTZ = 1.942E06 ! heat capacity of quartz (J kg-1 K-1)
    REAL, PARAMETER :: TCQUARTZ = 7.7      ! thermal conductivity of quartz (W m-1 K-1)
    REAL, PARAMETER :: TCH2O = 0.594       ! thermal conductivity of water (W m-1 K-1)
    REAL, PARAMETER :: TCORG = 0.25        ! thermal conductivity of organic matter (W m-1 K-1)
    REAL, PARAMETER :: SOILALBEDO = 0.15   ! Albedo of soil, without snow.
    REAL, PARAMETER :: DHEAT = 21.5e-6     ! molecular diffusivity for heat
    REAL, PARAMETER :: EMLEAF = 0.95       ! Emissivity of thermal radiation by leaf
    REAL, PARAMETER :: EMSOIL = 0.95       ! Emissivity of thermal radiation by soil
    REAL, PARAMETER :: H2OLV0 = 2.501e6    ! latent heat H2O (J/kg)
    REAL, PARAMETER :: H2OMW = 18.e-3      ! mol mass H2O (kg/mol)
    REAL, PARAMETER :: H2OVW = 18.05e-6    ! partial molal volume of water at 20C (m3 mol-1)
    REAL, PARAMETER :: RCONST = 8.314      ! universal gas constant (J/mol/K)
    REAL, PARAMETER :: SIGMA = 5.67e-8     ! Steffan Boltzman constant (W/m2/K4)
    REAL, PARAMETER :: GBHGBC = 1.32       ! Ratio of Gbh:Gbc
    REAL, PARAMETER :: GSVGSC = 1.57       ! Ratio of Gsw:Gsc
    REAL, PARAMETER :: GBVGBH = 1.075      ! Ratio of Gbw:Gbh
    REAL, PARAMETER :: ALPHAQ = 0.425      ! Quantum yield of RuBP regen (mol mol-1)
    REAL, PARAMETER :: SOLARC = 1370       ! Solar constant (J m-2 s-1)
    REAL, PARAMETER :: GCPERMOL = 12.0     ! Grams ! per mol !
    REAL, PARAMETER :: CPERDW = 0.5        ! fraction per DW
    REAL, PARAMETER :: VONKARMAN = 0.41    ! von Karman's constant
    REAL, PARAMETER :: GRAV = 9.8067       ! Gravitational acceleration
      
    ! Numbers of I/O units
    INTEGER, PARAMETER :: UCONTROL = 1        ! Confile.dat
    INTEGER, PARAMETER :: UTREES = 2          ! Trees.dat
    INTEGER, PARAMETER :: USTR = 3            ! Str.dat
    INTEGER, PARAMETER :: UPHY = 4            ! Phy.dat
    INTEGER, PARAMETER :: STDIN = 5           ! Standard in
    INTEGER, PARAMETER :: STDOUT = 6          ! Standard out
    INTEGER, PARAMETER :: UTEST = 7           ! Points.dat
    INTEGER, PARAMETER :: UDAILY = 10         ! Dayflx.dat
    INTEGER, PARAMETER :: UHRLY = 11          ! Hrflx.dat
    INTEGER, PARAMETER :: UERROR = 12         ! Maeserr.dat
    INTEGER, PARAMETER :: ULAY = 13           ! Layflx.dat
    INTEGER, PARAMETER :: UPT = 14            ! Ptflx.dat
    INTEGER, PARAMETER :: UHIST = 15          ! Histo.dat
    INTEGER, PARAMETER :: UTDAY = 16          ! Daytot.dat
    INTEGER, PARAMETER :: UTHR = 17           ! Hrtot.dat
    INTEGER, PARAMETER :: UTCON = 18          ! Sumcon.dat
    INTEGER, PARAMETER :: UTHIST = 19         ! Tothist.dat
    INTEGER, PARAMETER :: URESP = 20          ! Resp.dat
    INTEGER, PARAMETER :: UPOINTSI = 21       ! Points.dat
    INTEGER, PARAMETER :: UPOINTSO = 22       ! Testflx.dat
    INTEGER, PARAMETER :: USTOREYI = 23       ! Ustorey.dat
    INTEGER, PARAMETER :: USTOREYO = 24       ! Usout.dat
    INTEGER, PARAMETER :: USPTO = 25          ! Uspts.dat
    INTEGER, PARAMETER :: URESPHR = 26        ! Resphr.dat
    INTEGER, PARAMETER :: UWATBAL = 27        ! watbal.dat
    INTEGER, PARAMETER :: UWATPARS = 28       ! watpars.dat
    INTEGER, PARAMETER :: UWATLAY = 29        ! watlay.dat
    INTEGER, PARAMETER :: UWATUPT = 32        ! watupt.dat
    INTEGER, PARAMETER :: UWATTEST = 31       ! wattest.dat
    INTEGER, PARAMETER :: USOILT = 30         ! soilt.dat
    INTEGER, PARAMETER :: UPARUS = 33         ! uspar.dat
    INTEGER, PARAMETER :: UWATDAY = 34        ! watbalday.dat
    INTEGER, PARAMETER :: UDAYHDR = 35        ! Dayflx_hdr.asc
    INTEGER, PARAMETER :: UMET = 36            ! Met.dat
    INTEGER, PARAMETER :: UTUTD = 37           ! Tutd.dat
    INTEGER, PARAMETER :: UHRLYHDR = 38 
    INTEGER, PARAMETER :: ULAYHDR = 39 
    INTEGER, PARAMETER :: UHISTHDR = 40
    INTEGER, PARAMETER :: URESPHRHDR = 41
    INTEGER, PARAMETER :: UWATBALHDR = 42
    INTEGER, PARAMETER :: UWATLAYHDR = 43
    INTEGER, PARAMETER :: USOILTHDR = 44
    INTEGER, PARAMETER :: UWATTESTHDR = 45
    INTEGER, PARAMETER :: UWATUPTHDR = 46
    INTEGER, PARAMETER :: UWATDAYHDR = 47
    INTEGER, PARAMETER :: URESPHDR = 48
    INTEGER, PARAMETER :: USUNLA = 49    

    ! Flags passed to error handling subroutine
    INTEGER, PARAMETER :: IFATAL = 100        ! Error was fatal - stop program
    INTEGER, PARAMETER :: IWARN = 200         ! Error non-fatal: just print warning

    ! Flags for crown shape
    INTEGER, PARAMETER :: JCONE = 1           ! Conical shape
    INTEGER, PARAMETER :: JHELIP = 2          ! Half-ellipsoid
    INTEGER, PARAMETER :: JPARA = 3           ! Paraboloidal shape
    INTEGER, PARAMETER :: JFELIP = 4          ! Full ellipsoid
    INTEGER, PARAMETER :: JCYL = 5            ! Cylindrical
    INTEGER, PARAMETER :: JBOX = 6            ! Box-shaped

    ! Flags to indicate which program it is
    INTEGER, PARAMETER :: INORMAL = 0         ! Maestra, Maeshr
    INTEGER, PARAMETER :: ITEST = 1           ! Maestest

      
    !COMMON /HRS/ HHRS, KHRS, SPERHR ! Make KHRS and HHRS available throughout the program
    INTEGER  :: KHRS != 24                  ! Number of time intervals in a day
    REAL :: HHRS != (KHRS) / 2.0           ! Half a day length
    REAL :: SPERHR != 3600.0 * 24.0 / KHRS ! Seconds in one time interval

END MODULE maestcom
