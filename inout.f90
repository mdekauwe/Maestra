!**********************************************************************
! INOUT.FOR
! This file contains all the subroutines for reading and writing to
! input and output files.
! The main subroutines (called externally) are:
! OPENINPUTF - opens the input files
! OPENOUTPUTF - opens the output files
! CLOSEF - closes the files
! INPUTSTR - reads the canopy structure file str.dat
! INPUTPHY - reads the physiology file phy.dat
! INPUTTREE - reads the trees.dat file
! INPUTCON - reads the control file confile.dat
! INPUTSOIL - reads the soil/understorey file soil.dat
! OUTPUTDY - outputs daily fluxes
! OUTPUTHR - outputs hourly fluxes
! OUTPUTLAY - outputs to layer flux file
! OUTPUTHIST - outputs PAR histogram
! SORTTREES - sorts the trees into order of distance from the target tree
! INTERPOLATEP - calls the daily interpolation routines for physiology
! INTERPOLATET - calls the daily interpolation routines for tree dimensions
! NB The subroutines to do with reading the met file are in getmet.for.

! Subsidiary subroutines are:
! INPUTSTR
!   READCROWN - read crown shape 
!   READAERO - read aerodynamic properties (e.g. extinction of wind) 
!   READBETA - read leaf area density distribution
!   READALLOM - read parameters for biomass vs height & diameter
!   READLIA - read leaf incidence angle 
! READLIA uses the following
!   ANGLE - calculate the fraction of leaf area in each angle class
!   AVGLIA - calculate the average LIA for a given ellipsoidal parameter
!   INTEG - integrate leaf angle dist over specified range
!   FANG - ellipsoidal leaf angle density function
!   CALCELP - calculate ellipsoidal parameter from average LIA
! INPUTCON
!   READDATES - read in start & end dates of simulation
!   READMODEL - read in options for which model to use
!   READZEN - read in number of angle classes & layers to use
!   READHIST - read in details of PAR histogram required
!   READCCSCEN - read in climate change scenario
!   READOTC - read in effects of OTC on met data
! INPUTPHY
!   READAGEP - read number of age classes for physiology parameters
!   READPROP - read proportions of leaf area in each age class
!   READABSRP - get leaf absorptance, reflectance & transmittance arrays
!   READGS - read parameters of stomatal model
!   READLEAFN - read in leaf nitrogen concentrations
!   READJMAX - read in Jmax and Vcmax values
!   READRD - read in parameters of leaf respiration
!   READPHYARRAY - read array of physiological parameters 
!     (called by READLEAFN, READJMAX, READRD)
!   READRW - read in parameters for woody respiration
!   READRR - read in parameters for root respiration
!   READRB - read in parameters for branch respiration
!   RINTEG - convert physiology parameters spec. in layers to layers used
!     (called by READBASRP, READPHYARRAY)
! INPUTTREE
!   READPLOT - read dimensions of plot
!   READXYZ - read co-ordinates of trees
!   READTREEARRAY - read arrays of tree dimensions eg height, diameter
!   GETLEAFAREA - A subroutine to read in leaf area array.
!   CALCLAI - calculate plot LAI
!   READZPD - read canopy dimensions (roughness length etc)
!   READCONTREES - get list of target trees
!   INEDGES - determine whether given tree is in plot edges or not
! INTERPOLATEP
!   PHYINTERP - do interpolation of physiological parameters
! INTERPOLATET
!   PHENOL - calculate leaf area from phenology parameters
!   TREEINTERP - do interpolation of tree dimensions
! Unused: NOTREEGET, GETTREENO, READRSOIL
!**********************************************************************


!**********************************************************************
      SUBROUTINE OPENINPUTF(IODAILYI,IOHRLYI,IOTUTDI,IOHISTI,IORESPI, &
       CTITLE,TTITLE,PTITLE,STITLE,IPOINTSI,ISUNLAI,ISIMUS)
! This routine opens the input files.
! The filenames are defined in this routine.
!**********************************************************************

    USE maestcom
    
    IMPLICIT NONE
    INTEGER LEN1,IOERROR,IWATFILE
    INTEGER IOHRLY,IOTUTD,IOHIST,IORESP,IPOINTS,IOWATBAL,ISUNLA
    INTEGER IODAILY,IOTUTDI,IOHISTI,IORESPI,IODAILYI
    INTEGER IPOINTSI,ISUNLAI,IOHRLYI,ISIMUS,IOFORMAT
    
    CHARACTER(LEN=80) CTITLE, TTITLE, PTITLE, STITLE
    LOGICAL EXT

      NAMELIST /CONTROL/ IOHRLY,IOTUTD,IOHIST,IORESP,IPOINTS,IOWATBAL, &
           ISUNLA,IOFORMAT,ISIMUS

990   FORMAT (A60)     ! For reading titles in input files.

! Output file for errors and warnings
      OPEN (UERROR, FILE = 'Maeserr.dat', STATUS = 'UNKNOWN')

! Input file with control switches
      OPEN (UCONTROL, FILE = 'confile.dat', STATUS = 'OLD',IOSTAT=IOERROR)
      
      IF(IOERROR.NE.0)THEN
            CALL SUBERROR('ERROR: CONFILE.DAT DOES NOT EXIST', IFATAL,0)
      ENDIF
      
! Default values of control switches
      IOHRLY = 0   ! Controls daily, hourly, and/or layer output
      IOTUTD = 1   ! Controls transmittance file output
      IOHIST = 0   ! Controls histogram output
      IORESP = 0   ! Controls respiration output
      IODAILY = 1  ! Controls daily output: FIXED HERE 
      IPOINTS = 0  ! No maestest, by default.
      IOWATBAL = 0 ! For compatibility with Maespa style confile.dat
      ISUNLA = 0   ! No sunlit leaf areas
      IOFORMAT = 0 ! NOT USED (maespa compatability)
      ISIMUS = 0   ! No understorey simulation, by default.

! Read control file
      READ (UCONTROL, 990) CTITLE
      READ (UCONTROL, CONTROL, IOSTAT = IOERROR)
      IF (IOERROR.NE.0)  &
       CALL SUBERROR('WARNING: USING DEFAULT VALUES FOR CONTROL FILE',IWARN,IOERROR)

      IOHRLYI = IOHRLY
      IOTUTDI = IOTUTD
      IOHISTI = IOHIST
      IORESPI = IORESP
      IODAILYI = IODAILY
      IPOINTSI = IPOINTS
      ISUNLAI = ISUNLA

! Input file with data on tree position and size
      OPEN (UTREES, FILE = 'trees.dat', STATUS='OLD',  &
           IOSTAT=IOERROR)
      
      IF(IOERROR.NE.0)  &
            CALL SUBERROR('ERROR: TREES.DAT DOES NOT EXIST' , &
            IFATAL,0)
      
! Multi-species version does not open these files.
!! Input file with data on canopy structure
!      OPEN (USTR, FILE = 'str.dat', STATUS='OLD',
!     &      IOSTAT=IOERROR)
!      
!! Input file with data on physiology
!      OPEN (UPHY, FILE = 'phy.dat', STATUS='OLD')

! Input/output file with diffuse transmittances
      OPEN (UTUTD, FILE = 'tutd.dat', STATUS='UNKNOWN')
      
! Input file for understorey parameters.
      IF(ISIMUS.EQ.1)THEN
          OPEN(USTOREYI, FILE='ustorey.dat', STATUS='UNKNOWN',  &
            IOSTAT=IOERROR)

          IF(IOERROR.NE.0)THEN
          CALL SUBERROR(  &
           'USTOREY.DAT NOT FOUND. NO UNDERSTOREY SIMULATED.', &
           IWARN,IOERROR)
          ENDIF
      ENDIF

! Read titles from input files
      READ (UTREES, 990) TTITLE
      
      STITLE = ' '
      PTITLE = ' '
      ! Need to fix this for multi-species option,
      !READ (USTR, 990) STITLE
      !READ (UPHY, 990) PTITLE

      RETURN
      END !OpenInputF


!**********************************************************************
      SUBROUTINE OPENOUTPUTF(IODAILY,IOHRLY,IOHIST,IORESP, &
       ISIMUS,ISUNLA,CTITLE,TTITLE,PTITLES,STITLES,MTITLE,VTITLE, &
       NSPECIES,SPECIESNAMES)
! This routine opens the output files.
! The filenames are defined in this routine.
! It also writes initial comments to the output files.
!**********************************************************************

      USE maestcom
      IMPLICIT NONE
      INTEGER LENSTR(MAXSP)
      INTEGER IOHRLY,IOTUTD,IOHIST,IORESP,IPOINTS,IOWATBAL
      INTEGER ISUNLA,IODAILY,NSPECIES,I,J,ISIMUS
      INTEGER LEN1, LEN2, LEN3
      
      CHARACTER*80 CTITLE,TTITLE,MTITLE,VTITLE,STITLE,PTITLE
      CHARACTER*180 PTITLEM,STITLEM
      CHARACTER SPECIESNAMES(MAXSP)*30
      CHARACTER STITLES(MAXSP)*60
      CHARACTER PTITLES(MAXSP)*60

      
990   FORMAT (A80)
992   FORMAT (A120)
991   FORMAT (A12,A60) ! For writing comments to output files.

! Output file with daily fluxes
      IF (IODAILY.GT.0) THEN
        OPEN (UDAILY,FILE = 'Dayflx.dat', STATUS='UNKNOWN')
      END IF

      OPEN (UWATTEST, FILE='wattest.dat', status='unknown')

! Output file with hourly fluxes (if required).
      IF (IOHRLY.GT.0) THEN
        OPEN (UHRLY, FILE = 'hrflux.dat', STATUS='UNKNOWN')
      ENDIF
! Output file with layer fluxes (if required).
      IF (IOHRLY.GT.1) THEN
        OPEN (ULAY, FILE = 'layflx.dat', STATUS='UNKNOWN')
      ENDIF

! Output file with histogram (if required).
      IF (IOHIST.EQ.1) THEN
        OPEN (UHIST, FILE = 'histo.dat', STATUS = 'UNKNOWN')
      END IF

! Output file for respiration (if required).
      IF (IORESP.EQ.1) THEN
        OPEN (URESP, FILE = 'resp.dat', STATUS = 'UNKNOWN')
        OPEN (URESPHR, FILE = 'resphr.dat', STATUS = 'UNKNOWN')
      END IF

! Output file for understorey simulation (if required).
      IF(ISIMUS.EQ.1) THEN
        OPEN (UPARUS, FILE = 'uspar.dat', STATUS = 'UNKNOWN')
      ENDIF

998   FORMAT(I3,1X,A30)

      !!!!!
      IF(ISUNLA.EQ.1)THEN
      OPEN (USUNLA, FILE= 'sunla.dat', STATUS = 'UNKNOWN')
      ENDIF
      
! Write headings to daily flux files.
      IF (IODAILY.GT.0) THEN
        WRITE (UDAILY, 991) 'Program:    ', VTITLE
        WRITE (UDAILY, 991) 'Control:    ', CTITLE
        WRITE (UDAILY, 991) 'Trees:      ', TTITLE
        IF(NSPECIES.EQ.1)WRITE (UDAILY, 992) 'Structure:  ', STITLES(1)
        IF(NSPECIES.EQ.1)WRITE (UDAILY, 992) 'Physiology: ', PTITLES(1)
        WRITE (UDAILY, 991) 'Met data:   ', MTITLE
        WRITE (UDAILY, 990) ' '
        IF(NSPECIES.GT.1)THEN
            WRITE(UDAILY, *)'Species codes:'
            DO I = 1,NSPECIES
              LEN1 = LEN_TRIM(SPECIESNAMES(I))
              LEN2 = LEN_TRIM(STITLES(I))
              LEN3 = LEN_TRIM(PTITLES(I))
        WRITE (UDAILY, 998) I, ADJUSTL(" : " // SPECIESNAMES(I)(:LEN1))
              WRITE (UDAILY, *) "     STR FILE :" // STITLES(I)(:LEN2)
              WRITE (UDAILY, *) "     PHY FILE :" // PTITLES(I)(:LEN3)
            ENDDO            
            WRITE(UDAILY,990) ' '
        ENDIF
        WRITE (UDAILY,501)
        WRITE (UDAILY,502)
        WRITE (UDAILY,503)
        WRITE (UDAILY,504)
        WRITE (UDAILY,505)
        WRITE (UDAILY,506)
        WRITE (UDAILY,507)
        WRITE (UDAILY,508)
        WRITE (UDAILY,509)
        WRITE (UDAILY,510)
        WRITE (UDAILY,511)
        WRITE (UDAILY,512)
        WRITE (UDAILY, 990) ' '
        WRITE (UDAILY,513)
      END IF

! Comments to hourly output file (if required).
      if (IOHRLY.gt.0) then
        WRITE (UHRLY, 991) 'Program:    ', VTITLE
        WRITE (UHRLY, 991) 'Control:    ', CTITLE
        WRITE (UHRLY, 991) 'Trees:      ', TTITLE
        IF(NSPECIES.EQ.1)WRITE (UHRLY, 992) 'Structure:  ', STITLES(1)
        IF(NSPECIES.EQ.1)WRITE (UHRLY, 992) 'Physiology: ', PTITLES(1)
        WRITE (UHRLY, 991) 'Met data:   ', MTITLE
        WRITE (UHRLY, 990) ' '
        IF(NSPECIES.GT.1)THEN
            WRITE(UHRLY, *)'Species codes:'
            DO I = 1,NSPECIES
              LEN1 = LEN_TRIM(SPECIESNAMES(I))
              LEN2 = LEN_TRIM(STITLES(I))
              LEN3 = LEN_TRIM(PTITLES(I))
              WRITE (UHRLY, *) I,  " : " // SPECIESNAMES(I)(:LEN1) 
              WRITE (UHRLY, *) "     STR FILE :" // STITLES(I)(:LEN2)
              WRITE (UHRLY, *) "     PHY FILE :" // PTITLES(I)(:LEN3)
            ENDDO
            WRITE(UHRLY,990) ' '
        ENDIF
        WRITE (UHRLY,701)
        WRITE (UHRLY,702)
        WRITE (UHRLY,703)
        WRITE (UHRLY,704)
        WRITE (UHRLY,705)
        WRITE (UHRLY,706)
        WRITE (UHRLY,707)
        WRITE (UHRLY,708)
        WRITE (UHRLY,709)
        WRITE (UHRLY,710)
        WRITE (UHRLY,711)
        WRITE (UHRLY,712)
        WRITE (UHRLY,713)
        WRITE (UHRLY,714)
        WRITE (UHRLY,715)
        WRITE (UHRLY,716)
        WRITE (UHRLY, 990) ' '
        WRITE (UHRLY,717)
      END IF

! Comments to respiration output file (if required).
      if (IORESP.gt.0) then
        WRITE (URESP, 991) 'Program:    ', VTITLE
        WRITE (URESP, 991) 'Control:    ', CTITLE
        WRITE (URESP, 991) 'Trees:      ', TTITLE
        WRITE (URESP, 991) 'Structure:  ', STITLE
        WRITE (URESP, 991) 'Physiology: ', PTITLE
        WRITE (URESP, 991) 'Met data:   ', MTITLE
        WRITE (URESP, 990) ' '
        WRITE (URESP,601)
        WRITE (URESP,602)
        WRITE (URESP,603)
        WRITE (URESP,604)
        WRITE (URESP,605)
        WRITE (URESP,606)
        WRITE (URESP,607)
        WRITE (URESP,608)
        WRITE (URESP,609)
        WRITE (URESP,610)
        WRITE (URESP,611)
        WRITE (URESP, 990) ' '
        WRITE (URESP,612)

      WRITE (URESPHR, 991) 'Program:    ', VTITLE
        WRITE (URESPHR, 991) 'Control:    ', CTITLE
        WRITE (URESPHR, 991) 'Trees:      ', TTITLE
        WRITE (URESPHR, 991) 'Structure:  ', STITLE
        WRITE (URESPHR, 991) 'Physiology: ', PTITLE
        WRITE (URESPHR, 991) 'Met data:   ', MTITLE
        WRITE (URESPHR, 990) ' '
        WRITE (URESPHR,301)
        WRITE (URESPHR,302)
        WRITE (URESPHR,303)
        WRITE (URESPHR,304)
        WRITE (URESPHR,305)
        WRITE (URESPHR,306)
        WRITE (URESPHR, 990) ' '
        WRITE (URESPHR,307)

      END IF

! Write comments to layer output file (if required)
      IF (IOHRLY.GT.1) THEN
        WRITE (ULAY, 991) 'Program:    ', VTITLE
        WRITE (ULAY, 801)
        WRITE (ULAY, 802)
        WRITE (ULAY, 803)
        WRITE (ULAY, 804)
        WRITE (ULAY, *)
      END IF

501   format('DOY: simulation date')
502   FORMAT('Tree: tree number')
503   FORMAT('Spec: tree species number')
504   format('absPAR:   absorbed PAR              MJ tree-1 d-1')
505   format('absNIR:   absorbed NIR              MJ  tree-1 d-1')
506   format('absTherm: absorbed thermal          MJ tree-1 d-1')
507   format('totPs: gross photosynthesis         mol tree-1 d-1')
508   format('totRf: daily foliar respiration     mol tree-1 d-1')
509   format('netPs: photosyn. net of foliar resp   mol tree-1 d-1')
510   format('totLE1: daily transpiration         mol H2O tree-1 d-1')
511   format('totLE2: daily transpirn (CANOPY calc) mol H2O m-2 d-1')
512   format('totH:  daily sensible heat flux     MJ tree-1 d-1')
513   format('Columns: DOY Tree Spec absPAR absNIR absTherm', &
        ' totPs totRf netPs', &
        ' totLE1 totLE2 totH')

701   format('DOY: simulation date')
702   FORMAT('Tree: tree number')
703   FORMAT('Spec: tree species number')
704   format('Hour:  hour of the day')
705   format('hrPAR: absorbed PAR              umol tree-1 s-1')
706   format('hrNIR: absorbed NIR              W tree-1')
707   format('hrTHM: absorbed thermal          W tree-1')
708   format('hrPS: photosynthesis (net of leaf resp) umol tree-1 s-1')
709   format('hrRf:  hourly leaf respiration   umol tree-1 s-1')
710   format('hrRmW: hourly stem + branch Rm   umol tree-1 s-1')
711   format('hrLE:  hourly transpiration      mmol tree-1 s-1')
712   format('LECAN: hourly transpirn: CANOPY calc : mmol H2O m-2 s-1')
713   format('Gscan: canopy stomatal conductance : mol CO2 tree-1 s-1')
714   format('hrH:   hourly sensible heat flux:  MJ tree-1 s-1')
715    FORMAT('TCAN: Average foliage temperature (deg !)')
716   FORMAT('DECO: Decoupling coefficient (0-1)')
717   format('Columns: DOY Tree Spec HOUR hrPAR hrNIR hrTHM', &
        ' hrPs hrRf hrRmW hrLE', &
        ' LECAN Gscan hrH TCAN DECO')

801   FORMAT(' Fluxes for each layer on an hourly basis')
802   FORMAT(' Rows: absorbed PAR (umol m-2 leaf s-1) ')
803   FORMAT('       photosynthesis net of Rleaf (umol m-2 leaf s-1) ')
804   FORMAT('       transpiration (umol m-2 leaf s-1) ')

601   FORMAT('Daily maintenance and growth respiration components')
602   FORMAT('Rmf: Foliage maintenance resp.    mol m-2 d-1')
603   FORMAT('Rmw: Stem maintenance resp.       mol m-2 d-1')
604   FORMAT('RmB: Branch maintenance resp.     mol m-2 d-1')
605   FORMAT('Rmcr: Coarse root maintenance resp. mol m-2 d-1')
606   FORMAT('Rmfr: Fine root maintenance resp. mol m-2 d-1')
607   FORMAT('Rgf: Foliage growth resp.         mol m-2 d-1')
608   FORMAT('Rgw: Stem growth resp.            mol m-2 d-1')
609   FORMAT('Rgb: Branch growth resp.          mol m-2 d-1')
610   FORMAT('Rgcr: Coarse root growth resp.    mol m-2 d-1')
611   FORMAT('Rgfr: Fine root growth resp.      mol m-2 d-1')
612   FORMAT('Columns: DOY, Tree, Rmf, Rmw, Rmb, Rmcr, Rmfr, Rgf, Rgw,  &
        Rgb, Rgcr, Rgfr')

301   FORMAT('Hourly maintenance respiration components')
302   FORMAT('Rmf: Foliage maintenance resp.    umol m-2 s-1')
303   FORMAT('Rmw: Stem maintenance resp.       umol m-2 s-1')
304   FORMAT('RmB: Branch maintenance resp.     umol m-2 s-1')
305   FORMAT('Rmcr: Coarse root maintenance resp. umol m-2 s-1')
306   FORMAT('Rmfr: Fine root maintenance resp. umol m-2 s-1')
307   FORMAT('Columns: DOY, Hr, Rmf, Rmw, Rmb, Rmcr, Rmfr')

      RETURN
      END !OpenOutputF


!**********************************************************************
      SUBROUTINE CLOSEF(IODAILY,IOHRLY,IOHIST,IORESP)
! This routine closes the open files.
!**********************************************************************

      USE maestcom
      INTEGER IODAILY,IOHRLY,IOHIST,IORESP

      CLOSE(UCONTROL)
      CLOSE(UTREES)
      CLOSE(USTR)
      CLOSE(UPHY)
      CLOSE(UMET)
      CLOSE(UTUTD)
      CLOSE(UERROR)
      CLOSE(UWATTEST)
      IF (IODAILY.GT.0) CLOSE(UDAILY)
      IF (IOHRLY.GT.0) CLOSE(UHRLY)
      IF (IOHRLY.GT.1) CLOSE(ULAY)
      IF (IOHIST.EQ.1) CLOSE(UHIST)
      IF (IORESP.EQ.1) CLOSE(URESP)
      IF (IORESP.EQ.1) CLOSE(URESPHR)

      RETURN
      END !Closef


!**********************************************************************
      SUBROUTINE INPUTCON( &
       ISTART, IEND, NSTEP, &
       NUMPNT,NOLAY,PPLAY,NZEN,DIFZEN,NAZ, &
       MODELGS, MODELJM, MODELRD, MODELSS, MODELRW, ITERMAX, &
       IOHIST, BINSIZE, &
       ICC, CO2INC, TINC, &
       IOTC, TOTC, WINDOTC, PAROTC, FBEAMOTC, &
       NSPECIES, SPECIESNAMES, PHYFILES, STRFILES)
! Read in the information from the control file. 
!**********************************************************************

    USE maestcom
    IMPLICIT NONE

    REAL DIFZEN(MAXANG)
    INTEGER PPLAY,ISTART,IEND,NSTEP,NUMPNT,NOLAY,PPLY,NZEN,NAZ,MODELGS
    INTEGER MODELJM,MODELRD,MODELSS,MODELRW,ITERMAX,NSPECIES
    INTEGER IOTC,IOHIST,IWATFILE,ICC
    CHARACTER SPECIESNAMES(MAXSP)*30
    CHARACTER PHYFILES(MAXSP)*30
    CHARACTER STRFILES(MAXSP)*30
    !CHARACTER STITLE*60
    REAL BINSIZE,CO2INC,TINC,TOTC,WINDOTC,PAROTC,FBEAMOTC

      CALL READDATES(UCONTROL, ISTART, IEND, NSTEP)
            
      CALL READSPECIES(UCONTROL, NSPECIES, SPECIESNAMES, PHYFILES, STRFILES)

      CALL READZEN(UCONTROL, NUMPNT, NOLAY, PPLAY, NZEN, NAZ, DIFZEN)

      CALL READMODEL(UCONTROL, MODELGS, MODELJM, MODELRD, MODELSS, MODELRW, ITERMAX)

      CALL READHIST(UCONTROL,IOHIST,BINSIZE)
        
      CALL READCCSCEN(UCONTROL,ICC,CO2INC,TINC)

      CALL READOTC(UCONTROL,IOTC,TOTC,WINDOTC,PAROTC,FBEAMOTC)

      RETURN
      END  ! InputCon


!**********************************************************************
      


!**********************************************************************
      SUBROUTINE INPUTSTR( &
       NSPECIES,STRFILES, &
       JLEAF,BPT,RANDOM,NOAGEC, &
       JSHAPE,SHAPEC,EXTWIND, &
       NALPHA,ALPHA,FALPHA, &
       COEFFT,EXPONT,WINTERC, &
       BCOEFFT,BEXPONT,BINTERC, &
       RCOEFFT,REXPONT,RINTERC,FRFRAC, &
       STITLES)
! This routine reads in input data on canopy structure (from USTR)
! which is required to begin the simulation.
! This is the new multi-species version (RAD, March 2009).
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER NSPECIES,I,IOERROR
    REAL ALPHA(MAXANG,MAXSP),FALPHA(MAXANG,MAXSP)
    REAL BPT(8,MAXC,MAXSP)
    REAL SHAPEC(MAXSP),EXTWIND(MAXSP)
    REAL RANDOM(MAXSP)
    REAL COEFFT(MAXSP),EXPONT(MAXSP),WINTERC(MAXSP)
    REAL BCOEFFT(MAXSP),BEXPONT(MAXSP),BINTERC(MAXSP)
    REAL RCOEFFT(MAXSP),REXPONT(MAXSP)
    REAL RINTERC(MAXSP),FRFRAC(MAXSP)
    INTEGER NOAGEC(MAXSP),JLEAF(MAXSP),JSHAPE(MAXSP)
    INTEGER NALPHA(MAXSP)
    CHARACTER STRFILES(MAXSP)*30
    CHARACTER STITLES(MAXSP)*60

900   FORMAT(A60)

      DO I=1,NSPECIES

          OPEN(USTR, FILE=STRFILES(I), STATUS='OLD', IOSTAT=IOERROR)
          IF(IOERROR.NE.0)THEN
            CALL SUBERROR('ERROR: STR FILE ' // TRIM(STRFILES(I)) //  &
                         ' DOES NOT EXIST' , &
            IFATAL,0)
          ENDIF

          READ(USTR,900)STITLES(I)

          CALL READCROWN(USTR, JSHAPE(I), SHAPEC(I))
          
          CALL READAERO(USTR, EXTWIND(I))

          CALL READLIA(USTR, NALPHA(I), ALPHA(1:MAXANG,I),  &
            FALPHA(1:MAXANG,I))

          CALL READBETA(USTR, NOAGEC(I), JLEAF(I),  &
            BPT(1:8,1:MAXC,I), RANDOM(I))

          CALL READALLOM(USTR, COEFFT(I), EXPONT(I), WINTERC(I), &
              BCOEFFT(I), BEXPONT(I), BINTERC(I),  &
              RCOEFFT(I), REXPONT(I), RINTERC(I), FRFRAC(I))
            
          CLOSE(USTR)

      ENDDO

      RETURN
      END !InputStr


!**********************************************************************
      SUBROUTINE INPUTPHY( &
       NSPECIES,PHYFILES, &
       MODELJM,MODELRD,MODELGS,MODELRW, &
       NOLAY,NOAGEC,NOAGEP,PROPC,PROPP, &
       ABSRP,REFLEC,TRANS,RHOSOL, &
       JMAXTABLE,DATESJ,NOJDATES,IECO,EAVJ,EDVJ,DELSJ,THETA, &
       VCMAXTABLE,DATESV,NOVDATES,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
       SLATABLE,DATESSLA,NOSLADATES,NOADATES,DATESA,AJQTABLE, &
       RDTABLE,DATESRD,NORDATES,RTEMP,DAYRESP,TBELOW, &
       EFFYRW,RMW,RTEMPW,COLLA,COLLK,STEMSDW,RMWAREA,STEMFORM, &
       NOFQDATES,DATESFQ,Q10FTABLE,NOWQDATES,DATESWQ,Q10WTABLE, &
       RMFR,RMCR,Q10R,RTEMPR,EFFYRF, &
       RMB,Q10B,RTEMPB, &
       GSREF,GSMIN,PAR0,D0,VK1,VK2,VPD1,VPD2,VMFD0, &
       GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,WC1, WC2, SWPEXP, &
       G0TABLE,D0L,GAMMA,G1TABLE,GK,WLEAF,NSIDES,PTITLES)
! This routine reads in input data on physiology (from UPHY).
! Some parameters can change with time: JMAX25, VCMAX25, RD0
! - for these, arrays of dates and values are read in, and must
! be interpolated.
!**********************************************************************

      USE maestcom
      IMPLICIT NONE
      INTEGER NSPECIES,IOERROR,I,NOLAY,MODELGS,MODELJM,MODELRD,MODELRW
      REAL ABSRP(MAXLAY,3,MAXSP),REFLEC(MAXLAY,3,MAXSP)
      REAL TRANS(MAXLAY,3,MAXSP),RHOSOL(3,MAXSP)
      REAL PROPP(MAXC,MAXSP),PROPC(MAXC,MAXSP)
      REAL LEAFN(MAXDATE,MAXLAY,MAXC,MAXSP)
      REAL JMAXTABLE(MAXDATE,MAXLAY,MAXC,MAXSP)
      REAL VCMAXTABLE(MAXDATE,MAXLAY,MAXC,MAXSP)
      REAL RDTABLE(MAXDATE,MAXLAY,MAXC,MAXSP)
      REAL SLATABLE(MAXDATE,MAXLAY,MAXC,MAXSP)
      REAL AJQTABLE(MAXDATE,MAXLAY,MAXC,MAXSP)
      REAL Q10FTABLE(MAXDATE,MAXSP),Q10WTABLE(MAXDATE,MAXSP)
      INTEGER DATESN(MAXDATE,MAXSP),DATESJ(MAXDATE,MAXSP)
      INTEGER DATESV(MAXDATE,MAXSP),DATESRD(MAXDATE,MAXSP)
      INTEGER DATESSLA(MAXDATE,MAXSP),DATESA(MAXDATE,MAXSP)
      INTEGER DATESFQ(MAXDATE,MAXSP),DATESWQ(MAXDATE,MAXSP)
      INTEGER NOAGEC(MAXSP),NOAGEP(MAXSP),NSIDES(MAXSP)
      REAL GSREF(MAXSP), GSMIN(MAXSP), PAR0(MAXSP)
      REAL D0(MAXSP), VK1(MAXSP), VK2(MAXSP)
      REAL VPD1(MAXSP), VPD2(MAXSP), VMFD0(MAXSP) 
      REAL GSJA(MAXSP), GSJB(MAXSP), T0(MAXSP)
      REAL TREF(MAXSP), TMAX(MAXSP), SMD1(MAXSP)
      REAL SMD2(MAXSP), WC1(MAXSP), WC2(MAXSP)
      REAL SWPEXP(MAXSP), D0L(MAXSP)
      !REAL G0(MAXSP),G1(MAXSP)
      REAL G0TABLE(MAXSP),G1TABLE(MAXSP)
      REAL GAMMA(MAXSP),  GK(MAXSP), WLEAF(MAXSP)
      INTEGER NOJDATES(MAXSP),NOVDATES(MAXSP)
      INTEGER NOADATES(MAXSP),NOSLADATES(MAXSP)
      INTEGER NONDATES(MAXSP),NORDATES(MAXSP)
      INTEGER NOWQDATES(MAXSP),NOFQDATES(MAXSP)
      INTEGER IECO(MAXSP)
      REAL EAVJ(MAXSP), EDVJ(MAXSP)
      REAL DELSJ(MAXSP), EAVC(MAXSP)
      REAL EDVC(MAXSP), DELSC(MAXSP), TVJUP(MAXSP)
      REAL TVJDN(MAXSP), THETA(MAXSP)
      REAL RTEMP(MAXSP), DAYRESP(MAXSP), EFFYRF(MAXSP)
      REAL TBELOW(MAXSP),EFFYRW(MAXSP),RMW(MAXSP)
      REAL RTEMPW(MAXSP),COLLA(MAXSP),COLLK(MAXSP)
      REAL STEMSDW(MAXSP),RMWAREA(MAXSP),STEMFORM(MAXSP)
      REAL Q10R(MAXSP),RTEMPR(MAXSP),Q10B(MAXSP),RTEMPB(MAXSP)
      REAL RMCR(MAXSP),RMFR(MAXSP),RMB(MAXSP)

      CHARACTER PHYFILES(MAXSP)*30
      CHARACTER PTITLES(MAXSP)*60

900   FORMAT(A60)

      DO I = 1,NSPECIES
      
          OPEN(UPHY, FILE=PHYFILES(I), STATUS='OLD',IOSTAT=IOERROR)
          IF(IOERROR.NE.0)THEN
            CALL SUBERROR('ERROR: PHY FILE ' // TRIM(PHYFILES(I)) //  &
                         ' DOES NOT EXIST' , &
            IFATAL,0)
          ENDIF
          
          READ (UPHY, 900) PTITLES(I)           
          !LENSTR = LEN_TRIM(PTITLES(I))
  
          CALL READRR(UPHY,RMCR(I),RMFR(I),Q10R(I),RTEMPR(I))
          
          CALL READRB(UPHY,RMB(I),Q10B(I),RTEMPB(I))

          ! Note that NOAGEC is passed to, not read from.
          CALL READAGEP(UPHY, NOAGEC(I), NOAGEP(I))

          CALL READPROP(UPHY, NOAGEC(I), NOAGEP(I), PROPC(1:MAXC,I),  &
            PROPP(1:MAXC,I))

          CALL READABSRP(UPHY,NOLAY,ABSRP(1:MAXLAY,1:3,I), &
            REFLEC(1:MAXLAY,1:3,I),TRANS(1:MAXLAY,1:3,I),RHOSOL(1:3,I))

          CALL READGS(UPHY,MODELGS,  &
           GSREF(I),GSMIN(I),PAR0(I),D0(I),VK1(I),VK2(I),VPD1(I), &
           VPD2(I), VMFD0(I), GSJA(I), GSJB(I), T0(I), TREF(I),  & 
           TMAX(I), SMD1(I), SMD2(I), SWPEXP(I), &
           G0TABLE(I), D0L(I), GAMMA(I), G1TABLE(I), GK(I), WLEAF(I), NSIDES(I))

          CALL READLEAFN(UPHY, MODELJM, MODELRD, NOLAY, NOAGEP(I), &
           NONDATES(I), DATESN(1:MAXDATE,I),  &
           LEAFN(1:MAXDATE,1:MAXLAY,1:MAXC,I))

          CALL READJMAX(UPHY, MODELJM, NOLAY, NOAGEP(I), NONDATES(I), &
           DATESN(1:MAXDATE,I), LEAFN(1:MAXDATE,1:MAXLAY,1:MAXC,I), &
           NOJDATES(I), DATESJ(1:MAXDATE,I),  &
             JMAXTABLE(1:MAXDATE,1:MAXLAY,1:MAXC,I), &
           NOVDATES(I), DATESV(1:MAXDATE,I),  &
             VCMAXTABLE(1:MAXDATE,1:MAXLAY,1:MAXC,I), &
           NOADATES(I), DATESA(1:MAXDATE,I),  &
             AJQTABLE(1:MAXDATE,1:MAXLAY,1:MAXC,I), &
           IECO(I), EAVJ(I), EDVJ(I), DELSJ(I), EAVC(I), EDVC(I), &
           DELSC(I), TVJUP(I), TVJDN(I), THETA(I))

          CALL READPHYARRAY(UPHY,5,NOLAY,NOAGEP(I), &
             NOSLADATES(I),DATESSLA(1:MAXDATE,I), &
             SLATABLE(1:MAXDATE,1:MAXLAY,1:MAXC,I))

          CALL READRD(UPHY, MODELRD, NOLAY, NOAGEP(I), NONDATES(I), &
           DATESN(1:MAXDATE,I), LEAFN(1:MAXDATE,1:MAXLAY,1:MAXC,I), & 
           NORDATES(I), DATESRD(1:MAXDATE,I),  &
             RDTABLE(1:MAXDATE,1:MAXLAY,1:MAXC,I),  &
           NOFQDATES(I), DATESFQ(1:MAXDATE,I), &
            Q10FTABLE(1:MAXDATE,I), &
           RTEMP(I), DAYRESP(I), EFFYRF(I), TBELOW(I))

          CALL READRW(UPHY,MODELRW,EFFYRW(I),RMW(I),RTEMPW(I),  &
           NOWQDATES(I), DATESWQ(1:MAXDATE,I),  &
           Q10WTABLE(1:MAXDATE,I), &
           COLLA(I),COLLK(I),STEMSDW(I),RMWAREA(I),STEMFORM(I))

          CALL READRR(UPHY,RMCR(I),RMFR(I),Q10R(I),RTEMPR(I))
          CALL READRB(UPHY,RMB(I),Q10B(I),RTEMPB(I))
            
          CLOSE(UPHY)  

      ENDDO

      RETURN
      END !InputPhy


!**********************************************************************
      SUBROUTINE INPUTTREE(  &
       XSLOPE,YSLOPE,BEAR,X0,Y0,XMAX,YMAX,STOCKING, &
       ZHT,Z0HT,ZPD, &
       NOALLTREES,NOTREES,NOTARGETS,ITARGETS,SHADEHT, &
       NOXDATES,NOYDATES,NOZDATES,NOTDATES,NOLADATES,NODDATES, &
       DATESX,DATESY,DATESZ,DATEST,DATESLA,DATESD, &
       DX,DY,DZ,R1,R2,R3,TRUNK,FLT,TOTLAITABLE,DIAMA, &
       IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN, &
       NSPECIES,ISPECIES)
! This subroutine should read in the data from the trees.dat file on
! crown positions and dimensions. Some variables can change with time:
! radii, height, diameter & leaf area - for these, arrays of dates & values at
! those dates may be read in, & interpolated during the program.
! x- and y- co-ordinates of the crowns may be read in or calculated from
! plot size and stocking density.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IOERROR
    INTEGER DATESX(maxdate),DATESY(maxdate),DATESZ(maxdate)
    INTEGER DATEST(maxdate),DATESLA(maxdate),DATESD(maxdate)
    INTEGER ITARGETS(MAXT),ISPECIES(MAXT),IPLOTSHAPE,NOALLTREES
    INTEGER NOXDATES,NOYDATES,NOZDATES,NOTDATES,IFLUSH,NOLADATES
    INTEGER NODDATES,NOTREES,NOTARGETS,NSPECIES
    REAL DX(MAXT),DY(MAXT),DZ(MAXT),WEIGHTS(MAXT), EXPFACTORS(MAXT)
    REAL R1(maxdate,MAXT),R2(maxdate,MAXT),R3(maxdate,MAXT)
    REAL TRUNK(maxdate,MAXT),FLT(maxdate,MAXT),TOTLAITABLE(maxdate)
    REAL DIAMA(maxdate,MAXT)
    REAL X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE,BEAR,SHADEHT,STOCKING
    REAL ZHT,Z0HT,ZPD,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN

! Read in number of trees & number of target tree
      CALL READPLOT(UTREES, X0, Y0, XMAX, YMAX, NOALLTREES, &
       XSLOPE, YSLOPE, BEAR, SHADEHT, STOCKING, IPLOTSHAPE)

! Read in aerodynamic properties of canopy
      CALL READZPD(UTREES,ZHT,Z0HT,ZPD)   

! Get x, y, z co-ords of each tree
      CALL READXYZ(UTREES,NOALLTREES,X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE, &
       DX,DY,DZ)
     
! Get radii in x & y directions of each tree
      CALL READTREEARRAY(UTREES,1,NOALLTREES,NOXDATES,DATESX,R1)
      CALL READTREEARRAY(UTREES,2,NOALLTREES,NOYDATES,DATESY,R2)
      
! If radii in y direction missing, set to x (4/5/2009, RAD).
      IF(R2(1,1).EQ.0.0)THEN
        R2=R1
      ENDIF
      
! Get green crown height of each tree
      CALL READTREEARRAY(UTREES,3,NOALLTREES,NOZDATES,DATESZ,R3)
! Get trunk length of each tree
      CALL READTREEARRAY(UTREES,4,NOALLTREES,NOTDATES,DATEST,TRUNK)

! Get leaf area parameters 
      CALL GETLEAFAREA(UTREES,IFLUSH,DT1,DT2,DT3,DT4, &
       EXPTIME,APP,EXPAN,NOALLTREES,NOLADATES,DATESLA,FLT)

! Get diameter of each tree
      CALL READTREEARRAY(UTREES,6,NOALLTREES,NODDATES,DATESD,DIAMA)

! Calculate total LAI
      CALL CALCLAI(NOLADATES,FLT,NOALLTREES,XMAX,YMAX,XSLOPE,YSLOPE, &
       TOTLAITABLE)

! Read in how many of the trees form the subplot (from confile)
      CALL READCONTREES(UCONTROL,NOALLTREES,DX,DY,XMAX,YMAX, &
       NOTREES,NOTARGETS,ITARGETS,IPLOTSHAPE)

! Read species array, if provided.
      CALL READSPECLIST(UTREES, NSPECIES, ISPECIES)

      RETURN
      END !InputTree
        


!**********************************************************************
      SUBROUTINE SORTTREES( &
       NOALLTREES,NOTREES,NOTARGET, &
       DX,DY,DZ,R1,R2,R3,TRUNK,FLT,DIAMA, &
       DXT,DYT,DZT,RX,RY,RZ,FOLT,ZBC,DIAM, &
       ISPECIES,ISPECIEST)
! This routine selects the 'NOTREES' trees closest to the target tree.
! It sorts the information about each tree into order of distance from
! the target tree. Does this simply by calling SORTTREESP.
! Modified March 2009 (RAD): 
! - Now outputs IT (index of sorted tree numbers), sorted species index.
!**********************************************************************

      
    USE maestcom
    IMPLICIT NONE
    INTEGER IT(MAXT),ISPECIES(MAXT),ISPECIEST(MAXT)
    INTEGER NOALLTREES,NOTREES,NOTARGET
    REAL DX(MAXT),DY(MAXT),DZ(MAXT)
    REAL R1(maxdate,MAXT),R2(maxdate,MAXT),R3(maxdate,MAXT)
    REAL TRUNK(maxdate,MAXT),FLT(maxdate,MAXT)
    REAL DXT(MAXT),DYT(MAXT),DZT(MAXT)
    REAL RX(maxdate,MAXT),RY(maxdate,MAXT),RZ(maxdate,MAXT)
    REAL ZBC(maxdate,MAXT),FOLT(maxdate,MAXT)
    REAL DIAMA(maxdate,MAXT),DIAM(maxdate,MAXT)
    REAL X,Y
    
    
      X = DX(NOTARGET)
      Y = DY(NOTARGET)
      CALL SORTTREESP( &
       X,Y,NOALLTREES,NOTREES, &
       DX,DY,DZ,R1,R2,R3,TRUNK,FLT,DIAMA, &
       DXT,DYT,DZT,RX,RY,RZ,FOLT,ZBC,DIAM, &
       ISPECIES,ISPECIEST)

      RETURN
      END !SortTrees


!**********************************************************************
      SUBROUTINE SORTTREESP(  &
       X,Y,NOALLTREES,NOTREES, &
       DX,DY,DZ,R1,R2,R3,TRUNK,FLT,DIAMA, &
       DXT,DYT,DZT,RX,RY,RZ,FOLT,ZBC,DIAM, &
       ISPECIES,ISPECIEST)
! This routine selects the 'NOTREES' trees closest to the point (x,y,z).
! It sorts the information about each tree into order of distance from
! this point.
! Modified March 2009 (RAD): 
! - Now outputs IT (index of sorted tree numbers), sorted species index.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER N,NOALLTREES,I,J,ITEMP,NOTREES,MM,IDATE
    INTEGER IT(MAXT),ISPECIES(MAXT),ISPECIEST(MAXT),MT1,MT2
    REAL DX(MAXT),DY(MAXT),DZ(MAXT)
    REAL R1(maxdate,MAXT),R2(maxdate,MAXT),R3(maxdate,MAXT)
    REAL TRUNK(maxdate,MAXT),FLT(maxdate,MAXT)
    REAL DXT(MAXT),DYT(MAXT),DZT(MAXT)
    REAL RX(maxdate,MAXT),RY(maxdate,MAXT),RZ(maxdate,MAXT)
    REAL ZBC(maxdate,MAXT),FOLT(maxdate,MAXT)
    REAL DIAMA(maxdate,MAXT),DIAM(maxdate,MAXT)
    REAL DNT(MAXT),X,Y,TEMP

! We select NT trees which are closest to the tree (NOTARGET) we
! are concerned with.
      DO 10 N = 1,NOALLTREES
        DNT(N) = SQRT((X-DX(N))**2 + (Y-DY(N))**2)
        IT(N) = N
10    CONTINUE

! Perform a sort on the trees to find the NOTREES closest ones to NOTARGET.
      MT1 = NOALLTREES - 1
      DO 20 I = 1,MT1
         MT2 = I + 1
         DO 20 J = MT2,NOALLTREES
            IF (DNT(I).LE.DNT(J)) GO TO 20
            TEMP = DNT(I)
            DNT(I) = DNT(J)
            DNT(J) = TEMP
            ITEMP = IT(I)
            IT(I) = IT(J)
            IT(J) = ITEMP
20    CONTINUE

! Produce new arrays containing the information about the closest trees.
      DO 30 I = 1,NOALLTREES  ! RAD Change, Feb 2011. Return all trees sorted.
        MM = IT(I)
        DXT(I) = DX(MM)
        DYT(I) = DY(MM)
        DZT(I) = DZ(MM)
        ISPECIEST(I) = ISPECIES(MM)
        DO 30 IDATE = 1,MAXDATE
          RX(IDATE,I) =  R1(IDATE,MM)
          RY(IDATE,I) = R2(IDATE,MM)
          RZ(IDATE,I) = R3(IDATE,MM)
          ZBC(IDATE,I) = TRUNK(IDATE,MM)
          FOLT(IDATE,I) = FLT(IDATE,MM)
          DIAM(IDATE,I) = DIAMA(IDATE,MM)
30    CONTINUE

      RETURN
      END !SortTrees


!**********************************************************************
      SUBROUTINE ANGLE(ELP,NALPHA,FALPHA)
! This routine calculates the fraction of leaf area in each leaf
! angle class.
!**********************************************************************
    USE maestcom
    IMPLICIT NONE
    INTEGER NALPHA,N
    REAL FALPHA(MAXANG)
    REAL ELP,ECTR,TEMP,TOTAL,DALP,ALP1,ALP2,RSUM

      DALP = PID2/FLOAT(NALPHA)
      IF (ELP.NE.1) THEN
         IF (ELP.LT.1) THEN
            ECTR = SQRT(1.0-ELP*ELP)
            TEMP = 1.0 + ASIN(ECTR)/ (ECTR*ELP)
         END IF

         IF (ELP.GT.1) THEN
            ECTR = SQRT(1.0- (1.0/ELP)**2)
            TEMP = 1.0 + LOG((1.0+ECTR)/ (1.0-ECTR))/ (2.0*ELP*ELP*ECTR)
         END IF

         TOTAL = 0.000
         DO 100 N = 1,NALPHA
            ALP1 = (N-1)*DALP
            ALP2 = N*DALP
            RSUM = 0.000
            CALL INTEG(ALP1,ALP2,ELP,TEMP,RSUM)
            FALPHA(N) = RSUM
            TOTAL = TOTAL + RSUM
  100    CONTINUE
         DO 110 N = 1,NALPHA
            FALPHA(N) = FALPHA(N)/TOTAL
  110    CONTINUE
      ELSE
         DO 200 N = 1,NALPHA
            FALPHA(N) = COS((N-1)*DALP) - COS(N*DALP)
  200    CONTINUE
      END IF

      RETURN
      END ! Angle

!**********************************************************************
      SUBROUTINE INTEG(ALP1,ALP2,ELP,TEMP,RSUM)
!     this is subroutine to integrate the leaf area density
!     function over a specified angle range
!**********************************************************************
    IMPLICIT NONE
    INTEGER I
    REAL H,UP,ALP1,ALP2,DOWN,RSUM,ELP,TEMP,FANG
    
      H = (ALP2-ALP1)/50.0
      UP = ALP1 + H*0.42265
      DOWN = ALP1 + H*1.57735
      RSUM = FANG(ELP,TEMP,UP) + FANG(ELP,TEMP,DOWN)
      DO 100 I = 1,24
         UP = UP + 2.0*H
         DOWN = DOWN + 2.0*H
  100 RSUM = RSUM + FANG(ELP,TEMP,UP) + FANG(ELP,TEMP,DOWN)
      RSUM = H*RSUM

      RETURN
      END !Integ


!**********************************************************************
      REAL FUNCTION AVGLIA(ELP)
!     this function is the relation between the parameter ELP and the
!     average leaf angle. (see notes for details)
!**********************************************************************
    IMPLICIT NONE
    REAL ELP
    
      IF (ELP.GT.1.0) THEN
         AVGLIA = 1.0/ (0.5901*ELP+0.3037)
      ELSE
         AVGLIA = 1.0/ (0.3782*ELP+0.6131)
      END IF

      RETURN
      END  !AVGLIA


!**********************************************************************
      REAL FUNCTION CALCELP(AVGANG)
!     this function is the relation between the average leaf angle and
!     the parameter ELP. This is the inverse of function AVGLIA.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL AVGANG

      AVGANG = AVGANG*PID180    !Convert to radians
      IF (AVGANG.LT.1.0) THEN
         CALCELP = (1./AVGANG - 0.3037)/0.5901
      ELSE
         CALCELP = (1./AVGANG - 0.6131)/0.3782
      END IF

      RETURN
      END  !CalcElp


!**********************************************************************
      REAL FUNCTION FANG(ELP,TEMP,ALP)
!     this is the ellipsoidal leaf angular density function
!     reference: Campbell,G.S. (personal comm.)
!**********************************************************************
    IMPLICIT NONE
    REAL ELP,TEMP,ALP
    
      FANG = 2.0*ELP*ELP*SIN(ALP)/ (TEMP*  &
            (((COS(ALP))**2.0+ (ELP*SIN(ALP))**2)**2))

      RETURN
      END  !Fang



!**********************************************************************
      SUBROUTINE OUTPUTHR(IOHRLY,IDAY,IHOUR,ITREE,ISPEC,TCAN, &
       NOLAY,PPAR,PPS,PTRANSP,FOLLAY, &
       THRAB,FCO2,FRESPF,FRESPW,FRESPB, &
       FH2OT,GSCAN,FH2OCAN,FHEAT,DECOUPL)
! Output the hourly totals
!**********************************************************************
    USE maestcom
   
    IMPLICIT NONE
    INTEGER NOTARGETS,IDAY,IHOUR,ITAR,ITREE,ISPEC,I,NOLAY
    INTEGER ITARGETS(MAXT),ISPECIES(MAXT),IOHRLY
    
    REAL THRAB(MAXHRS,3),FCO2(MAXHRS),FRESPF(MAXHRS)
    REAL FRESPW(MAXHRS),FRESPB(MAXHRS)
    REAL GSCAN(MAXHRS),FH2OT(MAXHRS),FH2OCAN(MAXHRS)
    REAL FHEAT(MAXHRS)
    REAL PPAR(MAXLAY,MAXHRS),PPS(MAXLAY,MAXHRS)
    REAL PTRANSP(MAXLAY,MAXHRS)
    REAL FOLLAY(MAXLAY),TCAN(MAXHRS),VPD(MAXHRS),DECOUPL(MAXHRS)
    
      IF (IOHRLY.GE.1) THEN
        WRITE (UHRLY,500) IDAY,ITREE,ISPEC,IHOUR, &
         THRAB(IHOUR,1)*UMOLPERJ,THRAB(IHOUR,2),THRAB(IHOUR,3), &
         FCO2(IHOUR),FRESPF(IHOUR),FRESPW(IHOUR)+FRESPB(IHOUR), &
         FH2OT(IHOUR)*1e-3, &
         FH2OCAN(IHOUR)*1E-3,GSCAN(IHOUR), &
         FHEAT(IHOUR)*1E-3,TCAN(IHOUR),DECOUPL(IHOUR)
      END IF
500   FORMAT (I7,1X,3(I4,1X),3(F12.5,1X),9(F11.5,1X))

      IF (IOHRLY.GE.2) THEN
        WRITE (ULAY,610) 'DAY',IDAY,'HOUR',IHOUR
        IF (FOLLAY(1).GT.0.0) THEN 
          WRITE (ULAY,600) (PPAR(I,IHOUR)/FOLLAY(I),I=1,NOLAY)
          WRITE (ULAY,600) (PPS(I,IHOUR)/FOLLAY(I),I=1,NOLAY)
          WRITE (ULAY,600) (PTRANSP(I,IHOUR)/FOLLAY(I),I=1,NOLAY)
        ELSE
          WRITE (ULAY,*) 'NO FOLIAGE AT THIS TIME'
        END IF
      END IF
600   FORMAT (10(F10.2,1X))
610   FORMAT (A5,I5,A5,I5)

      RETURN
      END !OutputHr


!**********************************************************************
      SUBROUTINE OUTPUTDY(IDAY,ITREE,ISPEC,IODAILY,TDYAB,TOTCO2, &
       TOTRESPF,TOTRESPW,TOTRESPWG,TOTH2O,TOTH2OCAN,TOTHFX, &
       IORESP,TOTRESPCR,TOTRESPFR,TOTRESPFRG,TOTRESPCRG,TOTRESPFG, &
       TOTRESPB,TOTRESPBG)
! Output the daily totals
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER NOTARGETS,IDAY,ITAR,ITREE,ISPEC
    INTEGER ITARGETS(MAXT),ISPECIES(MAXT)
    INTEGER IODAILY,IORESP
    
    REAL TDYAB(3)
    REAL TOTCO2,TOTRESPF,TOTRESPWM
    REAL TOTRESPB,TOTRESPFR,TOTRESPCR
    REAL TOTH2O,TOTH2OCAN,TOTHFX
    REAL TOTRESPW
    
    REAL TOTRESPFG,TOTRESPWG,TOTRESPBG
    REAL TOTRESPCRG,TOTRESPFRG
   
      IF (IODAILY.EQ.1) THEN
        WRITE (UDAILY,500) IDAY,ITREE,ISPEC,TDYAB,TOTCO2+TOTRESPF, &
         TOTRESPF,TOTCO2,TOTH2O,TOTH2OCAN,TOTHFX
500   FORMAT (I7,1X,I4,1X,I4,1X,9(F12.4,1X))
      END IF

      IF (IORESP.EQ.1) THEN
        WRITE (URESP,510) IDAY,ITREE,ISPEC,TOTRESPF,TOTRESPW,TOTRESPB, &
         TOTRESPCR,TOTRESPFR,TOTRESPFG,TOTRESPWG,TOTRESPBG, &
         TOTRESPCRG,TOTRESPFRG 
      END IF
510   FORMAT (I7,1X,I4,1X,I4,1X,10(F12.5,1X))

      RETURN
      END !OutputDy


!**********************************************************************
      SUBROUTINE READDATES(UFILE, ISTART, IEND, NSTEPI)
! The routine must return start and end dates for the simulation,
! in days-since-1950 format. The function IDATE50 converts a DD/MM/YY
! date to this format. The routine must also return NSTEP where the
! program is to use met data for every NSTEP'th day (default 1).
!**********************************************************************

    IMPLICIT NONE
    INTEGER UFILE,ISTART,IEND,NSTEPI,NSTEP,IOERROR,IWARN
    INTEGER, EXTERNAL :: IDATE50
    CHARACTER*10 STARTDATE,ENDDATE
    NAMELIST /DATES/ STARTDATE,ENDDATE,NSTEP

      STARTDATE = '01/01/50'
      ENDDATE = '01/01/50'
      NSTEP = 1
      REWIND (UFILE)
      READ (UFILE, DATES, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: MISSING DATES: USING ALL MET DATA', &
         IWARN,IOERROR)
        ISTART = 0
        IEND = 0
      ELSE
        ISTART = IDATE50(STARTDATE)
        IEND = IDATE50(ENDDATE)
      END IF
      NSTEPI = NSTEP

      RETURN
      END !ReadDates


!**********************************************************************
      SUBROUTINE READPLOT(UFILE, X0I, Y0I, XMAXI, YMAXI, NOALLTREESI, &
        XSLOPEI, YSLOPEI, BEARI, SHADEHTI, STOCKING, IPLOTSHAPE)
! Read in plot details. Subroutine must return:
! X0,Y0 - if plot is offset from origin. All tree co-ords referenced to this.
! XMAX, YMAX - dimensions of plot, in m
! NOALLTREES - no of trees in plot
! XSLOPE, YSLOPE - slope of plot, in radians
! SHADEHT - height of shadecloth surrounding plot, if any
! STOCKING - no of stems per ha
!**********************************************************************


    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR,IPLOTSHAPE,PLOTSHAPE,NOTREES,NOALLTREESI
    
    REAL X0I, Y0I, XMAXI, YMAXI, XSLOPEI, YSLOPEI, BEARI
    REAL SHADEHTI, STOCKING, X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE,BEARING,SHADEHT
    
    NAMELIST /PLOT/ X0,Y0,XMAX,YMAX,NOTREES,XSLOPE,YSLOPE,BEARING,SHADEHT,PLOTSHAPE
    
      SHADEHT = 0.0
      X0 = 0.0
      Y0 = 0.0
      PLOTSHAPE = 0

      REWIND (UFILE)
      READ (UFILE, PLOT, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) &
       CALL SUBERROR('ERROR READING PLOT DETAILS',IFATAL,IOERROR)
      X0I = X0
      Y0I = Y0
      XMAXI = XMAX
      YMAXI = YMAX
      XSLOPEI = XSLOPE*PID180
      YSLOPEI = YSLOPE*PID180
      BEARI = BEARING*PID180
      SHADEHTI = SHADEHT
      IPLOTSHAPE = PLOTSHAPE

      IF (NOTREES.GT.MAXT) THEN
        CALL SUBERROR( &
       'WARNING: NO OF TREES IN TREES FILE EXCEEDED MAXIMUM', &
       IWARN,IOERROR)
        NOTREES = MAXT
      END IF
      NOALLTREESI = NOTREES

      IF (IPLOTSHAPE.EQ.0) THEN
        STOCKING = NOTREES/(XMAXI*YMAXI)
      ELSE
        STOCKING = NOTREES/(PI*XMAXI*YMAXI)
      END IF    !NB Round plot shape not fully supported - XY co-ords must be supplied

      RETURN
      END !ReadPlot


!**********************************************************************
      SUBROUTINE READSPECLIST(UFILE,NSPECIES,ISPECIESI)

! Read species array from the trees file. 
! March 2009, RAD.      
!**********************************************************************
      
    USE maestcom
    IMPLICIT NONE
    INTEGER ISPECIES(MAXT),ISPECIESI(MAXT)
    INTEGER UFILE,NSPECIES,IOERROR

      NAMELIST /SPECLIST/ ISPECIES
                
      ISPECIES = 1
      
      REWIND(UFILE)
      READ (UFILE, SPECLIST, IOSTAT = IOERROR)
      
      IF(MAXVAL(ISPECIES).GT.NSPECIES)THEN
        CALL SUBERROR('WARNING: NSPECIES DOES NOT MATCH SPECLIST. &
       ALL SET TO SPECIES 1', IWARN,IOERROR)
        ISPECIES = 1
      ENDIF
      
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING : SPECIES LIST NOT READ, &
       ALL SET TO SPECIES 1', IWARN,IOERROR)
      END IF
      
      ISPECIESI = ISPECIES
      
      RETURN
      END


!**********************************************************************
      SUBROUTINE READZPD(UFILE,ZHTI,Z0HTI,ZPDI)
! Read in z, z0 and d for the canopy boundary layer conductance calcn.
!**********************************************************************


    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR
    REAL ZHT,ZHTI,Z0HT,Z0HTI,ZPD,ZPDI
    
     NAMELIST /AERODYN/ ZHT,Z0HT,ZPD
     
      Z0HT = 0.0
      REWIND (UFILE)
      READ (UFILE, AERODYN, IOSTAT = IOERROR)
      IF ((IOERROR.NE.0).OR.(Z0HT.EQ.0.0)) THEN
        CALL SUBERROR('WARNING: NOT CALCULATING CANOPY BDRY LAYER COND', IWARN,IOERROR)
      END IF
      ZHTI = ZHT
      Z0HTI = Z0HT
      ZPDI = ZPD

      RETURN
      END !ReadZPD
      

!**********************************************************************
      SUBROUTINE READCROWN(UFILE, JSHAPE, SHAPE)
! Read in crown shape parameters.
! JSHAPE indicates the shape of the crowns:
!   JCONE = conical,
!   JHELIP = half-ellipsoid,  (default)
!   JPARA = paraboloid,
!   JFELIP = full ellipsoid,
!   JCYL = cylinder,
!   JBOX = box.
! SHAPE is the factor needed to calculate volume for that crown shape:
!   VOLUME = PI * R2 * H * SHAPE (1/3, 2/3, 1/2, 2/3, 1, 4/PI).
!**********************************************************************

    USE maestcom
    IMPLICIT NONE

    INTEGER UFILE,IOERROR,JSHAPE
    REAL SHAPE
    CHARACTER*5 CSHAPE
    NAMELIST /CANOPY/ CSHAPE
    LOGICAL SHAPEREAD

      CSHAPE = 'ELIP' !Default value
      REWIND (UFILE)
      READ (UFILE, CANOPY, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: USING DEFAULT VALUE (ELIP) FOR CROWN SHAPE',IWARN,IOERROR)
      END IF

      IF (cshape.EQ.'CONE') THEN
         JSHAPE = JCONE
         SHAPE = 1.0/3.0

      ELSE IF (cshape.EQ.'ELIP') THEN
         JSHAPE = JHELIP
         SHAPE = 2.0/3.0

      ELSE IF (cshape.EQ.'PARA') THEN
         JSHAPE = JPARA
         SHAPE = 0.5

      ELSE IF (cshape.EQ.'ROUND') THEN
         JSHAPE = JFELIP
         SHAPE = 2.0/3.0

      ELSE IF (cshape.EQ.'CYL') THEN
         JSHAPE = JCYL
         SHAPE = 1.0

      ELSE IF (cshape.EQ.'BOX') THEN
         JSHAPE = JBOX
         SHAPE = 4.0/PI

      END IF

      RETURN
      END !ReadCrown


!**********************************************************************
      SUBROUTINE READAGEP(UFILE, NOAGEC, NOAGEPI)
! Age classes:
! NOAGEC is the number of age classes of foliage for which distributions
! of leaf area density are provided (read from str.dat).
! NOAGEP is the number of age classes
! of foliage for which (some) physiological parameters are provided.
! If neither NOAGEC nor NOAGEP = 1, then NOAGEC must = NOAGEP.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,NOAGEC,NOAGEPI,NOAGEP,IOERROR
    NAMELIST /NOAGES/ NOAGEP

      NOAGEP = 1  !Default value

      REWIND (UFILE)
      READ (UFILE, NOAGES, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: DEFAULT VALUE: NOAGEP=1',IWARN,IOERROR)
      END IF

! Check value does not exceed maximum.
      IF (NOAGEP.GT.MAXC) THEN
        CALL SUBERROR('WARNING: NOAGEP EXCEEDED MAXIMUM AGE CLASSES',IWARN,0)
        NOAGEP = MAXC
      END IF

! Check same as no of structural age classes
      IF ((NOAGEC.NE.1).AND.(NOAGEP.NE.1).AND.(NOAGEC.NE.NOAGEP)) THEN
        CALL SUBERROR('ERROR IN SPECIFICATION OF NO OF AGE CLASSES',IFATAL,IOERROR)
      END IF

      NOAGEPI = NOAGEP

      RETURN
      END !ReadAGEP


!**********************************************************************
      SUBROUTINE READAERO(UFILE, EXTWINDI)
! Read in aerodynamic info about canopy
!**********************************************************************
    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR
    REAL EXTWIND,EXTWINDI
    NAMELIST /AERO/ EXTWIND
   
      EXTWIND = 0.0
      REWIND (UFILE)
      READ (UFILE,AERO,IOSTAT = IOERROR)
!      IF (IOERROR.NE.0) THEN
!        CALL SUBERROR(
!     &  'WARNING: USING DEFAULT VALUES FOR LEAF INCIDENCE ANGLE',
!     &  IWARN,IOERROR)
!      END IF

      EXTWINDI = EXTWIND

      RETURN
      END !ReadAero
      

!**********************************************************************
      SUBROUTINE READLIA(UFILE, NALPHAI, ALPHA, FALPHAI)
! Read in leaf incidence angles.
! Must return: NALPHA = number of angle classes,
!              ALPHA = angle of each angle class,
!              FALPHA = fraction of leaf area in each angle class.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,NALPHA,IOERROR,IALP,NALPHAI,IANG
    REAL FALPHA(MAXANG),FALPHAI(MAXANG),ALPHA(MAXANG)
    REAL ELP,AVGANG,AVGLIA,DALPHA
    REAL, EXTERNAL :: CALCELP
    NAMELIST /LIA/ ELP,NALPHA,FALPHA,AVGANG

! Alternatives: (a) AVGANG > 0
! The mean leaf inclination angle is given by AVGANG. This is used to
! generate the LIA distribution assuming an elliptical distribution.
!               (b) ELP > 0.0
!    NALPHA = 1: there is just one leaf angle class, average angle = AVGLIA(ELP)
!    NALPHA > 1 (max 9): there are NALPHA leaf angle classes, the distribution
! of angles is elliptical with parameter ELP.
!               (!) ELP < 0.0
! The proportion of leaf area in each angle class is read in. Number of angle
! classes is given by NALPHA.
! The distribution of angles is read into array FALPHA(MAXANG).

! Default values
      AVGANG = -1.0
      NALPHA = 1
      ELP = 1.0
      FALPHA(1) = 1.0

! Read file
      REWIND (UFILE)
      READ (UFILE,LIA,IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: USING DEFAULT VALUES FOR LEAF INCIDENCE ANGLE',IWARN,IOERROR)
      END IF

      IF (NALPHA.EQ.1) THEN
        IF (AVGANG.GT.0.0) THEN
          ALPHA(1) = AVGANG*PID180
        ELSE IF (ELP.GT.0.0) THEN
          ALPHA(1) = AVGLIA(ELP)
        END IF
        FALPHA(1) = 1.0

      ELSE
        IF (AVGANG.GT.0.0) ELP = CALCELP(AVGANG)
        IF (ELP.GT.0.00) THEN
          DALPHA = PID2/FLOAT(NALPHA)
          DO 60 IALP = 1,NALPHA
             ALPHA(IALP) = (IALP-0.5)*DALPHA
60        CONTINUE
          CALL ANGLE(ELP,NALPHA,FALPHA)
        END IF
      END IF

      NALPHAI = NALPHA
      DO 70 IANG = 1,MAXANG
        FALPHAI(IANG) = FALPHA(IANG)
70    CONTINUE

      RETURN
      END !ReadLIA


!**********************************************************************
      SUBROUTINE READSPECIES(UFILE,NSPECIES,SPECIESNAMES, &
                             PHYFILES,STRFILES)

! Read species namelist (RAD, March 2009).      
!**********************************************************************
 
    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR,NSPECIES
    NAMELIST /SPECIES/ NSPECIES, SPECIESNAMES, PHYFILES, STRFILES

    CHARACTER SPECIESNAMES(MAXSP)*30
    CHARACTER PHYFILES(MAXSP)*30
    CHARACTER STRFILES(MAXSP)*30
 
! Set defaults, used if namelist does not exist (as in old confile.dat).
      NSPECIES = 1
      PHYFILES(1) = 'phy.dat'
      STRFILES(1) = 'str.dat'
      
! Read file
      REWIND (UFILE)
      READ (UFILE, SPECIES, IOSTAT = IOERROR)

      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('NO MULTI SPECIES OPTIONS READ IN CONFILE.DAT',IWARN,0)
      END IF
      
      RETURN
      END


!**********************************************************************
      SUBROUTINE READZEN(UFILE,NUMPNT,NOLAYI,PPLAYI,NZENI,NAZI,DIFZEN)
! Read in number of layers and angles to consider.
! Vertical layers:
!   NOLAY is the number of vertical layers used to calculate radiation
! interception (NUMPNT = NOLAY * points per layer).
! Angles for diffuse radiation:
!   NZEN = no. of zenith angles (default 5)
!   NAZ = no. of azimuth angles (default 11) - no maximum enforced
!   DIFZEN = array containing the NZEN zenith angles (radians).
!**********************************************************************


    USE maestcom
    IMPLICIT NONE
    NAMELIST /DIFFANG/ NOLAY,PPLAY,NZEN,NAZ
    INTEGER UFILE,PPLAY,PPLAYI,NOLAY,NZEN,NAZ
    INTEGER IOERROR,NOLD,I,NUMPNT,NOLAYI,NZENI,NAZI
    INTEGER, EXTERNAL :: DIVFOURNEAR
    REAL DIFZEN(MAXANG),DFZ

! Default values
      NOLAY = 6
      NAZ = 11
      NZEN = 5
      PPLAY = 12

! Read file
      REWIND (UFILE)
      READ (UFILE, DIFFANG, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: USING DEFAULT VALUES FOR NOLAY,PPLAY,NAZ & NZEN', &
        IWARN,IOERROR)
      END IF

! Make sure PPLAY can be divided by 4.
      NOLD = PPLAY
! PPLAY = DIVFOURNEAR(PPLAY)

!      IF(NOLD.NE.PPLAY)THEN
!        CALL SUBERROR(
!     &  'WARNING: PPLAY SET TO NEAREST SMALLER VALUE DIVISIBLE BY 4.',
!     &  IWARN,0)      
!      ENDIF

! Check values do not exceed maxima.
      IF (NOLAY.GT.MAXLAY) THEN
        CALL SUBERROR('WARNING: EXCEEDED MAXIMUM NO OF LAYERS',IWARN,0)
        NOLAY = MAXLAY
      END IF
      IF (NZEN.GT.MAXANG) THEN
        CALL SUBERROR('WARNING: EXCEEDED MAXIMUM NO OF ZENITH ANGLES',IWARN,0)
        NZEN = MAXANG
      END IF

! Set zenith angle array
      DFZ = PID2/FLOAT(NZEN)
      DO 90 I = 1,NZEN
         DIFZEN(I) = (I-0.5)*DFZ
   90 CONTINUE

      PPLAYI = PPLAY
      NUMPNT = NOLAY * PPLAY
      NOLAYI = NOLAY
      NZENI = NZEN
      NAZI = NAZ

      RETURN
      END !ReadZen


!**********************************************************************
      SUBROUTINE READALLOM(UFILE, COEFFTI, EXPONTI, WINTERCI, &
       BCOEFFTI, BEXPONTI, BINTERCI, &
       RCOEFFTI, REXPONTI, RINTERCI, FRFRACI)
! Read in coefficients for allometric relation between stem
! height, diameter, and stem mass.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR
    REAL COEFFT,WINTERC,EXPONT,RCOEFFT,RINTERC,REXPONT
    REAL FRFRAC,BCOEFFT,BEXPONT,BINTERC,COEFFTI,EXPONTI
    REAL WINTERCI,BCOEFFTI,BEXPONTI,BINTERCI,RCOEFFTI
    REAL REXPONTI,RINTERCI,FRFRACI
    NAMELIST /ALLOM/ COEFFT, EXPONT, WINTERC
    NAMELIST /ALLOMB/ BCOEFFT, BEXPONT, BINTERC
    NAMELIST /ALLOMR/ RCOEFFT, REXPONT, RINTERC, FRFRAC
    
! Default values
      COEFFT = 0.0
      WINTERC = 0.0
      EXPONT = 2.
      RCOEFFT = 0.0
      RINTERC = 0.0
      REXPONT = 2.
      FRFRAC = 1.0
      BCOEFFT = 0.0
      BEXPONT = 2.

! Read file
      REWIND (UFILE)
      READ (UFILE, ALLOM, IOSTAT = IOERROR)
      IF (COEFFT.GT.0.0) THEN
        CALL SUBERROR('CALCULATING WOODY RESPIRATION',IWARN,IOERROR)
      END IF

      REWIND (UFILE)
      READ (UFILE, ALLOMB, IOSTAT = IOERROR)
      IF (RCOEFFT.GT.0.0) THEN
        CALL SUBERROR('CALCULATING BRANCH RESPIRATION',IWARN,IOERROR)
      END IF

      REWIND (UFILE)
      READ (UFILE, ALLOMR, IOSTAT = IOERROR)
      IF (RCOEFFT.GT.0.0) THEN
        CALL SUBERROR('CALCULATING ROOT RESPIRATION',IWARN,IOERROR)
      END IF

      COEFFTI = COEFFT
      EXPONTI = EXPONT
      WINTERCI = WINTERC
      BCOEFFTI = BCOEFFT
      BEXPONTI = BEXPONT
      BINTERCI = BINTERC
      RCOEFFTI = RCOEFFT
      REXPONTI = REXPONT
      RINTERCI = RINTERC
      FRFRACI = FRFRAC

      RETURN
      END !ReadAllom


!**********************************************************************
      SUBROUTINE READPROP(UFILE, NOAGEC, NOAGEP, PROPCI, PROPPI)
! Read in the proportion of leaf area in each age class. The number
! of age classes is given by NOAGEC (for beta distributions) and
! NOAGEP (for physiological parameters).
! BM 7/03 Getting into troble with PROP when NOAGEC = 1 and NOAGEP > 1:
! Fix by having different proportions
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IOERROR,UFILE,I,NAGE,NOAGEP,NOAGEC,IAGE
    NAMELIST /PHENOL/ PROP
    REAL PROP(MAXC),PROPCI(MAXC),PROPPI(MAXC)
    REAL TOTAL

      DO 10 I = 1,MAXC
        PROP(I) = 0.0
10    CONTINUE

! Find number of age classes
      NAGE = NOAGEC
      IF (NOAGEP.GT.NAGE) NAGE = NOAGEP

! Only need to read in proportions if there is > 1 age class.
      IF (NAGE.GT.1) THEN
        REWIND (UFILE)
        READ (UFILE, PHENOL, IOSTAT = IOERROR)
        IF (IOERROR.NE.0) THEN
          CALL SUBERROR('ERROR: NEED DATA ON AGE CLASS PROPORTIONS', &
          IFATAL,IOERROR)
        END IF
      ELSE
        PROP(1) = 1.0
      END IF

! Check proportions sum to one
      TOTAL = 0.0
      DO 20 I = 1,NAGE
        TOTAL = TOTAL + PROP(I)
20    CONTINUE
      IF ((TOTAL.LT.0.9999).OR.(TOTAL.GT.1.0001)) &
       CALL SUBERROR('ERROR: AGE CLASS PROPORTIONS DO NOT SUM TO ONE', &
       IFATAL,0)

      IF (NOAGEC.EQ.1) THEN
      PROPCI(1) = 1.0
    ELSE
        DO 30 IAGE = 1,MAXC
          PROPCI(IAGE) = PROP(IAGE)
30      CONTINUE
      END IF

      IF (NOAGEP.EQ.1) THEN
      PROPPI(1) = 1.0
    ELSE
        DO 40 IAGE = 1,MAXC
          PROPPI(IAGE) = PROP(IAGE)
40      CONTINUE
      END IF

      RETURN
      END !ReadProp


!**********************************************************************
      SUBROUTINE READBETA(UFILE, NOAGECI, JLEAFI, BPTI, RANDOMI)
! Read in beta distributions for leaf area. The number of age classes
! for which beta distributions are specified is given by NOAGEC.
! Function returns:
! switch JLEAF:
!   JLEAF = 0: Uniform leaf area density assumed
!   JLEAF = 1: Leaf area density is variable in vertical direction
!   JLEAF = 2: Leaf area density is variable in both horizontal & vertical
! array BPT: gives the coefficients of the beta distributions
! the clumping factor RANDOM (= shoot projected area: leaf projected area).
!**********************************************************************


    USE maestcom
    IMPLICIT NONE
    INTEGER JLEAF,JLEAFI,UFILE,NOAGECI
    INTEGER NOAGEC,I,J,IOERROR
    REAL RANDOM,RANDOMI
    
    NAMELIST /LADD/ NOAGEC,JLEAF,BPT,BPTEXT,RANDOM
    REAL BPT(6,MAXC),BPTEXT(2,MAXC),BPTI(8,MAXC)
    
! BM 11/99 added another parameter to BETA function - input via BPTEXT - default value 1

! Default values
      NOAGEC = 1
      JLEAF = 0
      RANDOM = 1.0
      DO 10 J = 1,MAXC
      BPTEXT(1,J) = 1.0
      BPTEXT(2,J) = 1.0
        DO 10 I = 1,6
          BPT(I,J) = 0.0
10    CONTINUE

! Read file
      REWIND (UFILE)
      READ (UFILE, LADD, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR('WARNING: DEFAULT VALUE: UNIFORM LEAF AREA DENSITY', &
        IWARN,IOERROR)
      END IF
      NOAGECI = NOAGEC
      JLEAFI = JLEAF
      RANDOMI = RANDOM
      DO 20 J = 1,MAXC
      DO 20 I = 0,1
      BPTI(I*4+1,J) = BPT(I*3+1,J)
      BPTI(I*4+2,J) = BPT(I*3+2,J)
      BPTI(I*4+3,J) = BPT(I*3+3,J)
      BPTI(I*4+4,J) = BPTEXT(I+1,J)
20    CONTINUE

! Elementary error checking:
! If not using uniform leaf area density
      IF (JLEAF.NE.0) THEN
! If the first coefficient of the first age class = 0, or
        IF ((BPT(1,1).EQ.0.0).OR. &
! If there's >1 age class and the first coefft of the second age class = 0
            ((NOAGEC.GT.1).AND.(BPT(1,2).EQ.0.0))) THEN
          CALL SUBERROR( &
          'ERROR: MISSING BETA FUNCTION COEFFICIENTS', &
          IFATAL,IOERROR)
        END IF
      END IF

      RETURN
      END !ReadBeta


!**********************************************************************
      SUBROUTINE READMODEL(UFILE, GSMODI, JMMODI, RDMODI,  &
        SSMODI, RWMODI, ITERMAXI)
! Read in flags which control the physiological models used.
! MODELGS - controls model of stomatal conductance
! MODELJM - whether JMAX,VCMAX read in (0) or calculated from leaf N (1)
! MODELRD - whether RD is a function of leaf area (0) or leaf N (1)
! MODELSS - whether photosynthesis is calculated for sunlit & shaded leaves
!   together or separately
! ITERMAX - no. of iterations to be used in solving for leaf temperature
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,GSMODI,JMMODI,RDMODI,SSMODI,RWMODI,ITERMAXI
    INTEGER MODELGS,MODELJM,MODELRD,MODELRW,MODELSS
    INTEGER ITERMAX,IOERROR
    
    NAMELIST /MODEL/ MODELGS,MODELJM,MODELRD,MODELRW, &
                        MODELSS,ITERMAX
    
! Default values
      MODELGS = 0     ! The Jarvis model of stom cond
      MODELJM = 0     ! Jmax & Vcmax read in directly
      MODELRD = 0     ! Rd0 read in directly
      MODELRW = 0     ! RW values read in directly
      MODELSS = 0     ! sunlit & shade calculations separate
      ITERMAX = 0     ! The leaf temperature is not calculated
     
! Read file
      REWIND (UFILE)
      READ (UFILE, MODEL, IOSTAT = IOERROR)
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR( &
        'WARNING: DEFAULT VALUES USED FOR PHYSIOLOGICAL CONTROL', &
        IWARN,IOERROR)
      END IF
      GSMODI = MODELGS
      RDMODI = MODELRD
      RWMODI = MODELRW
      JMMODI = MODELJM
      SSMODI = MODELSS
      ITERMAXI = ITERMAX

      RETURN
      END !ReadModel


!**********************************************************************
      SUBROUTINE READHIST(UFILE, IOHIST, BINSIZEI)
! Read in information for the PAR histogram if required.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOHIST,IOERROR
    REAL BINSIZE,BINSIZEI
    NAMELIST /HISTO/ BINSIZE

      IF (IOHIST.EQ.1) THEN

        REWIND (UFILE)
        READ (UFILE, HISTO, IOSTAT = IOERROR)
        IF (IOERROR.NE.0) THEN
          CALL SUBERROR( &
          'ERROR: MISSING INFO FOR HISTOGRAM', &
          IFATAL,IOERROR)
        END IF
        BINSIZEI = BINSIZE

      END IF
      RETURN
      END ! ReadHist


!**********************************************************************
      SUBROUTINE READCCSCEN(UFILE, ICC, CO2INCI, TINCI)
! Read in details of climate change scenario.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,ICC
    REAL CO2INC,TINC,CO2INCI,TINCI
    NAMELIST /CCSCEN/ CO2INC, TINC

! Default values
      ICC = 0
      CO2INC = 0.0
      TINC = 0.0
! Read file
      REWIND (UFILE)
      READ (UFILE, CCSCEN, IOSTAT = ICC)
      
      IF (ICC.EQ.0)  &
        CALL SUBERROR('APPLYING CLIMATE CHANGE SCENARIO',IWARN,0)
      CO2INCI = CO2INC
      TINCI = TINC

      RETURN
      END ! ReadCCScen


!**********************************************************************
      SUBROUTINE READOTC(UFILE, IOTC, TOTCI, WINDOTCI, &
        PAROTCI, FBEAMOTCI)
! Subroutine to read in parameters describing effect of OTC on met data.
!**********************************************************************
    
    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOTC
    REAL TOTCI,WINDOTCI,PAROTCI,FBEAMOTCI
    REAL TOTC,WINDOTC,PAROTC,FBEAMOTC
    NAMELIST /OTC/ TOTC,WINDOTC,PAROTC,FBEAMOTC
    

! Default values
      IOTC = 0
      TOTC = 0
      WINDOTC = -1
      PAROTC = 1.0
      FBEAMOTC = 1.0
! Read file
      REWIND (UFILE)
      READ (UFILE, OTC, IOSTAT = IOTC)
      IF (IOTC.EQ.0)  &
        CALL SUBERROR('APPLYING EFFECTS OF OTC ON MET DATA',IWARN,0)
      TOTCI = TOTC
      WINDOTCI = WINDOTC
      PAROTCI = PAROTC
      FBEAMOTCI = FBEAMOTC

      RETURN
      END ! ReadOTC

!**********************************************************************
      SUBROUTINE READABSRP(UFILE,NOLAY,ABSRP,REFLEC,TRANS,RHOSOLI)
! Read leaf absorptances and reflectances for 3 wavelengths.
! Required input: File unit (UFILE), No of layers (NOLAY).
!**********************************************************************
    
    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR,NOLAY,NOLAYERS,IP,J,I
    REAL ABSRP(MAXLAY,3),REFLEC(MAXLAY,3),TRANS(MAXLAY,3)
    REAL ARHO(MAXLAY*MAXC*maxdate),ATAU(MAXLAY*MAXC*maxdate)
    REAL RHOSOL(3),RHOSOLI(3),X(MAXLAY+1)
    REAL D1,D2
    REAL, EXTERNAL :: RINTEG
    NAMELIST /ABSORP/ NOLAYERS,RHOSOL,ATAU,ARHO
    

! Read file
      REWIND (UFILE)
      READ (UFILE,ABSORP,IOSTAT = IOERROR)
      IF (IOERROR.NE.0) &
       CALL SUBERROR( &
       'ERROR READING REFLECTANCE / TRANSMITTANCE', &
       IFATAL,IOERROR)

      DO 10 IP = 1,NOLAYERS+1
        X(IP) = (REAL(IP) - 1.0)/REAL(NOLAYERS)
10    CONTINUE

! Check values & set absorptance array.
      DO 30 J = 1,3
        RHOSOLI(J) = RHOSOL(J)
        DO 20 I = 1,NOLAY
          D1 = (REAL(I) - 1)/REAL(NOLAY)
          D2 = (REAL(I))/REAL(NOLAY)
          REFLEC(I,J) = RINTEG(D1,D2,X,ARHO,(J-1)*NOLAYERS,NOLAY)
          TRANS(I,J) = RINTEG(D1,D2,X,ATAU,(J-1)*NOLAYERS,NOLAY)
          ABSRP(I,J) = 1.0 - REFLEC(I,J) - TRANS(I,J)
20    CONTINUE

30    CONTINUE

100   FORMAT(10(1X,F8.3))

      RETURN
      END !ReadAbsrp


!**********************************************************************
      SUBROUTINE READGS(UFILE,MODELGS, &
       GSREFI,GSMINI,PAR0I,D0I,VK1I,VK2I,VPD1I,VPD2I,VMFD0I, &
       GSJAI,GSJBI,T0I,TREFI,TMAXI,SMD1I,SMD2I,SWPEXPI, &
       G0I,D0LI,GAMMAI,G1I,GKI,WLEAFI,NSIDESI)
! Get stomatal conductance parameters.
! Required input: File unit (UFILE), Stom cond model (MODELGS).
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    
    INTEGER UFILE,NODATES,MODELGS,NSIDES,NSIDESI,I
    INTEGER DATESGSI(maxdate),IOERROR,IDATE,NOGSDATES
    INTEGER, EXTERNAL :: IDATE50
    CHARACTER*10 DATES(maxdate)
    CHARACTER*3 CONDUNITS
    !REAL G0TABLEI(maxdate),G1TABLEI(maxdate)
    !REAL G0(maxdate),G1(maxdate)
    !REAL G0I(maxdate),G1I(maxdate)
    REAL G0,G1,G0I,G1I,I0
    REAL GSREFI,GSMINI,PAR0,D0,VPD1,VPD2,VMFD0,GSJA,GSJB
    REAL T0,TREF,TMAX,SMD1,SMD2,VK1,VK2,WLEAF,SWPEXP
    REAL GSREF,GSMIN,GAMMA,D0L,GAMMAI
    REAL SMD1I,SMD2I,SWPEXPI,D0LI
    REAL PAR0I,D0I,VPD1I,VPD2I,VK1I,VK2I,VMFD0I
    REAL GSJAI,GSJBI,T0I,TREFI,TMAXI,WLEAFI
    REAL GK,GKI
    
    NAMELIST /JARGS/ GSREF,GSMIN,PAR0,D0,VPD1,VPD2,VMFD0,GSJA,GSJB, &
                     T0,TREF,TMAX,SMD1,SMD2,VK1,VK2,WLEAF,NSIDES,SWPEXP
    NAMELIST /BBGS/ DATES,G0,G1,GAMMA,WLEAF,NSIDES, &
                    SMD1,SMD2,SWPEXP
    NAMELIST /BBLGS/DATES, G0,G1,D0L,GAMMA,WLEAF,NSIDES, &
                    SMD1,SMD2,SWPEXP
    NAMELIST /BBMGS/DATES, G0,G1,GAMMA,WLEAF,NSIDES, &
                    SMD1,SMD2,SWPEXP
    NAMELIST /BBTGS/DATES, G0,G1,GK,GAMMA,WLEAF,NSIDES
    NAMELIST /BBGSCON/ NODATES,CONDUNITS
    
! Defaults
      GSMIN = 0.001
      SMD1 = 0.0
      SMD2 = 0.0
      GAMMA = 0.0
      DATES(1) = '01/01/50'
      NODATES = 1

      REWIND(UFILE)

      ! Ball-Berry and BBL Dates, and optionally units of conductance.
      CONDUNITS = 'CO2'  ! DEFAULT
      IF (MODELGS.EQ.2.OR.MODELGS.EQ.3.OR.MODELGS.EQ.4.OR. &
         MODELGS.EQ.5) THEN   
        READ (UFILE,BBGSCON,IOSTAT = IOERROR)
      ENDIF

      IF (MODELGS.EQ.2) THEN    ! Ball-Berry model parameters
        READ (UFILE,BBGS,IOSTAT = IOERROR)
        
        ! If conductance pars given for water:
        IF(CONDUNITS.EQ.'H2O'.OR.CONDUNITS.EQ.'h2o')THEN
            G0 = G0 / GSVGSC
            G1 = G1 / GSVGSC
          CALL SUBERROR('GS PARAMETERS FOR H2O WERE CONVERTED TO CO2.', &
          IWARN,0)
        ELSE
          CALL SUBERROR('GS PARAMETERS ARE ASSUMED TO BE FOR CO2.', &
          IWARN,0)
        ENDIF
        
        G0I = G0
        G1I = G1
        GAMMAI = GAMMA
        SMD1I = SMD1
        SMD2I = SMD2
        SWPEXPI = SWPEXP

      ELSE IF (MODELGS.EQ.3) THEN   ! Ball-Berry Leuning model parameters
        READ (UFILE, BBLGS,IOSTAT = IOERROR)
        
        ! If conductance pars given for water:
        IF(CONDUNITS.EQ.'H2O'.OR.CONDUNITS.EQ.'h2o')THEN
            G0 = G0 / GSVGSC
            G1 = G1 / GSVGSC
          CALL SUBERROR('GS PARAMETERS FOR H2O WERE CONVERTED TO CO2.', &
          IWARN,0)
        ELSE
          CALL SUBERROR('GS PARAMETERS ARE ASSUMED TO BE FOR CO2.', &
          IWARN,0)
        ENDIF
        
        G0I = G0
        G1I = G1
        GAMMAI = GAMMA
        D0LI = D0L
        SMD1I = SMD1
        SMD2I = SMD2
        SWPEXPI = SWPEXP
        
        
        
        IF (D0L.LT.0.0) &
         CALL SUBERROR('ERROR IN GS PARAMETERS: D0L MUST BE > 0', &
         IFATAL,0)
      
      ELSE IF (MODELGS.EQ.4) THEN   ! Ball-Berry-Medlyn model parameters
        READ (UFILE, BBMGS,IOSTAT = IOERROR)
        
        DO IDATE = 1,NODATES
            DATESGSI(IDATE) = IDATE50(DATES(IDATE))
        ENDDO
        
        ! If conductance pars given for water:
        IF(CONDUNITS.EQ.'H2O'.OR.CONDUNITS.EQ.'h2o')THEN
            G0 = G0 / GSVGSC
            G1 = G1 / GSVGSC
          CALL SUBERROR('GS PARAMETERS FOR H2O WERE CONVERTED TO CO2.', &
          IWARN,0)
        ELSE
          CALL SUBERROR('GS PARAMETERS ARE ASSUMED TO BE FOR CO2.', &
          IWARN,0)
        ENDIF
        
        G0I = G0
        G1I = G1
        GAMMAI = GAMMA
        SMD1I = SMD1
        SMD2I = SMD2
        SWPEXPI = SWPEXP
        
        
      ! Three-parameter Ball-Berry
      ELSE IF (MODELGS.EQ.5) THEN
        READ (UFILE, BBTGS, IOSTAT = IOERROR)

        G0I = G0
        G1I = G1
        GAMMAI = GAMMA
        GKI = GK
        
        ! If conductance pars given for water:
        IF(CONDUNITS.EQ.'H2O'.OR.CONDUNITS.EQ.'h2o')THEN
            G0 = G0 / GSVGSC
            G1 = G1 / GSVGSC
          CALL SUBERROR('GS PARAMETERS FOR H2O WERE CONVERTED TO CO2.', &
          IWARN,0)
        ELSE
          CALL SUBERROR('GS PARAMETERS ARE ASSUMED TO BE FOR CO2.', &
          IWARN,0)
        ENDIF        

        G0I = G0
        G1I = G1
        GAMMAI = GAMMA


      ELSE    ! Jarvis model parameters
      I0 = 0.0
      TMAX = 0.0
      GSJA = 0.0
      GSJB = 0.0
      SMD1 = 0.0
        READ (UFILE,JARGS,IOSTAT = IOERROR)
        GSREFI = GSREF
        GSMINI = GSMIN
        PAR0I = PAR0
        D0I = D0
        VPD1I = VPD1
        VPD2I = VPD2
      VK1I = VK1
      VK2I = VK2
        VMFD0I = VMFD0
        GSJAI = GSJA
        GSJBI = GSJB
        T0I = T0
        TREFI = TREF
        TMAXI = TMAX
      SMD1I = SMD1
      SMD2I = SMD2
      SWPEXPI = SWPEXP

      END IF

      IF (IOERROR.NE.0) THEN
        CALL SUBERROR( &
        'ERROR READING STOMATAL CONDUCTANCE PARAMETERS', &
        IFATAL,IOERROR)
      END IF
      WLEAFI = WLEAF
      NSIDESI = NSIDES

      RETURN
      END !ReadGS


!**********************************************************************
      SUBROUTINE READLEAFN(UFILE,MODELJM,MODELRD,NOLAY,NOAGEP, &
        NONDATES,DATESN,VALUESN)
! If MODELJM or MODELRD = 1, then leaf N contents are required to
! calculate Jmax or Rd. This subroutine then reads in the leaf N.
! Leaf N may be specified for a maximum of MAXDATE dates -
! linear interpolation is used between those dates.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER DATESN(maxdate)
    INTEGER UFILE,MODELJM,MODELRD,NOLAY,NOAGEP,NONDATES
    REAL VALUESN(maxdate,MAXLAY,MAXC)

      IF ((MODELJM.EQ.1).OR.(MODELRD.EQ.1)) THEN  !Leaf n contents reqd.
        CALL READPHYARRAY(UFILE,1,NOLAY,NOAGEP, &
        NONDATES,DATESN,VALUESN)
      END IF

      RETURN
      END !ReadLeafN


!**********************************************************************
       SUBROUTINE READPHYARRAY(UFILE,NARRAY,NOLAY,NOAGEP, &
       NDATE,IDATES,VALUESI)
! Read in an array of physiological parameters from UFILE.
! NARRAY is the number of the array to be read (1 = LEAFN; 2 = JMAX;
! 3 = VCMAX; 4 = RD; 5 = SLA; 6 = ALPHAJ)
! Parameters are given for up to MAXDATE dates, up to MAXLAY layers,
! and up to MAXC age classes.
!**********************************************************************
    
    USE maestcom
    IMPLICIT NONE
    INTEGER IDATES(maxdate),NOLAY,NOAGEP,NDATE,NOLAYERS,NODATES,NOAGES
    INTEGER UFILE,I,NARRAY,IOERROR,IP,IDATE,IAGE,ILAY
    REAL VALUES(maxdate*MAXLAY*MAXC),VALUESI(maxdate,MAXLAY,MAXC)
    REAL X(MAXLAY+1),D1,D2
    CHARACTER*10 DATES(maxdate)
    INTEGER, EXTERNAL :: IDATE50
    REAL, EXTERNAL :: RINTEG

      NAMELIST /NFOLCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /NFOL/ DATES,VALUES
      NAMELIST /JMAXCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /JMAX/ DATES,VALUES
      NAMELIST /VCMAXCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /VCMAX/ DATES,VALUES
      NAMELIST /RDCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /RD/ DATES,VALUES
      NAMELIST /SLACON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /SLA/ DATES,VALUES
      NAMELIST /AJQCON/ NODATES,NOLAYERS,NOAGES
      NAMELIST /AJQ/ DATES,VALUES

! Default values
      NODATES = 1
      NOLAYERS = 1
      NOAGES = 1
      DATES(1) = '01/01/50'
      DO 10 I = 1,MAXDATE*MAXLAY*MAXC
        VALUES(I) = -1.0
10    CONTINUE

! Read array size from file
      REWIND(UFILE)
      IF (NARRAY.EQ.1) THEN
        READ(UFILE,NFOLCON,IOSTAT=IOERROR)
        CALL SUBERROR('LEAF N ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.2) THEN
        READ(UFILE,JMAXCON,IOSTAT=IOERROR)
        CALL SUBERROR('JMAX ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.3) THEN
        READ(UFILE,VCMAXCON,IOSTAT=IOERROR)
        CALL SUBERROR('VCMAX ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.4) THEN
        READ(UFILE,RDCON,IOSTAT=IOERROR)
        CALL SUBERROR('RD ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.5) THEN
        READ(UFILE,SLACON,IOSTAT=IOERROR)
        IF (IOERROR.EQ.0) THEN
          CALL SUBERROR('SLA ARRAY:',IWARN,0)
        ELSE
          CALL SUBERROR( &
            'NO SLA VALUES: FOLIAGE GROWTH RESP NOT CALCULATED', &
            IWARN,IOERROR)
          NDATE = 0
          RETURN ! No SLA values - end subroutine here
        END IF
      ELSE IF (NARRAY.EQ.6) THEN
        READ(UFILE,AJQCON,IOSTAT=IOERROR)
        IF (IOERROR.EQ.0) THEN
          CALL SUBERROR('AJQ ARRAY:',IWARN,0)
        ELSE
          CALL SUBERROR( &
            'NO AJQ ARRAY; DEFAULT VALUE USED', &
            IWARN,IOERROR)
          NODATES = 0
          RETURN ! No AJQ values - end subroutine here
        END IF
      END IF

! Error handling
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR( &
        'WARNING: PROBLEM READING ARRAY SIZE',IWARN,IOERROR)
      END IF
      IF ((NODATES.GT.MAXDATE).OR.(NOLAYERS.GT.MAXLAY).OR. &
          (NOAGES.GT.MAXC)) THEN
        CALL SUBERROR( &
       'WARNING: ARRAY EXCEEDED BOUNDS: SOME DATA LOST', &
       IWARN,IOERROR)
        IF (NODATES.GT.MAXDATE) NODATES = MAXDATE
        IF (NOLAYERS.GT.MAXLAY) NOLAYERS = MAXLAY
        IF (NOAGES.GT.MAXC) NOAGES = MAXC
      END IF
      IF ((NOAGES.GT.1).AND.(NOAGEP.NE.NOAGES)) THEN
        CALL SUBERROR( &
        'ERROR: PHYSIOLOGY CLASSES & PARAMETER CLASSES DO NOT COINCIDE', &
        IWARN,IOERROR)
      END IF

! Read arrays from file
      REWIND(UFILE)
      IF (NARRAY.EQ.1) THEN
        READ(UFILE,NFOL,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.2) THEN
        READ(UFILE,JMAX,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.3) THEN
        READ(UFILE,VCMAX,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.4) THEN
        READ(UFILE,RD,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.5) THEN
        READ(UFILE,SLA,IOSTAT=IOERROR)
      ELSE IF (NARRAY.EQ.6) THEN
        READ(UFILE,AJQ,IOSTAT=IOERROR)
      END IF

! Error handling
      IF (IOERROR.NE.0) THEN
        CALL SUBERROR( &
        'ERROR: PROBLEM READING ARRAY',IFATAL,IOERROR)
      END IF

! Check values
      DO 20 I = 1,NOAGES*NOLAYERS*NODATES
        IF (VALUES(I).LT.0.0) &
          CALL SUBERROR('ERROR: VALUES MISSING FROM ARRAY', &
          IFATAL,0)
20    CONTINUE

! Set up x array
      DO 30 IP = 1,NOLAYERS+1
        X(IP) = (REAL(IP) - 1.0)/REAL(NOLAYERS)
30    CONTINUE

! Assign arrays to data
      DO 40 IDATE = 1,NODATES
        IDATES(IDATE) = IDATE50(DATES(IDATE))
40    CONTINUE
      DO 50 IDATE = 1,NODATES
        DO 50 IAGE = 1,NOAGES
          DO 50 I = 1,NOLAY
            D1 = (REAL(I) - 1)/REAL(NOLAY)
            D2 = (REAL(I))/REAL(NOLAY)
            VALUESI(IDATE,I,IAGE) = &
        RINTEG(D1,D2,X,VALUES, &
        (IDATE-1)*NOAGES*NOLAYERS+(IAGE-1)*NOLAYERS,NOLAY)
50    CONTINUE

! Fill in array in case no of ages is less than NOAGEP
      IF (NOAGES.LT.NOAGEP) THEN 
        DO 60 IDATE = 1,NODATES
          DO 60 ILAY = 1,NOLAY
          DO 60 I = NOAGES+1,NOAGEP
            VALUESI(IDATE,ILAY,I) = VALUESI(IDATE,ILAY,1)
60        CONTINUE 
      END IF

      NDATE = NODATES

      RETURN
      END !ReadPhyArray


!**********************************************************************
      SUBROUTINE READJMAX(UFILE, MODELJM, NOLAY, NOAGEP,  &
       NONDATES, DATESN, LEAFN,  &
       NOJDATES, DATESJ, JMAXTABLE, &
       NOVDATES, DATESV, VCMAXTABLE, &
       NOADATES, DATESA, AJQTABLE, &
       IECOI, EAVJI, EDVJI, DELSJI, EAVCI, EDVCI, DELSCI, &
       TVJUPI, TVJDNI, THETAI)
! Read in all parameters to do with Jmax and Vcmax.
! If MODELJM = 1, use the leaf N contents to calculate Jmax and Vcmax;
! otherwise read in arrays directly.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER DATESN(maxdate), DATESJ(maxdate)
    INTEGER DATESV(maxdate), DATESA(maxdate)
    INTEGER UFILE,IECO,MODELJM,IOERROR,NOJDATES
    INTEGER NONDATES,NOVDATES,IDATE,ILAY,NOLAY,IAGE
    INTEGER NOAGEP,IECOI,NOADATES
    REAL LEAFN(maxdate,MAXLAY,MAXC)
    REAL JMAXTABLE(maxdate,MAXLAY,MAXC)
    REAL VCMAXTABLE(maxdate,MAXLAY,MAXC)
    REAL AJQTABLE(maxdate,MAXLAY,MAXC)
    REAL JMAXA,JMAXB,THETA,EAVJ,EDVJ,DELSJ,AJQ
    REAL EAVC,EDVC,DELSC,TVJUP,TVJDN,VCMAXA,VCMAXB
    REAL EAVCI,EDVCI,DELSCI,EAVJI,EDVJI,DELSJI,THETAI
    REAL TVJUPI,TVJDNI

      NAMELIST /JMAXPARS/ THETA,EAVJ,EDVJ,DELSJ,AJQ,IECO
      NAMELIST /VCMAXPARS/ EAVC,EDVC,DELSC,TVJUP,TVJDN
      NAMELIST /JMAXN/ JMAXA,JMAXB,VCMAXA,VCMAXB

      IF (MODELJM.EQ.1) THEN   ! Calculate Jmax, Vcmax from leaf N
        REWIND(UFILE)
        READ (UFILE,JMAXN,IOSTAT = IOERROR)
        IF (IOERROR.NE.0) &
         CALL SUBERROR('ERROR: MISSING CONSTANTS FOR JMAX/N', &
         IFATAL,IOERROR)
        NOJDATES = NONDATES
        NOVDATES = NONDATES
        DO 10 IDATE = 1,NONDATES
          DATESJ(IDATE) = DATESN(IDATE)
          DATESV(IDATE) = DATESN(IDATE)
          DO 10 ILAY = 1,NOLAY
            DO 10 IAGE = 1,NOAGEP
              JMAXTABLE(IDATE,ILAY,IAGE) = &
                JMAXA*LEAFN(IDATE,ILAY,IAGE) + JMAXB
              VCMAXTABLE(IDATE,ILAY,IAGE) = &
                VCMAXA*LEAFN(IDATE,ILAY,IAGE) + VCMAXB
10      CONTINUE

      ELSE                     ! Read in Jmax, Vcmax arrays
        CALL READPHYARRAY(UFILE,2,NOLAY,NOAGEP, &
          NOJDATES,DATESJ,JMAXTABLE)
        CALL READPHYARRAY(UFILE,3,NOLAY,NOAGEP, &
          NOVDATES,DATESV,VCMAXTABLE)
      END IF

! Read in additional parameters
      REWIND (UFILE)
      AJQ = ALPHAQ !Default values
      EDVC = 0.0
      DELSC = 0.0
      IECO = 1 ! Ecocraft formulation of T-deps of Km and Gamma. For Montpied formulation, put 0. 
      TVJUP = -100.0
      TVJDN = -100.0

      READ (UFILE,JMAXPARS,IOSTAT = IOERROR)
      IF (IOERROR.NE.0) &
        CALL SUBERROR('INPUT ERROR: MISSING JMAXPARS',IFATAL,IOERROR)
      REWIND (UFILE)
      READ (UFILE,VCMAXPARS,IOSTAT = IOERROR)
      IF (IOERROR.NE.0) &
        CALL SUBERROR('INPUT ERROR: MISSING VCMAXPARS',IFATAL,IOERROR)

      EAVCI = EAVC
      EDVCI = EDVC
      DELSCI = DELSC
      EAVJI = EAVJ
      EDVJI = EDVJ
      DELSJI = DELSJ
      THETAI = THETA
      IECOI = IECO
      TVJUPI = TVJUP
      TVJDNI = TVJDN

! Read in quantum yield array
      NOADATES = 0
      CALL READPHYARRAY(UFILE,6,NOLAY,NOAGEP, &
          NOADATES,DATESA,AJQTABLE)
      IF (NOADATES.EQ.0) THEN
      NOADATES = 1
        DO 20 ILAY = 1,NOLAY
          DO 20 IAGE = 1,NOAGEP
          AJQTABLE(1,ILAY,IAGE) = AJQ
20      CONTINUE
    END IF

      RETURN
      END !ReadJmax


!**********************************************************************
      SUBROUTINE READRD(UFILE, MODELRD, NOLAY, NOAGEP, NONDATES, &
       DATESN, LEAFN, NORDATES, DATESRD, RDTABLE, &
       NOFQDATES, DATESRFQ, Q10FTABLE, &
       RTEMPI, DAYRESPI, EFFYRFI, TBELOWI)
! Read in all parameters to do with leaf respiration rate.
! If MODELRD = 1, use the leaf N contents to calculate Rd0;
! otherwise read in array directly.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER DATESN(maxdate), DATESRD(maxdate), DATESRFQ(maxdate)
    INTEGER UFILE,MODELRD,NOLAY,NOAGEP,NONDATES,NOFQDATES
    INTEGER IOERROR,NORDATES,IDATE,ILAY,IAGE
    
    REAL LEAFN(maxdate,MAXLAY,MAXC)
    REAL RDTABLE(maxdate,MAXLAY,MAXC), Q10FTABLE(maxdate)
    REAL K10F,Q10F,RTEMP,TBELOW,DAYRESP,EFFYRF,RDA,RDB
    REAL K10FI,RTEMPI,DAYRESPI,EFFYRFI,TBELOWI

      NAMELIST /RDPARS/ Q10F, RTEMP, TBELOW, DAYRESP, EFFYRF
      NAMELIST /RDN/ RDA, RDB

      IF (MODELRD.EQ.1) THEN   ! Calculate Jmax, Vcmax from leaf N
        REWIND(UFILE)
        READ (UFILE,RDN,IOSTAT = IOERROR)
        IF (IOERROR.NE.0) &
         CALL SUBERROR('ERROR: MISSING CONSTANTS FOR RD/N', &
         IFATAL,IOERROR)
        NORDATES = NONDATES
        DO 10 IDATE = 1,NONDATES
          DATESRD(IDATE) = DATESN(IDATE)
          DO 10 ILAY = 1,NOLAY
            DO 10 IAGE = 1,NOAGEP
              RDTABLE(IDATE,ILAY,IAGE) = &
                RDA*LEAFN(IDATE,ILAY,IAGE) + RDB
10      CONTINUE

      ELSE                     ! Read in RD0 arrays
        CALL READPHYARRAY(UFILE,4,NOLAY,NOAGEP, &
          NORDATES,DATESRD,RDTABLE)
      END IF

      DAYRESP = 1.0            ! Default value - no effect of light
      RTEMP = 0.0              ! Default value - Rd at zero degrees
      EFFYRF = 0.0             ! Default value - don't calc growth respn
      Q10F = 0.0               ! Indicates missing value
      TBELOW = -100.0          ! Default value - Rd continues at all Ts
      NOFQDATES = 1
      Q10FTABLE(1) = 0.0
      DATESRFQ(1) = 0

! Read in additional parameters
      REWIND (UFILE)
      READ (UFILE,RDPARS,IOSTAT = IOERROR)
      Q10FTABLE(1) = Q10F
      RTEMPI = RTEMP
      DAYRESPI = DAYRESP
      EFFYRFI = EFFYRF
      TBELOWI = TBELOW

    CALL READARRAY(UFILE,4,NOFQDATES,DATESRFQ,Q10FTABLE)

      IF (Q10FTABLE(1).EQ.0.0) &
        CALL SUBERROR('INPUT ERROR: MISSING Q10F',IFATAL,IOERROR)

      RETURN
      END !ReadRD


!**********************************************************************
      SUBROUTINE READRW(UFILE,MODELRW,EFFYRWI,RMWI,RTEMPWI, &
       NOWQDATES, DATESRWQ, Q10WTABLE, &
       COLLA,COLLK,STEMSDWI,RMAI,STEMFORMI)
! Read in or calculate parameters to do with woody respiration rate.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE, DATESRWQ(maxdate), IOERROR, MODELRW
    INTEGER NOWQDATES
    
    REAL Q10WTABLE(maxdate)
    REAL COLLA, COLLK, STEMSDWI, RMAI, STEMFORMI
    REAL EFFY,RM,RMA,RTEMP,Q10W, STEMFORM,CA,CK
    REAL STEMSDW,RTEMPWI,EFFYRWI,RMWI
    
    NAMELIST /WRESP/ EFFY,RM,RMA,RTEMP,Q10W,STEMFORM
    NAMELIST /COLLWRESP/ CA,CK,STEMSDW

      EFFY = 0.0        ! Missing value - growth respn not calculated
      RM = 0.0          ! Missing value - maint respn not calculated
      RMA = 0.0         ! Missing value - maint respn not calculated
      RTEMP = 0.0       ! Default temp at which maint resp specified
      NOWQDATES = 1
      DATESRWQ(1) = 0
      Q10WTABLE(1) = 0.0 ! Missing value - will cause error if RM spec

      REWIND (UFILE)
      READ (UFILE,WRESP,IOSTAT = IOERROR)
      Q10WTABLE(1) = Q10W
    CALL READARRAY(UFILE,5,NOWQDATES,DATESRWQ,Q10WTABLE)

      IF (EFFY.EQ.0.0) &
       CALL SUBERROR('WARNING: WOODY GROWTH RESP NOT CALCULATED', &
       IWARN,IOERROR)
      IF ((MODELRW.NE.1).AND.RM.EQ.0.0.AND.RMA.EQ.0.0) THEN
        CALL SUBERROR('WARNING: WOODY MAINT RESP NOT CALCULATED',  &
        IWARN,IOERROR)
      ELSE IF (Q10WTABLE(1).EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10W',IFATAL,IOERROR)
      END IF

      RTEMPWI = RTEMP
      EFFYRWI = EFFY
      RMWI = RM
      RMAI = RMA
      STEMFORMI = STEMFORM

      IF (MODELRW.EQ.1) THEN ! Collelongo model
        REWIND(UFILE)        
        READ(UFILE,COLLWRESP,IOSTAT = IOERROR)
        IF (IOERROR.NE.0) &
          CALL SUBERROR('ERROR: MISSING INFO FOR COLLELONGO RW MODEL', &
          IFATAL,IOERROR)
        COLLA = CA
        COLLK = CK
        STEMSDWI = STEMSDW
      END IF

      IF (RM.EQ.0.0.AND.RMA.GT.0.0) MODELRW = 2

      RETURN
      END !ReadRW


!**********************************************************************
      SUBROUTINE READRR(UFILE,RMCRI,RMFRI,Q10RI,RTEMPRI)
! Read in or calculate parameters to do with root respiration rate.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR
    REAL RMCRI,RMCR,RMFRI,RMFR,Q10R,Q10RI
    REAL RTEMPR,RTEMPRI

    NAMELIST /RRESP/ RMCR,RMFR,Q10R,RTEMPR

      REWIND (UFILE)
      RMCR = 0.0        ! Missing value - maint respn not calculated
      RMFR = 0.0        ! Missing value - maint respn not calculated
      RTEMPR = 0.0      ! Default temp at which maint resp specified
      Q10R = 0.0        ! Missing value - will cause error if RM spec
      READ (UFILE,RRESP,IOSTAT = IOERROR)
      IF ((RMCR.EQ.0.0).AND.(RMFR.EQ.0.0)) THEN
        CALL SUBERROR('WARNING: ROOT MAINT RESP NOT CALCULATED', &
        IWARN,IOERROR)
      ELSE IF (Q10R.EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10R',IFATAL,IOERROR)
      END IF

      Q10RI = Q10R
      RTEMPRI = RTEMPR
      RMCRI = RMCR
      RMFRI = RMFR

      END !ReadRR


!**********************************************************************
      SUBROUTINE READRB(UFILE,RMBI,Q10BI,RTEMPBI)
! Read in or calculate parameters to do with branch respiration rate.
!**********************************************************************

     USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,IOERROR
    REAL RMBI,RMB,Q10BI,Q10B,RTEMPB,RTEMPBI
    
    NAMELIST /BRESP/ RMB,Q10B,RTEMPB

      REWIND (UFILE)
      RMB = 0.0        ! Missing value - maint respn not calculated
      RTEMPB = 0.0     ! Default temp at which maint resp specified
      Q10B = 0.0       ! Missing value - will cause error if RM spec
      READ (UFILE,BRESP,IOSTAT = IOERROR)
      IF (RMB.EQ.0.0) THEN
        CALL SUBERROR('WARNING: BRANCH MAINT RESP NOT CALCULATED', &
        IWARN,IOERROR)
      ELSE IF (Q10B.EQ.0.0) THEN
        CALL SUBERROR('INPUT ERROR: MISSING Q10B',IFATAL,IOERROR)
      END IF

      Q10BI = Q10B
      RTEMPBI = RTEMPB
      RMBI = RMB

      END !ReadRB

!**********************************************************************
    SUBROUTINE READARRAY(UFILE,INDEX,NDATESI,DATESI,RATES)
! Used to read in arrays of values specified by date only -
! only foliage and wood Q10's at present
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,DATESI(maxdate),IOERROR,INDEX,NDATES,I,NDATESI
    REAL RATES(maxdate),RATESI(maxdate)
    CHARACTER*10 DATES(maxdate)
    INTEGER, EXTERNAL :: IDATE50

      NAMELIST /WOODRESP/NDATES, DATES, RATES
      NAMELIST /BRANRESP/NDATES, DATES, RATES
      NAMELIST /SOILRESP/NDATES, DATES, RATES
      NAMELIST /FOLQ10/NDATES, DATES, RATES
      NAMELIST /WOODQ10/NDATES, DATES, RATES
      NAMELIST /BRANQ10/NDATES, DATES, RATES
      NAMELIST /SOILQ10/NDATES, DATES, RATES

      REWIND(UFILE)
    IF (INDEX.EQ.1) THEN
      READ(UFILE, WOODRESP, IOSTAT = IOERROR)
    ELSE IF (INDEX.EQ.2) THEN
      READ(UFILE, BRANRESP, IOSTAT = IOERROR)
    ELSE IF (INDEX.EQ.3) THEN
      READ(UFILE, SOILRESP, IOSTAT = IOERROR)
    ELSE IF (INDEX.EQ.4) THEN
      READ(UFILE, FOLQ10, IOSTAT = IOERROR)
    ELSE IF (INDEX.EQ.5) THEN
      READ(UFILE, WOODQ10, IOSTAT = IOERROR)
    ELSE IF (INDEX.EQ.6) THEN
      READ(UFILE, BRANQ10, IOSTAT = IOERROR)
    ELSE IF (INDEX.EQ.7) THEN
      READ(UFILE, SOILQ10, IOSTAT = IOERROR)
    END IF

    IF (IOERROR.EQ.0) THEN
      NDATESI = NDATES
      DO 10 I = 1,NDATES
        DATESI(I) = IDATE50(DATES(I))
        RATESI(I) = RATES(I)
10      CONTINUE
      DATESI(NDATES+1) = 100000 ! Limit Date
    END IF
    
    RETURN
    END


!!**********************************************************************
!! CURRENTLY UNUSED 3/6/98 BEM
!      SUBROUTINE READRSOIL(UFILE,RSI,RTEMPSI,Q10SI)
!! Read in or calculate parameters to do with soil respiration rate.
!! Initially assuming Q10 relationship - to be replaced by Hanson et al 1993.
!!**********************************************************************
!
!      USE maestcom
!      NAMELIST /SRESP/ RSOIL, RTEMPS, Q10S
!      INTEGER UFILE
!
!      REWIND (UFILE)
!      RSOIL = 0.0       ! Missing value - soil respn not calculated
!      RTEMPS = 0.0      ! Default temp at which soil resp specified
!      Q10S = 0.0        ! Missing value - will cause error if RS spec
!      READ (UFILE,SRESP,IOSTAT = IOERROR)
!      IF (RSOIL.EQ.0.0) THEN
!        CALL SUBERROR('WARNING: SOIL RESPIRATION NOT CALCULATED', &
!        IWARN,IOERROR)
!      ELSE IF (Q10S.EQ.0.0) THEN
!        CALL SUBERROR('INPUT ERROR: MISSING Q10S',IFATAL,IOERROR)
!      END IF
!
!      Q10SI = Q10S
!      RTEMPSI = RTEMPS
!      RSI = RSOIL
!
!      RETURN
!      END !ReadRSoil


!**********************************************************************
      SUBROUTINE PHYINTERP(IDATE,NODATES,IDATEARR,PARAMTABLE, &
        NOLAY,NOAGEP,PARAMS)
! Interpolate physiological parameters for this date from the
! date arrays read in from physiology file.
!**********************************************************************


    USE maestcom
    IMPLICIT NONE
    INTEGER IDATEARR(maxdate),IOERROR,NODATES,IDATE
    INTEGER NOLAY,NOAGEP,IAGE,ILAY,INDEX
    REAL PARAMTABLE(maxdate,MAXLAY,MAXC)
    REAL PARAMS(MAXLAY,MAXC),SLOPE,Y1,Y2

! If only one date, or before the starting date, take first value
      IF ((NODATES.EQ.1).OR.IDATE.LE.IDATEARR(1)) THEN
        DO 10 IAGE = 1,NOAGEP
          DO 10 ILAY = 1,NOLAY
            PARAMS(ILAY,IAGE) = PARAMTABLE(1,ILAY,IAGE)
10      CONTINUE

! If after the final date, take last value
      ELSE IF (IDATE.GT.IDATEARR(NODATES)) THEN
        DO 20 IAGE = 1,NOAGEP
          DO 20 ILAY = 1,NOLAY
            PARAMS(ILAY,IAGE) = PARAMTABLE(NODATES,ILAY,IAGE)
20      CONTINUE

! Otherwise have to interpolate
      ELSE
        INDEX = 1
        DO WHILE (IDATE.GT.IDATEARR(INDEX))
          INDEX = INDEX + 1
        END DO
        SLOPE = REAL((IDATE - IDATEARR(INDEX-1)))/ &
            REAL((IDATEARR(INDEX) - IDATEARR(INDEX-1)))
        DO 30 IAGE = 1,NOAGEP
          DO 30 ILAY = 1,NOLAY
            Y1 = PARAMTABLE(INDEX-1, ILAY, IAGE)
            Y2 = PARAMTABLE(INDEX, ILAY, IAGE)
            PARAMS(ILAY,IAGE) = Y1 + SLOPE*(Y2 - Y1)
30      CONTINUE

      END IF

      RETURN
      END !PhyInterp

!**********************************************************************
      SUBROUTINE PHYINTERP2(IDATE,NODATES,IDATEARR,PARAMTABLE, &
       PARAM)
! Same as PHYINTERP but for parameter which is not specified by age & layer
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATEARR(maxdate),IDATE,NODATES,INDEX
    REAL PARAMTABLE(maxdate),PARAM,SLOPE,Y1,Y2
    
! If only one date, or before the starting date, take first value
      IF ((NODATES.EQ.1).OR.IDATE.LE.IDATEARR(1)) THEN
        PARAM = PARAMTABLE(1)

! If after the final date, take last value
      ELSE IF (IDATE.GT.IDATEARR(NODATES)) THEN
        PARAM = PARAMTABLE(NODATES)

! Otherwise have to interpolate
      ELSE
        INDEX = 1
        DO WHILE (IDATE.GT.IDATEARR(INDEX))
          INDEX = INDEX + 1
        END DO
        SLOPE = REAL((IDATE - IDATEARR(INDEX-1)))/ &
            REAL((IDATEARR(INDEX) - IDATEARR(INDEX-1)))
        Y1 = PARAMTABLE(INDEX-1)
        Y2 = PARAMTABLE(INDEX)
        PARAM = Y1 + SLOPE*(Y2 - Y1)
      END IF

      RETURN
      END !PhyInterp2


!**********************************************************************
     REAL FUNCTION RINTEG(D1,D2,X,VALUES,IOFFSET,NOLAY)
! Estimate physiological parameters for each canopy layer by averaging.
! At the moment, this is done by integrating the physiological parameter
! over the height of the layer, assuming foliage is evenly distributed.
! Strictly speaking, should probably take the beta-distribution of
! foliage into account - but since the specification of physiology for
! different layers is woolly anyway, this is a good approx.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IP,IOFFSET,NOLAY
    REAL X(MAXLAY+1),VALUES(MAXLAY*MAXC*maxdate)
    REAL D1,D2
! NB The array VALUES can be of variable size - since this routine is used
! to handle reflectance and transmittance, and leaf N, Jmax, Vcmax and Rd
! arrays (as of 8/97). Here it should be set to the maximum size of
! any array to be passed.

      RINTEG = 0.0
      IP = 1

      DO WHILE (X(IP+1).LT.D1)
        IP = IP+1
      END DO

      DO WHILE (D2.GT.X(IP+1))
        RINTEG = RINTEG + NOLAY* &
            VALUES(IP+IOFFSET)*(X(IP+1)-AMAX1(D1,X(IP)))
        IP = IP+1
      END DO

      RINTEG = RINTEG + NOLAY* &
          VALUES(IP+IOFFSET)*(D2-AMAX1(D1,X(IP)))

100   FORMAT(4(1X,F8.3))

      RETURN
      END !RInteg

!**********************************************************************
      SUBROUTINE READCONTREES(UFILE,NOALLTREES,DX,DY,XMAX,YMAX, &
        NOTREESI,NOTARGETSI,ITARGETSI,IPLOTSHAPE)
! Read in controls about tree numbers: need 
! (1) no of trees to be included in the shading calculations
! (2) the numbers of the target trees. 
! (1) NOTREES is the number of trees to be considered in the calculation.
! (Default: NOTREES = all trees in plot).
! (2) ITARGETS is an array with the numbers of the target trees. 
! Just one can be specified with NOTARGET. 
! If NORANDOM is defined, then calculations are done for NORANDOM
! randomly-chosed target trees.
! If none of NORANDOM, NOTARGET and ITARGETS is defined, calculations 
! are done for all trees.
! Trees within EDGEDIST m of the plot edges are exempted. 
!**********************************************************************

      ! Comment when using g95!
    !USE IFPORT    ! For Intel Visual Fortran - allows use of random function
    USE maestcom
    IMPLICIT NONE
    INTEGER ITARGETS(MAXT),ITARGETSI(MAXT)
    INTEGER UFILE,NOALLTREES,IPLOTSHAPE,NOTREES,NOTARGET,NORANDOM
    INTEGER IOERROR,IRAN,ITREE,NOTARGETS,ITAR,NOTREESI,NOTARGETSI
    INTEGER(4) IFLAG
    INTEGER, EXTERNAL :: INEDGES
    REAL(4) RANVAL
    REAL DX(MAXT),DY(MAXT),WEIGHTSI(MAXT),WEIGHTS(MAXT)
    REAL XMAX,YMAX,EDGEDIST
    NAMELIST /TREESCON/ NOTREES,NOTARGET,ITARGETS,WEIGHTS,NORANDOM,EDGEDIST
    REAL :: RANNUM
    
! Default values
      NOTREES = 0
      NOTARGET = 0
      NORANDOM = 0
      EDGEDIST = 0
      ITARGETS = 0 ! ARRAY
      
! Read file
      REWIND (UFILE)
      READ (UFILE,TREESCON,IOSTAT = IOERROR)

! Check the namelist was there
      IF (IOERROR.NE.0) &
        CALL SUBERROR('INPUT ERROR: MISSING TREES IN CONTROL FILE', &
        IFATAL,IOERROR)

! Set number of trees used in calculations
      IF (NOTREES.EQ.0.OR.NOTREES.GT.NOALLTREES) NOTREES = NOALLTREES
      IF (NOTREES.GT.MAXT) THEN
        CALL SUBERROR( &
          'WARNING: NOTREES IN CONTROL FILE EXCEEDED MAXIMUM', &
          IWARN,0)
        NOTREES = MAXT
      END IF

! Set up the array of target tree numbers ITARGETS
! Case 1: One target tree specified
      IF (NOTARGET.GT.0) THEN
        IF (NOTARGET.GT.NOALLTREES) &
        CALL SUBERROR('INCORRECT TARGET TREE NUMBER', &
        IFATAL,0)
        NOTARGETS = 1
        ITARGETS(1) = NOTARGET  

! Case 2: Trees to be chosen randomly
      ELSE IF (NORANDOM.GT.0) THEN
        IFLAG = 0
        IF (NORANDOM.GT.NOALLTREES) &
          CALL SUBERROR('TOO MANY TARGET TREES SPECIFIED', &
          IFATAL,0)
        NOTARGETS = NORANDOM
        DO 20 IRAN = 1,NORANDOM
!        
!        ! When using g95
30      CALL RANDOM_NUMBER(RANNUM)
        RANNUM = RANNUM*REAL(NOALLTREES+1)
        ITREE = NINT(RANNUM)
        
!30        RANVAL = RANDOM(IFLAG)
!          RANVAL = RANVAL*REAL(NOALLTREES+1)
!          ITREE = NINT(RANVAL)
          IF (ITREE.EQ.0) ITREE = 1
          IF (INEDGES(DX(ITREE),DY(ITREE),XMAX,YMAX,EDGEDIST).LT.0) &
            GOTO 30          
          ITARGETS(IRAN) = ITREE
20      CONTINUE
          
! Case 3: All trees are to be used, except those in the edges
      ELSE IF (ITARGETS(1).EQ.0) THEN
        NOTARGETS = NOALLTREES
        DO 40 ITAR = 1,NOALLTREES
          IF (INEDGES(DX(ITAR),DY(ITAR),XMAX,YMAX, &
          EDGEDIST).LT.0.AND.IPLOTSHAPE.EQ.0) THEN
            NOTARGETS = NOTARGETS - 1
          ELSE          
            ITARGETS(ITAR + NOTARGETS - NOALLTREES) = ITAR
          END IF
40      CONTINUE

      ELSE
! Case 4: A series of target trees is given.
        NOTARGETS = 0
        DO WHILE (ITARGETS(NOTARGETS+1).GT.0)
          NOTARGETS = NOTARGETS + 1
        END DO
      END IF

      NOTREESI = NOTREES
      NOTARGETSI = NOTARGETS
      DO 50 ITAR = 1,MAXT
        ITARGETSI(ITAR) = ITARGETS(ITAR)
50    CONTINUE

      RETURN
      END !ReadConTrees


!**********************************************************************
      SUBROUTINE READXYZ(UFILE,NOALLTREES,X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE, &
        DX,DY,DZ)
! Read in the X, Y co-ordinates of each crown.
! Calculate the Z co-ordinates of each crown from slope.
! If they are not in the file, then the co-ordinates should be
! calculated from the stocking and the size of the plot.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,NOALLTREES,IOERROR,I
    REAL DX(MAXT),DY(MAXT),DZ(MAXT),XYCOORDS(MAXT*2)
    REAL X0,Y0,XMAX,YMAX,XSLOPE,YSLOPE,SPACING,NOX,ZADD
    NAMELIST /XY/ XYCOORDS

      REWIND (UFILE)
      READ (UFILE,XY,IOSTAT = IOERROR)

      IF (IOERROR.NE.0) THEN
! Must calculate x, y co-ordinates
        CALL SUBERROR('WARNING: CALCULATING X, Y CO-ORDINATES', &
        IWARN,IOERROR)
        SPACING = SQRT(XMAX*YMAX/NOALLTREES)
        NOX = XMAX/SPACING + 1
        DO 10 I = 1,NOALLTREES
          DX(I) = (I-1)/NOX * SPACING
          DY(I) = MOD(REAL(I)-1,NOX) * SPACING
10      CONTINUE

      ELSE
! Read in x, y co-ordinates
        DO 20 I = 1,NOALLTREES
          DX(I) = XYCOORDS(2*I - 1) - X0
          DY(I) = XYCOORDS(2*I) - Y0
20      CONTINUE
      END IF

! Move x0,y0 to origin
      XMAX = XMAX - X0
      YMAX = YMAX - Y0
      X0 = 0.0
      Y0 = 0.0

! Calculate the z co-ordinates from slopes
      ZADD=0.0
      IF (XSLOPE.LT.0.0) ZADD=XMAX*SIN(XSLOPE)
      IF (YSLOPE.LT.0.0) ZADD=ZADD+YMAX*SIN(YSLOPE)
      ZADD=ABS(ZADD)
!  X and Y distances are measured on the slope so the height is based
!  on the sin(slope).
      DO 30 I = 1,NOALLTREES
        DZ(I) = DX(I)*SIN(XSLOPE) + DY(I)*SIN(YSLOPE) + ZADD
30    CONTINUE

      RETURN
      END ! ReadXY


!**********************************************************************
      SUBROUTINE GETLEAFAREA(UFILE,IFLUSH,DT1I,DT2I,DT3I,DT4I, &
        EXPTIMEI,APP,EXPAN,NOALLTREES,NOLADATES,DATESLA,FLT)
! A subroutine to read in leaf area array.
! First tries to calculate it from phenology (Wang, Rey & Jarvis 1998).
! Otherwise reads in array directly.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,MAXLEAVES,IFLUSH,IOERROR,NOLADATES,NOALLTREES
    INTEGER, EXTERNAL :: IDATE50
    INTEGER DATESLA(maxdate)
    CHARACTER*10 FLUSHDATE
    REAL FLT(maxdate,MAXT)
    REAL DT1,DT1I,DT2,DT2I,DT3,DT3I,DT4,DT4I
    REAL EXPTIME,SIZELEAF,EXPTIMEI,APP,EXPAN

      NAMELIST /PHENOLOGY/ FLUSHDATE,DT1,DT2,DT3,DT4, &
        EXPTIME,MAXLEAVES,SIZELEAF

      REWIND(UFILE)
      IFLUSH = 0
      READ(UFILE,PHENOLOGY,IOSTAT=IOERROR)
      IF (IOERROR.EQ.0) THEN
        CALL SUBERROR('LEAF AREA FROM PHENOLOGY:',IWARN,0)
        NOLADATES = 0
        IFLUSH = IDATE50(FLUSHDATE)
        DT1I = DT1
        DT2I = DT2
        DT3I = DT3
        DT4I = DT4
        EXPTIMEI = EXPTIME
! Rate of leaf appearance (leaf per day)
        APP = 2. * MAXLEAVES / (DT1*(2*DT2-DT1))
! Rate of leaf area expansion (m2 leaf per day)
        EXPAN = SIZELEAF / EXPTIME
      ELSE
        CALL READTREEARRAY(UFILE,5,NOALLTREES,NOLADATES,DATESLA,FLT)
      END IF

      RETURN
      END ! GetLeafArea


!**********************************************************************
        SUBROUTINE PHENOL(IDAY,ISTART,IFLUSH,DT1,DT2,DT3,DT4, &
          EXPTIME,APP,EXPAN,STOCKING,NOTREES,FOLT,TOTLAI,NEWCANOPY)
! Implements phenology routine of Wang et al (1998) GCB to appear.
! Parameters:
!   IFLUSH = date of flushing (in days-since-1950 format)
!   DT1 = time from flushing to end of first flush (d)
!   DT2 = time from flushing to end of second flush (d)
!   DT3 = time from flushing to beginning of senescencs (d)
!   DT4 = time from flushing to leaf fall (d)
!   EXPTIME = time a leaf takes to expand (d)
!   APP = rate of leaf appearance (leaves / day)
!   EXPAN = rate of leaf expansion (m2 leaf / day)
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER I,IDAY,ISTART,IFLUSH,NEWCANOPY,ITREE,NOTREES
    REAL FOLT(MAXT)
    REAL EXPTIME, APP,EXPAN
    REAL DT1,DT2,DT3,DT4,TOTLAI,STOCKING
    REAL THRESH_FOLT

! Apply 7-part formula in Wang et al 1998
      I = IDAY + ISTART - IFLUSH
      IF (I.LE.0) THEN
        FOLT(1) = 0.0
      ELSE IF (I.LE.EXPTIME) THEN
        FOLT(1) = APP*EXPAN*(I**3)/6.
        NEWCANOPY = 1
      ELSE IF (I.LE.DT1) THEN
        FOLT(1) = APP*EXPAN*(I**3 - (I - EXPTIME)**3)/6.
        NEWCANOPY = 1
      ELSE IF (I.LE.(DT1+EXPTIME)) THEN
        FOLT(1) = APP*EXPAN*(DT1*(3.*I**2 + DT1**2 - 3.*DT1*I) -  &
          (I - EXPTIME)**3)/6.
        NEWCANOPY = 1
      ELSE IF (I.LE.DT2) THEN
        FOLT(1) = APP*EXPAN*DT1*EXPTIME*(I-0.5*(DT1+EXPTIME))
        NEWCANOPY = 1
      ELSE IF (I.LE.(DT2+EXPTIME)) THEN
        FOLT(1) = 0.5*APP*EXPAN*DT1*(2.*DT2*I-DT2**2-DT1*EXPTIME &
         - (EXPTIME - I)**2)
        NEWCANOPY = 1
      ELSE IF (I.LE.DT3) THEN
        FOLT(1) = APP*EXPAN*DT1*EXPTIME*(DT2-0.5*DT1)
      ELSE IF (I.LE.DT4) THEN
        FOLT(1) = (DT4-I)*(DT2-0.5*DT1)*APP*EXPAN*DT1*EXPTIME &
          / (DT4 - DT3)
        NEWCANOPY = 1
      ELSE 
        FOLT(1) = 0.0
      END IF

      DO 10 ITREE = 2,NOTREES
        FOLT(ITREE) = FOLT(ITREE)
10    CONTINUE
      TOTLAI = FOLT(1)*STOCKING

      RETURN
      END !Phenol

!**********************************************************************
      SUBROUTINE READTREEARRAY(UFILE,NARRAY,NOALLTREES, &
        NDATE,IDATES,VALUESI)
! Read in an array of tree parameters from UFILE.
! NARRAY is the number of the array to be read (1 = RADX; 2 = RADY;
! 3 = HTCROWN; 4 = HTTRUNK; 5 = AREALEAF; 6 = DIAM; 7 = LEAFN (for understorey))
! Either read in values for all trees (NOALLTREES, up to MAXT trees) or
! read in average values. All values can be given for a series of dates.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATES(maxdate)
    INTEGER UFILE, I, IDATE, ITREE,IOERROR,INDEX,ITREES,NDATE
    INTEGER NODATES,NARRAY,NOALLTREES
    INTEGER, EXTERNAL :: IDATE50
    CHARACTER*8 DATES(maxdate)
    REAL VALUES(maxdate*MAXT),VALUESI(maxdate,MAXT)

      NAMELIST /INDIVRADX/ NODATES, DATES, VALUES
      NAMELIST /ALLRADX/ NODATES, DATES, VALUES
      NAMELIST /INDIVRADY/ NODATES, DATES, VALUES
      NAMELIST /ALLRADY/ NODATES, DATES, VALUES
      NAMELIST /INDIVHTCROWN/ NODATES, DATES, VALUES
      NAMELIST /ALLHTCROWN/ NODATES, DATES, VALUES
      NAMELIST /INDIVHTTRUNK/ NODATES, DATES, VALUES
      NAMELIST /ALLHTTRUNK/ NODATES, DATES, VALUES
      NAMELIST /INDIVLAREA/ NODATES, DATES, VALUES
      NAMELIST /ALLLAREA/ NODATES, DATES, VALUES
      NAMELIST /INDIVDIAM/ NODATES, DATES, VALUES
      NAMELIST /ALLDIAM/ NODATES, DATES, VALUES
      NAMELIST /INDIVFOLN/ NODATES, DATES, VALUES
      NAMELIST /ALLFOLN/ NODATES, DATES, VALUES
      NAMELIST /INDIVUSLAI/ NODATES, DATES, VALUES ! not implemented
      NAMELIST /ALLUSLAI/ NODATES,DATES,VALUES

! Default values
      NODATES = 0
      DATES(1) = '01/01/50'  ! If only one date, doesn't matter what it is
      DO 10 I = 1,MAXDATE*MAXT
        VALUES(I) = -1.0
10    CONTINUE

! Try to read arrays for individual trees first
      REWIND(UFILE)
      IF (NARRAY.EQ.1) THEN
        READ(UFILE,INDIVRADX,IOSTAT=IOERROR)
        !CALL SUBERROR('X RADII ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.2) THEN
        READ(UFILE,INDIVRADY,IOSTAT=IOERROR)
        !CALL SUBERROR('Y RADII ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.3) THEN
        READ(UFILE,INDIVHTCROWN,IOSTAT=IOERROR)
        !CALL SUBERROR('CROWN HEIGHT ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.4) THEN
        READ(UFILE,INDIVHTTRUNK,IOSTAT=IOERROR)
        !CALL SUBERROR('TRUNK HEIGHT ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.5) THEN
        READ(UFILE,INDIVLAREA,IOSTAT=IOERROR)
        !CALL SUBERROR('LEAF AREA ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.6) THEN
        READ(UFILE,INDIVDIAM,IOSTAT=IOERROR)
        !CALL SUBERROR('DIAMETER ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.7) THEN
        READ(UFILE,INDIVFOLN,IOSTAT=IOERROR)
        !CALL SUBERROR('UNDERSTOREY N ARRAY:',IWARN,0)
      ELSE IF (NARRAY.EQ.8) THEN
        READ(UFILE,INDIVUSLAI,IOSTAT=IOERROR)
        !CALL SUBERROR('UNDERSTOREY LAI ARRAY:',IWARN,0)
      END IF

      IF (IOERROR.NE.-1) THEN
! Process arrays, if read in
        IF (NODATES.GT.MAXDATE) THEN
          CALL SUBERROR( &
          'WARNING: TOO MANY DATES: SOME DATA LOST', &
          IWARN,IOERROR)
          NODATES = MAXDATE
        END IF
        INDEX = 1
        DO 20 IDATE = 1,NODATES
          IDATES(IDATE) = IDATE50(DATES(IDATE))
20      CONTINUE
        DO 30 ITREE = 1,NOALLTREES
          DO 30 IDATE = 1,NODATES
            IF (VALUES(INDEX).LT.0) &
              CALL SUBERROR('MISSING DATA',IFATAL,0)
            VALUESI(IDATE,ITREE) = VALUES(INDEX)
            INDEX = INDEX+1
30      CONTINUE

      ELSE
! Read in values for all trees
        REWIND(UFILE)
        IF (NARRAY.EQ.1) THEN
          READ(UFILE,ALLRADX,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.2) THEN
          READ(UFILE,ALLRADY,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.3) THEN
          READ(UFILE,ALLHTCROWN,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.4) THEN
          READ(UFILE,ALLHTTRUNK,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.5) THEN
          READ(UFILE,ALLLAREA,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.6) THEN
          READ(UFILE,ALLDIAM,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.7) THEN
          READ(UFILE,ALLFOLN,IOSTAT=IOERROR)
        ELSE IF (NARRAY.EQ.8) THEN
          READ(UFILE,ALLUSLAI,IOSTAT=IOERROR)
        END IF

        IF ((IOERROR.NE.0).AND.(NARRAY.LT.6).AND.NARRAY.NE.2) &  ! Missing diam, rady, leaf N is OK
          CALL SUBERROR('ERROR: MISSING TREE DATA',IFATAL,IOERROR)
        IF (NODATES.GT.MAXDATE) THEN
          CALL SUBERROR( &
          'WARNING: TOO MANY DATES: SOME DATA LOST', &
          IWARN,IOERROR)
          NODATES = MAXDATE
        END IF
        
! Assign arrays to data
        DO 40 IDATE = 1,NODATES
          IDATES(IDATE) = IDATE50(DATES(IDATE))
40      CONTINUE
        DO 50 ITREES = 1,NOALLTREES
          DO 50 IDATE = 1,NODATES
            VALUESI(IDATE,ITREES) = VALUES(IDATE)
50      CONTINUE

      END IF

      NDATE = NODATES
      RETURN
      END !ReadTreeArray


!**********************************************************************
      SUBROUTINE CALCLAI(NOLADATES,FLT,NOALLTREES, &
        XMAX,YMAX,XSLOPE,YSLOPE,TOTLAITABLE)
! Calculate total LAI of the plot - based on horizontal plot area
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATE, ITREE
    INTEGER NOLADATES,NOALLTREES
    REAL FLT(maxdate,MAXT),TOTLAITABLE(maxdate)
    REAL XMAX,YMAX,XSLOPE,YSLOPE,GROUNDAREA
    REAL TOTLAI(MAXDATE)

      GROUNDAREA = XMAX*COS(XSLOPE)*YMAX*COS(YSLOPE)
      DO 10 IDATE = 1,NOLADATES
        TOTLAITABLE(IDATE) = 0.0
        DO 20 ITREE = 1,NOALLTREES
          TOTLAITABLE(IDATE) = TOTLAITABLE(IDATE) + FLT(IDATE,ITREE)
20      CONTINUE
        TOTLAITABLE(IDATE) = TOTLAITABLE(IDATE)/GROUNDAREA
10    CONTINUE

      RETURN
      END ! CalcLAI


!**********************************************************************
      SUBROUTINE TREEINTERP(IDAY,ISTART,NODATES,IDATEARR,PARAMTABLE, &
        NOTREES,PARAMS)
! Interpolate crown dimensions for this date from the
! date arrays read in from the trees file.
! Parameter NEWCANOPY indicates whether crown has changed - in which
! case grid points need to be reassigned.
! BM Jul07 NEWCANOPY is problematic. Change to (a) check change from previous values (not on a daily basis) and 
! (b) put thresholds for different dimensions. 
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER ITREE, INDEX, NOTREES,IDATE,IDAY,ISTART
    INTEGER IDATEARR(maxdate),NODATES
    
    REAL PARAMTABLE(maxdate,MAXT)
    REAL PARAMS(MAXT),SLOPE,Y1,Y2

      IDATE = IDAY + ISTART

! If no dates, take 0.0
      IF (NODATES.EQ.0) THEN
        DO 10 ITREE = 1,NOTREES
          PARAMS(ITREE) = 0.0
10      CONTINUE

! If only one date, or before the starting date, take first value
      ELSE IF ((NODATES.EQ.1).OR.IDATE.LE.IDATEARR(1)) THEN
        DO 20 ITREE = 1,NOTREES
          PARAMS(ITREE) = PARAMTABLE(1,ITREE)
20      CONTINUE

! If after the final date, take last value
      ELSE IF (IDATE.GT.IDATEARR(NODATES)) THEN
        DO 30 ITREE = 1,NOTREES
          PARAMS(ITREE) = PARAMTABLE(NODATES,ITREE)
30      CONTINUE

! Otherwise have to interpolate
      ELSE
        INDEX = 1
        DO WHILE (IDATE.GT.IDATEARR(INDEX))
          INDEX = INDEX + 1
        END DO
        SLOPE = REAL((IDATE - IDATEARR(INDEX-1)))/ &
            REAL((IDATEARR(INDEX) - IDATEARR(INDEX-1)))
        DO 40 ITREE = 1,NOTREES
            Y1 = PARAMTABLE(INDEX-1, ITREE)
            Y2 = PARAMTABLE(INDEX, ITREE)
            PARAMS(ITREE) = Y1 + SLOPE*(Y2 - Y1)
40      CONTINUE
      END IF

      RETURN
      END !InterpTree


!**********************************************************************
SUBROUTINE TREEINTERP2(IDAY,ISTART,NODATES,IDATEARR,PARAMTABLE,NOTREES,PARAMS)
! Like TREEINTERP, but for arrays that have dimension (maxdate), not
! (maxdate,MAXT), like most of them. 
! Addition to fix TOTLAI problems, Nov2010 (RAD).
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IDATEARR(maxdate),NODATES
    REAL PARAMTABLE(maxdate)
    REAL PARAMS,SLOPE,Y1,Y2
    INTEGER ITREE, INDEX, NOTREES,IDATE,IDAY,ISTART

    IDATE = IDAY + ISTART

    ! If no dates, take 0.0
    IF (NODATES.EQ.0) THEN
        PARAMS = 0.0
        
    ! If only one date, or before the starting date, take first value
    ELSE IF ((NODATES.EQ.1).OR.IDATE.LE.IDATEARR(1)) THEN
        PARAMS = PARAMTABLE(1)

    ! If after the final date, take last value
    ELSE IF (IDATE.GT.IDATEARR(NODATES)) THEN
            PARAMS = PARAMTABLE(NODATES)
    
    ! Otherwise have to interpolate
    ELSE
        INDEX = 1
        DO WHILE (IDATE.GT.IDATEARR(INDEX))
            INDEX = INDEX + 1
        END DO
        SLOPE = REAL((IDATE - IDATEARR(INDEX-1)))/ REAL((IDATEARR(INDEX) - IDATEARR(INDEX-1)))
            Y1 = PARAMTABLE(INDEX-1)
            Y2 = PARAMTABLE(INDEX)
            PARAMS = Y1 + SLOPE*(Y2 - Y1)
    END IF
    RETURN
END SUBROUTINE TREEINTERP2



!**********************************************************************
      SUBROUTINE INTERPOLATEP(IDAY, ISTART, &
       NOJDATES,DATESJ,JMAXTABLE, &
       NOVDATES,DATESV,VCMAXTABLE, &
       NORDATES,DATESRD,RDTABLE, &
       NOSLADATES,DATESSLA,SLATABLE, &
       NOADATES,DATESA,AJQTABLE, &
       NOFQDATES,DATESFQ,Q10FTABLE, &
       NOWQDATES,DATESWQ,Q10WTABLE, &
       NOLAY,NOAGEP, &
       JMAX25,VCMAX25,RD0,SLA,AJQ,Q10F,Q10W)
! Controls the calling of the interpolation routines to get daily values
! of Jmax, Vcmax, SLA, Rd.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER DATESJ(maxdate), DATESV(maxdate)
    INTEGER DATESRD(maxdate), DATESSLA(maxdate)
    INTEGER DATESA(maxdate), DATESGS(maxdate)
    INTEGER DATESFQ(maxdate),DATESWQ(maxdate)
    INTEGER IDAY,ISTART,NOJDATES,NOLAY,NOAGEP
    INTEGER NOVDATES,NORDATES,NOSLADATES,NOADATES,NOGSDATES
    INTEGER NOFQDATES,NOWQDATES
    
    REAL JMAXTABLE(maxdate,MAXLAY,MAXC)
    REAL VCMAXTABLE(maxdate,MAXLAY,MAXC)
    REAL RDTABLE(maxdate,MAXLAY,MAXC)
    REAL SLATABLE(maxdate,MAXLAY,MAXC)
    REAL AJQTABLE(maxdate,MAXLAY,MAXC)
    REAL Q10FTABLE(maxdate),Q10WTABLE(maxdate)
    REAL G0TABLE(maxdate),G1TABLE(maxdate)
    REAL JMAX25(MAXLAY,MAXC),VCMAX25(MAXLAY,MAXC)
    REAL RD0(MAXLAY,MAXC),SLA(MAXLAY,MAXC),AJQ(MAXLAY,MAXC)
    REAL G0,G1,Q10F,Q10W

! Interpolate to get daily values of physiological params
      CALL PHYINTERP(IDAY+ISTART,NOJDATES,DATESJ,JMAXTABLE, &
       NOLAY,NOAGEP,JMAX25)
      CALL PHYINTERP(IDAY+ISTART,NOVDATES,DATESV,VCMAXTABLE, &
       NOLAY,NOAGEP,VCMAX25)
      CALL PHYINTERP(IDAY+ISTART,NORDATES,DATESRD,RDTABLE, &
       NOLAY,NOAGEP,RD0)
      IF (NOSLADATES.GT.0) &
       CALL PHYINTERP(IDAY+ISTART,NOSLADATES,DATESSLA,SLATABLE, &
       NOLAY,NOAGEP,SLA)
      CALL PHYINTERP(IDAY+ISTART,NOADATES,DATESA,AJQTABLE, &
       NOLAY,NOAGEP,AJQ)
      CALL PHYINTERP2(IDAY+ISTART,NOFQDATES,DATESFQ,Q10FTABLE, &
       Q10F)
      CALL PHYINTERP2(IDAY+ISTART,NOWQDATES,DATESWQ,Q10WTABLE, &
       Q10W)

      RETURN
      END !InterpolateP


!**********************************************************************
      SUBROUTINE INTERPOLATET(IDAY, ISTART, &
       NOXDATES,DATESX,RXTABLE, &
       NOYDATES,DATESY,RYTABLE, &
       NOZDATES,DATESZ,RZTABLE, &
       NOTDATES,DATEST,ZBCTABLE, &
       NODDATES,DATESD,DIAMTABLE, &
       NOLADATES,DATESLA,FOLTABLE, &
       TOTLAITABLE,NOTREES, &
       RX,RY,RZ,ZBC,FOLT,TOTLAI,DIAM,STOCKING, &
       IFLUSH,DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN, &
       NEWCANOPY,CANOPYDIMS)
! Controls the calling of the interpolation routines to get daily values
! of crown heights and radii and leaf area.
! RXTABLE, RYTABLE, etc have all values of dimensions, for all trees and dates
! These are interpolated to give values for this date, for all trees, in RX, RY etc
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER DATESX(maxdate),DATESY(maxdate),DATESZ(maxdate)
    INTEGER DATEST(maxdate),DATESLA(maxdate),DATESD(maxdate)
    INTEGER IDAY,ISTART,NOXDATES,NOTREES,NOYDATES,NOZDATES,NOTDATES
    INTEGER NOLADATES,IFLUSH,NEWCANOPY,NODDATES,IHOUR

    REAL RX(MAXT),RY(MAXT),RZ(MAXT),ZBC(MAXT),FOLT(MAXT)
    REAL RXTABLE(maxdate,MAXT),RYTABLE(maxdate,MAXT)
    REAL RZTABLE(maxdate,MAXT),ZBCTABLE(maxdate,MAXT)
    REAL FOLTABLE(maxdate,MAXT),TOTLAITABLE(maxdate)  !!!
    REAL TOTLAITMP(MAXDATE)
    REAL DIAMTABLE(maxdate,MAXT),DIAM(MAXT), TOTLAI  !(MAXT)
    REAL CANOPYDIMS(6),DT1,DT2,DT3,DT4,EXPTIME,APP,EXPAN
    REAL STOCKING,THRESH_FOLT,THRESH_RX,THRESH_RY,THRESH_RZ
    REAL THRESH_ZBC,THRESH_LAI
         
! Interpolate to get daily values of crown dimensions
      CALL TREEINTERP(IDAY,ISTART,NOXDATES,DATESX,RXTABLE,NOTREES,RX)

      CALL TREEINTERP(IDAY,ISTART,NOYDATES,DATESY,RYTABLE,NOTREES,RY)

      CALL TREEINTERP(IDAY,ISTART,NOZDATES,DATESZ,RZTABLE,NOTREES,RZ)

      CALL TREEINTERP(IDAY,ISTART,NOTDATES,DATEST,ZBCTABLE,NOTREES,ZBC)

      IF (NOLADATES.GT.0) THEN
        CALL TREEINTERP(IDAY,ISTART,NOLADATES,DATESLA,FOLTABLE,NOTREES, &
          FOLT)
        CALL TREEINTERP2(IDAY,ISTART,NOLADATES,DATESLA,TOTLAITABLE,1,TOTLAI)
        !TOTLAI = TOTLAITMP(1)  ! Added, RAD Aug. 2010 to fix array dimensions.
      ELSE
        CALL PHENOL(IDAY,ISTART,IFLUSH,DT1,DT2,DT3,DT4, &
          EXPTIME,APP,EXPAN,STOCKING,NOTREES,FOLT,TOTLAI,NEWCANOPY)
      END IF

      IF (NODDATES.GT.0) &
        CALL TREEINTERP(IDAY,ISTART,NODDATES,DATESD,DIAMTABLE,NOTREES, &
        DIAM)

! BM 11/07/07 Fixed bug: we were not re-calculating points unless the DAILY change in the canopy
! was large. That allowed small changes to turn into a big change, without POINTS being recalculated.  
! Now, check the change in canopy size since POINTS were last calculated. 
! Thresholds for changes are set here: could potentially set elsewhere to make them easier to change? 
      NEWCANOPY = 0
      IF (IDAY.EQ.0) THEN
        NEWCANOPY = 1
      ELSE 
        THRESH_RX = RXTABLE(1,1)*0.01 ! 1% OF INITIAL CROWN RADIUS
        IF (RX(1)-CANOPYDIMS(1).GT.THRESH_RX) NEWCANOPY = 1
        THRESH_RY = RYTABLE(1,1)*0.01 
        IF (RY(1)-CANOPYDIMS(2).GT.THRESH_RY) NEWCANOPY = 1
        THRESH_RZ = RZTABLE(1,1)*0.01 
        IF (RZ(1)-CANOPYDIMS(3).GT.THRESH_RZ) NEWCANOPY = 1
        THRESH_ZBC = ZBCTABLE(1,1)*0.01 ! 
        IF (ZBC(1)-CANOPYDIMS(4).GT.THRESH_ZBC) NEWCANOPY = 1
        THRESH_FOLT = 0.05
        IF (FOLT(1)-CANOPYDIMS(5).GT.THRESH_FOLT) NEWCANOPY = 1
        IF (FOLT(1)-CANOPYDIMS(5).LT.-THRESH_FOLT) NEWCANOPY = 1
        THRESH_LAI = 0.1
        IF (TOTLAI-CANOPYDIMS(6).GT.THRESH_LAI) NEWCANOPY = 1
        IF (TOTLAI-CANOPYDIMS(6).LT.-THRESH_LAI) NEWCANOPY = 1
      END IF
      IF (NEWCANOPY.EQ.1) THEN
        CANOPYDIMS(1) = RX(1)
        CANOPYDIMS(2) = RY(1)
        CANOPYDIMS(3) = RZ(1)
        CANOPYDIMS(4) = ZBC(1)
        CANOPYDIMS(5) = FOLT(1)
        CANOPYDIMS(6) = TOTLAI
      END IF

      RETURN
      END !InterpolateT


!**********************************************************************
      SUBROUTINE OUTPUTLAY(UFILE,FOLLAY,JMAX25,VCMAX25,NOLAY)
! Daily output to layer flux file.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,I,NOLAY
    REAL FOLLAY(MAXLAY)
    REAL JMAX25(MAXLAY,MAXC),VCMAX25(MAXLAY,MAXC)

      WRITE(UFILE,*) 'LEAF AREA OF TARGET TREE IN EACH LAYER IN M2'
      WRITE(UFILE,500) (FOLLAY(I),I=1,NOLAY)
      WRITE(UFILE,*)
      WRITE(UFILE,*) 'JMAX (CURRENT) IN EACH LAYER IN UMOL M-2 S-1'
      WRITE(UFILE,500) (JMAX25(I,1),I=1,NOLAY)
      WRITE(UFILE,*)
      WRITE(UFILE,*) 'VCMAX (CURRENT) IN EACH LAYER IN UMOL M-2 S-1'
      WRITE(UFILE,500) (VCMAX25(I,1),I=1,NOLAY)
      WRITE(UFILE,*)

500   FORMAT(10(F10.5,1X))

      RETURN
      END !OutputLay


!**********************************************************************
      SUBROUTINE OUTPUTHIST(UFILE,ITREE,HISTO,BINSIZE)
! Write PAR histogram to file.
!**********************************************************************
    
    USE maestcom
    IMPLICIT NONE
    INTEGER UFILE,NOTARGETS,ITAR,IBIN,ITREE
    REAL HISTO(MAXHISTO)
    REAL BINSIZE

      WRITE (UFILE,500) ITREE
      WRITE (UFILE,510) 
      DO 10 IBIN = 1,MAXHISTO
        WRITE (UFILE,520) (IBIN-0.5)*BINSIZE,HISTO(IBIN)
10    CONTINUE
500   FORMAT ('TREE NO: ',I6)
510   FORMAT ('  PAR:         FREQUENCY (M^2.HR): ')
520   FORMAT (F8.2,1X,F12.6)

      RETURN
      END !OutputHist


!**********************************************************************
      INTEGER FUNCTION INEDGES(DX,DY,XMAX,YMAX,EDGEDIST)
! Checks to see if the chosen tree is within EDGEDIST m of the edge 
! of the plot. Returns -1 if yes and 1 if no.
!**********************************************************************
    IMPLICIT NONE
    REAL DX,DY,XMAX,YMAX,EDGEDIST
    
      IF ((DX.LT.EDGEDIST).OR.(DY.LT.EDGEDIST).OR. &
        (XMAX-DX.LT.EDGEDIST).OR.(YMAX-DY.LT.EDGEDIST)) THEN
        INEDGES = -1
      ELSE
        INEDGES = 1
      END IF

      RETURN
      END ! InEdges



!**********************************************************************
      SUBROUTINE GETPOINTSF(NUMTESTPNT,XL,YL,ZL,X0,Y0,XMAX,YMAX, &
        CTITLE,TTITLE,MTITLE,STITLE,VTITLE)
! Subroutine for testing radiation interception routines.
! Open input & output files and read information about sensor positions.
!**********************************************************************

      USE maestcom
      IMPLICIT NONE
      INTEGER NOPOINTS,INPUTTYPE,IOERROR,NUMTESTPNT,N,I
      
      CHARACTER*80 CTITLE, TTITLE, PTITLE, STITLE, MTITLE, VTITLE
      REAL XL(MAXP),YL(MAXP),ZL(MAXP),COORDS(MAXP*3)
      REAL X0,Y0,ANGLE,SPACING,ZHEIGHT,COSANG,SINANG,DIST
      REAL XMAX,YMAX

      NAMELIST /CONTROL/ NOPOINTS,INPUTTYPE
      NAMELIST /XYZ/ COORDS
      NAMELIST /TRANSECT/ ANGLE,SPACING,ZHEIGHT

! Open input file
      OPEN (UPOINTSI, FILE = 'points.dat', STATUS = 'OLD', &
         IOSTAT=IOERROR)
      
      IF (IOERROR.NE.0) &
        CALL SUBERROR('ERROR: POINTS INPUT FILE DOES NOT EXIST', &
        IFATAL,IOERROR)

! Read title from input file
990   FORMAT (A60)     ! For reading titles in input files.
      READ (UPOINTSI, 990) PTITLE

! Default values
      NOPOINTS = 0
      INPUTTYPE = 1

! Read control flags: no of points and type of input
      READ (UPOINTSI, CONTROL, IOSTAT = IOERROR)
      IF ((IOERROR.NE.0).OR.(NOPOINTS.EQ.0)) &
        CALL SUBERROR('ERROR: MISSING CONTROL INFO IN POINTS FILE', &
        IFATAL,IOERROR)
      IF (NOPOINTS.GT.MAXP) THEN
        CALL SUBERROR('WARNING: TOO MANY TEST POINTS SPECIFIED', &
        IWARN,IOERROR)
        NUMTESTPNT = MAXP
      ELSE 
        NUMTESTPNT = NOPOINTS
      END IF

! Read in list of points
      IF (INPUTTYPE.EQ.1) THEN
        READ (UPOINTSI, XYZ, IOSTAT = IOERROR)
        IF (IOERROR.NE.0) &
        CALL SUBERROR('ERROR READING GRID POINTS', &
        IFATAL,IOERROR)
        DO  N = 1,NUMTESTPNT
          XL(N) = COORDS((N-1)*3 + 1) - X0
          YL(N) = COORDS((N-1)*3 + 2) - Y0
          ZL(N) = COORDS(N*3)
        ENDDO

! Read in details of transect & construct points
      ELSE IF (INPUTTYPE.EQ.2) THEN
        READ (UPOINTSI, TRANSECT, IOSTAT = IOERROR)
        IF (IOERROR.NE.0) &
          CALL SUBERROR('ERROR READING TRANSECT DETAILS', &
          IFATAL,IOERROR)
        ANGLE = ANGLE*PID180
        COSANG = COS(ANGLE)
        SINANG = SIN(ANGLE)
        DIST = SPACING/2.0
        DO N = 1,NUMTESTPNT
          XL(N) = DIST*COSANG
          YL(N) = DIST*SINANG
          ZL(N) = ZHEIGHT
          DIST =  DIST + SPACING
        ENDDO
      END IF
  
      DO I = 1,NUMTESTPNT
         IF( (XL(I).LT.X0.OR.XL(I).GT.XMAX) .OR. (YL(I).LT.Y0.OR.YL(I).GT.YMAX) ) THEN
            CALL SUBERROR('WARNING: MAESTEST MAY CRASH WHEN POINTS ARE OUTSIDE PLOT BOUNDS. & 
                FIX POINTS.DAT IF INFINITE LOOP OCCURS!', IWARN, -1)
         ENDIF
      ENDDO
  
! Open output file
      OPEN (UPOINTSO, FILE = 'testflx.dat', STATUS = 'UNKNOWN')
! Write headings to output file
991   FORMAT (A12,A60) ! For writing comments to output files.
992   FORMAT (1X,3(A3,1X),9(A12,1X))
993   FORMAT (A60)
994   FORMAT (A60)

      WRITE (UPOINTSO, 991) 'Program:    ', VTITLE
      WRITE (UPOINTSO, 991) 'Control:    ', CTITLE
      WRITE (UPOINTSO, 991) 'Trees:      ', TTITLE
      WRITE (UPOINTSO, 991) 'Structure:  ', STITLE
      WRITE (UPOINTSO, 991) 'Points:     ', PTITLE
      WRITE (UPOINTSO, 991) 'Met data:   ', MTITLE
      WRITE (UPOINTSO, *)
      WRITE (UPOINTSO, 993) 'DAY: day number'
      WRITE (UPOINTSO, 993) 'HR: hour number'
      WRITE (UPOINTSO, 993) 'PT: point number'
      WRITE (UPOINTSO, 993) 'X,Y,Z, : coordinates of test point'
      WRITE (UPOINTSO, 993) 'PAR: incident PAR (umol m-2 s-1)'
      WRITE (UPOINTSO, 993) 'FBEAM: beam fraction of PAR'
      WRITE (UPOINTSO, 993)  &
        'SUNLA: sunlit leaf area at grid point (fraction)'
      WRITE (UPOINTSO, 993) &
         'TD: diffuse transmittance to grid point (fraction)'
      WRITE (UPOINTSO, 993) 'TSCAT: scattered radiation (umol m-2 s-1)'
      WRITE (UPOINTSO, 993) 'TTOT: total radiation (umol m-2 s-1)'
!      WRITE (UPOINTSO, 993) 'GSUN: gs for sunlit foliage (mol m-2 s-1)'
!      WRITE (UPOINTSO, 993) &
!            'GSHADE: gs for shaded foliage (mol m-2 s-1)'
!      WRITE (UPOINTSO, 993) &
!           'ASUN: A for sunlit foliage (mu mol m-2 s-1)'
!      WRITE (UPOINTSO, 993) &
!           'ASHADE: A for shaded foliage (mu mol m-2 s-1)'
!      WRITE (UPOINTSO, 993) 'RD: Dark respiration (mu mol m-2 s-1)'
!      WRITE (UPOINTSO, 993) 'ESUN: E for sunlit foliage (mmol m-2 s-1)'
!      WRITE (UPOINTSO, 993) &
!            'ESHADE: gs for sunlit foliage (mmol m-2 s-1)'

!      WRITE (UPOINTSO, 992) 'PT','BEAM','DIFF','SCAT','TOTAL'
!      WRITE (UPOINTSO,993) 'DAY','HR','PT','PAR','FBEAM','SUNLA','TD', &
!      'TSCAT','TTOT','GSUN','GSHADE','ASUN', &
!      'ASHADE','RD','ESUN','ESHADE'
       WRITE(UPOINTSO,993)' '
       WRITE (UPOINTSO,994)'DAY HR PT  X  Y  Z  PAR  FBEAM   SUNLA  TD   TSCAT  TTOT  ' 
!      'TSCAT','TTOT' 

!      WRITE (UPOINTSO, *)

      RETURN
      END ! GetPointsF


