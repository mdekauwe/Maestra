! PHYSIOL.FOR

!**********************************************************************
! This file contains all physiology calculations - photosynthesis,
! stomatal conductance, transpiration, respiration.
!
! The main subroutines (called externally) are:
! PSTRANSP - calls the other functions, does the leaf temp iteration calcn
! RESP - calculates maintenance respiration rate
! GRESP - calculates growth respiration rate
! CALCRMW - calculates stem maintenance respiration rate per unit mass
! CALCFBIOM - calculates foliar biomass from leaf area & SLA
! CALCWBIOM - calculates woody biomass from height & diameter
! For water balance: not finished
!   ETCAN - calculates canopy transpiration rate
!   PartitionPPT - partitions precip between drainage & canopy storage
!
! Subsidiary subroutines are
! 1. Photosynthesis
!   PHOTOSYN - calculates photosynthesis from the FvC model
!   GAMMAFN - calculates T dependence of Gamma*
!   KMFN - calculates T dependence of Km
!   ARRH - Arrhenius T dependence
!   JMAXTFN - calculates T dependence of Jmax
!   VCMAXTFN - calculates T dependence of Vcmax
! 2. Conductances  & transpiration
!   GBHFREE - calculates conductance to heat through free convection
!   GBHFORCED - calculates conductance to heat through forced convection
!   GRADIATION - calculates radiation conductance
!   GSJARVIS - calculates stomatal conductance using the Jarvis model
!   GBCAN - calculates boundary layer conductance of canopy
!   PENMON - implements the Penman-Monteith equation
!**********************************************************************


!**********************************************************************
      SUBROUTINE PSTRANSP( &
       RDFIPT,TUIPT,TDIPT,RNET,WIND, &
       PAR,TAIR,CA,RH,VPD,VMFD,PRESS,SOILMOIST, &
       JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
       THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
       MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0, &
       GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,SOILDATA,SWPEXP, &
       G0,D0L,GAMMA,G1,GK,WLEAF,NSIDES,ITERMAX, &
       GSC,ALEAF,RD,ET,FHEAT,TLEAF,GBH,DECOUP &
       )
! This subroutine calculates leaf photosynthesis and transpiration.
! These may be calculated by
! (1) assuming leaf temperature = air temperature, Cs = Ca and Ds = Da
! (2) using iterative scheme of Leuning et al (1995) (PC&E 18:1183-1200)
! to calculate leaf temp, Cs & Ca.
! Setting ITERMAX = 0 gives (1); ITERMAX > 0 (suggest 100) gives (2).
!**********************************************************************
    
    USE maestcom
    IMPLICIT NONE

    INTEGER MODELGS,SOILDATA,WSOILMETHOD,ITER
    INTEGER IECO,ITERMAX,NSIDES
    REAL JMAX25,I0,LHV,MINLEAFWP,KTOT,PSIL,K10F
    REAL TLEAF,TAIR,DLEAF,VPD,VMLEAF,VMFD,RHLEAF,RH,CS,CA
    REAL SLOPE,GRADN,RDFIPT,TUIPT,TDIPT,GBHU,PRESS,WIND
    REAL WLEAF,PAR,TMOVE,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC
    REAL DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP
    REAL TBELOW,GSREF,GSMIN,D0,VK1,VK2,VPD1,VPD2,VMDF0
    REAL GSJA,GSJB,T0,TREF,TMAX,SOILMOISTURE,EMAXLEAF
    REAL SMD1,SMD2,WC1,WC2,SWPEXP,FSOIL,G0,D0L,GAMMA,G1
    REAL GSC,ALEAF,RD,WEIGHTEDSWP,GBHF,GBH,GH,VMFD0,GBV,GSV,GV
    REAL ET,RNET,GBC,TDIFF,TLEAF1,FHEAT,ETEST
    REAL SOILMOIST,GK,EPSILON,DECOUP,GAMMAC
    
    REAL, EXTERNAL :: SATUR
    REAL, EXTERNAL :: GRADIATION
    REAL, EXTERNAL :: GBHFORCED
    REAL, EXTERNAL :: GBHFREE
    REAL, EXTERNAL :: PENMON
    REAL, EXTERNAL :: CALCRMW,GRESP
    
    !f2py  intent(in,out) :: ALEAF
    !f2py  intent(in,out) :: ET
    !f2py  intent(in,out) :: GBH
    !f2py  intent(in,out) :: DECOUP
    !f2py  intent(in,out) :: TLEAF
    
! Set initial values of leaf temp and surface CO2 & VPD
      TLEAF = TAIR
      DLEAF = VPD
      VMLEAF = VMFD
      RHLEAF = RH
      CS = CA

! Following calculations do not depend on TLEAF
! Latent heat of water vapour at air temperature (J mol-1)
      LHV = (H2OLV0 - 2.365E3 * TAIR) * H2OMW
! Const s in Penman-Monteith equation  (Pa K-1)
      SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1
! Radiation conductance (mol m-2 s-1)
      GRADN = GRADIATION(TAIR,RDFIPT,TUIPT,TDIPT)
! Boundary layer conductance for heat - single sided, forced convection
      GBHU = GBHFORCED(TAIR,PRESS,WIND,WLEAF)
      

!**********************************************************************
      ITER = 0  ! Counter for iterations - finding leaf temperature
100   CONTINUE  ! Return point for iterations

       CALL PHOTOSYN( &
       PAR,TLEAF,CS,RHLEAF,DLEAF,VMLEAF,SOILMOIST, &
       JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
       THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
       MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0, &
       GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,SOILDATA,SWPEXP, &
       G0,D0L,GAMMA,G1,GK, &
       GSC,ALEAF,RD)
      
! Boundary layer conductance for heat - single sided, free convection
      GBHF = GBHFREE(TAIR,TLEAF,PRESS,WLEAF)
! Total boundary layer conductance for heat
      GBH = GBHU + GBHF

! Total conductance for heat - two-sided
      GH = 2.*(GBH + GRADN)
! Total conductance for water vapour
      GBV = GBVGBH*GBH
      GSV = GSVGSC*GSC
!      GV = NSIDES*(GBV*GSV)/(GBV+GSV) ! already one-sided value
      GV = (GBV*GSV)/(GBV+GSV)

!  Call Penman-Monteith equation
      ET = PENMON(PRESS,SLOPE,LHV,RNET,VPD,GH,GV)
      !ET2 = (VPD/PRESS)*GSV
      !ET3 = 1E-3*VPD*GSV / (115.8 + 0.423*TAIR) 
      
! Calculate decoupling coefficient (McNaughton and Jarvis 1986)
      GAMMAC = CPAIR*AIRMA*PRESS/LHV
      EPSILON = SLOPE/GAMMAC
      DECOUP = (1+EPSILON) / (1 + EPSILON + GBV/GSV)

! End of subroutine if no iterations wanted.
      IF (ITERMAX.EQ.0.OR.ALEAF.LE.0.0) GOTO 200

! Otherwise, calculate new TLEAF, DLEAF, RHLEAF & CS
      GBC = GBH/GBHGBC
      CS = CA - ALEAF/GBC
      TDIFF = (RNET - ET*LHV) / (CPAIR * AIRMA * GH)
      TLEAF1 = TAIR + TDIFF
      DLEAF = ET * PRESS / GV
      RHLEAF = 1. - DLEAF/SATUR(TLEAF1)
      VMLEAF = DLEAF/PRESS*1E-3

! Check to see whether convergence achieved or failed
      IF (ABS(TLEAF - TLEAF1).LT.TOL) GOTO 200
      IF (ITER.GT.ITERMAX) THEN
        CALL SUBERROR('FAILED CONVERGENCE IN PSTRANSP',IWARN,0)
      GOTO 200
      END IF

! Update temperature & do another iteration
      TLEAF = TLEAF1
      ITER = ITER + 1
      GOTO 100

200   FHEAT = RNET - LHV*ET
!      FHEAT = (TLEAF - TAIR)*2.*GBH*CPAIR*AIRMA  !BM 12/05 Not correct - use energy bal instead
      ET = ET*1E6  ! Return ET,EI in umol m-2 s-1

      !WRITE(99,*)et,et2*1E6,tdiff


      RETURN
      END ! PsTransp


!**********************************************************************
      SUBROUTINE PHOTOSYN( &
       PAR,TLEAF,CS,RH,VPD,VMFD,SOILMOIST, &
       JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
       THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
       MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0, &
       GSJA,GSJB,T0,TREF,TMAX,SMD1,SMD2,SOILDATA,SWPEXP, &
       G0,D0L,GAMMA,G1,GK, &
       GS,ALEAF,RD)
! This subroutine calculates photosynthesis according to the ECOCRAFT
! agreed formulation of the Farquhar & von Caemmerer (1982) equations.
! Stomatal conductance may be calculated according to the Jarvis,
! Ball-Berry or BB-Leuning models.
! NB ALEAF is NET leaf photosynthesis. 
! NB The effect of soil water content is currently ignored - needs
! to be added as part of soil water balance routines.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER MODELGS,IQERROR,IECO,SOILDATA, WSOILMETHOD
    
    REAL PAR,TLEAF,CS,RH,VPD,VMFD,PSIL,KTOT
    REAL JMAX25,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN
    REAL THETA,AJQ,RD0,Q10F,K10F,RTEMP,DAYRESP,TBELOW
    REAL GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0
    REAL GSJA,GSJB,T0,TREF,TMAX
    REAL G0,D0L,GAMMA,G1
    REAL GS,ALEAF,RD,MINLEAFWP,RD0ACC
    REAL TMOVE,RESP,FSOIL,SOILMOISTURE,SMD1,SMD2,WC1,WC2,SWPEXP
    REAL VPDG,ETEST,WEIGHTEDSWP,EMAXLEAF,GSV

    REAL GAMMASTAR,KM,JMAX,VCMAX,J,VJ
    REAL A,B,C,AC,AJ,GSDIVA,CIC,CIJ
    REAL KMFN,JMAXTFN,SOILMOIST,GK
    REAL, EXTERNAL :: GAMMAFN
    REAL, EXTERNAL :: VCMAXTFN
    REAL, EXTERNAL :: QUADM
    REAL, EXTERNAL :: CALCFSOIL
    REAL, EXTERNAL :: QUADP
    REAL, EXTERNAL :: GSJARVIS
    
    !f2py  intent(in,out) :: ALEAF
   
! Calculate photosynthetic parameters from leaf temperature.
      GAMMASTAR = GAMMAFN(TLEAF,IECO)                   ! CO2 compensation point, umol mol-1
      KM = KMFN(TLEAF,IECO)                             ! Michaelis-Menten for Rubisco, umol mol-1
      JMAX = JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Potential electron transport rate, umol m-2 s-1
      VCMAX = VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Maximum Rubisco activity, umol m-2 s-1
      RD = RESP(RD0,TLEAF,Q10F,RTEMP,DAYRESP,TBELOW)    ! Day leaf respiration, umol m-2 s-1
      J = QUADM(THETA,-(AJQ*PAR+JMAX),AJQ*PAR*JMAX,IQERROR) ! Actual e- transport rate, umol m-2 s-1
      VJ = J/4.0                                        ! RuBP-regen rate, umol m-2 s-1

! Deal with extreme cases
    IF ((JMAX.LE.0.0).OR.(VCMAX.LE.0.0)) THEN
      ALEAF = -RD
        IF ((MODELGS.NE.2).AND.(MODELGS.NE.3).AND.(MODELGS.NE.4).AND. &
             (MODELGS.NE.5)) THEN          ! Minimum Gs value, mol m-2 s-1    
        GS = GSMIN
      ELSE
        GS = G0
      END IF
      RETURN
    END IF

! Calculation varies according to the stomatal model chosen.

! Jarvis model
      IF ((MODELGS.NE.2).AND.(MODELGS.NE.3).AND.(MODELGS.NE.4).AND. &
         (MODELGS.NE.5)) THEN
        GS = GSJARVIS(MODELGS,PAR,I0,VPD,D0,VK1,VK2,VPD1,VPD2,VMFD, &
           VMFD0,CS,GSJA,GSJB,SOILMOIST,SMD1,SMD2,SOILDATA,SWPEXP, &
           TLEAF,T0,TREF,TMAX,GSREF,GSMIN)

        IF (GS.LE.GSMIN) THEN
          ALEAF = -RD
          GS = GSMIN
          RETURN
        END IF

! Photosynthesis when Rubisco is limiting
        A = 1./GS
        B = (RD - VCMAX)/GS - CS - KM
        C = VCMAX * (CS - GAMMASTAR) - RD * (CS + KM)
        AC = QUADM(A,B,C,IQERROR)

        IF (IQERROR.EQ.1) THEN
          GS = GSMIN
          AC = - RD
        END IF

! Photosynthesis when electron transport is limiting
        B = (RD - VJ)/GS - CS - 2*GAMMASTAR
        C = VJ * (CS - GAMMASTAR) - RD * (CS + 2*GAMMASTAR)
        AJ = QUADM(A,B,C,IQERROR)

        IF (IQERROR.EQ.1) THEN
          GS = GSMIN
          AJ = - RD
        END IF

        ALEAF = AMIN1(AC,AJ)        ! Jarvis model solution

      ELSE
! Calculate soil moisture modifying factor
      FSOIL = CALCFSOIL(SOILMOIST,SOILDATA,SMD1,SMD2,SWPEXP)

      !write(uwattest,'(11F15.5)')km,jmax,vcmax,rd,j,gammastar,fsoil,rh,cs,g0,g1

        IF (MODELGS.EQ.2) THEN
! Ball-Berry model
          GSDIVA = G1 * RH / (CS - GAMMA) * FSOIL
        ELSE IF (MODELGS.EQ.3) THEN
! Ball-Berry-Leuning model
          GSDIVA = G1 / (CS - GAMMA) / (1 + VPD/D0L) * FSOIL
        ELSE IF (MODELGS.EQ.4) THEN
! Ball-Berry-Medlyn model
        IF(VPD.LT.50)THEN
            VPDG = 50.0/1000.0
        ELSE
            VPDG = VPD/1000.0
        ENDIF     
        GSDIVA = G1 / (CS - GAMMA) / SQRT(VPDG) * FSOIL
        
        ELSE IF (MODELGS.EQ.5) THEN
! Three-parameter Ball-Berry
        IF(VPD.LT.50)THEN
            VPDG = 50.0/1000.0
        ELSE
            VPDG =VPD/1000.0
        ENDIF     
        GSDIVA = G1 / (CS - GAMMA) / VPDG**GK
        
        END IF

! Following calculations are used for both BB & BBL models.
! Solution when Rubisco activity is limiting

        A = G0 + GSDIVA * (VCMAX - RD)
        B = (1. - CS*GSDIVA) * (VCMAX - RD) + G0 * (KM - CS) &
            - GSDIVA * (VCMAX*GAMMASTAR + KM*RD)
        C = -(1. - CS*GSDIVA) * (VCMAX*GAMMASTAR + KM*RD) - G0*KM*CS
        CIC = QUADP(A,B,C,IQERROR)

        !write(uwattest,*)CS
        !write(uwattest,'(13F15.5)')gsdiva,cic,fsoil,g0,g1,vcmax,km,rd,gammastar,gamma,vpdg,cs

        IF ((IQERROR.EQ.1).OR.(CIC.LE.0.0).OR.(CIC.GT.CS)) THEN
          AC = 0.0
        ELSE
          AC = VCMAX * (CIC - GAMMASTAR) / (CIC + KM)
        END IF

! Solution when electron transport rate is limiting
        A = G0 + GSDIVA * (VJ - RD)
        B = (1. - CS*GSDIVA) * (VJ - RD) + G0 * (2.*GAMMASTAR - CS) &
            - GSDIVA * (VJ*GAMMASTAR + 2.*GAMMASTAR*RD)
        C = -(1. - CS*GSDIVA) * GAMMASTAR * (VJ + 2.*RD) &
            - G0*2.*GAMMASTAR*CS
        CIJ = QUADP(A,B,C,IQERROR)

        AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
      IF (AJ-RD.LT.1E-6) THEN      ! Below light compensation point
        CIJ = CS
          AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        END IF

        ALEAF = AMIN1(AC,AJ) - RD  ! Solution for Ball-Berry model
        GS = G0 + GSDIVA*ALEAF
        IF (GS.LT.G0) GS = G0

      END IF

      RETURN
      END !Photosyn

!**********************************************************************
      REAL FUNCTION GAMMAFN(TLEAF,IECO)
! This subroutine calculates Gamma(star), or the CO2 compensation point
! in the absence of non-photorespiratory respiration.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IECO
    REAL TLEAF
    REAL, EXTERNAL :: ARRH

      IF (IECO.EQ.1) THEN 
! Ecocraft fomulation; based on Brooks & Farquhar and von Caemmerer et al. 
! If TLEAF < -1.0 then calculate Gamma for T = -1 (quadratic not applicable)
        IF (TLEAF.LT.-1.0) THEN
          GAMMAFN = 36.9 + 1.88*(-26.0) + 0.036*(-26.0)*(-26.0)
        ELSE
          GAMMAFN = 36.9 + 1.88*(TLEAF-25) + 0.036*(TLEAF-25)*(TLEAF-25)
        END IF

    ELSE    ! Bernacchi et al 2001 PCE 24: 253-260
      GAMMAFN = ARRH(42.75,37830.0,TLEAF,25.0)
    ENDIF

      RETURN
      END !Gammafn

!**********************************************************************
      REAL FUNCTION KMFN(TLEAF,IECO)
! This subroutine calculates Km, or the effective Michaelis-Menten
! coefficient of Rubisco activity.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IECO
    REAL OI,KC25,KO25,KCEA,KOEA,KC,KO, TLEAF
    REAL, EXTERNAL :: ARRH
    
      OI = 205000         ! Oxygen partial pressure (umol mol-1)
      IF (IECO.EQ.1) THEN
! Physiological constants - values agreed by Ecocraft - Badger & Collatz values
        KC25 = 404          ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 248000       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 59400        ! Temp. response of Kc (J mol-1)
        KOEA = 36000        ! Temp. response of Ko (J mol-1)
    ELSE !  Bernacchi et al 2001 PCE 24: 253-260
        KC25 = 404.9        ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 278400       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 79430        ! Temp. response of Kc (J mol-1)
        KOEA = 36380        ! Temp. response of Ko (J mol-1)
    END IF

! This function is well-behaved for TLEAF < 0.0
      KC = ARRH(KC25,KCEA,TLEAF,25.0)
      KO = ARRH(KO25,KOEA,TLEAF,25.0)
        KMFN = KC * (1. + OI/KO)

      RETURN
      END !KmFn


!**********************************************************************
      REAL FUNCTION JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)
! This subroutine calculates the potential electron transport rate
! (Jmax) at the leaf temperature.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL JMAX25,TLEAF,EAVJ,EDVJ,DELSJ, TVJUP, TVJDN
    REAL, EXTERNAL :: TK
    
! This function is well-behaved for TLEAF < 0.0
       JMAXTFN = JMAX25 * EXP((TLEAF-25)*EAVJ/(RCONST*TK(TLEAF)*TK(25.))) &
       * (1.+EXP((DELSJ*TK(25.)-EDVJ)/(RCONST*TK(25.)))) &
       / (1.+EXP((DELSJ*TK(TLEAF)-EDVJ)/(RCONST*TK(TLEAF))))

! Function allowing Vcmax to be forced linearly to zero at low T - 
! introduced for Duke data
    IF (TLEAF.LT.TVJDN) THEN
      JMAXTFN = 0.0 
      ELSE IF (TLEAF.LT.TVJUP) THEN
      JMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*JMAXTFN
    END IF
     
      RETURN
      END !JmaxTFn

!**********************************************************************
      REAL FUNCTION VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)
! This subroutine calculates the maximum Rubisco activity
! (Vcmax) at the leaf temperature.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL VCMAX25,TLEAF,EAVC,EDVC,DELSC
    REAL TVJDN,TVJUP
    REAL, EXTERNAL :: TK
    
! There is still disagreement as to whether this function has an
! optimum or not. Both forms are available here. If EDVC <= 0 (default 0)
! then the no-optimum form is used.
! Both functions are well-behaved for TLEAF < 0.0

      IF (EDVC.LE.0.0) THEN
         VCMAXTFN = VCMAX25 * EXP(EAVC*(TLEAF - 25) &
        / (TK(25.)*RCONST*TK(TLEAF)))
      ELSE
        VCMAXTFN = VCMAX25 &
       * EXP((TLEAF-25.)*EAVC/(RCONST*TK(TLEAF)*TK(25.))) &
       * (1.+EXP((DELSC*TK(25.)-EDVC)/(RCONST*TK(25.)))) &
       / (1.+EXP((DELSC*TK(TLEAF)-EDVC)/(RCONST*TK(TLEAF))))
      END IF

! Function allowing Vcmax to be forced linearly to zero at low T - 
! introduced for Duke data
    IF (TLEAF.LT.TVJDN) THEN
      VCMAXTFN = 0.0 
      ELSE IF (TLEAF.LT.TVJUP) THEN
      VCMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*VCMAXTFN
    END IF
     
      RETURN
      END !VcmaxTFn


!**********************************************************************
      REAL FUNCTION RESP(RD0,TLEAF,Q10F,RTEMP,DAYRESP,TBELOW)
! This function calculates respiration from temperature
! using a Q10 (exponential) formulation.
!**********************************************************************

      IMPLICIT NONE      
      REAL RD0,TLEAF,Q10F,RTEMP,DAYRESP,K10F,RD0ACC,TMOVE,TBELOW

      IF (TLEAF.GE.TBELOW) THEN
        RESP = RD0 * EXP(Q10F * (TLEAF-RTEMP)) * DAYRESP
      ELSE
        RESP = 0.0
      END IF

    
      RETURN
      END !Resp

!**********************************************************************
      REAL FUNCTION ARRH(KT,EA,T,TREF)
! The Arrhenius function.
! KT is the value at Tref deg !; Ea the activation energy (J mol-1) and T the temp (deg !). 
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL KT, EA, T, TREF
    
         ARRH = KT*EXP(EA*(T-TREF)/(RCONST*(T-ABSZERO)*(TREF-ABSZERO)))
      RETURN
      END !Arrh


!**********************************************************************
       REAL FUNCTION GSJARVIS( &
       MODELGS,PAR,I0,VPD,D0,VK1,VK2,VPD1,VPD2,VMFD,VMFD0,CS,GSJA,GSJB, &
       SOILMOIST,SMD1,SMD2,SOILDATA,SWPEXP, &
       TLEAF,T0,TREF,TMAX,GSREF,GSMIN &
       )
! Calculate stomatal conductance according to the Jarvis model.
! This model calculates gs by multiplying together functions of several
! environmental variables: light, VPD, CO2 & temperature.
!**********************************************************************
      
      USE maestcom
      IMPLICIT NONE
      INTEGER MODELGS,SOILDATA
      REAL PAR,I0,VPD,D0,VK1,VK2,VMFD,VMFD0,CS,GSJA,GSJB
      REAL TLEAF,T0,TREF,TMAX,GSREF,GSMIN,P,SOILMOIST,SMD1,SMD2
      REAL FLIGHT,FVPD,FCO2,FTEMP,FSOIL,VPD1,VPD2,SWPEXP
      REAL, EXTERNAL :: CALCFSOIL

! Defaults for the different factors
      FLIGHT = 1.0
      FVPD = 1.0
      FCO2 = 1.0
      FTEMP = 1.0
      FSOIL = 1.0

! Response to incident radiation - PAR & I0 in umol m-2 s-1
      IF (I0.GT.0.0) FLIGHT = PAR / (PAR + I0)

! Hyperbolic decline with VPD (VPD in Pa) * BRAY (ALEX BOSC)
! or Lohammer response to VPD (VPD & D0 in Pa) * ECOCRAFT
! or mole fraction deficit (VMFD in mmol mol-1) * MARK RAYMENT
! or linear decline with VPD (VPD1, VPD2 in Pa) * GROMIT (TIM RANDLE)
      IF (MOD(MODELGS,100).GE.30) THEN
      IF (VPD.GT.0.0) FVPD = 1./(VK1*VPD**VK2)
    ELSE IF (MOD(MODELGS,100).GE.20) THEN
        IF (VPD.GE.VPD1) FVPD = 1. - (VPD - VPD1)/(VPD2 - VPD1)
      ELSE IF (MOD(MODELGS,100).GE.10) THEN
        IF (VMFD0.GT.0.0) FVPD = 1.0 - VMFD/VMFD0
      ELSE
        IF (D0.GT.0.0) FVPD = 1. / (1. + VPD/D0) 
      END IF
    IF (FVPD.GT.1.0) FVPD = 1.0
      IF (FVPD.LT.0.0) FVPD = 0.0

! Two possible responses to CO2 - linear or non-linear. CO2 in umol mol-1.
      IF (MOD(MODELGS,10).EQ.0) THEN
        IF (GSJA.NE.0.0) FCO2 = 1 - GSJA * (CS - 350.0)
      ELSE
        IF (GSJB.NE.0.0) FCO2 = (GSJB + 350.0) / (GSJB + CS)
      END IF
      IF (FCO2.LT.0.0) FCO2 = 0.0

! Temperature function is optional. To not use it, set TMAX <= 0.0. (Default).
      IF (TMAX.LE.0.0) THEN
        FTEMP = 1.0
      ELSE IF (MOD(MODELGS,1000).GE.100) THEN
        P = (TMAX - TREF)/(TREF - T0)
        FTEMP = (TLEAF - T0)*((TMAX - TLEAF)**P) /  &
              ((TREF - T0)*((TMAX - TREF)**P))
      ELSE
        FTEMP = (TLEAF - T0) * (2*TMAX - T0 - TLEAF) / &
                ((TREF - T0) * (2*TMAX - T0 - TREF))
      END IF      
      IF (FTEMP.LT.0.0) FTEMP = 0.0

! Soil moisture function 
      FSOIL = CALCFSOIL(SOILMOIST,SOILDATA,SMD1,SMD2,SWPEXP)

! Multiply factors together.
! Gsref is in mol m-2 s-1; Gsjarvis required in mol m-2 s-1.
      GSJARVIS = (GSREF-GSMIN)*FLIGHT*FVPD*FCO2*FTEMP*FSOIL + GSMIN

      RETURN
      END !GsJarvis


!**********************************************************************
      REAL FUNCTION ETCAN( &
       WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,RNET,VPD,GSCAN,STOCKING &
     )
! Calculate transpiration by applying Penman-Monteith to whole canopy.
! Returns umol m-2 s-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL LHV,WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,RNET,VPD,GSCAN,STOCKING
    REAL GB,GSV,RNETM2,SLOPE,GH,GV
    REAL, EXTERNAL :: GBCAN
    REAL, EXTERNAL :: HEATEVAP
    REAL, EXTERNAL :: SATUR
    REAL, EXTERNAL :: PENMON
    
! Get boundary layer conductance
      GB = GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR)
      GSV = GSCAN*GSVGSC*STOCKING !convert mol CO2/tree/s to mol H2O/m2/s
      RNETM2 = RNET*STOCKING 
      
      IF (GB*GSV.GT.0.0) THEN 
! Latent heat of water vapour at air temperature (J mol-1)
        LHV = (H2OLV0 - 2.365E3 * TAIR) * H2OMW
       
! Const s in Penman-Monteith equation  (Pa K-1)
        SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1
! Call Penman-Monteith
        GH = GB
        GV = 1./(1./GSV + 1./GB)
        ETCAN = PENMON(PRESS,SLOPE,LHV,RNETM2,VPD,GH,GV)*1E6
      ELSE
        ETCAN = 0.0
      END IF

      RETURN
      END !ETCan


!**********************************************************************
      REAL FUNCTION PENMON( &
       PRESS,SLOPE,LHV,RNET,VPD,GH,GV &
     )
! This subroutine calculates evapotranspiration by leaves using the Penman-Monteith equation. 
! Inputs:    PRESS atmospheric pressure, Pa
!        SLOPE slope of VPD/T curve, Pa K-1
!        LHV latent heat of water at air T, J mol-1
!        RNET net radiation, J m-2 s-1
!        VPD vapour pressure deficit of air, Pa
!        GH boundary layer conductance to heat (free & forced & radiative components), mol m-2 s-1
!        GV conductance to water vapour (stomatal & bdry layer components), mol m-2 s-1
! Result in mol H2O m-2 s-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL PRESS,SLOPE,LHV,RNET,VPD,GH,GV,GAMMA
    REAL ET
    
      GAMMA = CPAIR*AIRMA*PRESS/LHV

      IF (GV.GT.0.0) THEN
        ET = (SLOPE * RNET + VPD * GH * CPAIR * AIRMA) / &
              (SLOPE + GAMMA * GH/GV)
      ELSE
        ET = 0.0
      END IF
      PENMON = ET / LHV
!    IF (PENMON.LT.0.0) PENMON = 0.0        ! BM 12/05 Should not be negative

      RETURN
      END !PenMon


!**********************************************************************
      REAL FUNCTION GRADIATION(TAIR,RDFIPT,TUIPT,TDIPT)
! Returns the 'radiation conductance' at given temperature.
! Formula from Ying-Ping's version of Maestro.
! See also Jones (1992) p. 108.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL TAIR,RDFIPT,TDIPT,TUIPT
    REAL, EXTERNAL :: TK

      GRADIATION = 4.*SIGMA*(TK(TAIR)**3.) &
        * RDFIPT/TDIPT * EMLEAF * (TDIPT + TUIPT) / (CPAIR * AIRMA)

      RETURN
      END !GRadiation


!**********************************************************************
      REAL FUNCTION GBHFORCED(TAIR,PRESS,WIND,WLEAF)
! Boundary layer conductance for heat - single sided, forced convection
! in mol m-2 s-1
! See Leuning et al (1995) PC&E 18:1183-1200 Eqn E1
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL TAIR,PRESS,WIND,WLEAF,CMOLAR
    REAL, EXTERNAL :: TK

      CMOLAR = PRESS / (RCONST * TK(TAIR))
      GBHFORCED = 0.003 * SQRT(WIND/WLEAF) * CMOLAR

      RETURN
      END !GBHForced


!**********************************************************************
      REAL FUNCTION GBHFREE(TAIR,TLEAF,PRESS,WLEAF)
! Boundary layer conductance for heat - single sided, free convection
! in mol m-2 s-1
! See Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4
!**********************************************************************
    
    USE maestcom
    IMPLICIT NONE
    REAL CMOLAR,PRESS,TAIR,TLEAF,GRASHOF,WLEAF
    REAL, EXTERNAL :: TK

      CMOLAR = PRESS / (RCONST * TK(TAIR))
      IF ((TLEAF-TAIR).NE.0.0) THEN
        GRASHOF = 1.6E8 * ABS(TLEAF-TAIR) * (WLEAF**3.) ! Grashof number
        GBHFREE = 0.5 * DHEAT * (GRASHOF**0.25) / WLEAF * CMOLAR
      ELSE
        GBHFREE = 0.0
      END IF

      RETURN
      END !GBHFree


!**********************************************************************
      REAL FUNCTION GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR)
! Canopy boundary layer conductance (from Jones 1992 p 68) 
! in mol m-2 s-1
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,CMOLAR
    REAL, EXTERNAL :: TK
    
          IF (Z0HT.GT.0.0) THEN

! Formula from Jones 1992 p 68
        GBCAN = WIND*(VONKARMAN**2)/(LOG((ZHT - ZPD)/Z0HT))**2
        
! Convert from mm s-1 to mol m-2 s-1
        CMOLAR = PRESS / (RCONST * TK(TAIR))
        GBCAN = GBCAN * CMOLAR 
        
      ELSE 
        GBCAN = 0.0
      END IF

      RETURN
      END !GBCan


!**********************************************************************
      SUBROUTINE CALCWBIOM(IDAY,HT,DIAM,COEFFT,EXPONT,WINTERC, &
        WBIOM,WBINC)
! Calculate the woody biomass (kg DW) on the given day from the height
! (m) and diameter (m). Also calculate the increment in woody biomass
! since previous day (g DW). Needed to calculate woody respiration.
!**********************************************************************
    IMPLICIT NONE
    INTEGER IDAY
    REAL PREVWBIOM,HT,DIAM,COEFFT,EXPONT,WINTERC,WBIOM,WBINC
    

      PREVWBIOM = WBIOM
      WBIOM = COEFFT*HT*(DIAM**EXPONT) + WINTERC
      IF (IDAY.EQ.0) PREVWBIOM = WBIOM
      WBINC = (WBIOM - PREVWBIOM)*1E3

      RETURN
      END ! CalcWBiom

!**********************************************************************
      SUBROUTINE CALCFBIOM(IDAY,NOSLADATES,FOLLAY,SLA,PROP,NOLAY,NOAGEP, &
        FBIOM,FBINC)
! Calculate foliage biomass from SLA and leaf area - done in layers.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER NOSLADATES,IDAY,NOLAY,NOAGEP,I,J
    REAL FOLLAY(MAXLAY)
    REAL SLA(MAXLAY,MAXC)
    REAL PROP(MAXC)
    REAL FBIOM,FBINC,PREVFBIOM

      IF (NOSLADATES.GT.0) THEN
        PREVFBIOM = FBIOM
        FBIOM = 0.0
        DO 10 I = 1,NOLAY
          DO 10 J = 1,NOAGEP
            FBIOM = FBIOM + FOLLAY(I)*PROP(J)/SLA(I,J)
10      CONTINUE
        IF (IDAY.EQ.0) PREVFBIOM = FBIOM
        FBINC = (FBIOM - PREVFBIOM)*1E3
      ELSE
        FBIOM = 0.0
        FBINC = 0.0
      END IF

      RETURN
      END !CalcFBiom


!**********************************************************************
    REAL  FUNCTION CALCRMW(MODELRW,COLLA,COLLK,STEMSDW, &
        DIAM,HT,STEMFORM,RMWAREA,WBIOM,RMW)
! Calculate stem respiration rate per unit biomass if necessary
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER MODELRW
    REAL COLLA,COLLK,STEMSDW,DIAM,HT,STEMFORM,RMWAREA,WBIOM,RMW
    REAL STEMAREA

      IF (MODELRW.EQ.1) THEN
        RMW = COLLA*EXP(COLLK*DIAM*100.0)*STEMSDW
      ELSE IF (MODELRW.EQ.2) THEN
      STEMAREA = STEMFORM*PI*(DIAM**2)*HT
      RMW = RMWAREA*STEMAREA/WBIOM
      END IF

      CALCRMW = RMW
    RETURN
    END !CALCRMW


!**********************************************************************
    REAL  FUNCTION GRESP(BINC,EFFY)
! Calculate the growth respiration. Use increment in biomass
! (g DW tree-1 d-1) and the efficiency of growth (g g-1 !).
! Returns a value in mol CO2 tree-1 d-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL BINC,EFFY

      IF (BINC.GT.0.0) THEN
        GRESP = BINC * EFFY / GCPERMOL * CPERDW
      ELSE
        GRESP = 0.0
      ENDIF

      RETURN
      END ! GResp

!**********************************************************************
     REAL FUNCTION CALCFSOIL(SOILMOIST,SOILDATA,SMD1,SMD2,SWPEXP)
! Calculate the effect of soil water content on stomatal conductance
! Two alternative forms:
! Granier & Loustau 1994 Fs = 1-a exp(b SMD) where SMD is soil moisture deficit, dimnless
! Otherwise negative exponential function of soil water potential
!**********************************************************************

    USE maestcom
    USE metcom
    IMPLICIT NONE
    INTEGER SOILDATA, WSOILMETHOD, SETFSOIL, IOERROR
    REAL SMD1,SMD2,WC1,WC2,SWPEXP,SOILMOISTURE,SOILMOIST

      CALCFSOIL = 1.0

! Exponential relationship with potential: parameter = SWPEXP
      IF (SOILDATA.EQ.POTENTIAL) THEN
      IF (SWPEXP.GT.0.0) THEN
        CALCFSOIL = EXP(SWPEXP*SOILMOIST)
      END IF
! Exponential relationship with deficit: params SMD1,SMD2 
      ELSE IF (SOILDATA.EQ.DEFICIT) THEN
        IF (SMD1.GT.0.0) THEN
        CALCFSOIL = 1.0 - SMD1*EXP(SMD2*SOILMOIST)
! Linear decline with increasing deficit: pUT SMD1 < 0, param SMD2
      ELSE IF (SMD2.GT.0.0) THEN
        IF ((1.0-SOILMOIST).LT.SMD2) THEN
          CALCFSOIL = (1.0-SOILMOIST)/SMD2
        END IF
      END IF
    END IF
      
    IF (CALCFSOIL.LT.0.0) CALCFSOIL = 0.0
     
    RETURN
    END ! Fsoil


