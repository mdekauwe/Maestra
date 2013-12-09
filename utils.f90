!**********************************************************************
! UTILS.FOR
!
! This file contains all the 'utility' functions which are called throughout
! the other files. I have set the files up so that, when writing a test program
! which calls one of the program subroutines in file xxx.for, the program
! can be compiled by compiling test.for, xxx.for, and utils.for. 
!
! The functions are:
! IDATE50 - translates string-format date in days-since-1950 format
! IDATEYR - given string-format date, finds Julian day
! JDATE - given days-since-1950-format date, finds Julian day
! SUBERROR - error-handling subroutine
! BETA - implements beta function
! QUADM, QUADP - solve quadratic
! AVERAGEVAL - find mean of array of values
! NIGHT - whether an hour qualifies as night-time or not
!**********************************************************************

!**********************************************************************
      FUNCTION IDATE50(STRDATE)
! This function translates a string-format date (DD/MM/YY) into the number 
! of days since 1/1/1950. Provides a numerical way of handling dates.
!**********************************************************************

    IMPLICIT NONE
    CHARACTER*8 STRDATE
    INTEGER IDAY,IYEAR,IYRD,IDATE50
    INTEGER, EXTERNAL :: IDATEYR

! Get the day of year and the year number from function IDATEYR.
      IDAY = IDATEYR(STRDATE,IYEAR)

! Calculate how many days in full years have passed since 1950.
      IYRD = 365*(IYEAR - 50)
      IYRD = IYRD + (IYEAR - 49)/4
      IDATE50 = IYRD + IDAY

      RETURN
      END !Idate50


!**********************************************************************
     INTEGER FUNCTION IDATEYR(STRDATE,IYEAR)
! Given a date strdate, this function calculates the number of the day 
! from 1 on 1st Jan to 365/366 on 31st Dec.
!**********************************************************************

    IMPLICIT NONE
    CHARACTER*8 STRDATE
    LOGICAL LEAPYR
    INTEGER IFD(12),IDAY,IMON,IYEAR
    DATA IFD/0,31,59,90,120,151,181,212,243,273,304,334/

      READ (STRDATE,10) IDAY,IMON,IYEAR
10    FORMAT(T1,I2,T4,I2,T7,I2)
! Code to handle passage to 2000: note this will only work until 2050
! I hope this model has been superseded by then!
    IF (IYEAR.LT.50) IYEAR = IYEAR+100

      LEAPYR = .FALSE.
      IF (4.0* (IYEAR/4).EQ.IYEAR) LEAPYR = .TRUE.

      IDATEYR = IFD(IMON) + IDAY
      IF (LEAPYR .AND. (IMON.GE.3)) IDATEYR = IDATEYR + 1

      RETURN
      END !Idateyr


!**********************************************************************
      INTEGER FUNCTION MONTH(IDATE50)
! This function returns the month given a date in days-since-1950 format.
!**********************************************************************
    IMPLICIT NONE
    INTEGER JD,IDATE50
    INTEGER, EXTERNAL :: JDATE
    
    JD = JDATE(IDATE50)
    MONTH = JD/30 + 1
    IF (MONTH.GT.12) MONTH = 12

    RETURN
    END !Mon


!**********************************************************************
      INTEGER FUNCTION JDATE(IDATE50)
! This function returns the julian day given a date in days-since-1950
! format.
!**********************************************************************
    IMPLICIT NONE
    INTEGER IYR,IYRD,IDATE50
    
      IYR = IDATE50 / 365
      IYRD = 365*IYR
      IYRD = IYRD + (IYR - 1)/4 
      JDATE = IDATE50 - IYRD + 1

      RETURN
      END !Jdate


!**********************************************************************
      SUBROUTINE SUBERROR(MESSAGE,IFLAG,IOERROR)
! The error-handling subroutine. When called, writes MESSAGE to the
! error file (UERROR). Uses IFLAG to determine whether to terminate program
! or not. IOERROR is a FORTRAN error number, reported in the error file. 
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    CHARACTER MESSAGE *(*)
    INTEGER IFLAG,IOERROR

      WRITE(UERROR,*) MESSAGE
      IF (IOERROR.GT.0) WRITE(UERROR,20) IOERROR
!10    FORMAT (A80)
20    FORMAT ('FORTRAN ERROR CODE NO: ',I10)

      IF (IFLAG.EQ.IFATAL) THEN
        STOP 'FATAL ERROR - SEE MAESERR.DAT FOR DETAILS.'
      END IF
      
      RETURN
      END !SubError


!**********************************************************************
     REAL FUNCTION BETA(BC1,BC2,BC3,BC4,RP)
! This a function to define the vertical and horizontal leaf area
! density distribution.
! BM 11/99 added fourth parameter BC4
!**********************************************************************
    IMPLICIT NONE
    REAL BC1,BC2,BC3,BC4,RP
    
      IF (RP.EQ.0.0) RP=0.0000001
      IF (RP.EQ.1.0) RP=0.9999999
    IF (RP.GE.BC4) THEN
      BETA = 0.0
    ELSE
        BETA = BC1* (RP**BC2)* ((BC4-RP)**BC3)
    END IF

      RETURN
      END !Beta


!**********************************************************************
      REAL FUNCTION QUADM(A,B,C,IQERROR)
! Solves the quadratic equation - finds smaller root.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL A, B, C
    INTEGER IQERROR
    
      IQERROR = 0

      IF ((B*B - 4.*A*C).LT.0.0) THEN
        CALL SUBERROR('WARNING:IMAGINARY ROOTS IN QUADRATIC', &
       IWARN,0)
        IQERROR = 1
        QUADM = 0.0
      ELSE
        IF (A.EQ.0.0) THEN
          IF (B.EQ.0.0) THEN
            QUADM = 0.0
            IF (C.NE.0.0) & 
             CALL SUBERROR('ERROR: CANT SOLVE QUADRATIC',IWARN,0)
          ELSE
            QUADM = -C/B
          END IF
        ELSE
          QUADM = (- B - SQRT(B*B - 4*A*C)) / (2.*A)
        END IF
      END IF

      RETURN
      END !QuadM


!**********************************************************************
      REAL FUNCTION QUADP(A,B,C,IQERROR)
! Solves the quadratic equation - finds larger root.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL A,B,C
    INTEGER IQERROR
   
      IQERROR = 0

      IF ((B*B - 4.*A*C).LT.0.0) THEN
        CALL SUBERROR('WARNING:IMAGINARY ROOTS IN QUADRATIC',IWARN,0)
        IQERROR = 1
        QUADP = 0.0
      ELSE
        IF (A.EQ.0.0) THEN
          IF (B.EQ.0.0) THEN
            QUADP = 0.0
            IF (C.NE.0.0) &
             CALL SUBERROR('ERROR: CANT SOLVE QUADRATIC',IWARN,0)
          ELSE
            QUADP = -C/B
          END IF
        ELSE
          QUADP = (- B + SQRT(B*B - 4*A*C)) / (2.*A)
        END IF
      END IF

      RETURN
      END !QuadP

         
!**********************************************************************
      REAL FUNCTION AVERAGEVAL(ARR,NUMVAL)
! Finds mean value of an array of NUMVAL values.
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
    INTEGER NUMVAL,I
    INTEGER, PARAMETER :: MAXP = 10000
    REAL ARR(MAXP)!,ARRSUM

      AVERAGEVAL = 0.0
    DO 10 I = 1,NUMVAL
      AVERAGEVAL = AVERAGEVAL + ARR(I)
10    CONTINUE
      AVERAGEVAL = AVERAGEVAL/REAL(NUMVAL)

    RETURN
    END ! AVERAGEVal


!**********************************************************************
      INTEGER FUNCTION NIGHT(ZEN,PAR)
! Function to determine whether a particular hour is at night or not
!**********************************************************************
    IMPLICIT NONE
    REAL ZEN,PAR
 
      IF ((ABS(ZEN).LE.1.57).AND.(PAR.GT.0.1)) THEN
        NIGHT = 0
      ELSE
        NIGHT = 1    
      END IF
    
      RETURN
      END !Night

!**********************************************************************
      REAL FUNCTION SATUR(TAC)
! Calculate saturated water vapour pressure (Pa) at temperature TAC (Celsius)
! from Jones 1992 p 110 (note error in a - wrong units)
!**********************************************************************
      IMPLICIT NONE
      REAL TAC
      
      SATUR = 613.75*EXP(17.502*TAC/(240.97+TAC))

      RETURN
      END !Satur

!**********************************************************************
    REAL FUNCTION TK(TCELSIUS)
! Converts Celsius temperature to Kelvin.
!**********************************************************************

      USE maestcom
    IMPLICIT NONE
    REAL TCELSIUS
    
      TK = TCELSIUS - ABSZERO
      RETURN
      END !TCelsius

!**********************************************************************
REAL FUNCTION STDEV(ARR,N)
!**********************************************************************
    IMPLICIT NONE      
    INTEGER N,I
    REAL ARR(N)
    REAL ARRMEAN2,ARR2MEAN

    ARRMEAN2 = (SUM(ARR(1:N))/REAL(N))**2
    ARR2MEAN = 0.0
    DO I = 1,N
        ARR2MEAN = ARR2MEAN + ARR(I)**2 / REAL(N)
    END DO

    ! May happen due to floating point inaccuracies.
    IF(ARRMEAN2.GT.ARR2MEAN)THEN
        STDEV = 0.0
    ELSE
        STDEV = SQRT(ARR2MEAN - ARRMEAN2)
    ENDIF

    RETURN
END FUNCTION STDEV

!**********************************************************************
INTEGER FUNCTION DIVFOURNEAR(NUM)
! Finds the nearest smaller number that can be divided by 4.
! RAD, March 2009.      
!********************************************************************
    
    IMPLICIT NONE
    INTEGER NUM,NEWNUM
    
    NEWNUM = NUM
    DO WHILE (MOD(NEWNUM,4).GT.0)
        NEWNUM = NEWNUM - 1
    END DO

    DIVFOURNEAR = NEWNUM

    RETURN
END FUNCTION DIVFOURNEAR
      


