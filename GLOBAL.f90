!****************************************************************************                                                              
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *  
!**************************************************************************** 
MODULE GLOBAL
  IMPLICIT NONE
           INTEGER , PARAMETER :: N=400
           DOUBLE PRECISION, PARAMETER :: H =0.1D0
           DOUBLE PRECISION :: Y(5,N), DY(5,N)
           DOUBLE PRECISION :: MU(N), T(N), DMU(N)
           DOUBLE PRECISION :: P, Q, R
           DOUBLE PRECISION :: S, V
           INTEGER :: I, J, K
           DOUBLE PRECISION, PARAMETER ::   PR=0.7D0, M =15.0D0 , PINF =101325.0D0
           DOUBLE PRECISION, PARAMETER :: GAMMAN =1.4D0, RCONST = 287.0D0, CP =1005.2D0
           DOUBLE PRECISION :: WMU, FMU 
           DOUBLE PRECISION :: TTOT, THETA0, TWALL, TINST , U0, TINF
           DOUBLE PRECISION :: SOS, UINF, H0E, H0, GWALL, RHOINF, RHO
           DOUBLE PRECISION :: TMU1, TTK1, UDASH, FUNCG , DDMU
           INTEGER :: ITER
           INTEGER :: MA
END MODULE GLOBAL

