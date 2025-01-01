!****************************************************************************                                                              
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *  
!****************************************************************************
      PROGRAM COMPRESSIBLE_BL
      USE GLOBAL
      IMPLICIT NONE
      
      DOUBLE PRECISION :: C1, C2, LD, P1, Q1, C11, C22
      INTEGER :: MA1, NA
      
      
      


         TINF = 700.0D0                                           ! FREESTREAM TEMPERATURE (GIVEN)
         SOS =SQRT(GAMMAN * RCONST * TINF)                        ! SPEED OF SOUND
         UINF = SOS*M                                             ! FREESTREAM VELOCITY


         IF (TINF .GT. 110.4D0) THEN
            FMU =(1.458D0*TINF**1.5)/(TINF+110.4D0)*1.0D-         ! FREESTREAM VISCOSITY
         ELSE
            FMU = (0.693873D0*1.0D-06)*TINF
         END IF

         DO I = 1, N
           T(I) = 1.0D0                                           ! NONDIMENSIONAL TEMPERATURE
           MU(I) = 1.0D0                                          ! NONDIMENSIONAL VISCOSITY
           DMU(I) = 0.0D0                                         ! DERIVATIVE OF NONDIMENISONAL VISCOSITY
         END DO



!====================================================================================
!====================================================================================

            DO NA =1, 50

!-----------SOLVING BL MOMENTUM EQUATION---------------------------------------------
!               BC
                Y(1,1) = 0.0D0                                    ! U AT SURFACE
                Y(3,1) = 0.0D0                                    ! G AT SUARFCE

                P = 0.01D0                                        ! TWO INITIAL GUESSES FOR U' AT SURFACE TO USE SECANT ITERATION
                Q =0.020D0

                Y(2,1) =P
                CALL RK41 
                C1 = Y(1,N)-1.0D0                                 ! Y(1,N) IS U IN FREESTREAM THAT MUST SATISFY 1. 
                                                                  ! HENCE, FOR PERFECT INITIAL GUESS C1 OR C2 MUST BE ZERO.
                Y(2,1) = Q
                CALL RK41
                C2 = Y(1,N)-1.0D0

             
             DO  MA = 1, 10                                       ! SECANT ITERATION STARTS WITH INITIAL GUESS

                R = P - (C1 *(P-Q))/(C1- C2)
                C2 = C1
                Q = P
                Y(2,1) = R
                CALL RK41
                C1 = Y(1,N)-1.0D0
                P = R

                IF (DABS(C1) .LT. 1.0D-10) THEN
                 EXIT
                END IF

            END DO

!-----------SOLVING BL ENERGY EQUATION---------------------------------------------
!             BC

             Y(4,1) = 0.01666667D0                                ! FIXED WALL TEMP -THETA (SCALED TEMPERATURE) AT SURFACE -(GIVEN)

             P1 = 0.4D0                                           ! THETA' (DERIVATIVE OF SCALED TEMPERATURE) AT SUFACE AS INTIAL GUESS
             Q1 = 0.5D0

             Y(5,1) = P1
             CALL RK42
             C11 = Y(4,N)                                         ! Y(4,N) SCALED TEMPRATURE IN FREESTREAM MUST SATISFY ZERO
             Y(5,1) = Q1
             CALL RK42
             C22 = Y(4,N)
          
             DO MA1 = 1, 10                                       ! SECANT ITERATION STARTS WITH INITIAL GUESS

                  R = P1-(C11*(P1-Q1))/(C11-C22)
                  C22 = C11
                  Q1 = P1
                  Y(5,1) = R
                  CALL RK42
                  C11 = Y(4,N)
                  P1 = R

                  IF (DABS(C11) .LT. 1.0D-10) THEN
                     EXIT
                  END IF 
 
              END DO

            DO  J =1, N                                           ! UPDATE VISCOITY USING UPDATED TEMPERATURE
              T(J) =  1+0.5D0*(GAMMAN-1.0D0)*M*M*Y(4,J)
              TINST = T(J)*TINF
              IF (TINST .GT. 110.4D0) THEN
                 WMU =(1.458D0*TINST**1.5)/(TINST+110.4D0)*1.0D-05
              ELSE
                 WMU = (0.693873D0*1.0D-06)*TINST
              END IF
              MU(J) = WMU/FMU    
            END DO

            DO J = 1, N 
                                                                  ! UPDATE VISCOITY DERIVATIVE
             IF ( J .EQ. 1) THEN
               DMU(J)= (-3.0D0*MU(J) + 4.0D0*MU(J+1) - MU(J+2))/(2.0D0*H)
             ELSE IF (J .EQ. N) THEN
               DMU(J) = (3.0D0*MU(J)-4.0D0*MU(J-1)+MU(J-2))/(2.0D0*H)
             ELSE
               DMU(J) =  (MU(J+1)-MU(J-1))/(2.0D0*H)
             END IF
              
           END DO

             PRINT*, 'ITERATION NO.' , NA

          END DO

!=================================================================================
!=================================================================================
         OPEN(2,FILE="DATA.DAT")
         WRITE(2,*),'NON-DIMENSIONAL U,', 'NON-DIMENSIONAL T'
         PRINT*, 'NON-DIMENSIONAL U,', 'NON-DIMENSIONAL T'
         DO  J =1, N
             T(J) =  1+0.5D0*(GAMMAN-1.0D0)*M*M*Y(4,J)
             PRINT*, Y(1,J), T(J)
             WRITE(2,*),Y(1,J), T(J)
         END DO


       END PROGRAM COMPRESSIBLE_BL

