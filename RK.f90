!****************************************************************************                                                              
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *  
!****************************************************************************

      SUBROUTINE RK41
      USE global
      DOUBLE PRECISION :: K1(3), K2(3), K3(3), K4(3)
      DOUBLE PRECISION :: F1, A, B, C, D, E


      

       DO I=1, N-1
   
         TMU1 = MU(I)
         TTK1 = T(I)
         DDMU = DMU(I)

      

        A=Y(1,I)
        B=Y(2,I)
        C=Y(3,I)
        

        DO J =1,3
           K1(J) = F1(A, B, C, J)
        END DO
     

        A=Y(1,I)+0.5*K1(1)*H
        B=Y(2,I)+0.5*K1(2)*H
        C=Y(3,I)+0.5*K1(3)*H
        

        DO J =1,3
           K2(J) = F1(A, B, C, J)
        END DO

        A=Y(1,I)+0.5*K2(1)*H
        B=Y(2,I)+0.5*K2(2)*H
        C=Y(3,I)+0.5*K2(3)*H
      

        DO J =1,3
           K3(J) = F1(A, B, C, J)
        END DO

        A=Y(1,I)+K3(1)*H
        B=Y(2,I)+K3(2)*H
        C=Y(3,I)+K3(3)*H
        

        DO J =1,3
           K4(J) = F1(A, B, C, J)          
        END DO
      
        DO J=1,3
           Y(J,I+1)=Y(J,I) + H*(K1(J) + 2.0D0*K2(J) +2.0D0*K3(J) +K4(J))/6.0D0
        END DO
      
      
      END DO
    

      RETURN
      END SUBROUTINE
 
!__________________________________________________________________________________________________
!________________________________________FUNCTION_________________________________________________

      FUNCTION F1(Y1,Y2,Y3,W)
      USE global
      DOUBLE PRECISION :: Y1,Y2,Y3, F1
      INTEGER :: W
      
      IF (W .EQ. 1) THEN
      F1= Y2
      ELSE IF (W .EQ. 2) THEN
      F1 = (-Y3*Y2- Y2*DDMU)/TMU1 
      ELSE 
      F1 = Y1/(2.0D0*TTK1)
      END IF
 
      RETURN
      END FUNCTION

!================================================================

      SUBROUTINE RK42
      USE global
      DOUBLE PRECISION :: K11(5), K22(5), K33(5), K44(5)
      DOUBLE PRECISION :: F2, A2, B2


      

       DO I=1, N-1
   
          UDASH = Y(2,I)
          FUNCG = Y(3,I)
          TMU1 = MU(I)
          DDMU = DMU(I)
    
 
        A2=Y(4,I)
        B2=Y(5,I)
       
        

        DO J =4,5
           K11(J) = F2(A2, B2,J)
        END DO
     

        A2=Y(4,I)+0.5*K11(4)*H
        B2=Y(5,I)+0.5*K11(5)*H
        
        

        DO J =4,5
           K22(J) = F2(A2, B2, J)
        END DO

        A2=Y(4,I)+0.5*K22(4)*H
        B2=Y(5,I)+0.5*K22(5)*H
        
      

        DO J =4,5
           K33(J) = F2(A2, B2, J)
        END DO

        A2=Y(4,I)+K33(4)*H
        B2=Y(5,I)+K33(5)*H
        
        

        DO J =4,5
           K44(J) = F2(A2, B2, J)          
        END DO
      
        DO J=4,5
           Y(J,I+1)=Y(J,I) + H*(K11(J) + 2.0D0*K22(J) +2.0D0*K33(J) +K44(J))/6.0D0
        END DO
      
      
      END DO
    

      RETURN
      END SUBROUTINE
 
!__________________________________________________________________________________________________
!________________________________________FUNCTION_________________________________________________

      FUNCTION F2(Y4,Y5,W2)
      USE global
      DOUBLE PRECISION :: Y4, Y5, F2
      INTEGER :: W2
      
      IF (W2 .EQ. 4) THEN
      F2 = Y5
      ELSE 
      F2 = (-FUNCG*Y5 - 2.0D0*TMU1*UDASH**2 - Y5*DDMU/PR)*PR/TMU1
      END IF
 
      RETURN
      END FUNCTION


