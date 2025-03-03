       PROGRAM IC_FC_HT_TEMP

               USE DeterminantModule_IC  ! IMPORT DETERMINANT MODULE 

               USE MPI

       IMPLICIT NONE

       REAL*8, DIMENSION(:,:), ALLOCATABLE :: RJ, RJT, OMEGAI, &
                               OMEGAF, AS1, AK_DOUBLE_BAR, &
                               BK_DOUBLE_BAR, NACME,P, P2

       REAL*8, DIMENSION(:,:), ALLOCATABLE :: WI, WF, D, NAC

       REAL*8, DIMENSION(:,:), ALLOCATABLE :: DT, NAC_TRANS

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: D_COMPLEX, V1, H1, H2

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: DT_COMPLEX, Y6

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: A, B, BB, B_INV, &
            AA, B2, AA_INV, E2, E1, AB_INV, RJ_COMPLEX, &
            RJT_COMPLEX, X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5, &
            E, B_A_INV, B_A_DOUBLE_BAR_COMPLEX, B_A, &
            B_A_DOUBLE_BAR, BK_DOUBLE_BAR_COMPLEX, &
            AK_DOUBLE_BAR_COMPLEX, BS2, AS2, &
            OMEGAI_COMPLEX, AK_BAR, BK_BAR, NACME_COMPLEX, &
            G11, G12, G21, G22, BBB, P2_COMPLEX, BBB_COPY

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: AK, AK_INV, G, GK_INV

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: F, AKINV_F, H

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: FT, AKINV_F_TRANS, &
                  H_TRANS, V2, V3, AKINV_T_G

       REAL*8 :: TMIN, TMAX, HH, CM_AU, S_AU, T, AMU_AU, AU_HZ, &
               BETA, KB, KBT, TEMP, EV_AU, PF, AMP, DET_I, DET_R, &
               DELE_F_I, NUMR, NUMI, &
               T1, ETA, THETA, SCAL, NUM1R, NUM2R, DET_P2, &
               NUM1I, NUM2I, TRACE_GK_INV_REAL, TRACE_GK_INV_IMAG

       COMPLEX*16 :: DET, NUM, IC_FC, IC_FC_HT_TOTAL, IC_HT, &
               TRACE, TRACE_GK_INV, IC_HT_TEMPO, IC_HT_2

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: NUM1, NUM2, NUM3, NUM4

       INTEGER :: I, J, K, M, N, K1, N1, M1, INFO, M2, NP, &
                  IERROR, MYID, NP1

       INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV

       COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WORK

       COMPLEX*16, DIMENSION(:), ALLOCATABLE :: IC_HT_1

       !Character for reading input file data

        CHARACTER(LEN = 100) :: NAM1, NAM2, NAM3, NAM4, NAM5, NAM6, &
                               NAM7, NAM8, NAM0, NAM


!MPI INITILIALIZATION 

       CALL MPI_INIT(IERROR)

!INITIALIZE THE NUMBER OF PROCESSOR (NP) 
!MYID =  PROCESSORS ID FOR EACH JOB 

       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, IERROR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERROR)      


        OPEN(30,FILE='input_ic.txt')


      ! N       NORMAL MODES OF THE INVESTIGATED MOLECULE
      ! M1      M1 = NUMBER OF INTERVAL FOR INTEGRATION
      ! ETA     ETA = DAMPING PARAMETER IN CM^(-1)
      ! TEMP    TEMP = TEMPERATURE IN K
      ! TMIN    TMIN = LOWER LIMIT OF TIME IN SEC
      ! TMAX    TMAX = UPPER LIMIT OF TIME IN SEC
      ! NP      NP = NUMBER OF PROCESSORS USED


       READ (30, *) NAM

       READ (30, *) NAM0
       READ (30, *) N

       READ (30, *) NAM1
       READ (30, *) M1

       READ (30, *) NAM2
       READ (30, *) TMAX

       READ (30, *) NAM3
       READ (30, *) TMIN

       READ (30, *) NAM4
       READ (30, *) TEMP

       READ (30, *) NAM5
       READ (30, *) ETA

       READ (30, *) NAM6
       READ (30, *) DELE_F_I

       READ(30,*)NAM7
       READ(30,*)NP

       READ(30,*)NAM8

       N1 = 2 * N


      OPEN(66, FILE = 'J.TXT')   ! J is the Duschinsky rotation matrix
 
      OPEN(77, FILE = 'WI.TXT')  ! frequencies of the initial state

      OPEN(88, FILE = 'WF.TXT')  ! frequencies of the final state

      OPEN(99, FILE = 'D.TXT')   ! displacement vector 

      OPEN(109, FILE = 'NAC.TXT') ! nonadiabatic coupling vector

      OPEN(129, FILE = 'IC_FC_HT.TXT') ! real part of the time-dependent correlation function 


  

       !INPUT THE MATRIX AND READ THE MATRIX


       ALLOCATE (RJ(N,N), WI(N,1), WF (N,1), D(N,1), NAC(N,1))


       READ(66, *) (( RJ (I, J), J = 1, N), I = 1, N)       !RJ = DUSCHINSKY ROTATION MATRIX

       READ(77, *) ( WI (I, 1), I = 1, N)                   !WI = FREQUENCIES OF THE INITAL STATE

       READ(88, * ) ( WF (I, 1), I = 1, N )                 !WF = FREQUENCIES OF THE FINAL STATE

       READ(99, *) ( D( I, 1), I = 1, N )                   !D = DISPLACEMENT VECTOR

       READ(109, *) ( NAC (I, 1), I = 1, N )                !NAC = NONADIABATIC COUPLING VECTOR



       !TRANSFER THE DATA TO ATOMIC UNIT

       AU_HZ = 6.579D0 * (10.0D0 ** 15.0D0)              ! ATOMIC UNIT TO HERZ

       CM_AU = 4.55634D0 * (10.0D0 ** ( - 6.0D0))        ! CM^(-1) TO ATOMIC UNIT 

       S_AU = 0.4134D0 * (10.0D0 ** 17.0D0)              ! SECOND TO ATOMIC UNIT 

       AMU_AU = ( 1.82289D0 * (10.0D0 ** 3.0D0 ))        ! ATOMIC MASS UNIT TO ATOMIC UNIT 

       AMU_AU = (( AMU_AU ) ** 0.5D0 )

       EV_AU = 0.036749844D0                             ! EV TO ATOMIC UNIT

       KB = 1.3806452D0 * (10.0D0 ** ( - 23.0D0))        ! BOLTZMANN CONSTANT 

       KBT = KB * TEMP                                   ! THERMAL ENERGY 

       KBT = KBT * (6.242D0 * (10.0D0 ** (18.0D0)))      !JOULE To EV

       KBT = KBT * EV_AU                                 !EV TO ATOMIC UNIT CONVERSION OF THERMAL ENERGY 

       BETA = (1.0D0 / KBT)

     

       DELE_F_I = DELE_F_I * EV_AU         


       DO I = 1, N

          WI(I,1) = WI(I,1) * CM_AU

          WF(I,1) = WF(I,1) * CM_AU

       ENDDO


       TMAX = TMAX * S_AU
       TMIN = TMIN * S_AU
       ETA = ETA * CM_AU


       
      !GENERATE THE DIAGONAL FREQUENCY MATRIX

      ALLOCATE (OMEGAI (N,N), OMEGAF(N,N), OMEGAI_COMPLEX(N,N))

       DO I = 1, N

         DO J = 1, N

         IF (I.EQ.J) THEN

         OMEGAI(I,J) = WI(I,1)

         OMEGAF(I,J) = WF(I,1)

         OMEGAI_COMPLEX(I,J) = COMPLEX(WI(I,1), 0.0D0)

         ELSE

         OMEGAI(I,J) = 0.0D0

         OMEGAF(I,J) = 0.0D0

         OMEGAI_COMPLEX(I,J) = (0.0D0, 0.0D0)

         ENDIF

         ENDDO

       ENDDO

!GENERATE THE TRANSPOSE OF RJ AND D

       ALLOCATE ( RJT(N,N), RJ_COMPLEX(N,N), RJT_COMPLEX(N,N), &
                DT(1,N), NAC_TRANS(1,N), D_COMPLEX(N,1), &
                DT_COMPLEX(1,N), NACME(N,N), NACME_COMPLEX(N,N), &
                P(N,N), P2(N,N), P2_COMPLEX(N,N) )

       DO I = 1, N

          DO J = 1, N

          RJT(J,I) = RJ(I,J)           !RJTRANSPOSE OF DUSCHINSKY ROTATION MATRIX RJ

          RJ_COMPLEX(I,J) = COMPLEX(RJ(I,J),0.0D0)

          RJT_COMPLEX(J,I) = COMPLEX(RJT(J,I),0.0D0)

          ENDDO

       ENDDO


!GENERATE THE TRANSPOSE OF THE DISPLACEMENT VECTOR
       DO I = 1, N

         DT(1,I) = D(I,1)                                            !DT=TRANSPOSE OF DISPLACEMENT VECTOR

         NAC_TRANS(1,I) = NAC(I,1)

         D_COMPLEX(I,1) = COMPLEX(D(I,1), 0.0D0)

         DT_COMPLEX(1,I) = COMPLEX(DT(1,I), 0.0D0)

       ENDDO

!NACME = NONADIABATIC COUPLING MATRIX ELEMENTS

         CALL DGEMM('N','N',N, N, 1, 1.0D0, NAC, N, NAC_TRANS, 1, &
                0.0D0, NACME, N)


        DO I = 1, N

          DO J = 1, N

          NACME_COMPLEX(I, J) = COMPLEX(NACME(I, J), 0.0D0)

          ENDDO

        ENDDO

        DO I = 1, N

          DO J = 1, N

          IF (I.EQ.J) THEN

                  P(I,J) = 2.0D0 * SINH(OMEGAI(I,J) * BETA / 2.0D0)

          ELSE

                  P(I,J) = 0.0D0

          ENDIF

          ENDDO

       ENDDO

       CALL DGEMM('N','N',N, N, N, 1.0D0, P, N, P, N, &
                0.0D0, P2, N)

        DET_P2 = 1.0D0

        DO I = 1, N

         DET_P2 = DET_P2 * P2(I,I)

         ENDDO


       DO I = 1, N

         DO J = 1, N

         P2_COMPLEX(I, J) = COMPLEX(P2(I, J), 0.0D0)

         ENDDO

       ENDDO



!........................................................................................................................................................\\
       !GENERATE THE REQUIRED MATRIX AND MATRIX MULTIPLICATION WITHIN
       !THE TIME LOOP FOR DETERMINANT CALCULATION IN THE FC REGION

      


       DO K1 = 1, M1 + 1

        HH = ( TMAX - TMIN ) / FLOAT( M1 )  !HH=GAP BETWEEN THE TWO SUCCESSIVE POINTS

        M2 = (M1 / 2) + 1

        IF (K1.NE.M2) THEN

        T = TMIN + (HH * (K1 - 1))

        ELSE

        T = 1.0D0 * (10.0D0 ** ( - 24.0D0))

        T = T * S_AU

        ENDIF

        ALLOCATE ( AK_DOUBLE_BAR (N,N), AK_DOUBLE_BAR_COMPLEX(N,N), &
                  BK_DOUBLE_BAR(N,N), BK_DOUBLE_BAR_COMPLEX(N,N), &
                  AS1(N,N), AK_BAR(N,N), BK_BAR(N,N) )


          DO I = 1, N

            DO J = 1, N

            IF (I.EQ.J) THEN


!BS=OMEGAS*TAN(OMEGAS*T/2)

          AK_DOUBLE_BAR(I,J) = OMEGAF(I,J) / SIN(T * OMEGAF(I,J))

          AK_DOUBLE_BAR_COMPLEX(I,J) = COMPLEX(AK_DOUBLE_BAR(I,J), &
                                      0.0D0)



!AS=OMEGAS*SIN(OMEGAS*T/2)

          BK_DOUBLE_BAR(I, J) = OMEGAF(I,J) / TAN( T * OMEGAF(I,J))

          BK_DOUBLE_BAR_COMPLEX(I, J) = COMPLEX(BK_DOUBLE_BAR(I, J), 0.0D0)


!FORMATION OF AT

          AS1(I,J) = (COSH( - 2.0D0 * BETA * OMEGAI(I,J))) - (COS( - 2.0D0 * T &
          *OMEGAI(I,J)))
          AK_BAR(I,J) = COMPLEX((( 2.0D0 * SIN(- OMEGAI(I,J) * T) * &
          COSH( - BETA * OMEGAI(I,J))) / AS1(I,J)), (( - 2.0D0 * &
          COS( - T * OMEGAI(I,J)) * SINH( - BETA * OMEGAI(I,J))) &
          / AS1(I,J)))

!FORMATION OF BT

          BK_BAR(I,J) = COMPLEX(((SIN( - 2.0D0 * OMEGAI(I,J) * T)) &
          / AS1(I,J)), (( - SINH( - 2.0D0 * BETA * OMEGAI(I,J))) &
          / AS1(I,J)))

          ELSE

          AK_BAR(I,J) = (0.0D0, 0.0D0)

          AK_DOUBLE_BAR(I,J) = 0.0D0

          BK_BAR(I,J) = (0.0D0,0.0D0)

          BK_DOUBLE_BAR(I,J) = 0.0D0

          AS1(I,J) = 0.0D0

          AK_DOUBLE_BAR_COMPLEX(I,J) = (0.0D0, 0.0D0)

          BK_DOUBLE_BAR_COMPLEX(I,J) = (0.0D0, 0.0D0)

          ENDIF

          ENDDO
          
       ENDDO

 

      !FORMATION OF K MATRIX AND ITS INVERSE. HERE, WE DENOTES IT AS AK

      !FORMATION OF A AND B MATRIX

      ALLOCATE (AS2(N,N), BS2(N,N), X1(N,N), X2(N,N), X3(N,N), &
               X4(N,N), A(N,N), B(N,N), B_INV(N,N), IPIV (N1), &
               WORK (N1) )

      CALL ZGEMM('N','N', N, N, N, (1.0D0, 0.0D0), OMEGAI_COMPLEX, N, &
                AK_BAR, N, (0.0D0, 0.0D0), AS2, N)


      CALL ZGEMM('N','N', N, N, N, (1.0D0,0.0D0), OMEGAI_COMPLEX, N, &
              BK_BAR, N, (0.0D0,0.0D0), BS2, N)


      CALL ZGEMM('N','N', N, N, N, (1.0D0,0.0D0), RJT_COMPLEX, N, &
                AS2, N, (0.0D0, 0.0D0), X1, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), X1, N, RJ_COMPLEX, &
                 N, (0.0D0, 0.0D0), X2, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, N, &
                BS2, N, (0.0D0, 0.0D0), X3, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), X3, N, RJ_COMPLEX, &
                 N, (0.0D0, 0.0D0), X4, N)


      DO I = 1, N

       DO J = 1, N

        A(I,J) = (AK_DOUBLE_BAR_COMPLEX(I, J) + X2(I, J))

        B(I,J) = (BK_DOUBLE_BAR_COMPLEX(I, J) + X4(I, J))

        ENDDO

      ENDDO



      !INVERSE OF B 

      DO I = 1, N

       DO J = 1, N

       B_INV(I,J) = B(I,J)

       ENDDO

      ENDDO

      CALL ZGETRF(N, N, B_INV, N, IPIV, INFO)

      CALL ZGETRI(N, B_INV, N, IPIV, WORK, N, INFO)

!FORMATION OF THE DENOMINATOR MATRIX WITHIN THE SQUARE ROOT OG THE
!DETERMINANT WHICH IS DENOTED HERE AS AA


      ALLOCATE ( AB_INV(N,N), E1(N,N), E2(N,N), B2(N,N), AA(N,N), &
               AA_INV (N,N) )
      
      CALL ZGEMM('N','N',N, N, N, (1.0D0,0.0D0), A, N, B_INV, N, &
              (0.0D0, 0.0D0), AB_INV, N)
 

      CALL ZGEMM('N','N',N, N, N, (1.0D0,0.0D0), AB_INV, N, A, N, &
              (0.0D0, 0.0D0), E1, N)
      

      CALL ZGEMM('N','N',N, N, N, (1.0D0,0.0D0), B, N, E1, N, &
              (0.0D0, 0.0D0), E2, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), B, N, B, N, &
              (0.0D0, 0.0D0), B2, N)


      DO I = 1, N

       DO J = 1, N

       AA(I,J) = (B2(I,J) - E2(I,J))

       ENDDO

      ENDDO

      !INVERSE OF AA MATRIX (INVERSE OF THE DENOMINATOR WITHIN THE SQUARE
      !ROOT OF THE DETERMINANT)


      DO I = 1, N

       DO J = 1, N 

       AA_INV(I, J) = AA(I, J)

       ENDDO

      ENDDO

      CALL ZGETRF(N, N, AA_INV, N, IPIV, INFO)

      CALL ZGETRI(N, AA_INV, N, IPIV, WORK, N, INFO)

     

     !FORMATION OF K MATRIX

      ALLOCATE ( AK(N1,N1), AK_INV(N1,N1), X5(N,N), BB(N,N), &
               BBB(N,N), BBB_COPY(N,N) )


      AK(1 : N, 1 : N) = B

      AK(1 : N, N + 1 : N1) = - A

      AK(N + 1 : N1, 1 : N) = - A

      AK(N + 1 : N1, N + 1 : N1) = B


      !INVERSE OF AK 
      !STORED THE AK MATRIX AS AK_INV FOR INVERSION CALCULATION USING BLAS/LAPACK

      DO I = 1, N1

       DO J = 1, N1

        AK_INV(I, J) = AK(I, J)

       ENDDO
       
      ENDDO


      CALL ZGETRF(N1, N1, AK_INV, N1, IPIV, INFO)

      CALL ZGETRI(N1, AK_INV, N1, IPIV, WORK, N1, INFO)

        

!XX=NUMERATOR MATRIX WITHIN THE SQUARE ROOT OF THE DETERMINANT

      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0),  &
                AK_DOUBLE_BAR_COMPLEX, N, AS2, N, &
                (0.0D0, 0.0D0), X5,N)


      !FORMATION OF THE COMPLETE DETERMINANT DENOTED BY BB


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), X5, N, AA_INV, N, &
              (0.0D0, 0.0D0), BB, N)

      
      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), BB, N, P2_COMPLEX, &
                 N, (0.0D0, 0.0D0), BBB, N)



       !COPY BBB TO BBB_COPY 

       DO I = 1, N
          
          DO J = 1, N

          BBB_COPY (I, J) = BBB(I, J)

          ENDDO

       ENDDO
!DETERMINANAT OF BB WHICH IS N BY N MATRIX

      CALL DETERMINANT(N, BBB_COPY, DET)


!AMPLITUDE AND PHASE OF THE DETERMINANT 
      AMP = ABS( DET )

      AMP = SQRT( AMP )

      DET_I = AIMAG( DET )

      DET_R = REAL( DET )

      THETA = DET_I / DET_R

      THETA = THETA / 2.0D0

      WRITE( *, *) T, DET

!........................................................................................................................................................

!........................................................................................................................................................

      !SIMPLIFICATION OF THE EXPONENTIAL PART 

      !FORMATION OF E MATRIX OF N BY N 



      ALLOCATE (E(N,N), B_A(N,N), B_A_INV(N,N), Y1(N,N), V1(N,1), &
                F(N1,1), FT(1,N1), V2(1,N1), NUM1(1,1), V3(1,N), &
                NUM2(1,1))

      DO I = 1, N

       DO J = 1, N

       E(I, J) = (BS2(I, J) - AS2(I, J))

       B_A(I, J) = B(I, J) - A(I, J)

       ENDDO

      ENDDO


       DO I = 1, N

       DO J = 1, N

        B_A_INV(I, J) = B_A(I, J)

       ENDDO

      ENDDO

      CALL ZGETRF(N, N, B_A_INV, N, IPIV, INFO)

      CALL ZGETRI(N, B_A_INV, N, IPIV, WORK, N, INFO)




      !FORMATION OF F VECTOR AND ITS TRANSPOSE FT CONTAINING 2N ELEMENTS


      CALL ZGEMM('N','N',N, N, N,(1.0D0, 0.0D0), RJT_COMPLEX, N, E, N, &
                (0.0D0, 0.0D0), Y1, N)

      CALL ZGEMM('N','N',N, 1, N, (1.0D0, 0.0D0), Y1, N, D_COMPLEX, N, &
                (0.0D0,0.0D0), V1, N)


      F(1 : N, 1) = V1(:,1)

      F(N + 1 : N1, 1 ) = V1(:,1)

      DO I = 1, N1

       FT(1, I) = F(I, 1)

      ENDDO

!MATRIX MULTIPLICATION FOR THE EXPONENTIAL PART


      CALL ZGEMM('N','N',1, N1, N1, (1.0D0,0.0D0), FT, 1, AK_INV, N1, &
              (0.0D0,0.0D0), V2, 1)


      CALL ZGEMM('N','N',1, 1, N1, (1.0D0, 0.0D0), V2, 1, F, N1, &
              (0.0D0, 0.0D0), NUM1, 1)


      CALL ZGEMM('N','N',1, N, N, (1.0D0, 0.0D0), DT_COMPLEX, 1, E, N, &
                (0.0D0, 0.0D0), V3, 1)


      CALL ZGEMM('N','N',1, 1, N, (1.0D0, 0.0D0), V3, 1, D_COMPLEX, N, &
            (0.0D0, 0.0D0), NUM2, 1)


      NUM1R = REAL(NUM1(1, 1)) / 2.0D0

      NUM1I = AIMAG(NUM1(1, 1)) / 2.0D0

      NUM2R = REAL(NUM2(1, 1))

      NUM2I = AIMAG(NUM2(1, 1))

      NUMI = - NUM1R + NUM2R + THETA - ( T * DELE_F_I ) 

      NUMR = NUM1I - NUM2I


      !TOTAL PART WITHIN THE EXPONENTIAL 

      NUM = COMPLEX( NUMR, NUMI )

      ! CORRELATION FUNCTION IN THE FRANCK-CONDON REGION OR X(T) TERM

      IC_FC = AMP * EXP( NUM )
      
    

!CONSTRUCTION OF THE HERZBERG-TELLER(HT) PART OR X_K,M(T) TERM

!PARALLELIZATION OF THE MOST TIME CONSUMING LOOP BY USING ALL PROCESSORS
!THROUGH  DIVIDING THE LOOP ITERATIONS AMONG ALL PROCESSSORS


       

       IC_HT = (0.0D0, 0.0D0)

       IC_HT_2 = (0.0D0, 0.0D0)

       NP1 = MYID + 1

       DO K = NP1, N, NP

          DO  M = 1, N

            ALLOCATE (G11(N,N), G12(N,N), G21(N,N), G22(N,N), &
                     H1(N,1), H2(N,1))

             DO I = 1, N

               DO J = 1, N

                IF (I.EQ.K) THEN

                G11(I, J) = - BK_DOUBLE_BAR_COMPLEX(K, K) * X2(M, J)

                G12(I, J) = BK_DOUBLE_BAR_COMPLEX(K, K) * X4(M, J)

                G21(I, J) = AK_DOUBLE_BAR_COMPLEX(K, K) * X2(M, J)

                G22(I, J) = - AK_DOUBLE_BAR_COMPLEX(K, K) * X4(M, J)

                ELSE

                G11(I, J) = (0.0D0, 0.0D0)
                G12(I, J) = (0.0D0, 0.0D0)
                G21(I, J) = (0.0D0, 0.0D0)
                G22(I, J) = (0.0D0, 0.0D0)

                ENDIF

               ENDDO

               IF (I.EQ.K) THEN

               H1(I, 1) = BK_DOUBLE_BAR_COMPLEX(K, K) * V1(M, 1)
               H2(I, 1) = - AK_DOUBLE_BAR_COMPLEX(K, K) * V1(M, 1)

               ELSE

               H1(I, 1) = (0.0D0, 0.0D0)

               H2(I, 1) = (0.0D0, 0.0D0)

               ENDIF

              ENDDO


          ALLOCATE (G(N1,N1), H(N1,1), H_TRANS(1,N1), GK_INV(N1,N1), &
                AKINV_F(N1,1), AKINV_F_TRANS(1,N1), AKINV_T_G(1,N1), &
                NUM3(1,1), NUM4(1,1))
     

         G(1 : N, 1 : N ) = G11

         G( 1 : N, N + 1 : N1 ) = G12

         G( N + 1 : N1, 1 : N )= G21

         G(N + 1 : N1, N + 1 : N1) = G22

         H(1 : N, 1 ) = H1(:,1)

         H(N + 1 : N1, 1 ) = H2(:,1)

!TRANSPOSE OF THE H VECTOR

         DO I = 1, N1

           H_TRANS(1, I) = H(I, 1)

        ENDDO

!MATRIX MULTIPLICATION FOR THE SPIN-VIBRONIC PART

        !THIS IS FOR THE Tr(GK^(-1))+ [(K^(-1)F)^(T)G(K^(-1)F)]-H^(T)K^(-1)F


     CALL ZGEMM('N','N', N1, N1, N1, (1.0D0, 0.0D0), G, N1, AK_INV, &
               N1, (0.0D0, 0.0D0), GK_INV, N1)


      CALL ZGEMM('N','N',N1, 1, N1, (1.0D0, 0.0D0), AK_INV, N1, F, N1, &
              (0.0D0, 0.0D0), AKINV_F, N1)



      DO I = 1, N1

       AKINV_F_TRANS(1, I) = AKINV_F(I, 1)   !TRANSPOSE OF AKINV_F i.e. K^(-1)F

      ENDDO


      CALL ZGEMM('N','N',1, N1, N1, (1.0D0, 0.0D0), AKINV_F_TRANS, 1, &
                G, N1, (0.0D0, 0.0D0), AKINV_T_G, 1)


      CALL ZGEMM('N','N',1, 1, N1, (1.0D0, 0.0D0), AKINV_T_G, 1, &
                AKINV_F, N1, (0.0D0, 0.0D0), NUM3, 1)


      CALL ZGEMM('N','N', 1, 1, N1, (1.0D0, 0.0D0), H_TRANS, 1, &
               AKINV_F, N1, (0.0D0, 0.0D0), NUM4, 1)


      TRACE_GK_INV = (0.0D0, 0.0D0)


      DO I = 1, N1

        DO J = 1, N1

        IF (I.EQ.J) THEN

        TRACE_GK_INV = TRACE_GK_INV + (GK_INV(I, J))

        ELSE

        TRACE_GK_INV = TRACE_GK_INV + (0.0D0, 0.0D0)

        ENDIF

        ENDDO

      ENDDO

      TRACE_GK_INV_REAL = REAL( TRACE_GK_INV )

      TRACE_GK_INV_IMAG = AIMAG( TRACE_GK_INV )


      TRACE = COMPLEX( - TRACE_GK_INV_IMAG, TRACE_GK_INV_REAL )


      !TEMPORARY HT PART FOR EACH PAIR OF NORMAL MODES
      
      
       IC_HT_TEMPO = ( TRACE + NUM3(1, 1) - NUM4(1, 1)) * & 
                     NACME_COMPLEX(K, M)
    
       IC_HT = IC_HT + ( IC_HT_TEMPO )

      DEALLOCATE (G, H, H_TRANS, GK_INV, &
              AKINV_F, AKINV_F_TRANS, AKINV_T_G, NUM3, NUM4)

      DEALLOCATE (G11, G12, G21, G22, H1, H2 )

      
      ENDDO   !END LOOP FOR K 

      ENDDO   !END LOOP FOR M

      CALL MPI_BARRIER( MPI_COMM_WORLD, IERROR )


!CALLING ALL LOOP ITERATIONS RESULTS INTO THE MASTER ID I.E. NP=0 BY
!CALLING THE MPI_GATHER FUNCTION


       ALLOCATE ( IC_HT_1(NP) )

       CALL MPI_GATHER(IC_HT, 1, MPI_DOUBLE_COMPLEX, IC_HT_1, 1 , &
       MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, IERROR)


!CALCULATION OF TOTAL SV FOR EACH TIME


      IF (MYID == 0 ) THEN

      DO I = 1, NP

      IC_HT_2 = IC_HT_2 + IC_HT_1(I)

      ENDDO


      T1 = ABS( T )

      IC_FC_HT_TOTAL = AMP * EXP( NUM ) * IC_HT_2 * EXP( - ETA * T1 )


      WRITE( 129, * ) T, REAL( IC_FC_HT_TOTAL )



      ENDIF         !ENDIF MYID 


      DEALLOCATE ( AK_DOUBLE_BAR, AK_DOUBLE_BAR_COMPLEX, &
                  BK_DOUBLE_BAR, BK_DOUBLE_BAR_COMPLEX, AS1, &
                  AK_BAR, BK_BAR )

      DEALLOCATE (AS2, BS2, X1, X2, X3, &
             X4, A, B, B_INV, IPIV, WORK )

      DEALLOCATE ( AB_INV, E1, E2, B2, AA, &
             AA_INV )

      DEALLOCATE ( AK, AK_INV, X5, BB, &
             BBB, BBB_COPY )

      DEALLOCATE (E, B_A, B_A_INV, Y1, V1, &
              F, FT, V2, NUM1, V3, NUM2 )



     DEALLOCATE ( IC_HT_1 )

 
      ENDDO         !END LOOP FOR T


      DEALLOCATE (RJ, WI, WF, D, NAC)


      DEALLOCATE (OMEGAI, OMEGAF, OMEGAI_COMPLEX)
     

      DEALLOCATE ( RJT, RJ_COMPLEX, RJT_COMPLEX, &
              DT, NAC_TRANS, D_COMPLEX, DT_COMPLEX, NACME, &
              NACME_COMPLEX, P, P2, P2_COMPLEX )


       !DESTRUCTION OF THE MPI 

       CALL MPI_FINALIZE( IERROR )


       END PROGRAM IC_FC_HT_TEMP


