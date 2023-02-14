       PROGRAM ISC_DSO_CORRELATION_TEMP
       IMPLICIT NONE
       REAL*8,DIMENSION(30,30)::RJ,RJT,OMEGAS,OMEGAT,OMEGAS2, &
       OMEGA2,OMEGA,ST,STI,BT,BTI,RJ1,RJ2,B1,B2,B3,B4,C1,C2,C3, &
       C4,R,R1,OMEGAT2,P,P1,P2,A1,P5IIM,P5IR,BSR,BSIR,BSIM,BSIIM, &
       SSIR,SSIIM,B5,B6,B7,C5,C6,E1,E2,P3,P4,X4,X5,X1,X2,X3,X3I
       REAL*8,DIMENSION(30,1)::WS,WT,A3,A4,A5,A6,F4,F5,D,A2,F3
       REAL*8,DIMENSION(1,1)::A7,F7,F6
       REAL*8,DIMENSION(1,30)::DT
       COMPLEX*16,DIMENSION(30,30)::RT,RTI,RT1,P5I,P5,E,BS,SS,BSI,SSI
       REAL*8::TMIN,TMAX,H,CM_AU,S_AU,T,AMP,SQ_DET,K1,TEMP_AU, &
       K2,DELE_AB,SOC,GEN_FUNC,KISC_DIR,AMU_AU,CONST,AU_HZ,KBT, &
       ETA,U2,U1,GEN_FN,THETA,U,BETA,TEMP,SUM_Z,EV_AU,KB, &
       SCAL,SCAL1
       COMPLEX*16::DET
       INTEGER::I,J,N,K,M,INFO
       INTEGER,DIMENSION(30)::IPIV
       COMPLEX*16,DIMENSION(30)::WORK

! THIS INPUT SECTION IS FOR INTERSYSTEM CROSSING RATE CALCULATION OF
! URACIL MOLECULE 

       OPEN(30,FILE='JMAT-URACIL-S1-T1.TXT')
       OPEN(3,FILE='WS1-URACIL.TXT')
       OPEN(2,FILE='WT1-URACIL.TXT')
       OPEN(10,FILE='SHIFT-VECTOR-URACIL-S1-T1.TXT')
       OPEN(34,FILE='GEN_FN.TXT')
      

!       WRITE(*,*)'GIVE N,TMIN,TMAX,M,ETA,TEMP'    !INPUT PARAMETER FOR NUMERICAL INTEGRATION
!       READ(*,*)N,TMIN,TMAX,M,ETA,TEMP            !TEMP=TEMPERATURE IN KELVIN

       N=30
       TMIN=1.0D0*(10.0D0**(-24.0D0))
       TMAX=20.0D0*(10.0D0**(-12.0D0))
       M=20000
       ETA=2.0D0
       TEMP=300.0D0

       !READ THE DUSCHINSKY ROTATION MATRIX, FREQUENCIES AND SHIFT
       !VECTOR

       READ(30,*)((RJ(I,J),J=1,N),I=1,N)       !RJ=DUSCHINSKY ROTATION MATRIX     
       READ(3,*)(WS(I,1),I=1,N)                !WS=FREQUENCY VECTOR OF SINGLET STATE
       READ(2,*)(WT(I,1),I=1,N)                !WT=FREQUENCY VECTOR OF TRIPLET STATE
       READ(10,*)(D(I,1),I=1,N)                !D=DISPLACEMENT VECTOR 

       !TRANSFER THE DATA TO ATOMIC UNIT

       AU_HZ=6.579D0*(10.0D0**15.0D0)
       CM_AU=4.55634D0*(10.0D0**(-6.0D0))
       S_AU=0.4134D0*(10.0D0**17.0D0)
       AMU_AU=(1.82289D0*(10.0D0**3.0D0))
       AMU_AU=((AMU_AU)**0.5D0)
       EV_AU=0.036749844D0
       DELE_AB=0.810D0       !ENERGY GAP BETWEEN S1 AND T1 IN EV
       DELE_AB=DELE_AB*EV_AU
       SOC=53.10D0           !SPIN-ORBIT COUPLING BETWEEN S1 AND T1 IN CM^(-1)
       SOC=SOC*CM_AU       
       KB=1.3806452D0*(10.0D0**(-23.0D0))
       KBT=KB*TEMP
       KBT=KBT*(6.242D0*(10.0D0**(18.0D0)))    !JOULE To EV    
       KBT=KBT*EV_AU                           !EV TO AU
       BETA=(1.0D0/KBT)

       DO I=1,N
          WS(I,1)=WS(I,1)*CM_AU
          WT(I,1)=WT(I,1)*CM_AU
       ENDDO

       TMAX=TMAX*S_AU
       TMIN=TMIN*S_AU
       ETA=ETA*CM_AU

       !GENERATE THE DIAGONAL FREQUENCY MATRIX

       DO I=1,N
         DO J=1,N
         IF(I.EQ.J)THEN
         OMEGAS(I,J)=WS(I,1)
         OMEGAT(I,J)=WT(I,1)
         ELSE
         OMEGAS(I,J)=0.0D0
         OMEGAT(I,J)=0.0D0
         ENDIF
         ENDDO
       ENDDO

       !GENERATE THE TRANSPOSE OF RJ AND D

       DO I=1,N
          DO J=1,N
          RJT(I,J)=RJ(J,I)                                              !RJTRANSPOSE OF DUSCHINSKY ROTATION MATRIX RJ
          ENDDO
       ENDDO   

       !GENERATE THE TRANSPOSE OF THE DISPLACEMENT VECTOR

       DO I=1,N
         DT(1,I)=D(I,1)                                                 !DT=TRANSPOSE OF DISPLACEMENT VECTOR
       ENDDO

       !GENERATE THE SQUARE OF THE FREQUENCY OMEGAT AND OMEGAS MATRIX

       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAS,N,OMEGAS,N,0.0D0,OMEGAS2,N)
       
       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAT,N,OMEGAT,N,0.0D0,OMEGAT2,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAS,N,OMEGAT,N,0.0D0,OMEGA2,N)

       !CALCULATION OF THE CANONICAL PARTITION FUNCTION

       SUM_Z=0.0D0
       SCAL=1.0D0
       DO I=1,N
         DO J=1,N
         IF(I.EQ.J)THEN
         SUM_Z=SUM_Z+(EXP(-OMEGAS(I,J)*BETA))                        !SUM_OMEGAS=SUMMATION OF DIAGONAL ELEMENT OF WS ELEMENT
         SCAL=SCAL*(4.0D0*(((COSH(BETA*OMEGAS(I,J)/2.0D0))**2.0D0)-1))
         ELSE
         ENDIF
         ENDDO
       ENDDO
       WRITE(*,*)SUM_Z,SCAL,SQRT(SCAL)
     
!........................................................................................................................................................\\
       !GENERATE THE REQUIRED MATRIX AND MATRIX MULTIPLICATION WITHIN
       !THE TIME LOOP FOR DETERMINANT CALCULATION
       KISC_DIR=0.0D0
       H=(TMAX-TMIN)/FLOAT(M)
       DO K=1,M
       T=TMIN+(H*(K-1))
          DO I=1,N
          DO J=1,N
          IF(I.EQ.J)THEN
          ST(I,J)=SIN(OMEGAT(I,J)*T)                                   !FORM OF ST DIAGONAL MATRIX
          STI(I,J)=(1.0D0)/ST(I,J)                                      !INVERSE OF ST MATRIX
          BT(I,J)=TAN((OMEGAT(I,J)*T)/2.0D0)                            !FORM OF DIAGONAL BT MATRIX
          BTI(I,J)=(1.0D0)/(BT(I,J))                                    !INVERSE OF ST MATRIX

!FORM THE DIAGONAL SS AND BS MATRIX

          X1(I,J)=(OMEGAS(I,J)*BETA)/2.0D0   
          X2(I,J)=(OMEGAS(I,J)*T)/2.0D0
          X4(I,J)=(EXP(X1(I,J)*2.0D0)+EXP(-2.0D0*X1(I,J)))/2.0D0        !DEFINITION OF COSH(X)
          X5(I,J)=(EXP(X1(I,J)*2.0D0)-EXP(-2.0D0*X1(I,J)))/2.0D0        !DEFINITION OF SINH(X)
          X3(I,J)=X4(I,J)+COS(X2(I,J)*2.0D0)
          X3I(I,J)=1.0D0/X3(I,J)
          BS(I,J)=CMPLX((X5(I,J)*X3I(I,J)),((-SIN(X2(I,J)*2.0D0) &
          *X3I(I,J))))                                                  !FORM OF BS DIAGONAL MATRIX
          BSI(I,J)=(1.0D0/(BS(I,J)))                                    !INVERSE OF BS MATRIX IS BSI
          SS(I,J)=CMPLX((X5(I,J)*COS(2.0D0*X2(I,J))),(-X4(I,J) &
          *SIN(2.0D0*X2(I,J))))                                         !FORM OF SS DIAGONAL MATRIX
          SSI(I,J)=(1.0D0/(SS(I,J)))                                    !INVERSE OF SS MATRIX IS SSI
          ELSE
          ST(I,J)=0.0D0
          STI(I,J)=0.0D0
          BT(I,J)=0.0D0
          BTI(I,J)=0.0D0
          BS(I,J)=(0.0D0,0.0D0)
          BSI(I,J)=(0.0D0,0.0D0)
          SS(I,J)=(0.0D0,0.0D0)
          SSI(I,J)=(0.0D0,0.0D0)
          ENDIF
          ENDDO
       ENDDO

       DO I=1,N
          DO J=1,N
          BSR(I,J)=REAL(BS(I,J))                    !REAL PART OF BS COMPLEX MATRIX
          BSIM(I,J)=AIMAG(BS(I,J))                  !IMAGINARY PART OF BS COMPLEX MATRIX
          BSIR(I,J)=REAL(BSI(I,J))                  !REAL PART OF INVERSE OF BS MATRIX
          BSIIM(I,J)=AIMAG(BSI(I,J))                !IMAGINARY PART OF INVERSE OF BS MATRIX
          SSIR(I,J)=REAL(SSI(I,J))                  !REAL PART OF INVERSE OF SS MATRIX
          SSIIM(I,J)=AIMAG(SSI(I,J))                !IMAGINARY PART OF INVERSE OF SS MATRIX
          ENDDO
       ENDDO

       !GENERATE THE MATRIX RJT*OMRGAT2*RJ

       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAT2,N,RJ,N,0.0D0,RJ1,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,RJT,N,RJ1,N,0.0D0,RJ2,N)


       !GENERATE THE MATRIX OMEGAS*BSR*RJT*OMEGAT*BTI*RJ

       CALL DGEMM('N','N',N,N,N,1.0D0,BTI,N,RJ,N,0.0D0,B1,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAT,N,B1,N,0.0D0,B2,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,RJT,N,B2,N,0.0D0,B3,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,BSR,N,B3,N,0.0D0,B4,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAS,N,B4,N,0.0D0,B5,N)

       !GENERATE THE TERM OMEGAS*BSIM*RJT*OMEGAT*BTI*RJ

       CALL DGEMM('N','N',N,N,N,1.0D0,BSIM,N,B3,N,0.0D0,B6,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAS,N,B6,N,0.0D0,B7,N)


       !GENERATE THE MATRIX RJT*OMEGAT*BT*RJ*OMEGAS

       CALL DGEMM('N','N',N,N,N,1.0D0,RJ,N,OMEGAS,N,0.0D0,C1,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,BT,N,C1,N,0.0D0,C2,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAT,N,C2,N,0.0D0,C3,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,RJT,N,C3,N,0.0D0,C4,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,C4,N,BSIR,N,0.0D0,C5,N)


       CALL DGEMM('N','N',N,N,N,1.0D0,C4,N,BSIIM,N,0.0D0,C6,N)


       DO I=1,N
          DO J=1,N
          R(I,J)=RJ2(I,J)+OMEGAS2(I,J)+B7(I,J)-C6(I,J)
          R1(I,J)=C5(I,J)-B5(I,J)
          RT(I,J)=CMPLX(R(I,J),R1(I,J))
          ENDDO
       ENDDO

       !GENERATE THE INVERSE OF DENOMINATOR MATRIX RT  RT=FORM OF
       !DENOMINATOR MATRIX WITHIN SQUARE ROOT


       DO I=1,N
         DO J=1,N
         RTI(I,J)=RT(I,J)
         ENDDO
       ENDDO

       CALL ZGETRF(N,N,RTI,N,IPIV,INFO)
       CALL ZGETRI(N,RTI,N,IPIV,WORK,N,INFO)

       !GENERATE THE MATRIX SSI*STI*OMEGAT*OMEGAS        !MATRIX FORM OF NUMERATOR WITHIN SQUARE ROOT

       CALL DGEMM('N','N',N,N,N,1.0D0,STI,N,OMEGA2,N,0.0D0,OMEGA,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,SSIR,N,OMEGA,N,0.0D0,E1,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,SSIIM,N,OMEGA,N,0.0D0,E2,N)

       DO I=1,N
         DO J=1,N
          E(I,J)=CMPLX(-E2(I,J),E1(I,J))
       ENDDO
       ENDDO

       !GENERATE THE TOTAL MATRIX FORM FOR DETERMINANT CALCULATION =OMEGAI*RTI
       
       CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),E,N,RTI, &
       N,(0.0D0,0.0D0),RT1,N)

       !CALCULATE THE DETERMINANT OF THE MATRIX RT1

       CALL DETERMINANT(N,RT1,DET)                        !RT1=COMPLETE DETERMINANT WITHIN SQUARE ROOT

       U1=AIMAG(DET)
       U2=REAL(DET)
       AMP=ABS(DET)                                      !AMP=AMPLITUDE OF DETERMINANT
       THETA=(ATAN(U1/U2))                                !THETA=PHASE FACTOR OF THE INDIVIDUAL DETERMINANT AT EACH TIME
       THETA=THETA/2.0D0
       SQ_DET=SQRT(AMP)                                    !SQ_RT=SQUARE ROOT OF THE ABSLOUTE VALUE OF DETERMINANT
!........................................................................................................................................//
       !GENERATE THE REAL EXPONENTIAL PART

       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAT,N,BT,N,0.0D0,A1,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,A1,N,RJ,N,0.0D0,P1,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,RJT,N,P1,N,0.0D0,P,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAS,N,BSR,N,0.0D0,P2,N)

       CALL DGEMM('N','N',N,N,N,1.0D0,OMEGAS,N,BSIM,N,0.0D0,P3,N)

      DO I=1,N
          DO J=1,N
          P4(I,J)=(P(I,J)+P3(I,J))
          P5(I,J)=CMPLX(P2(I,J),P4(I,J))    
          ENDDO
       ENDDO
    
       DO I=1,N
         DO J=1,N
         P5I(I,J)=P5(I,J)
         ENDDO
       ENDDO

       CALL ZGETRF(N,N,P5I,N,IPIV,INFO)
       CALL ZGETRI(N,P5I,N,IPIV,WORK,N,INFO)
      
       DO I=1,N
         DO J=1,N
         P5IR(I,J)=REAL(P5I(I,J))
         P5IIM(I,J)=AIMAG(P5I(I,J))
         ENDDO
       ENDDO

       CALL DGEMM('N','N',N,1,N,1.0D0,A1,N,D,N,0.0D0,A2,N)

       CALL DGEMM('N','N',N,1,N,1.0D0,RJT,N,A2,N,0.0D0,A3,N)

       CALL DGEMM('N','N',N,1,N,1.0D0,P5IR,N,A3,N,0.0D0,A4,N)

       CALL DGEMM('N','N',N,1,N,1.0D0,RJ,N,A4,N,0.0D0,A5,N)

       CALL DGEMM('N','N',N,1,N,1.0D0,A1,N,A5,N,0.0D0,A6,N)

       CALL DGEMM('N','N',1,1,N,1.0D0,DT,1,A6,N,0.0D0,A7,1)
       K1=-A7(1,1)                           !K1=TERM WITHIN THE EXPONENTIAL TERM

!..................................................................................................................................................//
!GENERATE THE TERM WITHIN COSINE

       CALL DGEMM('N','N',N,1,N,1.0D0,P5IIM,N,A3,N,0.0D0,F3,N)

       CALL DGEMM('N','N',N,1,N,1.0D0,RJ,N,F3,N,0.0D0,F4,N)

       CALL DGEMM('N','N',N,1,N,1.0D0,A1,N,F4,N,0.0D0,F5,N)
       
       CALL DGEMM('N','N',1,1,N,1.0D0,DT,1,F5,N,0.0D0,F6,1)
 
       CALL DGEMM('N','N',1,1,N,1.0D0,DT,1,A2,N,0.0D0,F7,1)

       K2=(-F6(1,1)-F7(1,1)+(T*DELE_AB))        !K2=TERM WITHIN THE COSINE TERM


!FINAL EXPRESSION OF THE REAL PART

       CONST=(SOC**2.0D0)*(1.0D0/SUM_Z)*2.0D0
       GEN_FUNC=SQ_DET*EXP(K1)*COS(K2+THETA) &
       *EXP(-ETA*(T**2.0D0))

       WRITE(34,*)T,GEN_FUNC

!......................................................................................................................................................//
!INTEGRATION USING SIMPSIN'S ONE-THIRD RULE WITHIN TIME DOMAIN
      IF((K.EQ.1).OR.(K.EQ.(M)))THEN
          KISC_DIR=KISC_DIR+GEN_FUNC
          ELSE
              IF(MOD(K,2).EQ.0.0D0)THEN
              KISC_DIR=KISC_DIR+(GEN_FUNC*4.0D0)
              ELSE
              KISC_DIR=KISC_DIR+(2.0D0*GEN_FUNC)
              ENDIF
      ENDIF

      ENDDO
       KISC_DIR=(KISC_DIR*H*CONST)/3.0D0
       KISC_DIR=KISC_DIR*AU_HZ

       WRITE(*,*)'KISC_DIR IS',KISC_DIR,CONST

       END PROGRAM ISC_DSO_CORRELATION_TEMP

 
!..................................................................................................................................................................\\
!SUBROUTINE FOR DETERMINANT CALCULATION USING LAPACK-BLAS LIBRARY

      SUBROUTINE DETERMINANT(N,A,DET)
      COMPLEX*16::A(N,N)
      COMPLEX*16::DET,SGN
      INTEGER::I,INFO,N
      INTEGER::IPIV(N)

      CALL ZGETRF(N, N, A, N, IPIV,INFO)
      DET =(1.0D0,0.0D0)
       DO I =1,N
        DET = DET*A(I,I)
       ENDDO
      SGN =(1.0D0,0.0D0)
      DO I = 1, N
      IF(IPIV(I) /= I) THEN
        SGN = -SGN
      END IF
      ENDDO
      DET= SGN*DET
      RETURN
      END

