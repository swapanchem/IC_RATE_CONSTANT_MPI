    MODULE DeterminantModule_IC

       IMPLICIT NONE

      CONTAINS

        SUBROUTINE DETERMINANT(N, E, DET)

        INTEGER, INTENT(IN) :: N

        COMPLEX*16, INTENT(INOUT) :: E(N, N)  

        COMPLEX*16, INTENT(OUT) :: DET

        INTEGER :: I, INFO, SWAPS

        INTEGER, ALLOCATABLE :: IPIV(:)

        REAL*8 :: SGN

        COMPLEX*16, ALLOCATABLE :: E_copy(:,:)

        ! ALLOCATE LOCAL MEMORY 

        ALLOCATE(E_copy(N, N), IPIV(N))

        E_copy = E  ! COPY ORIGINAL MATRIX

        ! LU DECOMPOSITION

        CALL ZGETRF(N, N, E_copy, N, IPIV, INFO)

        IF (INFO /= 0) THEN

            PRINT*, "ERROR: Singular matrix, determinant undefined. INFO =", INFO

            DET = (0.0D0, 0.0D0)

            DEALLOCATE(E_copy, IPIV)

            RETURN

        END IF

        ! EVALUATE DETERMINANT

        DET = (1.0D0, 0.0D0)

        DO I = 1, N

            DET = DET * E_copy(I, I)

        END DO

        ! COMPUTE THE SIGN CHNAGE FOR SWAPPING

        SGN = 1.0D0

        SWAPS = 0

        DO I = 1, N

            IF (IPIV(I) /= I) SWAPS = SWAPS + 1

        END DO

        SGN = SGN * (-1.0D0) ** SWAPS

        DET = SGN * DET

        ! DEALLOCATE MEMORY 

        DEALLOCATE(E_copy, IPIV)

        END SUBROUTINE DETERMINANT

        END MODULE DeterminantModule_IC

