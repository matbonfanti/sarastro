!***************************************************************************************
!*                              MODULE MatrixInversion
!***************************************************************************************
!
!>  \brief     Inversion, regularization and analysis
!>  \details   This module defines the procedures which are necessary to diagnose \n
!>             and handle regularization problems in matrix inversion. The module \n
!>             is based on P.E.'s Qoupè code.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             1 August 2017
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          N.A.
!
!***************************************************************************************

#if !defined(WITH_LAPACK)
#error "MatrixInversion: the module cannot be compiled without BLAS/LAPACK libraries."
#endif

MODULE MatrixInversion
#include "preprocessoptions.cpp"
   IMPLICIT NONE

   PRIVATE

   ! Logical variable to check the setup status of the module
   ! Module is ready to use when the initialisation subroutine has been called, defining the necessary threshold and stuff
   LOGICAL :: ModuleIsSetup = .FALSE.


   ! Definitions of the various type of regularization implemented
   INTEGER, PUBLIC, PARAMETER :: EIGENVALUES = 1,  & ! small eigenvalues are increased
                                 TIKHONOV    = 2,  & ! tikhonov pseudo-inverse is computed
                                 PSEUDOINV   = 3     ! pseudo-inverse


   ! Arrays with the definitions of the regularization setup for each type of matrix which is defined in the calling program
   CHARACTER(20), DIMENSION(:), ALLOCATABLE :: MatrixLabels
   INTEGER, DIMENSION(:), ALLOCATABLE       :: MatrixRegulType
   REAL, DIMENSION(:), ALLOCATABLE          :: MatrixSmallEps


   ! Public subroutines of the module
   PUBLIC :: MatrixInversionSetup           ! Module setup subroutine
   PUBLIC :: MatrixInversionDispose         ! Module memory deallocation subroutine
   PUBLIC :: MatrixInversionDo              ! Compute inverse matrix, with regularization if needed (based on condition number)


   CONTAINS

!===========================================================================================================
!                               PUBLIC SUBROUTINES AND FUNCTIONS
!===========================================================================================================


!*******************************************************************************
!          MatrixInversionSetup
!*******************************************************************************
!> Setup all the necessary parameters and stuff to diagnose regularization
!> problems and compute inverse of (in case) regularized matrices.
!>
!> @param NrMatrixType   Nr of types of matrices, for which different setup will be defined
!> @param Labels         Vector of strings with the labels of the different type of matrices
!> @param RegulType      Vector of integers with the type of regularization adopted per each matrix
!> @param SmallEps       Vector of real numbers with the relevant threshold of the regularization method
!*******************************************************************************
   SUBROUTINE MatrixInversionSetup( NrMatrixType, InpLabels, InpRegulType, InpSmallEps )
      IMPLICIT NONE
      INTEGER, INTENT(IN)                               :: NrMatrixType
      CHARACTER(*), DIMENSION(NrMatrixType), INTENT(IN) :: InpLabels
      INTEGER, DIMENSION(NrMatrixType), INTENT(IN)      :: InpRegulType
      REAL, DIMENSION(NrMatrixType), INTENT(IN)         :: InpSmallEps
      INTEGER :: i

      ! give error if module has been already set up
      CALL ERROR( ModuleIsSetup, " MatrixInversionSetup: module is already initialized ", ERR_MODULE_SETUP )

      ! Allocate internal array to store the setup variables
      ALLOCATE( MatrixLabels(NrMatrixType), MatrixRegulType(NrMatrixType), MatrixSmallEps(NrMatrixType) )

      ! Store the names of the types of inverted matrices
      DO i = 1, NrMatrixType
         MatrixLabels(i) = TRIM(ADJUSTL( InpLabels(i) ))
      END DO

      ! Define the type of regularization to be used per each matrix type
      DO i = 1, NrMatrixType
         CALL ERROR( InpRegulType(i) /= EIGENVALUES .AND. InpRegulType(i) /= TIKHONOV .AND. InpRegulType(i) /= PSEUDOINV, &
                     " MatrixInversionSetup: wrong input type of regularization ", ERR_WRONG_INP )
         MatrixRegulType(i) = InpRegulType(i)
      END DO

      ! Store the regularization threshold
      DO i = 1, NrMatrixType
         MatrixSmallEps(i) = InpSmallEps(i)
      END DO

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "MatrixInversionSetup: module has been setup with: "
      __WRITE_TO_LOG " ",NrMatrixType," types of matrix to be inverted "
      DO i = 1, NrMatrixType
         SELECT CASE( MatrixRegulType(i) )
            CASE( EIGENVALUES )
               __WRITE_TO_LOG MatrixLabels(i), "  eigenvalues regularization with eps = ", MatrixSmallEps(i)
            CASE( TIKHONOV )
               __WRITE_TO_LOG MatrixLabels(i), "  Tikhonov regularization with eps = ", MatrixSmallEps(i)
         END SELECT
      END DO
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif

      ! Now the module is ready to be used
      ModuleIsSetup = .TRUE.

   END SUBROUTINE MatrixInversionSetup


!*******************************************************************************
!          MatrixInversionDispose
!*******************************************************************************
!> Deallocate memory allocated by MatrixInversionSetup.
!>
!*******************************************************************************
   SUBROUTINE MatrixInversionDispose(  )
      IMPLICIT NONE

      ! Deallocate internal arrays
      IF ( ALLOCATED(MatrixLabels) )    DEALLOCATE( MatrixLabels )
      IF ( ALLOCATED(MatrixRegulType) ) DEALLOCATE( MatrixRegulType )
      IF ( ALLOCATED(MatrixSmallEps) )  DEALLOCATE( MatrixSmallEps )

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "MatrixInversionDispose: module internal memory was deallocated. "
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif

      ! Now the memory of the module is not allocated
      ModuleIsSetup = .FALSE.

   END SUBROUTINE MatrixInversionDispose


!*******************************************************************************
!          MatrixInversionDo
!*******************************************************************************
!> Compute the inverse matrix of an input matrix, in case using a regularization
!> scheme defined with the initialization call to MatrixInversionSetup.
!> The matrix can optionally give back to the caller the conditioning number
!> of the matrix, but not that the condition number is not computed in all
!> the regularization schemes. When it is not available, a dummy 0.0 value is returned.
!> When the integer array InvertMask is given, the matrix which is actually
!> inverted is the square submatrix defined by the indices of InvertMask.
!> In the subroutine returned by MatrixInversionDo, in this case,
!> the elements which are not computed by inversion are set equal to a unit matrix.
!> IMPORTANT!!!!! only the upper half of the matrix is computed by MatrixInversionDo!!!!
!>
!> @param   InpMatrix  Input matrix
!> @param   SetupNr    Int nr defining the setup to be used for regularization
!> @param   InvCondNr  optional, on output gives the inverse condition number of the matrix
!> @param   InvertMask optional int array with the indices of the sub-matrix to actually invert
!> @returns SetupNr    Int nr defining the setup to be used for regularization
!*******************************************************************************
   SUBROUTINE MatrixInversionDo( InpFullMatrix, FullInverse, SetupNr, InvCondNr, InvertMask )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN)         :: InpFullMatrix
      COMPLEX , DIMENSION(:,:), INTENT(OUT)       :: FullInverse
      INTEGER, INTENT(IN)                         :: SetupNr
      REAL, INTENT(OUT), OPTIONAL                 :: InvCondNr
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: InvertMask

      COMPLEX, DIMENSION(:,:), ALLOCATABLE                   :: InpMatrix, Inverse
      INTEGER(SHORT_INTEGER_KIND), DIMENSION(:), ALLOCATABLE :: Pivot
      COMPLEX , DIMENSION(:), ALLOCATABLE                    :: Workspace
      INTEGER( SHORT_INTEGER_KIND )                          :: NShort, Info
      LOGICAL :: DoRegularization
      REAL    :: MaxDiag, CondNr
      INTEGER :: i, j, N

      ! give error if module is not setup
      CALL ERROR( .NOT. ModuleIsSetup, " MatrixInversionDo: module has not been initialized yet ", ERR_MODULE_SETUP )
      ! Check that the input SetupNr actually corresponds to a defined regularization setup
      CALL ERROR( SetupNr > SIZE(MatrixRegulType) .OR. SetupNr <= 0," MatrixInversionDo: incorrect input SetupNr ", ERR_SUBS_INPUT )

      N = SIZE(InpFullMatrix,1)
      CALL ERROR( N /= SIZE(InpFullMatrix,2), " MatrixInversionDo: input matrix is not square ", ERR_SUBS_INPUT )
      CALL ERROR( N /= SIZE(FullInverse,1),   " MatrixInversionDo: array dimension mismatch ", ERR_SUBS_INPUT )
      CALL ERROR( N /= SIZE(FullInverse,2),   " MatrixInversionDo: array dimension mismatch ", ERR_SUBS_INPUT )

      IF ( N == 1 ) THEN
         FullInverse(1,1)  = CMPLX(1.0,0.0) / InpFullMatrix(1,1)
         RETURN
      ENDIF

      ! start clock
      CALL StartTimer(MatrixInversionClock)

      ! Allocate memory
      IF ( PRESENT(InvertMask) ) THEN
         N = SIZE(InvertMask)
      ELSE
         N = SIZE(InpFullMatrix,1)
      END IF
      ALLOCATE( InpMatrix(N,N), Inverse(N,N), Pivot(N), Workspace(N) )

      ! Copy sub-matrix which needs to be inverted
      IF ( PRESENT(InvertMask) ) THEN
         InpMatrix = InpFullMatrix(InvertMask,InvertMask)
      ELSE
         InpMatrix = InpFullMatrix
      END IF

      ! Decide how to operate based on the value of MatrixRegulType(SetupNr), which gives the type of regularization
      SELECT CASE( MatrixRegulType(SetupNr) )

         CASE( TIKHONOV )          ! tikhonov regularized inverse is used regardless of the invertibility of the input matrix

            ! Call to the pseudo-inverse tikhonov computing routine
            CALL RegularInversion_Tikhonov( InpMatrix, Inverse, MatrixSmallEps(SetupNr) )
            ! no need to compute condition number, which is set to 0.0 by default
            CondNr = 1.E+99

         CASE( EIGENVALUES, PSEUDOINV )               ! in this case first perform a matrix analysis,
                                                      ! if it is bad conditioned, use regularization
            ! Initialize regularization check
            DoRegularization = .FALSE.

            ! First check the diagonal elements, if they are all below the threshold, perform regularization
            MaxDiag = 0.0
            DO i = 1, size(InpMatrix, 1)
               IF ( ABS(InpMatrix(i,i)) > MaxDiag )  MaxDiag = ABS(InpMatrix(i,i))
            END DO
            IF ( MaxDiag < MatrixSmallEps(SetupNr) ) DoRegularization = .TRUE.

            IF ( .NOT. DoRegularization ) THEN
               ! First compute block factorization of the matrix, if factorization is succesful compute condition number
               IF ( ComputeBlockFactorization( InpMatrix, Inverse, Pivot ) ) THEN
                  CondNr = ComputeInverseConditionNumber( InpMatrix, Inverse, Pivot )
               ELSE  ! otherwise, assume that the matrix is bad conditioned
                  CondNr = 1.E+99
               END IF

               ! If the condition number of the matrix is below the threshold, perform regularization
               IF ( CondNr < MatrixSmallEps(SetupNr) ) DoRegularization = .TRUE.
            END IF

            ! If everything went smoothly up to now, compute the inverse with ZHETRI
            IF ( .NOT. DoRegularization ) THEN
               NShort = size(InpMatrix, 1)
               CALL ZHETRI( "Upper", NShort, Inverse, NShort, Pivot, Workspace, Info )
               ! again, if ZHETRI failed, enforce regularization
               IF ( Info /= 0 ) DoRegularization = .TRUE.
            END IF

            ! Now, if DoRegularization is still .FALSE., the upper triangle of Inverse matrix contains the inverse
            ! otherwise, we need to compute the inverse with RegularInversion_Eigenvalues
            IF ( DoRegularization .AND. MatrixRegulType(SetupNr) == EIGENVALUES ) THEN
                  CALL RegularInversion_Eigenvalues( InpMatrix, Inverse, MatrixSmallEps(SetupNr) )
            ELSE IF ( DoRegularization .AND. MatrixRegulType(SetupNr) == PSEUDOINV ) THEN
                  CALL RegularInversion_PseudoInverse( InpMatrix, Inverse, MatrixSmallEps(SetupNr) )
            END IF

      END SELECT

      ! Copy inverse matrix to the output array, including 1.0 and 0.0 where needed
      IF ( PRESENT(InvertMask) ) THEN
         ! inverse is initialized as Identity matrix
         ForAll(i = 1:SIZE(FullInverse,1), j = 1:SIZE(FullInverse,2)) FullInverse(i,j) = CMPLX(0.0,0.0)
         ForAll(i = 1:SIZE(FullInverse,1)) FullInverse(i,i) = CMPLX(1.0,0.0)
         ! and the rest is given by Inverse
         ForAll(i = 1:N, j = 1:N) FullInverse(InvertMask(i),InvertMask(j)) = Inverse(i,j)
      ELSE
         FullInverse = Inverse
      END IF

      IF ( PRESENT(InvCondNr) ) InvCondNr = CondNr

      ! deallocate memory
      DEALLOCATE( InpMatrix, Inverse, Pivot, Workspace )

      ! stop clock
      CALL StopTimer(MatrixInversionClock)

  END SUBROUTINE MatrixInversionDo


!===========================================================================================================
!                               PRIVATE SUBROUTINES AND FUNCTIONS
!===========================================================================================================


!*******************************************************************************
!          RegularInversion_Tikhonov
!*******************************************************************************
!> Computes the tikhonov pseudo-inverse for a hermitian matrix, according
!> to the formula:
!>  tilde{A}^-1 = (A**H * A + G**H * G)^-1 * A**H
!>
!> @param InpMatrix            Input matrix to invert (full matrix)
!> @param RegFactor            Diagonal element of the G**H*G matrix
!> @returns RegInverse         Tikhonov regularized inverse tilde{A}^-1
!*******************************************************************************
   SUBROUTINE RegularInversion_Tikhonov( InpMatrix, RegInverse, RegFactor )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN)  :: InpMatrix
      COMPLEX, DIMENSION(:,:), INTENT(OUT) :: RegInverse
      REAL , INTENT(IN) :: RegFactor

      COMPLEX , DIMENSION(SIZE(InpMatrix,1),SIZE(InpMatrix,1)) :: TempMat
      INTEGER :: N, i, j, k
      INTEGER(SHORT_INTEGER_KIND) :: Info, NShort
      
      ! dimension of the input matrix
      N = SIZE(InpMatrix,1); NShort = N

  
      ! computing A**2 + G
      DO j = 1, N
         ! diagonal element
         TempMat(j,j) = CMPLX(RegFactor, 0.0)
         DO k = 1, N
            IF ( k <= j ) THEN
               TempMat(j,j) = TempMat(j,j) + ABS(InpMatrix(k,j))**2
            ELSE
               TempMat(j,j) = TempMat(j,j) + ABS(InpMatrix(j,k))**2
            END IF         
         END DO
         ! lower upper triangle
         DO i = 1, j-1
            TempMat(i,j) = CMPLX(0.0, 0.0)
            DO k = 1, N
               IF ( k <= j .AND. k <= i ) THEN
                  TempMat(i,j) = TempMat(i,j) + CONJG(InpMatrix(k,i))*InpMatrix(k,j)
               ELSE IF ( k > j .AND. k <= i ) THEN
                  TempMat(i,j) = TempMat(i,j) + CONJG(InpMatrix(k,i))*CONJG(InpMatrix(j,k))
               ELSE IF ( k <= j .AND. k > i ) THEN
                  TempMat(i,j) = TempMat(i,j) + InpMatrix(i,k)*InpMatrix(k,j)
               ELSE
                  TempMat(i,j) = TempMat(i,j) + InpMatrix(i,k)*CONJG(InpMatrix(j,k))
               END IF
            END DO
         END DO
      END DO
      
!       CALL ZHERK('U','C', NShort,NShort, 1.0,FullMat,NShort, 1.0,TempMat,NShort )

      ! computing the pseudo inverse using cholesky factorization
      ! to solve the complex system:  ( A^H A + G )* X = A^H
      ! on exit of ZPOSV, RegInverse contains the solution X = tilde{A}^-1

      RegInverse = InpMatrix
      DO j = 1, N
         DO i = 1, j-1
            RegInverse(j,i) = CONJG(InpMatrix(i,j))
         END DO
      END DO
      CALL ZPOSV('U', NShort,NShort, TempMat,NShort, RegInverse,NShort, Info)

      ! give error if the linear system solution subroutine failed
      CALL ERROR( Info /= 0, " RegularInversion_Tikhonov: could not solve linear system with ZPOSV ", ERR_BLAS_LAPACK, Info )

   END SUBROUTINE RegularInversion_Tikhonov


!*******************************************************************************
!          RegularInversion_Eigenvalues
!*******************************************************************************
!> Returns the regularized pseudo-inverse for a hermitian matrix, computing
!> the eigen-decomposition of the matrix and adjusting the small eigenvalues
!> ( EigenVal < 64. * RegFactor ) with the formula:
!> EigenVal = EigenVal + RegFactor * EXP(-EigenVal/RegFactor)
!>
!> @param InpMatrix            Input matrix to invert (full matrix)
!> @param RegFactor            RegFactor, threshold for regularizing eigenvalues
!> @param RegInverse           regularized inverse matrix
!*******************************************************************************
   SUBROUTINE RegularInversion_Eigenvalues( InpMatrix, RegInverse, RegFactor )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN)    :: InpMatrix
      COMPLEX, DIMENSION(:,:), INTENT(OUT)   :: RegInverse
      REAL , INTENT(IN)                      :: RegFactor

      INTEGER :: N, i, k
      COMPLEX, DIMENSION(SIZE(InpMatrix,1),SIZE(InpMatrix,1)) :: EigenVectors
      REAL, DIMENSION(SIZE(InpMatrix,1))                      :: EigenValues

      COMPLEX, DIMENSION(:), ALLOCATABLE     :: Workspace
      COMPLEX, DIMENSION(1)                  :: OptDim
      INTEGER( SHORT_INTEGER_KIND )          :: NShort, Info, LWork
      CHARACTER(300)                         :: ErrMsg
      REAL, DIMENSION(SIZE(InpMatrix,1)*3-2) :: RealWork

      ! dimension of the input matrix
      N = SIZE(InpMatrix,1)

      ! Compute Eigenvectors and Eigenvalues
      EigenVectors = InpMatrix
      ! first find optimal dimension of Workspace
      LWork = -1; NShort = N
      CALL ZHEEV( 'Vectors', 'Upper', NShort, EigenVectors, NShort, EigenValues, OptDim, LWork, RealWork, Info )
      ! Now allocate workspace  and compute diagonalization
      LWork = INT( REAL(OptDim(1)) ); ALLOCATE( Workspace( LWork ) )
      CALL ZHEEV( 'Vectors', 'Upper', NShort, EigenVectors, NShort, EigenValues, Workspace, LWork, RealWork, Info )
      DEALLOCATE( Workspace )

      ! Error message and warning in case diagonalization did not work
      IF ( Info < 0 ) THEN
         WRITE(ErrMsg, *) " RegularInversion_Eigenvalues: the ",-Info,"-th argument had an illegal value."
         CALL AbortWithError( ErrMsg )
      END IF
      IF ( Info > 0 ) THEN
         WRITE(ErrMsg, "(A,I7,A)") " RegularInversion_Eigenvalues: the algorithm failed to converge; ",Info, &
                " off-diagonal elements of an intermediate tridiagonal form did not converge to zero."
         CALL ShowWarning( ErrMsg )
         STOP
      END IF

      ! Now substitute the small eigenvalues with a regularized function of them
      DO i = 1, N
         IF ( ABS(EigenValues(i)) < 64. * RegFactor )  &
            EigenValues(i) = EigenValues(i) + RegFactor * EXP(-EigenValues(i)/RegFactor)
      END DO

      ! Reconstruct inverse matrix by Eigenvectors * Diag_Eigenvalues * Eigenvectors^dagger
      DO k = 1, N
         EigenVectors(:,k) = EigenVectors(:,k) / SQRT(EigenValues(k))
      END DO
      CALL ZHERK('U','N', NShort, NShort, 1.0, EigenVectors, NShort, 0.0, RegInverse, NShort )
!       DO j = 1, N
!          DO i = j+1, N
!             RegInverse(i,j) = CONJG( RegInverse(j,i) )
!          END DO
!       END DO

   END SUBROUTINE RegularInversion_Eigenvalues


!*******************************************************************************
!          RegularInversion_PseudoInverse
!*******************************************************************************
!>
!>
!> @param InpMatrix            Input matrix to invert (full matrix)
!> @param RegFactor            RegFactor, threshold for regularizing eigenvalues
!> @param RegInverse           regularized inverse matrix
!*******************************************************************************
   SUBROUTINE RegularInversion_PseudoInverse( InpMatrix, RegInverse, RegFactor )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN)    :: InpMatrix
      COMPLEX, DIMENSION(:,:), INTENT(OUT)   :: RegInverse
      REAL , INTENT(IN)                      :: RegFactor

      INTEGER :: N, i, k
      COMPLEX, DIMENSION(SIZE(InpMatrix,1),SIZE(InpMatrix,1)) :: EigenVectors
      REAL, DIMENSION(SIZE(InpMatrix,1))                      :: EigenValues

      COMPLEX, DIMENSION(:), ALLOCATABLE     :: Workspace
      COMPLEX, DIMENSION(1)                  :: OptDim
      INTEGER( SHORT_INTEGER_KIND )          :: NShort, Info, LWork
      CHARACTER(300)                         :: ErrMsg
      REAL, DIMENSION(SIZE(InpMatrix,1)*3-2) :: RealWork

      ! dimension of the input matrix
      N = SIZE(InpMatrix,1)

      ! Compute Eigenvectors and Eigenvalues
      EigenVectors = InpMatrix
      ! first find optimal dimension of Workspace
      LWork = -1; NShort = N
      CALL ZHEEV( 'Vectors', 'Upper', NShort, EigenVectors, NShort, EigenValues, OptDim, LWork, RealWork, Info )
      ! Now allocate workspace  and compute diagonalization
      LWork = INT( REAL(OptDim(1)) ); ALLOCATE( Workspace( LWork ) )
      CALL ZHEEV( 'Vectors', 'Upper', NShort, EigenVectors, NShort, EigenValues, Workspace, LWork, RealWork, Info )
      DEALLOCATE( Workspace )

      ! Error message and warning in case diagonalization did not work
      IF ( Info < 0 ) THEN
         WRITE(ErrMsg, *) " RegularInversion_PseudoInverse: the ",-Info,"-th argument had an illegal value."
         CALL AbortWithError( ErrMsg )
      END IF
      IF ( Info > 0 ) THEN
         WRITE(ErrMsg, "(A,I7,A)") " RegularInversion_PseudoInverse: the algorithm failed to converge; ",Info, &
                " off-diagonal elements of an intermediate tridiagonal form did not converge to zero."
         CALL ShowWarning( ErrMsg )
      END IF

      ! Now substitute the small eigenvalues with a regularized function of them
      DO i = 1, N
         IF ( ABS(EigenValues(i)) < RegFactor ) THEN
            EigenValues(i) = 0.0
         ELSE
            EigenValues(i) = 1./EigenValues(i)
         END IF
      END DO

      ! Reconstruct inverse matrix by Eigenvectors * Diag_Eigenvalues * Eigenvectors^dagger
      DO k = 1, N
         EigenVectors(:,k) = EigenVectors(:,k) * SQRT(EigenValues(k))
      END DO
      CALL ZHERK('U','N', NShort, NShort, 1.0, EigenVectors, NShort, 0.0, RegInverse, NShort )


   END SUBROUTINE RegularInversion_PseudoInverse


!*******************************************************************************
!          ComputeBlockFactorization
!*******************************************************************************
!> Computes the factorization of a complex hermitian matrix returned by ZHETRF.
!> If any of the steps fail, the function returns FALSE (which is
!> meant to enforce regularization).
!>
!> @param   InpMatrix      Input matrix
!> @param   Factorization  on output give the factorization of InpMatrix
!> @param   Pivot          details on the interchanges of the factorization
!> @returns Check          tells whether factorization was successful
!*******************************************************************************
   LOGICAL FUNCTION ComputeBlockFactorization( InpMatrix, Factorization, Pivot )  RESULT( Check )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN)                    :: InpMatrix
      COMPLEX, DIMENSION(:,:), INTENT(OUT)                   :: Factorization
      INTEGER(SHORT_INTEGER_KIND), DIMENSION(:), INTENT(OUT) :: Pivot

      COMPLEX, DIMENSION(:), ALLOCATABLE     :: Workspace
      COMPLEX, DIMENSION(1)                  :: OptDim
      INTEGER( SHORT_INTEGER_KIND )          :: NShort, Info, LWork

      INTEGER :: N

      ! NShort is a short integer copy of N, to avoid troubles with BLAS/LAPACK version using short integers
      N = SIZE(InpMatrix,1); NShort = N

      ! Initialize Check
      Check = .TRUE.

      ! Work on a copy of the input matrix
      Factorization = InpMatrix

      ! first find optimal dimension of Workspace, by calling ZHETRF with LWork = -1
      LWork = -1;  CALL ZHETRF( 'Upper', NShort, Factorization, NShort, Pivot, OptDim, LWork, Info )
      IF ( Info /= 0 ) Check = .FALSE.

      ! Now allocate workspace  and compute diagonalization
      IF ( Check ) THEN
         LWork = INT( REAL(OptDim(1))+0.5 ); ALLOCATE( Workspace( LWork ) )
         CALL ZHETRF( 'Upper', NShort, Factorization, NShort, Pivot, Workspace, LWork, Info )
         DEALLOCATE( Workspace )
         IF ( Info /= 0 ) Check = .FALSE.      ! if the ZHETRF factorization failed, exit with Check = .FALSE.
      END IF

   END FUNCTION ComputeBlockFactorization


!*******************************************************************************
!          ComputeInverseConditionNumber
!*******************************************************************************
!> Computes the inverse of the condition number of an input matrix,
!> using the LAPACK routine ZHECON. As input needs the factorization
!> of a complex hermitian matrix returned by ZHETRF.
!> If any of the steps fail, the function returns zero (and this is
!> meant to enforce regularization).
!>
!> @param   InpMatrix      Input matrix
!> @param   Factorization  input factorization of InpMatrix, returned by ZHETRF
!> @param   Pivot          details on the interchanges of the factorization, returned by ZHETRF
!> @returns CondNr         Condition number of the input matrix
!*******************************************************************************
   REAL FUNCTION ComputeInverseConditionNumber( InpMatrix, Factorization, Pivot )  RESULT( CondNr )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN)                   :: InpMatrix
      COMPLEX, DIMENSION(:,:), INTENT(IN)                   :: Factorization
      INTEGER(SHORT_INTEGER_KIND), DIMENSION(:), INTENT(IN) :: Pivot

      COMPLEX, DIMENSION(:), ALLOCATABLE     :: Workspace
      INTEGER( SHORT_INTEGER_KIND )          :: NShort, Info
      REAL, DIMENSION(:), ALLOCATABLE        :: RealWork

      INTEGER :: N
      REAL :: Norm, ZLANHE

      ! NShort is a short integer copy of N, to avoid troubles with BLAS/LAPACK version using short integers
      N = SIZE(InpMatrix,1); NShort = N

      ! Use ZLANHE to compute the matrix 1-norm of a complex hermitian matrix, which is needed by ZHECON
      ALLOCATE( RealWork(N) )
      Norm = ZLANHE('1', 'Upper', NShort, InpMatrix, NShort, RealWork)
      DEALLOCATE( RealWork )
      
      ! Now finally use ZHECON to calculate the reciprocal of the condition number
      ALLOCATE( Workspace(2*N) )
      CALL ZHECON( 'Upper', NShort, Factorization, NShort, Pivot, Norm, CondNr, Workspace, Info)
      DEALLOCATE(Workspace)

      IF ( Info /= 0 ) CondNr = 0.0     ! if the ZHECON call failed, exit with Condition Number = 0.0

   END FUNCTION ComputeInverseConditionNumber

!
!
!
!
!   SUBROUTINE count_regularizations(reason,msg)
!
!     IMPLICIT NONE
!
!     INTEGER, INTENT(IN) :: reason
!     INTEGER :: mattyp
!     CHARACTER(len=*), INTENT(IN) :: msg
!     LOGICAL, SAVE :: init=.FALSE.
!
!     IF (.NOT.init) THEN
!       reg_count = 0_ilong
!       init=.true.
!     END IF
!     IF (reason .EQ. -1) THEN
!       WRITE(6,*) "reg reason is -1"
!       STOP
!     ELSE IF (reason .LT. -1 .AND. reason .GT. 6) THEN
!       STOP "reason for regularization count not in range!"
!     END IF
!
!     !which matrix
!     SELECT CASE(msg)
!       CASE('Overlap-Matrix')
!         mattyp = 1
!       CASE('C-Matrix')
!         mattyp = 2
!       CASE DEFAULT
!         mattyp = 3
!     END SELECT
!
!     reg_count(mattyp,reason) = reg_count(mattyp,reason) + iOne
!     ! reason: meaning
!     ! 0: did not regularize
!     ! 1: too small maxdiag
!     ! 2: error getting correct lwork for zhetrf
!     ! 3: error in zhetrf
!     ! 4: error in zhetri
!     ! 5: reziprocal conditioning number < threshold
!     ! 6: error in zhecon
!
!   END SUBROUTINE count_regularizations
!
!   !---------------------------------------------------------------------------------------------------------
!
!   SUBROUTINE print_regul()
!
!     IMPLICIT NONE
!     INTEGER :: i,j, io_unit, f_status
!
!     io_unit = freeunit()
!     OPEN(UNIT=io_unit,FILE='reg_stat.log',STATUS='REPLACE',ACTION='WRITE',IOSTAT=f_status)
!     IF(f_status.NE.0)THEN
!       WRITE(6,'("WARNING: could not write file reg_stat.log")')
!     ELSE
!       WRITE(io_unit,'(5X,"NO REG/REASON",3X,"|",5X,"OVL",5X,"|",4X,"C-MAT",4X,"|",5X,"ANY",5X,"|")')
!       WRITE(io_unit,'(64("="))')
!       WRITE(io_unit,'(2X," No regularization |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,0)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(32(" -"))')
!       WRITE(io_unit,'(2X," too small maxdiag |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,1)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(32(" -"))')
!       WRITE(io_unit,'(2X,"  lwork for zhetrf |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,2)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(32(" -"))')
!       WRITE(io_unit,'(2X,"   error in zhetrf |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,3)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(32(" -"))')
!       WRITE(io_unit,'(2X,"   error in zhetri |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,4)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(32(" -"))')
!       WRITE(io_unit,'(2X,"   rcond too small |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,5)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(32(" -"))')
!       WRITE(io_unit,'(2X,"   error in zhecon |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') reg_count(i,6)
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!
!       WRITE(io_unit,'(64("="))')
!       WRITE(io_unit,'(2X," total # inversion |")',ADVANCE='NO')
!       DO i=1,3
!         WRITE(io_unit,'(2X,I10,1X,"|")',ADVANCE='NO') SUM(reg_count(i,:))
!       END DO
!       WRITE(io_unit,'(1A)',ADVANCE='YES') " "
!     END IF
!
!     CLOSE(io_unit)
!
!   END SUBROUTINE print_regul
!
!  !##################################################################################
!  !#                          Subroutine log_condnumbers                            #
!  !#                                                                                #
!  !# Conditioning Numbers for a Matrix. Not Used at the moment...                   #
!  !#                                                                                #
!  !##################################################################################
!  ! Calls in subroutine LOG_CONDNUMBERS:
!  ! => ZHEEV (on line <3236>)
!  ! ==> ZHEEV (on line <3241>)
!  ! ===> ZHEEV (on line <3276>)
!  ! ====> ZHEEV (on line <3281>)
!  SUBROUTINE log_condnumbers(C_MAT, overlap, time)
!    IMPLICIT NONE
!    COMPLEX , DIMENSION(:,:), INTENT(IN) :: C_MAT
!    COMPLEX , DIMENSION(:,:), INTENT(IN) :: overlap
!    REAL, INTENT(IN) :: time
!
!    INTEGER, PARAMETER :: CMatrixUnit  = 400
!    INTEGER, PARAMETER :: OvMatrixUnit = 401
!    LOGICAL :: OpenedUnit
!
!    INTEGER , SAVE :: step=0
!    INTEGER  :: dim1
!    COMPLEX , DIMENSION(:,:), ALLOCATABLE :: matrix
!    INTEGER  :: e, info
!    REAL  :: epss=10d-14
!
!    REAL , DIMENSION(:), ALLOCATABLE :: eigenval
!    INTEGER      :: lwork
!    COMPLEX , DIMENSION(:), ALLOCATABLE   :: work
!    REAL  , DIMENSION(:), ALLOCATABLE :: rwork
!
!    INQUIRE( UNIT=CMatrixUnit , OPENED=OpenedUnit )
!    IF ( .NOT. OpenedUnit ) THEN
!       OPEN( UNIT=CMatrixUnit, FILE="eigen_cmatrix.log" )
!       WRITE(CMatrixUnit,*) "# conditioning of C-Matrix"
!       WRITE(CMatrixUnit,*) "# step, time, min(eigval), max(eigval), condn "
!    END IF
!    INQUIRE( UNIT=OvMatrixUnit , OPENED=OpenedUnit )
!    IF ( .NOT. OpenedUnit ) THEN
!       OPEN( UNIT=OvMatrixUnit, FILE="eigen_ovmatrix.log" )
!       WRITE(OvMatrixUnit,*) "# conditioning of S-Matrix (overlap)"
!       WRITE(OvMatrixUnit,*) "# step, time, min(eigval), max(eigval), condn "
!    END IF
!
!    step=step+1
!    info = 0
!
!    dim1=SIZE(C_MAT,1)
!
!    ALLOCATE(matrix(1:dim1,1:dim1))
!    matrix=cmplx(dZero,dZero,KIND=dop)
!    matrix=C_MAT(:,:)
!
!    ALLOCATE(rwork(1:dim1*3 - 2))
!    rwork = dZero
!
!    ALLOCATE(eigenval(1:dim1))
!    eigenval=  dZero
!
!    ALLOCATE(work(10))
!    CALL ZHEEV('V','U',dim1,matrix,dim1,eigenval,work,-1,rwork,info)
!    lwork = int(real(work(1)))
!    DEALLOCATE(work)
!    ALLOCATE(work(lwork))
!
!    CALL ZHEEV('V','U',dim1,matrix,dim1,eigenval,work,lwork,rwork,info)
!    IF ( info.NE.0) THEN
!      WRITE(*,*) " Error while computing eigenvectors of C-Matrix"
!    ENDIF
!
! !    !from gaussian.f
! !    DO e=1,dim1
! !      IF( eigenval(e) .LT. 24.d0*epss) THEN
! !        eigenval(e)=eigenval(e)+epss*exp(-eigenval(e)/epss)
! !      ENDIF
! !    ENDDO
!
!    WRITE(CMatrixUnit,'(F16.9,x,E11.4,E11.4,E11.4,E14.6)') &
!      & time, minval(eigenval),maxval(eigenval), SUM(eigenval)/dim1, maxval(eigenval)/minval(eigenval)
!
!    DEALLOCATE(matrix)
!    DEALLOCATE(rwork)
!    DEALLOCATE(work)
!    DEALLOCATE(eigenval)
!
!    info = 0
!
!    dim1=SIZE(overlap,1)
!
!    ALLOCATE(matrix(1:dim1,1:dim1))
!    matrix=cmplx(dZero,dZero,KIND=dop)
!    matrix=overlap(:,:)
!
!    ALLOCATE(rwork(1:dim1*3 - 2))
!    rwork = dZero
!
!    ALLOCATE(eigenval(1:dim1))
!    eigenval=  dZero
!
!    ALLOCATE(work(10))
!    CALL ZHEEV('V','U',dim1,matrix,dim1,eigenval,work,-1,rwork,info)
!    lwork = int(real(work(1)))
!    DEALLOCATE(work)
!    ALLOCATE(work(lwork))
!
!    CALL ZHEEV('V','U',dim1,matrix,dim1,eigenval,work,lwork,rwork,info)
!    IF ( info.NE.0) THEN
!      WRITE(*,*) " Error while computing eigenvectors of S-Matrix"
!    ENDIF
!
! !    !from gaussian.f
! !    DO e=1,dim1
! !      IF( eigenval(e) .LT. 24.d0*epss) THEN
! !        eigenval(e)=eigenval(e)+epss*exp(-eigenval(e)/epss)
! !      ENDIF
! !    ENDDO
!
!    WRITE(OvMatrixUnit,'(F16.9,x,E11.4,E11.4,E11.4,E14.6)') &
!      & time, minval(eigenval),maxval(eigenval), SUM(eigenval)/dim1, maxval(eigenval)/minval(eigenval)
!
!    DEALLOCATE(matrix)
!    DEALLOCATE(rwork)
!    DEALLOCATE(work)
!    DEALLOCATE(eigenval)
!
! END SUBROUTINE log_condnumbers

END MODULE MatrixInversion
