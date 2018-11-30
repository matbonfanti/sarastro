!***************************************************************************************
!*                              MODULE EmptyOptimize
!***************************************************************************************
!
!>  \brief     Optimize the initial empty configurations of a gaussian wavefunction
!>  \details   This module contains the subroutines that take as input an initial \n
!>             gaussian wavefunction and optimize the unoccupied gaussian configurations \n
!>             according to a 2nd order estimate of the short time propogation error.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             28 September 2017
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          .............
!
!***************************************************************************************
MODULE EmptyOptimize
#include "preprocessoptions.cpp"
   USE PsiObject
   USE GauConf
   USE OperatorDefine
   USE MatrixInversion

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OptimizeEmptyGaussians

   REAL, PARAMETER :: GuessOverlap = 0.8


   CONTAINS

!===========================================================================================================
!                               PUBLIC SUBROUTINES AND FUNCTIONS
!===========================================================================================================


   SUBROUTINE OptimizeEmptyGaussians( Psi, Hamiltonian )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      TYPE(OperatorData), INTENT(IN)   :: Hamiltonian

      COMPLEX, ALLOCATABLE, DIMENSION(:,:)   :: OptBVector   !< array to store the coefficients of the optimized configurations
      COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: OptGaussPar  !< array to store the gaussian parameters of the optimized configurations
      TYPE(WaveFunct)                        :: RedPsi       !< temporary psi object to store the occupied configs + the optimized unocc ones
      COMPLEX, DIMENSION(3,Psi%GDim)       :: SglConfig

      INTEGER :: iCfg, iOrdCfg, NOccCfg

      REAL :: ErrorEstimate, Lambda, DeltaPlus, DeltaMinus
      COMPLEX, DIMENSION(Psi%GDim) :: Derivatives, OldXi, OldDerivs
      COMPLEX, DIMENSION(:), ALLOCATABLE :: Derivs

      ! FOR THE TIME BEING A SINGLE STATE CALCULATION IS ASSUMED

      ! allocate temporary arrays to store optimized values of the wavefunction gaussian parameters
      ALLOCATE( OptBVector(Psi%NrCfg, Psi%NrStates), OptGaussPar(3, Psi%NrPrimGau, Psi%NrGauSets) )

      ! EXTRACT OCCUPIED CONFIGURATIONS
      iOrdCfg = 0
      DO iCfg = 1,Psi%NrCfg
         IF ( ABS(Psi%Bvector(iCfg, 1)) > 1.E-8 ) THEN
            iOrdCfg = iOrdCfg + 1
            OptBVector(iOrdCfg, 1) = Psi%Bvector(iCfg, 1)
            OptGaussPar(:, (iOrdCfg-1)*Psi%GDim+1:iOrdCfg*Psi%GDim, 1) = Psi%GaussPar(:, (iCfg-1)*Psi%GDim+1:iCfg*Psi%GDim, 1)
         END IF
      END DO

      ! Store the number of the occupied configurations
      NOccCfg = iOrdCfg

      ! Append the unoccupied ones
      DO iCfg = 1,Psi%NrCfg
         IF ( ABS(Psi%Bvector(iCfg, 1)) <= 1.E-8 ) THEN
            iOrdCfg = iOrdCfg + 1
            OptBVector(iOrdCfg, 1) = Psi%Bvector(iCfg, 1)
            OptGaussPar(:, (iOrdCfg-1)*Psi%GDim+1:iOrdCfg*Psi%GDim, 1) = Psi%GaussPar(:, (iCfg-1)*Psi%GDim+1:iCfg*Psi%GDim, 1)
         END IF
      END DO

      ! wavefunction with "core" state which are not optimized
      CALL SetPsi( RedPsi, 1, NOccCfg, Psi%GDim, Psi%WFType )
      RedPsi%Bvector(:,1)  = OptBVector(1:NOccCfg, 1)
      RedPsi%GaussPar(:,:,1) = OptGaussPar(:, 1:NOccCfg*Psi%GDim, 1)
      CALL UpdateWaveFunctInstance( RedPsi )

      ! OPTIMIZE THE CURRENT CONFIGURATION
      CALL OptimizeEmptySet( RedPsi, Psi%NrCfg-NOccCfg, OptGaussPar(:, NOccCfg*Psi%GDim+1:Psi%NrCfg*Psi%GDim, 1), Hamiltonian )

      ! REDEFINE BVEC AND GAUSSPAR OF INITPSI WITH THE OPTIMIZED VALUES
      Psi%Bvector(:,1) =  OptBVector(:, 1)
      Psi%GaussPar(:, :, 1) = OptGaussPar(:, :, 1)
      ! UPDATE OTHER ARRAYS OF INITPSI
      CALL UpdateWaveFunctInstance( Psi  )

      ! DEALLOCATE MEMORY
      DEALLOCATE( OptBVector, OptGaussPar )
      CALL DisposePsi( RedPsi )

   END SUBROUTINE OptimizeEmptyGaussians


!===========================================================================================================
!                               PRIVATE SUBROUTINES AND FUNCTIONS
!===========================================================================================================

   SUBROUTINE OptimizeEmptySet( Psi, NEmpty, EmptyConfig, Hamiltonian )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)                          :: Psi
      INTEGER, INTENT(IN)                                     :: NEmpty
      COMPLEX, DIMENSION(3,Psi%GDim*NEmpty), INTENT(INOUT)  :: EmptyConfig
      TYPE(OperatorData), INTENT(IN)                          :: Hamiltonian

      COMPLEX, DIMENSION(2,Psi%GDim*NEmpty) :: Derivs, OldXi, OldDerivs
      REAL :: ErrorEstimate, OldErrorEstimate, Lambda
      INTEGER :: i,j

!       DO iOrdCfg = -50.,50
!          OptGaussPar(2, NOccCfg*Psi%GDim+1, 1) = CMPLX( REAL(iOrdCfg)*0.01, 0.0 )
!          CALL GauConf_Normalize(OptGaussPar(:, NOccCfg*Psi%GDim+1:, 1))
!          CALL EmptySetErrorEstimate( RedPsi, Psi%NrCfg-NOccCfg, OptGaussPar(:, NOccCfg*Psi%GDim+1:, 1), Hamiltonian, ErrorEstimate )
!          CALL ErrorEstimateDerivative( RedPsi, Psi%NrCfg-NOccCfg, OptGaussPar(:, NOccCfg*Psi%GDim+1:, 1), Hamiltonian, Derivs )
!          WRITE(900,"(100F20.8)") REAL(iOrdCfg), ErrorEstimate, REAL(Derivs), AIMAG(Derivs)
!       END DO

      CALL EmptySetErrorEstimate( Psi, NEmpty, EmptyConfig, Hamiltonian, ErrorEstimate, .FALSE. )
      CALL ErrorEstimateDerivative( Psi, NEmpty, EmptyConfig, Hamiltonian, Derivs )

      WRITE(901,"(100F30.16)") REAL(0), ErrorEstimate, MAXVAL(ABS(Derivs(:,:))), 0.0

      ! Store the previous coordinates and the previous derivatives
      OldXi(:,:) = EmptyConfig(1:2, :)
      OldDerivs(:,:) = Derivs(:,:)
      OldErrorEstimate = ErrorEstimate

      ! First step with lambda = 100.0
      Lambda = 100.0

      DO i = 1,200

         DO
            ! Move along gradient
            EmptyConfig(1:2,:) = OldXi + Lambda * OldDerivs
            CALL GauConf_Normalize(EmptyConfig)
            ! compute new ErrorEstimate and derivatives
            CALL EmptySetErrorEstimate( Psi, NEmpty, EmptyConfig, Hamiltonian, ErrorEstimate, .FALSE. )

            IF ( ErrorEstimate > OldErrorEstimate ) THEN
               EXIT
            ELSE
               Lambda = Lambda * 0.7
            END IF
            IF ( Lambda < 1.E-8 ) EXIT
         END DO
         CALL ErrorEstimateDerivative( Psi, NEmpty, EmptyConfig, Hamiltonian, Derivs )

         WRITE(901,"(100F30.16)") REAL(i), ErrorEstimate, MAXVAL(ABS(Derivs(:,:))), Lambda

         IF ( MAXVAL(ABS(Derivs(:,:))) < 1.E-8 .OR. Lambda < 1.E-8 ) EXIT

         ! We now have new and old derivatives ( Derivs, OldDerivs ) and new and old coordinates ( EmptyConfig, OldXi ), compute lambda
         Lambda = ComputeLambda( EmptyConfig(1:2,:), OldXi, Derivs, OldDerivs )
         IF ( Lambda < 0.0 ) Lambda = 100.0

         ! Store the previous coordinates and the previous derivatives
         OldXi(:,:) = EmptyConfig(1:2, :)
         OldDerivs(:,:) = Derivs(:,:)
         OldErrorEstimate = ErrorEstimate

!          DO j = 1, Psi%GDim*NEmpty
!             WRITE(900,"(100F20.8)") REAL(EmptyConfig(2,j)) / 2.0 / REAL(EmptyConfig(1,j)), AIMAG(EmptyConfig(2,j))
!             WRITE(900,"(100F20.8)") REAL(EmptyConfig(2,j)), AIMAG(EmptyConfig(2,j))
!          END DO
!          WRITE(900,*) "  "

      END DO

   END SUBROUTINE OptimizeEmptySet


   REAL FUNCTION ComputeLambda( NewXi, OldXi, NewDerivs, OldDerivs  )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN) ::  NewXi, OldXi, NewDerivs, OldDerivs
      REAL :: GradTimesX, GradTimesGrad
      COMPLEX :: DeltaX, DeltaG
      INTEGER :: j, i

      GradTimesX = 0.0
      GradTimesGrad = 0.0

      DO i = 1, SIZE(NewXi,2)
         DO j = 1, SIZE(NewXi,1)
            DeltaX = NewXi(j,i) - OldXi(j,i)
            DeltaG = NewDerivs(j,i) - OldDerivs(j,i)
            GradTimesX = GradTimesX + REAL(DeltaG)*REAL(DeltaX) + AIMAG(DeltaG)*AIMAG(DeltaX)
            GradTimesGrad = GradTimesGrad + REAL(DeltaG)*REAL(DeltaG) + AIMAG(DeltaG)*AIMAG(DeltaG)
         END DO
      END DO
      ComputeLambda = - GradTimesX / GradTimesGrad

   END FUNCTION ComputeLambda


   SUBROUTINE ErrorEstimateDerivative( Psi, NEmpty, EmptyConfig, Hamiltonian, Derivative )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)                          :: Psi
      INTEGER, INTENT(IN)                                     :: NEmpty
      COMPLEX, DIMENSION(3,Psi%GDim*NEmpty), INTENT(IN)     :: EmptyConfig
      TYPE(OperatorData), INTENT(IN)                          :: Hamiltonian
      COMPLEX, DIMENSION(2,Psi%GDim*NEmpty), INTENT(OUT)    :: Derivative

      COMPLEX, DIMENSION(3,Psi%GDim*NEmpty)       :: DisplConfig
      REAL :: ErrorEstimate
      INTEGER :: i, j, k

      ! Finite difference - 4 points formula for first derivative
      !> Displacements in units of delta for 4pts first derivative formula
      REAL, DIMENSION(4), PARAMETER :: DeltasI = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
      !> Coefficients for 4pts first derivative formula
      REAL, DIMENSION(4), PARAMETER :: CoeffsI = (/ +1./12., -8./12., +8./12., -1./12. /)
      ! Small displacement delta
      REAL, DIMENSION(2), PARAMETER :: SmallDelta = (/ 0.00001, 0.001 /)

      ! Initialize gradient to 0
      Derivative(:,:) = CMPLX(0.0, 0.0)

      DO i = 1, Psi%GDim*NEmpty          ! Cycle over the primitive gaussians to displace
         DO j = 1,2                        ! loop over the paramters to optimize
            DO k = 1, size(DeltasI)         !  Cycle over the finite displacements

               ! REAL PART OF THE DERIVATIVE

               ! Define small displacement from the point where compute the derivative
               DisplConfig(:, :) = EmptyConfig
               DisplConfig(j, i) = DisplConfig(j, i) + CMPLX( DeltasI(k)*SmallDelta(j), 0.0 )
               CALL GauPrim_Normalize( DisplConfig(:, i) )

               ! Compute the error estimate in the displaced coordinate
               CALL EmptySetErrorEstimate( Psi, NEmpty, DisplConfig, Hamiltonian, ErrorEstimate, .FALSE. )
               ! Increment numerical derivative of the analytical derivative
               Derivative(j,i) = Derivative(j,i) + CMPLX( CoeffsI(k)*ErrorEstimate, 0.0 )

               ! IMAGINARY PART OF THE DERIVATIVE

               ! Define small displacement from the point where compute the derivative
               DisplConfig(:, :) = EmptyConfig
               DisplConfig(j, i) = DisplConfig(j, i) + CMPLX( 0.0, DeltasI(k)*SmallDelta(j) )
               CALL GauPrim_Normalize( DisplConfig(:, i) )

               ! Compute the error estimate in the displaced coordinate
               CALL EmptySetErrorEstimate( Psi, NEmpty, DisplConfig, Hamiltonian, ErrorEstimate, .FALSE. )
               ! Increment numerical derivative of the analytical derivative
               Derivative(j,i) = Derivative(j,i) + CMPLX( 0.0, CoeffsI(k)*ErrorEstimate )

            END DO
         END DO
      END DO

      DO j = 1,2                        ! loop over the paramters to optimize
         Derivative(j,:) = Derivative(j,:)/SmallDelta(j)
      END DO

   END SUBROUTINE ErrorEstimateDerivative


   SUBROUTINE EmptySetErrorEstimate( Psi, NEmpty, EmptyConfig, Hamiltonian, ErrorEstimate, Check )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)                          :: Psi
      INTEGER, INTENT(IN)                                     :: NEmpty
      COMPLEX, DIMENSION(3,Psi%GDim*NEmpty), INTENT(IN)     :: EmptyConfig
      TYPE(OperatorData), INTENT(IN)                          :: Hamiltonian
      REAL, INTENT(OUT)                                       :: ErrorEstimate
      LOGICAL, INTENT(IN)                                     :: Check

      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrGauSets)              ::  OverlapInv
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrStates,Psi%NrStates)  ::  HMatrix
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg)                            ::  InvOverlapTimesHamilt
      COMPLEX, DIMENSION(NEmpty,Psi%NrCfg)                               ::  UnoccOccOverlap
      COMPLEX, DIMENSION(NEmpty,NEmpty)                                  ::  ProjOverlap, ProjOverlapInv
      COMPLEX, DIMENSION(NEmpty,Psi%NrCfg)                               ::  UnoccOccHamiltonian
      COMPLEX, DIMENSION(NEmpty,Psi%NrCfg)                               ::  TmpMatrix
      COMPLEX, DIMENSION(NEmpty)                                         ::  TmpVector

      REAL    :: InvCondNr
      COMPLEX :: Delta
      INTEGER :: i, j, l, m

      REAL :: Cross, OvlpMax
      INTEGER :: imax, jmax

      ! In the notes below, we will indicate:
      !  | g_i >       i = 1,2, ...    configurations of the "starting" wavefunction Psi
      !  | gtilde_l >  k = 1,2,...     additional configuration included, we want to estimate the effect of their inclusion
      !                                on the short time propagation error
      !   P                            the projector on the space spanned by the | g_alpha >

      ! Compute the matrix inverse of the overlap between the configurations | g_i >
      CALL MatrixInversionDo( Psi%Overlap(:,:,1,1), OverlapInv(:,:,1), 1, InvCondNr )
      ! Only upper diagonal has been computed, complete the inverse overlap
      DO j = 1, Psi%NrCfg
         DO i = j+1, Psi%NrCfg
            OverlapInv(i,j,1) = CONJG(OverlapInv(j,i,1))
         END DO
      END DO

      ! Calculate all the Hamiltonian and Overlap matrix elements over the Psi configurations + the additional one
      ! HMatrix = < g_i | H | g_j >                 InvOverlapTimesHamilt = SUM_j [ S^{-1} ]_i,j  < g_j | H | g_k >
      HMatrix = GetOperatorMatrix( Hamiltonian, Psi )
      CALL TheOneWithMatrixMultiplication(InvOverlapTimesHamilt, OverlapInv(:,:,1), HMatrix(:,:,1,1), "N","N")

      ! UnoccOccOverlap = < gtilde_l | g_i > and UnoccOccHamiltonian = < gtilde_l | H | g_i >
      DO i = 1, Psi%NrCfg
         DO l = 1, NEmpty
            UnoccOccOverlap(l,i) = GauConf_Overlap( EmptyConfig(:,  (l-1)*Psi%GDim+1:l*Psi%GDim), &
                                                    Psi%GaussPar(:, (i-1)*Psi%GDim+1:i*Psi%GDim, 1) )
            UnoccOccHamiltonian(l,i) =  OperatorMatrixElement( Hamiltonian, EmptyConfig(:,  (l-1)*Psi%GDim+1:l*Psi%GDim), &
                                                                            Psi%GaussPar(:, (i-1)*Psi%GDim+1:i*Psi%GDim, 1), 1, 1 )
         END DO
      END DO

      ! ProjOverlap = < gtilde_l | 1-P | gtilde_m > and its inverse ProjOverlapInv
      CALL TheOneWithMatrixMultiplication( TmpMatrix, UnoccOccOverlap(:,:), OverlapInv(:,:,1), "N","N")
      CALL TheOneWithMatrixMultiplication( ProjOverlap, TmpMatrix(:,:), UnoccOccOverlap(:,:),  "N","C")
      DO m = 1, NEmpty
         DO l = 1, NEmpty
            ProjOverlap(l,m) = GauConf_Overlap( EmptyConfig(:,  (l-1)*Psi%GDim+1:l*Psi%GDim), &
                                                EmptyConfig(:,  (m-1)*Psi%GDim+1:m*Psi%GDim) ) - ProjOverlap(l,m)
         END DO
      END DO
      IF ( Check ) THEN
         OvlpMax = 0.0
         DO j = 1, NEmpty
            DO i = j+1, NEmpty
               Cross = ABS(ProjOverlap(i,j)) / SQRT(ABS(ProjOverlap(i,i))) / SQRT(ABS(ProjOverlap(j,j)))
               IF (Cross > OvlpMax)  THEN
                  OvlpMax = Cross; imax = i; jmax = j
               END IF
            END DO
         END DO
         WRITE(*,"(2I5,100F9.5)") imax, jmax, OvlpMax, ABS(ProjOverlap(imax,imax)), ABS(ProjOverlap(jmax,jmax))
      END IF

!       CALL  TheOneWithMatrixPrintedLineAfterLine( ABS(ProjOverlap) )
      CALL MatrixInversionDo( ProjOverlap(:,:), ProjOverlapInv(:,:), 3, InvCondNr )

      ! Now construct the integral < gtilde_l | (1-P) H | g_i > and update UnoccOccHamiltonian with it
      CALL TheOneWithMatrixMultiplication( TmpMatrix, UnoccOccOverlap(:,:), InvOverlapTimesHamilt(:,:),  "N","N")
      UnoccOccHamiltonian = UnoccOccHamiltonian - TmpMatrix
      ! Sum with coefficients to get < gtilde_l | (1-P) H | Psi >
      TmpVector = TheOneWithMatrixVectorProduct( UnoccOccHamiltonian, Psi%Bvector(:, 1) )

      ! Compute Delta
      Delta = CMPLX(0.0,0.0)
      DO m = 1, NEmpty
         DO l = 1, NEmpty
            Delta = Delta + CONJG(TmpVector(l)) * ProjOverlapInv(l,m) * TmpVector(m)
         END DO
      END DO

      ErrorEstimate = REAL(Delta)

   END SUBROUTINE EmptySetErrorEstimate


END MODULE EmptyOptimize
