!***************************************************************************************
!*                              MODULE AdaptGBasis
!***************************************************************************************
!
!>  \brief     Adapt the gaussian basis set, by adding or removing gaussian functions    \n
!>  \details   This module contains the subroutines that take as a wavefunction at given \n
!>             step, and based on some given thresholds removes or add gaussian.         \n
!>             Removal is done when linear dependence problems are identified, and it is \n
!>             done by projecting the actual wavefunction on the N-1 gaussian set.       \n
!>             Addition is done when it is seen that the addition of an empty gaussian   \n
!>             significantly reduce the second-order short time estimated error of the   \n
!>             propagation.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             16 October 2017
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          .............
!
!***************************************************************************************
MODULE AdaptGBasis
#include "preprocessoptions.cpp"
   USE PsiObject
   USE GauConf
   USE OperatorDefine
   USE MatrixInversion

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: AddOrRemoveGaussianFunctions

   REAL, PARAMETER :: OverlapThreshold = 1.0
   REAL, PARAMETER :: NormErrorThreshold = 1.E-10

   REAL, PARAMETER :: GuessOverlap = 0.8
   REAL, PARAMETER :: OptimizationThreshold = 1.E+99
   REAL, PARAMETER :: SecondOrderThreshold = 1.E+99

!    REAL, PARAMETER :: OptimizationThreshold = 1.5E-05
!     REAL, PARAMETER :: SecondOrderThreshold = 2.8E-05
   INTEGER, PARAMETER :: NGaussMax = 40

   TYPE(WaveFunct), SAVE :: PsiReduced

   INTEGER, SAVE :: TimeOut

   CONTAINS

!===========================================================================================================
!                               PUBLIC SUBROUTINES AND FUNCTIONS
!===========================================================================================================




   SUBROUTINE AddOrRemoveGaussianFunctions( Psi, Hamiltonian, ActualTime, LogUnit )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      TYPE(OperatorData), INTENT(IN)   :: Hamiltonian
      REAL, INTENT(IN)                 :: ActualTime
      INTEGER, INTENT(IN)              :: LogUnit

      COMPLEX, DIMENSION(3,Psi%GDim) :: GuessConfig, BestGuessConfig
      REAL, DIMENSION(Psi%GDim)      :: Center, Displ
      INTEGER :: iRightCfg, iLeftCfg, j, k, l
      REAL :: NormError, NormReduced, ErrorEstimate, MaxErrorEstimate
      LOGICAL :: Check

      REAL :: BestNorm
      INTEGER :: iRemove

      IF ( TimeOut > 0 ) THEN
         TimeOut = TimeOut-1
         RETURN
      END IF

      CALL ERROR( Psi%WFType == vMCG_MULTISET, " AddOrRemoveGaussianFunctions not implemented for multi-set wavefunctions ")

      !     *******************   ADD GAUSSIANS  ***************************

      !      in some guessed position, check the effect that a gaussian
      !      might have at short time with the second order error estimate,
      !      keep those gaussians that are above a certain threshold,
      !      optimize their position and then include them in the basis set
      !      as unoccupied states


      IF ( Psi%NrCfg < NGaussMax ) THEN

      ! Initialize a variable to store the maximum error reduction
      MaxErrorEstimate = 0.0

      ! Loop over the configurations included in the Psi wavefunction
      DO j = 1,Psi%NrCfg

         ! compute the center of this configuration
         Center = REAL(Psi%GaussPar(2, (j-1)*Psi%GDim+1:j*Psi%GDim, 1)) / 2.0 / REAL(Psi%GaussPar(1, (j-1)*Psi%GDim+1:j*Psi%GDim, 1))
         ! compute the displacements of Re[Xi]
         Displ = - 2.0 * REAL(Psi%GaussPar(1, (j-1)*Psi%GDim+1:j*Psi%GDim, 1)) * &
               SQRT(2.0*LOG(GuessOverlap)/REAL(Psi%GaussPar(1, (j-1)*Psi%GDim+1:j*Psi%GDim, 1)))

         ! Loop over the degrees of freedom of the wavefunction and over positive and negative displacements
         DO k = 1, Psi%GDim
            DO l = -1, 1, 2

               ! translate the j-th configuration along the k-th degree of freedom and normalize it
               GuessConfig = Psi%GaussPar(:, (j-1)*Psi%GDim+1:j*Psi%GDim, 1)
               GuessConfig(2,k) = GuessConfig(2,k) + REAL(l) * CMPLX(Displ(k),0.0)
               CALL GauPrim_Normalize(GuessConfig(:,k))

               ! compute errror decrease estimate
               CALL EmptyGaussianErrorEstimate( Psi, GuessConfig, Hamiltonian, ErrorEstimate )

               ! it the configuration has the maximum error reduction so far, store it
               IF ( ErrorEstimate > MaxErrorEstimate ) THEN
                  BestGuessConfig = GuessConfig
                  MaxErrorEstimate = ErrorEstimate
               END IF

            END DO
         END DO
      END DO

      ! now, if the maximum error decrese estimate is greater than an optimization threshold ...
      IF ( MaxErrorEstimate > OptimizationThreshold ) THEN

         ! ... first optimize the empty state ...
         GuessConfig = BestGuessConfig
         CALL OptimizeEmptyGaussian( Psi, BestGuessConfig, Hamiltonian )
         CALL EmptyGaussianErrorEstimate( Psi, BestGuessConfig, Hamiltonian, ErrorEstimate )

         ! and if the optimized error is greater than the threshold value ...
         IF ( ErrorEstimate > SecondOrderThreshold ) THEN

            WRITE(*,"(A,F9.3,A)") "Time = ",ActualTime/MyConsts_fs2AU, &
                  " : empty gaussian added to the basis set "
            WRITE(LogUnit,"(A,F9.3,A)") "# Time = ",ActualTime/MyConsts_fs2AU, &
                  " : empty gaussian added to the basis set "

            ! ... include it in the gaussian basis set
            CALL AddEmptyGaussianToPsi( Psi, GuessConfig )

            TimeOut = 100

         END IF
      END IF

      END IF

      !     *******************   REMOVE GAUSSIANS  ***************************

      !      check those gaussians with overlap greater than a given threshold
      !      and then remove gaussians only when the norm of the projected
      !      wavefunction is reduced only by a small amount (another given
      !      threshold). The condition of equal norm is equivalent to
      !      overlap < Psi^prime | Psi > = 1.0


      Check = .FALSE.
      outer: DO iRightCfg = 1, Psi%NrCfg
         DO iLeftCfg = iRightCfg+1, Psi%NrCfg

            ! check if overlap is greater then the threshold, in case proceed with removing one gaussian from the basis set
            IF ( ABS(Psi%Overlap(iLeftCfg,iRightCfg,1,1)) > OverlapThreshold ) THEN

               Check = .TRUE.



               BestNorm = 0.0
               DO k = 1, Psi%NrCfg
                  CALL RemoveSingleGaussian(Psi, k, PsiReduced, NormReduced)
                  IF ( NormReduced > BestNorm ) THEN
                     BestNorm = NormReduced
                     iRemove = k
                  END IF
                  CALL DisposePsi( PsiReduced )
               END DO

               NormError = ABS(WFNorm( Psi ) - NormReduced)
               WRITE(800,*) ActualTime/MyConsts_fs2AU, NormError, iRemove
               IF ( NormError < NormErrorThreshold ) THEN
                  WRITE(*,"(A,F9.3,A,F14.6)") "Time = ",ActualTime/MyConsts_fs2AU, &
                        " : two gaussians have been merged, new norm = ", NormReduced
                  WRITE(LogUnit,"(A,F9.3,A,F14.6)") "# Time = ",ActualTime/MyConsts_fs2AU, &
                        " : two gaussians have been merged, new norm = ", NormReduced

                  CALL RemoveSingleGaussian(Psi, iRemove, PsiReduced, NormReduced)
                  CALL DisposePsi( Psi )
                  CALL SetPsi( Psi, PsiReduced%NrStates, PsiReduced%NrCfg, PsiReduced%GDim, PsiReduced%WFType )
                  Psi%Bvector = PsiReduced%Bvector
                  Psi%GaussPar = PsiReduced%GaussPar
                  CALL UpdateWaveFunctInstance( Psi )
                  CALL NormalizePsi( Psi )
                  CALL DisposePsi( PsiReduced )

                  TimeOut = 100

                  EXIT outer

               END IF

!                ! Compute the projection of Psi on the reduced basis set in which nCfg1 and nCfg2 are merged
! !                CALL RemoveOverlappingGaussians( Psi, iLeftCfg, iRightCfg, PsiReduced, NormReduced )
!                IF ( ABS(Psi%BVector(iLeftCfg,1)) > ABS(Psi%BVector(iRightCfg,1)) ) THEN
!                   CALL RemoveSingleGaussian(Psi, iRightCfg, PsiReduced, NormReduced)
!                ELSE
!                   CALL RemoveSingleGaussian(Psi, iLeftCfg, PsiReduced, NormReduced)
!                END IF
!                ! and estimate error of the gaussian merging as ABS(Norm(Psi)-Norm(PsiReduced))
!                NormError = ABS(WFNorm( Psi ) - NormReduced)
!
! !                WRITE(230,"(E18.8,2I10,100E18.8)") ActualTime/MyConsts_fs2AU, iLeftCfg, iRightCfg, NormError, NormReduced
!
!                IF ( NormError < NormErrorThreshold ) THEN
!                   WRITE(*,"(A,F9.3,A,F14.6)") "Time = ",ActualTime/MyConsts_fs2AU, &
!                         " : two gaussians have been merged, new norm = ", NormReduced
!                   WRITE(LogUnit,"(A,F9.3,A,F14.6)") "# Time = ",ActualTime/MyConsts_fs2AU, &
!                         " : two gaussians have been merged, new norm = ", NormReduced
!
!                   CALL DisposePsi( Psi )
!                   CALL SetPsi( Psi, PsiReduced%NrStates, PsiReduced%NrCfg, PsiReduced%GDim, PsiReduced%WFType )
!                   Psi%Bvector = PsiReduced%Bvector
!                   Psi%GaussPar = PsiReduced%GaussPar
!                   CALL UpdateWaveFunctInstance( Psi )
!                   ! CALL NormalizePsi( Psi )
!                   CALL DisposePsi( PsiReduced )
!                   EXIT outer
!
!                END IF
!
!                CALL DisposePsi( PsiReduced )
            END IF

         END DO
      END DO outer




     ! ******************

   END SUBROUTINE AddOrRemoveGaussianFunctions



!===========================================================================================================
!                               PRIVATE SUBROUTINES AND FUNCTIONS
!===========================================================================================================


!*******************************************************************************
!          RemoveOverlappingGaussians
!*******************************************************************************
!> Given an input Psi and the indices of two gaussian configurations,
!> computes the projection of the original wavefunction on a basis
!> of N-1 configurations, in which the given two configurations have
!> been replaced by a gaussian centered in the middle.
!>
!> @param Psi            WaveFunct object
!> @param nCfg1, nCfg2   indices of the two configurations to remove
!> @param PsiReduced     on output gives the projection of Psi on the reduced basis
!> @param NormReduced    on output gives the norm of PsiReduced
!*******************************************************************************
   SUBROUTINE RemoveOverlappingGaussians(Psi, nCfg1, nCfg2, PsiReduced, NormReduced)
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      INTEGER, INTENT(IN)              :: nCfg1, nCfg2
      TYPE(WaveFunct), INTENT(INOUT)   :: PsiReduced
      REAL, INTENT(OUT)                :: NormReduced

      INTEGER :: jCfg, iCfg, n
      COMPLEX, DIMENSION(Psi%NrCfg-1,Psi%NrCfg)       :: NewOldOverlap, Projector
      COMPLEX, DIMENSION(Psi%NrCfg-1,Psi%NrCfg-1)     :: NewNewInverse

      ! Allocate memory for the new reduced wavefunction
      CALL SetPsi( PsiReduced, Psi%NrStates, Psi%NrCfg-1, Psi%GDim, Psi%WFType )

      ! Define new set of gaussians removing the two gaussians with large overlap
      ! at the same time, define the matrix with the overlap between the old and the new basis
      n = 1
      DO jCfg = 1, Psi%NrCfg
         IF ( jCfg == nCfg1 .OR. jCfg == nCfg2 ) CYCLE
         ! Copy the gaussian parameters
         PsiReduced%GaussPar(:, (n-1)*Psi%GDim+1:n*Psi%GDim, 1) = Psi%GaussPar(:, (jCfg-1)*Psi%GDim+1:jCfg*Psi%GDim, 1)
         ! Copy the overlap with the old gaussians (which obviously remains the same)
         NewOldOverlap(n,:) = Psi%Overlap(jCfg,:,1,1)
         n = n+1
      END DO

      ! Add a new gaussian halfway between the two overlapping gaussians
      n = Psi%NrCfg-1
      PsiReduced%GaussPar(:, (n-1)*Psi%GDim+1:n*Psi%GDim, 1) = &
            0.5 * Psi%GaussPar(:, (nCfg1-1)*Psi%GDim+1:nCfg1*Psi%GDim, 1) + &
            0.5 * Psi%GaussPar(:, (nCfg2-1)*Psi%GDim+1:nCfg2*Psi%GDim, 1)
      CALL GauConf_Normalize( PsiReduced%GaussPar(:, (n-1)*Psi%GDim+1:n*Psi%GDim, 1) )

      ! and compute the overlap of the new gaussian with the old ones
      DO jCfg = 1, Psi%NrCfg
         NewOldOverlap(n,jCfg) = GauConf_Overlap( PsiReduced%GaussPar(:, (n-1)*Psi%GDim+1:n*Psi%GDim,       1),  &
                                                         Psi%GaussPar(:, (jCfg-1)*Psi%GDim+1:jCfg*Psi%GDim, 1) )
      END DO

      ! Now starting from the NewOldOverlap matrix, construct the overlap matrix of the new set of gaussians
      n = 1
      DO jCfg = 1, Psi%NrCfg
         IF ( jCfg == nCfg1 .OR. jCfg == nCfg2 ) CYCLE
         PsiReduced%Overlap(:,n,1,1) = NewOldOverlap(:,jCfg)
         n = n+1
      END DO
      ! only last column of NewNewOverlap is missing, and this can be recostructed from the conjugate values
      n = Psi%NrCfg-1
      DO jCfg = 1, Psi%NrCfg-2
         PsiReduced%Overlap(jCfg,n,1,1) = CONJG(PsiReduced%Overlap(n,jCfg,1,1))
      END DO
      PsiReduced%Overlap(n,n,1,1) = CMPLX(1.0,0.0)

      ! Now invert the matrix
      CALL MatrixInversionDo( PsiReduced%Overlap(:,:,1,1), NewNewInverse, 1 )
      ! Only upper diagonal has been computed, complete the inverse overlap
      DO jCfg = 1, Psi%NrCfg
         DO iCfg = jCfg+1, Psi%NrCfg
            NewNewInverse(iCfg,jCfg) = CONJG(NewNewInverse(jCfg,iCfg))
         END DO
      END DO
      ! And compute projector
      CALL TheOneWithMatrixMultiplication(Projector, NewNewInverse, NewOldOverlap)
      ! Compute new Bvector
      PsiReduced%BVector(:,1) = TheOneWithMatrixVectorProduct( Projector, Psi%BVector(:,1) )

      ! Now update all the internal elements of PsiReduced
!       CALL UpdateWaveFunctInstance( PsiReduced )
      ! and compute Norm
      NormReduced = WFNorm( PsiReduced )

   END SUBROUTINE RemoveOverlappingGaussians



!*******************************************************************************
!          RemoveSingleGaussian
!*******************************************************************************
!> Given an input Psi and the indices of a single gaussian configurations,
!> computes the projection of the original wavefunction on a basis
!> of N-1 configurations, in which the given configuration has been removed
!>
!> @param Psi            WaveFunct object
!> @param RemoveCfgNr    index of the two configuration to remove
!> @param PsiReduced     on output gives the projection of Psi on the reduced basis
!> @param NormReduced    on output gives the norm of PsiReduced
!*******************************************************************************
   SUBROUTINE RemoveSingleGaussian(Psi, RemoveCfgNr, PsiReduced, NormReduced)
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      INTEGER, INTENT(IN)              :: RemoveCfgNr
      TYPE(WaveFunct), INTENT(INOUT)   :: PsiReduced
      REAL, INTENT(OUT)                :: NormReduced

      INTEGER :: iCfg, jCfg, nRed, mRed
      COMPLEX, DIMENSION(Psi%NrCfg-1,Psi%NrCfg)       :: NewOldOverlap, Projector
      COMPLEX, DIMENSION(Psi%NrCfg-1,Psi%NrCfg-1)     :: NewNewInverse

      ! Allocate memory for the new reduced wavefunction
      CALL SetPsi( PsiReduced, Psi%NrStates, Psi%NrCfg-1, Psi%GDim, Psi%WFType )

      ! Define new set of gaussian removing the given gaussian configuration,
      ! at the same time, construct overlap matrix which are necessary (new-old config, new-new config)
      nRed = 1
      DO jCfg = 1, Psi%NrCfg
         ! skip the old configuration if it is the one to remove
         IF ( jCfg == RemoveCfgNr ) CYCLE

         ! Copy the gaussian parameters
         PsiReduced%GaussPar(:, (nRed-1)*Psi%GDim+1:nRed*Psi%GDim, 1) = Psi%GaussPar(:, (jCfg-1)*Psi%GDim+1:jCfg*Psi%GDim, 1)
         ! Copy the overlap of the new gaussian with the old gaussians (which obviously remains the same, except for one row which is removed)
         NewOldOverlap(nRed,:) = Psi%Overlap(jCfg,:,1,1)

         mRed = 1
         DO iCfg = 1, Psi%NrCfg
            ! skip the old configuration if it is the one to remove
            IF ( iCfg == RemoveCfgNr ) CYCLE

            ! Copy the overlap of the new configurations extracting from the old overlap
            PsiReduced%Overlap(nRed, mRed, 1,1) = Psi%Overlap(jCfg,iCfg,1,1)

            mRed = mRed+1
         END DO

         nRed = nRed+1
      END DO

      ! Now invert the overlap matrix for the new configurations
      CALL MatrixInversionDo( PsiReduced%Overlap(:, :, 1,1), NewNewInverse, 1 )
      ! Only upper diagonal has been computed, complete the inverse overlap
      DO jCfg = 1, Psi%NrCfg-1
         DO iCfg = jCfg+1, Psi%NrCfg-1
            NewNewInverse(iCfg,jCfg) = CONJG(NewNewInverse(jCfg,iCfg))
         END DO
      END DO
      ! And compute projector
      CALL TheOneWithMatrixMultiplication(Projector, NewNewInverse, NewOldOverlap)
      ! Compute new Bvector
      PsiReduced%BVector(:,1) = TheOneWithMatrixVectorProduct( Projector, Psi%BVector(:,1) )

      ! Now update all the internal elements of PsiReduced
!       CALL UpdateWaveFunctInstance( PsiReduced )
      ! and compute Norm
      NormReduced = WFNorm( PsiReduced )

   END SUBROUTINE RemoveSingleGaussian



   SUBROUTINE EmptyGaussianErrorEstimate( Psi, EmptyConfig, Hamiltonian, ErrorEstimate, Derivatives )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)                          :: Psi
      COMPLEX, DIMENSION(3,Psi%GDim), INTENT(IN)            :: EmptyConfig
      TYPE(OperatorData), INTENT(IN)                          :: Hamiltonian
      REAL, INTENT(OUT)                                       :: ErrorEstimate
      COMPLEX, DIMENSION(Psi%GDim), INTENT(OUT), OPTIONAL   :: Derivatives

      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrGauSets)              ::  OverlapInv
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrStates,Psi%NrStates)  ::  HMatrix

      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrGauSets)                        ::  UnoccOccOverlap
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrGauSets)                        ::  UnoccOccHamiltonian
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrGauSets,Psi%GDim)             ::  UnoccOccFirstMom
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrGauSets,Psi%GDim)             ::  UnoccOccHamilFirstMom
      COMPLEX, DIMENSION(Psi%GDim)                                     ::  UnoccUnoccFirstMom
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg)                            ::  InvOverlapTimesHamilt

      REAL    :: InvCondNr, ProjNorm
      COMPLEX :: SumOv, g_1minP_H_Psi, gPHPsi
      COMPLEX, DIMENSION(Psi%GDim) :: ProjFirstMom, g_x_1minP_H_Psi
      INTEGER :: i, j, k

      ! In the notes below, we will indicate:
      !  | g_alpha >  alpha = 1,2, ...    configurations of the "starting" wavefunction Psi
      !  | Gtilde >                       additional configuration included, we want to estimate the effect of the inclusion of |Gtilde> in the wavefunction
      !     x_k  k = 1,2...               one of the degrees of freedom of the space in which the configurations are defined
      !   P                               the projector on the non-orthogonal configurations  sum_alpha,beta | g_alpha > [ S^{-1} ]_alpha,beta < g_beta |

      ! Compute the matrix inverse of the overlap between the configurations | g_alpha >

      CALL MatrixInversionDo( Psi%Overlap(:,:,1,1), OverlapInv(:,:,1), 1, InvCondNr )
      ! Only upper diagonal has been computed, complete the inverse overlap
      DO j = 1, Psi%NrCfg
         DO i = j+1, Psi%NrCfg
            OverlapInv(i,j,1) = CONJG(OverlapInv(j,i,1))
         END DO
      END DO

      ! Calculate all the Hamiltonian and Overlap matrix elements over the Psi configurations + the additional one
      ! HMatrix = < g_alpha | H | g_beta >                    UnoccOccHamiltonian = < Gtilde | H | g_alpha >
      ! UnoccOccOverlap = < Gtilde | g_alpha >                UnoccOccFirstMom = < Gtilde | x_k | g_alpha >
      ! UnoccOccHamilFirstMom = < Gtilde | x_k H | g_beta >   UnoccUnoccFirstMom = < Gtilde | x_k | Gtilde >
      ! InvOverlapTimesHamilt = SUM_beta [ S^{-1} ]_alpha,beta  < g_beta | H | g_gamma >

      HMatrix = GetOperatorMatrix( Hamiltonian, Psi )
      CALL TheOneWithMatrixMultiplication(InvOverlapTimesHamilt, OverlapInv(:,:,1), HMatrix(:,:,1,1), "N","N")

      DO j = 1, Psi%NrCfg
         UnoccOccHamiltonian(j,1) = OperatorMatrixElement( Hamiltonian, EmptyConfig(:,:), &
                      Psi%GaussPar(:, (j-1)*Psi%GDim+1:j*Psi%GDim, 1), 1, 1 )
         UnoccOccOverlap(j,1) = GauConf_Overlap( EmptyConfig(:,:), Psi%GaussPar(:, (j-1)*Psi%GDim+1:j*Psi%GDim, 1) )
         IF ( PRESENT(Derivatives)) THEN
            DO k = 1, Psi%GDim
               UnoccOccFirstMom(j,1,k) = UnoccOccOverlap(j,1) * ( &
                     GauPrim_Moment( EmptyConfig(:,k), Psi%GaussPar(:, (j-1)*Psi%GDim+k, 1), 1 )  )
               UnoccOccHamilFirstMom(j,1,k) = OperatorMatrixElement( Hamiltonian, EmptyConfig(:,:), &
                  Psi%GaussPar(:, (j-1)*Psi%GDim+1:j*Psi%GDim, 1), 1, 1, AugDim = k, AugF = 1)
            END DO
         END IF
      END DO
      IF ( PRESENT(Derivatives)) THEN
         DO k = 1, Psi%GDim
            UnoccUnoccFirstMom(k) = GauPrim_Moment( EmptyConfig(:,k), EmptyConfig(:,k), 1 )
         END DO
      END IF

      ! COMPUTE  ProjNorm = < Gtilde | 1-P | Gtilde > and ProjFirstMom = < Gtilde | x_k (1-P) | Gtilde >

      SumOv = CMPLX(0.0,0.0)
      DO j = 1, Psi%NrCfg
         DO i = 1, Psi%NrCfg
            SumOv = SumOv + UnoccOccOverlap(i,1) * OverlapInv(i,j,1) * CONJG(UnoccOccOverlap(j,1))
         END DO
      END DO
      ProjNorm = 1 - SumOv

      IF ( PRESENT(Derivatives)) THEN
         ProjFirstMom = CMPLX(0.0,0.0)
         DO k = 1, Psi%GDim
            DO j = 1, Psi%NrCfg
               DO i = 1, Psi%NrCfg
                  ProjFirstMom(k) = ProjFirstMom(k) + UnoccOccFirstMom(i,1,k) * OverlapInv(i,j,1) * CONJG(UnoccOccOverlap(j,1))
               END DO
            END DO
            ProjFirstMom(k) = UnoccUnoccFirstMom(k) - ProjFirstMom(k)
         END DO
      END IF

      g_1minP_H_Psi = CMPLX(0.0,0.0)
      DO j = 1, Psi%NrCfg
         gPHPsi = CMPLX(0.0,0.0)
         DO i = 1, Psi%NrCfg
            gPHPsi = gPHPsi + UnoccOccOverlap(i,1) * InvOverlapTimesHamilt(i,j)
         END DO
         g_1minP_H_Psi = g_1minP_H_Psi + (UnoccOccHamiltonian(j,1) - gPHPsi) * Psi%Bvector(j, 1)
      END DO

      IF ( PRESENT(Derivatives)) THEN
         DO k = 1, Psi%GDim
            g_x_1minP_H_Psi(k) = CMPLX(0.0,0.0)
            DO j = 1, Psi%NrCfg
               gPHPsi = CMPLX(0.0,0.0)
               DO i = 1, Psi%NrCfg
                  gPHPsi = gPHPsi + UnoccOccFirstMom(i,1,k) * InvOverlapTimesHamilt(i,j)
               END DO
               g_x_1minP_H_Psi(k) = g_x_1minP_H_Psi(k) + (UnoccOccHamilFirstMom(j,1,k) - gPHPsi) * Psi%Bvector(j, 1)
            END DO
            g_x_1minP_H_Psi(k) = g_x_1minP_H_Psi(k)
         END DO
      END IF

      ErrorEstimate = ABS(g_1minP_H_Psi)**2 / ProjNorm
      IF ( PRESENT(Derivatives)) THEN
         DO k = 1, Psi%GDim
            Derivatives(k) = 2.0 * ( CONJG(g_1minP_H_Psi) * g_x_1minP_H_Psi(k) - ErrorEstimate * ProjFirstMom(k) ) / ProjNorm
         END DO
      END IF

   END SUBROUTINE EmptyGaussianErrorEstimate


   SUBROUTINE OptimizeEmptyGaussian( Psi, SglConfig, Hamiltonian )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)                          :: Psi
      COMPLEX, DIMENSION(3,Psi%GDim), INTENT(INOUT)            :: SglConfig
      TYPE(OperatorData), INTENT(IN)                          :: Hamiltonian

      INTEGER  :: iter
      COMPLEX, DIMENSION(Psi%GDim) :: Derivatives, OldXi, OldDerivs
      REAL :: ErrorEstimate, Lambda, OldErrorEstimate

      ! Compute new derivatives
      CALL EmptyGaussianErrorEstimate( Psi, SglConfig, Hamiltonian, ErrorEstimate, Derivatives )
      WRITE(901,"(I5,100F30.16)") 0, ErrorEstimate, MAXVAL(ABS(Derivatives(:))), 0.0

      ! Store the previous coordinates and the previous derivatives
      OldXi = SglConfig(2,:)
      OldDerivs = Derivatives(:)
      OldErrorEstimate = ErrorEstimate

      ! First step with lambda = 100.0
      Lambda = 100.0

      ! ITERATIONS for the optimization steepest descent
      DO iter = 1,100

         ! gradually reduce Lambda to get an actual increase of ErrorEstimate
         DO
            ! Move along gradient
            SglConfig(2,:) = OldXi + Lambda * OldDerivs(:)
            CALL GauConf_Normalize(SglConfig)
            ! compute new ErrorEstimate and derivatives
            CALL EmptyGaussianErrorEstimate( Psi, SglConfig, Hamiltonian, ErrorEstimate )

            IF ( ErrorEstimate > OldErrorEstimate ) THEN
               EXIT
            ELSE
               Lambda = Lambda * 0.7
            END IF
            IF ( Lambda < 1.E-8 ) EXIT
         END DO

         ! compute new derivatives
         CALL EmptyGaussianErrorEstimate( Psi, SglConfig, Hamiltonian, ErrorEstimate, Derivatives )
         WRITE(901,"(I5,100F30.16)") iter, ErrorEstimate, MAXVAL(ABS(Derivatives(:))), Lambda

         ! check convergence criteria
         IF ( MAXVAL(ABS(Derivatives(:))) < 1.E-6 .OR. Lambda < 1.E-8 ) EXIT

         ! We now have new and old derivatives ( Derivatives, OldDerivs ) and new and old coordianates ( SglConfig, OldXi ), compute lambda
         Lambda = ComputeLambda( SglConfig(2,:), OldXi, Derivatives, OldDerivs )
         IF ( Lambda < 0.0 ) Lambda = 100.0

         ! Store the previous coordinates and the previous derivatives
         OldXi = SglConfig(2,:)
         OldDerivs = Derivatives(:)

         ! update the value of the configurations
         SglConfig(2,:) = SglConfig(2,:) + Lambda * Derivatives(:)
         CALL GauConf_Normalize(SglConfig)

      END DO

      WRITE(901,*) " "

   END SUBROUTINE OptimizeEmptyGaussian


   REAL FUNCTION ComputeLambda( NewXi, OldXi, NewDerivs, OldDerivs  )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:), INTENT(IN) ::  NewXi, OldXi, NewDerivs, OldDerivs
      REAL :: GradTimesX, GradTimesGrad
      COMPLEX :: DeltaX, DeltaG
      INTEGER :: j

      GradTimesX = 0.0
      GradTimesGrad = 0.0

      DO j = 1, SIZE(NewXi)
         DeltaX = NewXi(j) - OldXi(j)
         DeltaG = NewDerivs(j) - OldDerivs(j)
         GradTimesX = GradTimesX + REAL(DeltaG)*REAL(DeltaX) + AIMAG(DeltaG)*AIMAG(DeltaX)
         GradTimesGrad = GradTimesGrad + REAL(DeltaG)*REAL(DeltaG) + AIMAG(DeltaG)*AIMAG(DeltaG)
      END DO
      ComputeLambda = - GradTimesX / GradTimesGrad

   END FUNCTION ComputeLambda



   SUBROUTINE AddEmptyGaussianToPsi( Psi, SglConfig )
      TYPE(WaveFunct), INTENT(INOUT)                   :: Psi
      COMPLEX, DIMENSION(3,Psi%GDim), INTENT(IN)     :: SglConfig
      TYPE(WaveFunct) :: PsiTmp

      ! make a temporary copy of Psi
      PsiTmp = Psi

      ! reallocate memory of Psi
      CALL DisposePsi( Psi )
      CALL SetPsi( Psi, PsiTmp%NrStates, PsiTmp%NrCfg + 1, PsiTmp%GDim, PsiTmp%WFType )

      ! copy to psi the previous gaussian basis
      Psi%GaussPar(:,1:PsiTmp%NrPrimGau,1) = PsiTmp%GaussPar(:,:,1)
      ! and add the one given from input
      Psi%GaussPar(:,PsiTmp%NrPrimGau+1:,1) = SglConfig

      ! copy to psi the previous coefficients and add an unoccupied one
      Psi%Bvector(1:PsiTmp%NrCfg,1) = PsiTmp%Bvector(:,1)
      Psi%Bvector(PsiTmp%NrCfg+1,1) = CMPLX(0.0,0.0)

      ! now recompute the rest of the internal matrices and stuff
      CALL UpdateWaveFunctInstance( Psi )

      ! deallocate memory for the temporary wavefunction
      CALL DisposePsi( PsiTmp )

   END SUBROUTINE AddEmptyGaussianToPsi

END MODULE AdaptGBasis


!    SUBROUTINE ReduceAChangeVectors( RemoveCfgNr )
!       IMPLICIT NONE
!       INTEGER :: RemoveCfgNr
!       REAL, DIMENSION(SIZE(AChangeInitTime))     :: CopyAChangeInitTime
!       INTEGER, DIMENSION(SIZE(AChangeDirection)) :: CopyAChangeDirection
!       INTEGER :: i, n, RedSize
!
!       ! Save copy of the two arrays
!       CopyAChangeInitTime = AChangeInitTime
!       CopyAChangeDirection = AChangeDirection
!
!       RedSize = SIZE(CopyAChangeInitTime)-1
!       ! Reallocate the vectors with one element less
!       DEALLOCATE( AChangeInitTime, AChangeDirection )
!       ALLOCATE( AChangeInitTime(RedSize), AChangeDirection(RedSize) )
!
!       ! Cycle over the dimension of the copies and save right elements to the reallocated arrays
!       n = 1
!       DO i = 1, SIZE(CopyAChangeInitTime)
!          IF ( i == RemoveCfgNr ) CYCLE
!          AChangeInitTime(n) = CopyAChangeInitTime(i)
!          AChangeDirection(n) = CopyAChangeDirection(i)
!          n = n+1
!       END DO
!
!       PRINT*, AChangeDirection
!
!
!    END SUBROUTINE ReduceAChangeVectors


!    LOGICAL FUNCTION ActivateAImpulsiveChange( Psi, OverlapThreshold, ActualTime ) RESULT(Check)
!       TYPE(WaveFunct), INTENT(INOUT)   :: Psi
!       REAL, INTENT(IN)                 :: OverlapThreshold
!       REAL, INTENT(IN)                 :: ActualTime
!
!       INTEGER :: iLeftCfg, iRightCfg
!       REAL, PARAMETER :: DeltaT = 30.0*MyConsts_fs2AU
!
!       ! Allocate memory when needed: AChangeDirection defines whether shrinking is active or not
!       ! when gaussias is shrinking, AChangeInitTime define the initial time
!       IF ( .NOT. ALLOCATED(AChangeDirection) ) THEN
!          ALLOCATE( AChangeInitTime(Psi%NrCfg) )
!          AChangeInitTime = 0.0
!          ALLOCATE( AChangeDirection(Psi%NrCfg) )
!          AChangeDirection = 0
!       END IF
!
!       Check = .FALSE.
!
!       DO iRightCfg = 1, Psi%NrCfg
!          DO iLeftCfg = iRightCfg+1, Psi%NrCfg
!
!             ! check if overlap is greater then the threshold, in case proceed with removing one gaussian from the basis set
!             IF ( ABS(Psi%Overlap(iLeftCfg,iRightCfg,1,1)) > OverlapThreshold ) THEN
!
!                IF ( AChangeDirection(iLeftCfg) /= 0 .OR. AChangeDirection(iRightCfg) /= 0 ) CYCLE
!
!                AChangeDirection(iLeftCfg) = +1
!                AChangeDirection(iRightCfg) = 0
!
!                AChangeInitTime(iLeftCfg) = ActualTime
!                AChangeInitTime(iRightCfg) = 0.0
!
!                Check = .TRUE.
!
!                RETURN
!
!             END IF
!
!          END DO
!       END DO
!
!    END FUNCTION ActivateAImpulsiveChange
