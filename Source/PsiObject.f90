!***************************************************************************************
!*                              MODULE PsiObject
!***************************************************************************************
!
!>  \brief     Definition of Psi in vMCG format
!>  \details   This module contains an object to define the Psi in vMCG format \n
!>             and all the necessary methods to operate on the object.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             12 July 2017
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          Define begin and end of gaussian configurations in the combined
!>                 index that label the primitive functions in the WaveFunct object
!>                 this map could be useful because in most case the not primitive
!>                 gaussians but configurations needs to be accessed in the GaussPar
!>                 array.
!
!***************************************************************************************
MODULE PsiObject
#include "preprocessoptions.cpp"
   USE GauConf
   USE OperatorDefine
   USE MatrixInversion
   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WaveFunct                                ! New data type definying the Psi object instance
   PUBLIC :: Derivative                               ! New data type to store derivative of finite difference of the wavefunction

   PUBLIC :: SetPsi, DisposePsi                       ! constructor and destructor of a WaveFunct instance
   PUBLIC :: ReadPsiFromFile, WritePsiToFile          ! read and write wavefunction from/to single file
   PUBLIC :: NormalizePsi                             ! normalize wavefunction
   PUBLIC :: UpdateWaveFunctInstance                  ! after new coeff.s and param.s are defined, update WaveFunct instance

   PUBLIC :: EquationsOfMotionSetup, FixGaussianConfigurationsWithSmallPopulations
   PUBLIC :: ComputeDerivative, DisposeDerivative     ! constructor and desctructur of a Derivative instance

   PUBLIC :: WFDifference                             ! Compute the difference between two wfs accoring to some given metrics
   PUBLIC :: WFExpectation, WFNorm
   PUBLIC :: WFOverlap, GetOverlapMax
   PUBLIC :: WFGauCenters, WFGauWidth                 ! Get centers and width of the gaussians along a certain degree of freedom and electronic set
   PUBLIC :: WFGauPopulations
   PUBLIC :: GetOperatorMatrix

   PUBLIC :: PsiPlusDeltaPsi
!    PUBLIC :: ASSIGNMENT(=), OPERATOR(+), OPERATOR(*)  ! Overloaded operators
!
!    INTERFACE ASSIGNMENT (=)           ! Assigment between two WaveFunct or two Derivatives   Psi1 = Psi2
!       MODULE PROCEDURE CopyPsi, CopyDerivative
!    END INTERFACE

!    INTERFACE OPERATOR (+)             ! Sum of two derivatives                   DeltaPsi1 + DeltaPsi2
!       MODULE PROCEDURE SumDerivatives
!    END INTERFACE
!
!    INTERFACE OPERATOR (*)             ! Multiplication of derivative by scalar   alpha * DeltaPsi1
!       MODULE PROCEDURE MultiplyDerivativeByScalar
!    END INTERFACE

   !===========================================================================================================

   ! Format of the b-vector file
   CHARACTER(100), PARAMETER :: BVectorFormat = '(I4,1X,E23.16,1X,E23.16)'
   ! Format of the gaussian parameters file
   CHARACTER(100), PARAMETER :: GaussFormat = '(I4,3(1X,E23.16,1X,E23.16),1X,E23.16)'


   ! Definitions of the various type of wavefunction implemented
   INTEGER, PUBLIC, PARAMETER :: vMCG_SINGLESET   = 1,  & ! vMCG ansatz, single-set wavefunction
                                 vMCG_MULTISET    = 2     ! vMCG ansatz, multi-set wavefunction

   ! Definitions of the various type of wavefunction distances implemented
   INTEGER, PUBLIC, PARAMETER :: WFDISTANCE_H_SCALARPROD = 1, &  ! usual distance, induced by wf scalar product in Hilbert space
                                 WFDISTANCE_L1_NORM      = 2, &  ! l1-norm sum|D|
                                 WFDISTANCE_L2_NORM      = 3, &  ! l2-norm SQRT(sum|D|^2)
                                 WFDISTANCE_LINF_NORM    = 4     ! linf-norm sup(|D|)
                                    ! with D difference between WFs seen as a sequence of numbers, for both coeffs and parameters

   ! Definitions of the type of equations of motion implemented
   INTEGER, PUBLIC, PARAMETER :: EOM_STANDARD            = 1, &  ! standard EOMs: (1-P) in the parameter equations
                                 EOM_INVERSE             = 2, &  ! alternative set of EOMs: (1-P) in the coefficient equations
                                 EOM_SINGLEINV           = 3     ! alternative set of EOMs: coeffs and pars are obtained by single inversion

   !===========================================================================================================

   !> Data type to store the information on the wavefunction at given time.
   TYPE WaveFunct
                     ! a "public" definition of the components of WaveFunct is a possible source of error
                     ! (since the access to the WaveFunct components is now much less controllable)
                     ! however, in current FORTRAN it is not yet possible to extend module with other modules
                     ! and if WaveFunct components are not public, this module ends up having ALL the lines of
                     ! code that manipulate the wavefunction, too much stuff to keep the code readable
                     ! Other modules which access the WaveFunct components: EmptyOptimize

      ! General dimensions of the wavefunction
      INTEGER :: WFType                            !< Variable to define which kind of wavefunction is used
      INTEGER :: NrStates                          !< number of electronic states
      INTEGER :: NrCfg                             !< number of vMCG configurations
      INTEGER :: GDim                              !< number of dimensions of the gaussians
      INTEGER :: NrPrimGau                         !< number of primitive gaussians (GauDim*NrCfg)
      INTEGER :: NrGauSets                         !< number of sets of gaussians (1 for single-set, NrStates for multi-set)

      ! Storage area of the wavefunction
      COMPLEX, ALLOCATABLE, DIMENSION(:,:)   :: BVector   !< linear coefficients of the vMCG expansion
      COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: GaussPar  !< gaussian parameters of the vMCG expansion

      ! The following arrays store the overlap and the moments of the gaussian functions and needs to be updated any time
      ! the value of the gaussian parameters are changed. Overlap store the integral <gi|gj> and needs to be fully computed
      ! any time the gaussian are displaced. qCoord store the actual values of the mean q of the primitive gaussians.
      COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)   :: Overlap        !< store the actual overlap between the gaussian functions
      REAL, ALLOCATABLE, DIMENSION(:,:)          :: qCoord         !< store the actual value of the mean q of the prim gaussian

      LOGICAL :: WaveFunctIsSetup = .FALSE.             !< Logical flag to check whether the datatype is ready or not
   END TYPE WaveFunct

   !===========================================================================================================

   !> Data type to store the variation or the derivative of the wavefunction at given time.
   TYPE Derivative
      PRIVATE
      INTEGER :: WFType                                  ! Type of the wavefunction (same definitions as the PsiObject type

      ! Storage area for each of the components of the wavefunction
      COMPLEX, ALLOCATABLE, DIMENSION(:,:)   :: BVector   !< linear coefficients of the vMCG expansion
      COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: GaussPar  !< gaussian parameters of the vMCG expansion

      LOGICAL :: DerivativeIsSetup = .FALSE.             !< Logical flag to check whether the datatype is ready or not

!       CONTAINS
!          FINAL :: DisposeDerivative
   END TYPE Derivative

   !===========================================================================================================

!    ! The following arrays are used by the module AdaptGBasis and the subroutine
!    ! DerivativeAImpulsiveChange (here in this module) to shrink gaussians during propagation
!    PUBLIC :: AChangeInitTime, AChangeDirection
!    REAL, ALLOCATABLE, DIMENSION(:), SAVE    :: AChangeInitTime
!    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: AChangeDirection

   !===========================================================================================================

   ! ==============================================================
   ! FINE SETUP FOR THE DEFINITION OF THE EQUATIONS OF MOTION
   ! ==============================================================

   ! Type of equation of motions used for the propagation
   INTEGER :: EquationType          = EOM_STANDARD
   ! Flag to use the form of the equations with the simplification of the coefficients from the C matrix
   LOGICAL :: SimplifyCoefficients  = .FALSE.
   !> threshold that defines the regularization of the inversion of the coefficients vector
   REAL    :: SimplifyThreshold

   ! Freeze gaussians with low populations
   !> logical variable to activate the freezing of the low-populated gaussians
   LOGICAL :: FreezeLowPopGaussians  = .FALSE.
   !> thresholds that define which gaussian to freeze
   REAL    ::  PopMinThreshold, PopMaxThreshold

   ! Freeze gaussians with low abs(coefficients)
   !> logical variable to activate the freezing of the low-coefficient gaussians
   LOGICAL :: FreezeLowCoeffGaussians  = .FALSE.
   !> thresholds that define which gaussian to freeze
   REAL    ::  CoeffMinThreshold

   ! Following arrays define which gaussian configurations are frozen during the dynamics
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC  ::  CMask
   INTEGER, SAVE, PUBLIC :: NInv

   CONTAINS

!===========================================================================================================
!                               PUBLIC SUBROUTINES AND FUNCTIONS
!===========================================================================================================


   INTEGER FUNCTION Prim2Cfg( PrimitiveId, GauDim )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: PrimitiveId, GauDim
      Prim2Cfg = (PrimitiveId-1)/GauDim + 1
   END FUNCTION Prim2Cfg

   INTEGER FUNCTION Prim2Dim( PrimitiveId, GauDim )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: PrimitiveId, GauDim
      Prim2Dim = MOD(PrimitiveId-1, GauDim) + 1
   END FUNCTION Prim2Dim

   INTEGER FUNCTION Cfg2Prim( Cfg, iDimen, GauDim )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Cfg, iDimen, GauDim
      Cfg2Prim = (Cfg-1)*GauDim + iDimen
   END FUNCTION Cfg2Prim

   FUNCTION CfgSec( Cfg, GauDim )  RESULT( Indices )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Cfg, GauDim
      INTEGER, DIMENSION(GauDim) :: Indices
      INTEGER :: i
      Indices = (/ ( (Cfg-1)*GauDim+i, i = 1, GauDim  ) /)
   END FUNCTION CfgSec

   INTEGER FUNCTION CfgBeg( Cfg, GauDim )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Cfg, GauDim
      CfgBeg = (Cfg-1)*GauDim + 1
   END FUNCTION CfgBeg

   INTEGER FUNCTION CfgEnd( Cfg, GauDim )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Cfg, GauDim
      CfgEnd = Cfg*GauDim
   END FUNCTION CfgEnd


!*******************************************************************************
!          SetPsi
!*******************************************************************************
!> WaveFunct constructor: set the parameters defining the wavefunction
!> Ansatz and allocate memory to store the wavefunction coefficients and
!> parameters.
!>
!> @param Psi            WaveFunct istance to costruct
!> @param InpNrStates    Nr of electronic states of the wavefunction
!> @param InpNrCfg       Nr of gaussian configurations
!> @param InpGauDim      Nr of degrees of freedom of the gaussians
!> @param InpType        Type of wavefunction ansatz
!*******************************************************************************
   SUBROUTINE SetPsi( Psi, InpNrStates, InpNrCfg, InpGauDim, InpType )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT) :: Psi
      INTEGER, INTENT(IN) :: InpNrStates, InpNrCfg, InpGauDim, InpType

      ! Nr states different from 1 is not yet implemented
      CALL ERROR( InpNrStates /= 1, " PsiObject: wavefunction with nstates > 1 is not implemented yet" )

      ! Store the values of the wavefunction type and dimensional parameters
      Psi%WFType = InpType
      Psi%NrStates = InpNrStates
      Psi%NrCfg = InpNrCfg
      Psi%GDim = InpGauDim
      SELECT CASE( Psi%WFType )
         CASE( vMCG_SINGLESET )
            Psi%NrGauSets = 1
         CASE( vMCG_MULTISET )
            Psi%NrGauSets = InpNrStates
      END SELECT
      Psi%NrPrimGau = Psi%NrCfg*Psi%GDim

      ! In case memory is already allocated, deallocate (printing warning to be safe)
      IF ( Psi%WaveFunctIsSetup ) THEN
         CALL ShowWarning( " PsiObject.SetPsi: disposing memory that was previosly allocated for a WaveFunct object" )
         CALL DisposePsi( Psi )
      END IF

      ! Allocate memory according to the defined dimensions
      ALLOCATE( Psi%Bvector(Psi%NrCfg, Psi%NrStates) )
      ALLOCATE( Psi%GaussPar(3, Psi%NrPrimGau, Psi%NrGauSets) )
      ALLOCATE( Psi%Overlap(Psi%NrCfg,Psi%NrCfg,Psi%NrGauSets,Psi%NrGauSets) )
      ALLOCATE( Psi%qCoord(Psi%NrPrimGau, Psi%NrGauSets) )

      ! Initialize allocated memory
      Psi%Bvector  = CMPLX(0.0,0.0)
      Psi%GaussPar = CMPLX(0.0,0.0)
      Psi%Overlap  = CMPLX(0.0,0.0)
      Psi%qCoord   = 0.0

      ! Wavefunction is now setup
      Psi%WaveFunctIsSetup = .TRUE.

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "SetPsi: an instance of WaveFunct has been constructed with: "
      SELECT CASE( Psi%WFType )
         CASE( vMCG_SINGLESET )
            __WRITE_TO_LOG " Single-Set type wavefunction "
         CASE( vMCG_MULTISET )
            __WRITE_TO_LOG " Multi-Set type wavefunction "
      END SELECT
      __WRITE_TO_LOG "* NrStates = ", NumberToString(Psi%NrStates)
      __WRITE_TO_LOG "* NrCfg = ", NumberToString(Psi%NrCfg)
      __WRITE_TO_LOG "* GauDim = ", NumberToString(Psi%GDim)
      __WRITE_TO_LOG "* Bvector dimension ", NumberToString(Psi%NrCfg*Psi%NrStates)
      __WRITE_TO_LOG "* GaussPar dimension ", NumberToString(3*Psi%NrPrimGau*Psi%NrStates)
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif
   END SUBROUTINE SetPsi


!*******************************************************************************
!          DisposePsi
!*******************************************************************************
!> WaveFunct destructor: deallocate memory of an instance of WaveFunct
!>
!> @param Psi            WaveFunct istance to destruct
!*******************************************************************************
   SUBROUTINE DisposePsi( Psi )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT) :: Psi

      ! In case memory is already allocated, deallocate
      IF ( ALLOCATED(Psi%Bvector) ) DEALLOCATE( Psi%Bvector )
      IF ( ALLOCATED(Psi%GaussPar) ) DEALLOCATE( Psi%GaussPar )
      IF ( ALLOCATED(Psi%Overlap) ) DEALLOCATE( Psi%Overlap )
      IF ( ALLOCATED(Psi%qCoord) ) DEALLOCATE( Psi%qCoord )

      Psi%WaveFunctIsSetup = .FALSE.

! #if defined(LOG_FILE)
!       __OPEN_LOG_FILE;
!       __WRITE_TO_LOG "DisposePsi: an instance of WaveFunct has been desctructed"
!       __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
! #endif

   END SUBROUTINE DisposePsi


!*******************************************************************************
!          ReadPsiFromFile
!*******************************************************************************
!> Set the coefficients and parameters of a WaveFunct object with the numbers
!> written in an input file formatted as usual.
!> For vMCG, two files need to be provided, 1) the coefficients in this format:
!> %%% nr of line, real part, imaginary part
!> the gaussians parameters in this format:
!> %%% nr of line, a_real, a_imag, xi_real, xi_imag, eta_real, eta_imag
!>
!> @param Psi            WaveFunct object
!> @param FileNameList   List of strings, with the necessary filenames to setup Psi
!*******************************************************************************
   SUBROUTINE ReadPsiFromFile( Psi, FileNameList )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)          :: Psi
      CHARACTER(*), DIMENSION(2), INTENT(IN)  :: FileNameList

      INTEGER :: InpUnit, nGau, nAconfig, IOStatus, nLine
      REAL    :: RealPart, ImagPart, AReal, AImg, XiReal, XiImg, EtaReal, EtaImg
      LOGICAL :: FileExists

      ! Check that the WaveFunct object has been constructed already
      CALL ERROR( .NOT. Psi%WaveFunctIsSetup, " ReadPsiFromFile: WaveFunct instance is not setup ", ERR_OBJ_MISUSE )

      ! Check the existence of the file
      INQUIRE(FILE=FileNameList(1), EXIST=FileExists)
      CALL ERROR( .NOT. FileExists, " ReadPsiFromFile: input file " //trim(FileNameList(1))// " does not exist", ERR_FILE_MISSING )

      ! Open input
      InpUnit = LookForFreeUnit(); OPEN(UNIT=InpUnit, FILE=TRIM(ADJUSTL(FileNameList(1))))

      ! Loop over the configurations and set the appropriate value of the Bvector
      DO nAconfig = 1,Psi%NrCfg
         READ(InpUnit,*,IOSTAT=IOStatus) nLine, RealPart, ImagPart
         CALL ERROR( IOStatus /= 0, "ReadPsiFromFile: Error reading file "//TRIM(FileNameList(1)), ERR_INP_READ )
         CALL ERROR( nLine /= nAconfig, "ReadPsiFromFile: wrong file "//TRIM(FileNameList(1))//" at line "//NumberToString(nLine),&
                     ERR_INP_READ )
         ! Store the linear coefficients
         Psi%Bvector(nAconfig, 1) = CMPLX(RealPart, ImagPart)
      END DO

      ! Close file
      CLOSE(InpUnit)

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "ReadPsiFromFile: the Bvector of a WaveFunct instance has been set reading from file ", FileNameList(1)
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif

      ! Check the existence of the file
      INQUIRE(FILE=FileNameList(2), EXIST=FileExists)
      CALL ERROR( .NOT. FileExists, " ReadPsiFromFile: input file "//trim(FileNameList(2))//" does not exist",ERR_FILE_MISSING)

      ! Open input
      OPEN(UNIT=InpUnit, FILE=TRIM(ADJUSTL(FileNameList(2))))

      ! Loop over the primitive gaussians and set the appropriate value of the gaussian parameters array
      DO nGau = 1,Psi%NrPrimGau
            READ(InpUnit,*,IOSTAT=IOStatus) nLine, AReal, AImg, XiReal, XiImg, EtaReal, EtaImg
            CALL ERROR( IOStatus /= 0, "ReadPsiFromFile: Error reading file "//TRIM(FileNameList(2)), ERR_INP_READ )
            CALL ERROR( nLine /= nGau, &
                "ReadPsiFromFile: wrong file "//TRIM(FileNameList(2))//" at line "//NumberToString(nLine),ERR_INP_READ )
            ! Store the gaussian parameters
            Psi%GaussPar(1, nGau, 1) = CMPLX(AReal, AImg)
            Psi%GaussPar(2, nGau, 1) = CMPLX(XiReal, XiImg)
            Psi%GaussPar(3, nGau, 1) = CMPLX(EtaReal, EtaImg)
      END DO

      ! Close file
      CLOSE(InpUnit)

      ! Update the components of the WaveFunct instance
      CALL  UpdateWaveFunctInstance( Psi  )

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "ReadPsiFromFile: the gaussian param.s of a WaveFunct instance have been set from file ", FileNameList(2)
      __WRITE_TO_LOG "ReadPsiFromFile: normalization of the primitive gaussians have been imposed "
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif

   END SUBROUTINE ReadPsiFromFile


!*******************************************************************************
!          WritePsiToFile
!*******************************************************************************
!> Write the coefficients vector of a WaveFunct object to an output file
!> In this case a single Bvector file is written, the format is
!> %%% nr of line, real part, imaginary part     -      1 line for each nconfig
!> Write the gaussian parameters array of a WaveFunct object to an output file
!> In this case the array of a single wavefunction is written, the format is
!> %%% nr of line, real part, imaginary part for A, Xi, Eta - 1 line for each nconfig,ndim
!>
!> @param Psi            WaveFunct object
!> @param FileNameList   List of strings, with the necessary filenames to setup Psi
!*******************************************************************************
   SUBROUTINE WritePsiToFile( Psi, FileNameList )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      CHARACTER(*), DIMENSION(2), INTENT(IN)  :: FileNameList
      INTEGER :: InpUnit, nGau, nAconfig
      LOGICAL :: FileExists

      ! Check that the WaveFunct object has been constructed already
      CALL ERROR( .NOT. Psi%WaveFunctIsSetup, " WritePsiToFile: WaveFunct instance is not setup ", ERR_OBJ_MISUSE )

      ! Check the existence of the file
      INQUIRE(FILE=FileNameList(1), EXIST=FileExists)
      CALL WARN( FileExists, " WritePsiToFile: input file "//trim(FileNameList(1))//" will be overwritten" )

      ! Open input, loop over the configurations and write coefficient, then close file
      InpUnit = LookForFreeUnit()
      OPEN(UNIT=InpUnit, FILE=TRIM(ADJUSTL(FileNameList(1))))
      DO nAconfig = 1,Psi%NrCfg
         WRITE(InpUnit,BVectorFormat) nAconfig, REAL(Psi%Bvector(nAconfig, 1)), AIMAG(Psi%Bvector(nAconfig, 1))
      END DO
      CLOSE(InpUnit)

      ! Check the existence of the file
      INQUIRE(FILE=FileNameList(2), EXIST=FileExists)
      CALL WARN( FileExists, " WriteGaussParToFile: input file "//trim(FileNameList(2))//" will be overwritten" )

      ! Open input, loop over the configurations and write coefficient, then close file
      InpUnit = LookForFreeUnit()
      OPEN(UNIT=InpUnit, FILE=TRIM(ADJUSTL(FileNameList(2))))
      DO nGau = 1,Psi%NrCfg*Psi%GDim
         WRITE(InpUnit,GaussFormat) nGau, REAL(Psi%GaussPar(1, nGau, 1)), AIMAG(Psi%GaussPar(1, nGau, 1)), &
                          &      REAL(Psi%GaussPar(2, nGau, 1)), AIMAG(Psi%GaussPar(2, nGau, 1)), &
                          &      REAL(Psi%GaussPar(3, nGau, 1)), AIMAG(Psi%GaussPar(3, nGau, 1))
      END DO
      CLOSE(InpUnit)

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "WritePsiToFile: the Bvector of a WaveFunct instance has been written to the file ", FileNameList(1)
      __WHITELINE_TO_LOG;
      __WRITE_TO_LOG "WritePsiToFile: the gaussian param.s of a WaveFunct instance have been written to the file ", FileNameList(2)
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif
   END SUBROUTINE WritePsiToFile


!*******************************************************************************
!          UpdateWaveFunctInstance
!*******************************************************************************
!> After the definition of new coefficients and parameters, update the
!> components of the WaveFunct instance to reflect the new values.
!> More precisely: 1) normalize gaussian functions, 2) update overlap
!> 3) initialize gaussian moments storage, 4) compute gaussian centers
!>
!> @param Psi            WaveFunct object
!*******************************************************************************
   SUBROUTINE UpdateWaveFunctInstance( Psi  )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      INTEGER :: iEl, nGau

      ! Check that the WaveFunct object has been constructed already
      CALL ERROR( .NOT. Psi%WaveFunctIsSetup, " UpdateWaveFunctInstance: WaveFunct instance is not setup ", ERR_OBJ_MISUSE )

      ! By default, renormalize the gaussian functions
      CALL NormalizeGaussians(Psi)
      ! Update the value of the overlap matrix ...
      Psi%Overlap = UpdateOverlapMatrix(Psi)
      ! Update the value of the average q's for the primitive gaussians
      DO iEl = 1, Psi%NrGauSets
         DO nGau = 1, Psi%NrPrimGau
            Psi%qCoord(nGau,iEl) =  - REAL(Psi%GaussPar(2, nGau, iEl)) / 2.0 / REAL(Psi%GaussPar(1, nGau, iEl))   ! REAL(Eta) / 2.0 / A
         END DO
      END DO

   END SUBROUTINE UpdateWaveFunctInstance


!*******************************************************************************
!          NormalizePsi
!*******************************************************************************
!> Normalize wavefunction by scaling its linear coefficients.
!>
!> @param Psi            WaveFunct to normalize
!*******************************************************************************
   SUBROUTINE NormalizePsi( Psi)
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT) :: Psi
      REAL :: Norm

      ! Compute norm
      Norm = WFNorm( Psi )

      ! Scale B vector by the square root of the norm
      Psi%Bvector(:, :) = Psi%Bvector(:, :)/ SQRT(Norm)

   END SUBROUTINE NormalizePsi


!*******************************************************************************
!          EquationsOfMotionSetup
!*******************************************************************************
!> Compute derivative of an input wavefunction and store it in a new instance
!> of the Derivative object.
!>
!> @param    InpEquationType
!> @param    InpPopThreshold
!> @param    InpPopThreshold
!*******************************************************************************
   SUBROUTINE EquationsOfMotionSetup( InpEquationType, InpSimplifyCoeffs, InpPopThreshold, InpCoeffThreshold )
      IMPLICIT NONE
      INTEGER, INTENT(IN)                        ::  InpEquationType
      LOGICAL, INTENT(IN)                        ::  InpSimplifyCoeffs
      REAL, DIMENSION(2), INTENT(IN), OPTIONAL   ::  InpPopThreshold
      REAL, INTENT(IN), OPTIONAL                 ::  InpCoeffThreshold

      ! Define the type of equations
      EquationType = InpEquationType
      
      ! Set flag to simplify the coefficients from the expression of the C matrix
      SimplifyCoefficients = InpSimplifyCoeffs
      SimplifyThreshold    = 1.E-8

      ! Set flag to freeze the Gaussian depending on the populations
      IF ( PRESENT(InpPopThreshold) ) THEN
         FreezeLowPopGaussians = .TRUE.
         PopMinThreshold = InpPopThreshold(1)
         PopMaxThreshold = InpPopThreshold(2)
      ELSE
         FreezeLowPopGaussians = .FALSE.
      END IF

      ! Set flag to freeze the Gaussian depending on the coefficients
      IF ( PRESENT(InpCoeffThreshold) ) THEN
         FreezeLowCoeffGaussians = .TRUE.
         CoeffMinThreshold = InpCoeffThreshold
      END IF

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "EquationsOfMotionSetup: setup of the equations of motion done "
      SELECT CASE( EquationType )
         CASE(EOM_STANDARD)
            __WRITE_TO_LOG "EquationsOfMotionSetup: using standard equations of motion "
         CASE(EOM_INVERSE)
            __WRITE_TO_LOG "EquationsOfMotionSetup: using inverse equations of motion (complementary proj. appearing in the coeff. equation)"
      END SELECT
      IF ( SimplifyCoefficients ) __WRITE_TO_LOG "WritePsiToFile: coeff.s are simplified from the parameter equations"
      IF ( FreezeLowPopGaussians ) THEN
         __WRITE_TO_LOG "WritePsiToFile: fixing GWPs with Mulliken populations in the interval", PopMinThreshold, " - ", PopMaxThreshold
      END IF
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif
   END SUBROUTINE EquationsOfMotionSetup   


!*******************************************************************************
!          ComputeDerivative
!*******************************************************************************
!> Compute derivative of an input wavefunction and store it in a new instance
!> of the Derivative object.
!>
!> @param    Psi          input wavefunction
!> @param    Hamiltonian  Hamiltonian of the system
!> @param    Forces       Operators to compute Ehrenfest forces in case of problematic gaussians
!> @param    InvMatCondNr array with the inverse condition number of the inverted matrices (optional)
!> @returns  PsiDeriv     on output derivative of the input wavefunction Psi
!*******************************************************************************
    SUBROUTINE ComputeDerivative( Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      TYPE(Derivative), INTENT(INOUT)  :: PsiDeriv
      TYPE(OperatorData), INTENT(IN)   :: Hamiltonian
      REAL, INTENT(IN)                 :: ActualTime
      REAL, DIMENSION(2), OPTIONAL, INTENT(OUT)  :: InvMatCondNr

      ! Set initial time
      CALL StartTimer(DerivativesClock)

      ! Check incompatibility between the EOM scheme chosen and the type of WF
      IF ( SimplifyCoefficients .AND. PsiDeriv%WFType == vMCG_SINGLESET ) THEN
         CALL AbortWithError( " ComputeDerivative: coeff. simplification cannot be used with single-set WF", ERR_MODULE_SETUP )
      END IF

      ! In case derivative has not been setup, allocate memory with the correct dimensionality of the wavefunction
      IF ( .NOT. PsiDeriv%DerivativeIsSetup ) THEN

         ! Store the type of wavefunction
         PsiDeriv%WFType = Psi%WFType

         ! ALLOCATE MEMORY FOR DERIVATIVES
         SELECT CASE( PsiDeriv%WFType )
            CASE( vMCG_SINGLESET, vMCG_MULTISET)
               ALLOCATE( PsiDeriv%BVector( SIZE(Psi%BVector,1), SIZE(Psi%BVector,2) ) )
               ALLOCATE( PsiDeriv%GaussPar( SIZE(Psi%GaussPar,1), SIZE(Psi%GaussPar,2), SIZE(Psi%GaussPar,3) ) )
            CASE DEFAULT
               CALL AbortWithError( " ComputeDerivative: do not recognize input wavefunction type", ERR_OBJ_MISUSE )
         END SELECT

         ! Derivative is now well defined
         PsiDeriv%DerivativeIsSetup = .TRUE.

      ELSE  ! Only checks if the type of wavefunction corresponds

         CALL ERROR( PsiDeriv%WFType /= Psi%WFType, " ComputeDerivative: wavefunction type mismatch", ERR_OBJ_MISUSE )

         ! ****************************************************
         ! ADD CHECK ON THE ACTUAL DIMENSIONS
         ! *****************************************************
      END IF

      ! Evaluate derivatives using the appropriate subroutine depending on the EOM type chosen
      IF ( EquationType == 1 ) THEN
          CALL ComputeDerivative_OldScheme(  Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      ELSE IF ( EquationType == 2 ) THEN
          CALL ComputeDerivative_NovelScheme(  Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      ELSE IF ( EquationType == 3 ) THEN
          CALL ComputeDerivative_SingleInversion(  Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      END IF

      ! Set final time
      CALL StopTimer(DerivativesClock)

   END SUBROUTINE ComputeDerivative
      
   SUBROUTINE ComputeDerivative_OldScheme( Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)      :: Psi
      TYPE(Derivative), INTENT(INOUT)  :: PsiDeriv
      TYPE(OperatorData), INTENT(IN)   :: Hamiltonian
      REAL, INTENT(IN)                 :: ActualTime
      REAL, DIMENSION(2), OPTIONAL, INTENT(OUT)  :: InvMatCondNr

      ! TEMPORARY MEMORY FOR DERIVATIVE CONSTRUCTION: General matrices
      COMPLEX, DIMENSION(:,:), ALLOCATABLE     :: FirstMom      ! matrix of the GWP config first moments
      COMPLEX, DIMENSION(:,:), ALLOCATABLE     :: SbetaDotSinv  ! product S^{beta,0} * S^{-1}
      COMPLEX, DIMENSION(:,:), ALLOCATABLE     :: ProjGDeriv    ! projection onto the GWP basis of the GWP derivatives
      COMPLEX, DIMENSION(:,:), ALLOCATABLE     :: ProjHamilt    ! projection onto the GWP basis of the Hamiltonian
      COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE :: HMatrix       ! representation of the Hamiltonian with the GWP basis

      ! TEMPORARY MEMORY FOR DERIVATIVE CONSTRUCTION: Matrices for the gaussian coefficients derivatives
      COMPLEX, DIMENSION(:,:,:), ALLOCATABLE   ::  OverlapInv   ! inverse  of the overlap matrix
      COMPLEX, DIMENSION(:), ALLOCATABLE       ::  CoeffRHS     ! R.H.S of the coefficient equation

      ! TEMPORARY MEMORY FOR DERIVATIVE CONSTRUCTION: Matrices for the gaussian parameters derivatives
      COMPLEX, DIMENSION(:,:), ALLOCATABLE     ::  CMatrix      ! C matrix
      COMPLEX, DIMENSION(:,:), ALLOCATABLE     ::  CMatrixInv   ! inverse of the  C matrix
      COMPLEX, DIMENSION(:), ALLOCATABLE       ::  YVector      ! Y vector

      ! integer indices
      INTEGER :: iElL,   iElR            ! Indices for loops over the electronic states
      INTEGER :: iCfgL,  iCfgR           ! Indices for loops over the GWPs
      INTEGER :: iDR,    iDL,    iD      ! Indices for loops over the dimensions of the system
      INTEGER :: iPrimR, iPrimL          ! Indices for loops over the primitive Gaussian functions
      INTEGER :: GauSet                  ! Index which defines the Gaussian set considered in the current iteration
      INTEGER :: iPtr,   jPtr            ! Integer index pointers
      INTEGER( SHORT_INTEGER_KIND )   :: NShort  ! Integer index for LAPACK calls

      ! other complex and real variables
      COMPLEX :: OverIFactor, Rho, HAlpha, TauMatrix
      REAL    :: InvCondNr, CMatrixCondMin, OvMatrixCondMin

      ! Compute 1/i (=-i) factor which is in the equations of motion
      OverIFactor = CMPLX(0.0,-1.0)

      ! ***************************************************************************************
      ! Preliminary construction of arrays that are needed to evaluate the derivatives
      ! ***************************************************************************************

      ALLOCATE( HMatrix(Psi%NrCfg,Psi%NrCfg,Psi%NrStates,Psi%NrStates), OverlapInv(Psi%NrCfg,Psi%NrCfg,Psi%NrGauSets) )

      ! 1) Compute Inverse of the Gaussians Overlap
      OvMatrixCondMin = 1.E+99
      DO GauSet = 1, Psi%NrGauSets
         CALL MatrixInversionDo( Psi%Overlap(:,:,GauSet,GauSet), OverlapInv(:,:,GauSet), 1, InvCondNr )
         IF ( InvCondNr < OvMatrixCondMin ) OvMatrixCondMin = InvCondNr
         ! Only upper diagonal has been computed, complete the inverse overlap
         DO iCfgR = 1, Psi%NrCfg
            DO iCfgL = iCfgR+1, Psi%NrCfg
               OverlapInv(iCfgL,iCfgR,GauSet) = CONJG(OverlapInv(iCfgR,iCfgL,GauSet))
            END DO
         END DO
      END DO

      ! 2) Calculate Hamiltonian matrix over configurations
      HMatrix = GetOperatorMatrix( Hamiltonian, Psi )

      ! ***************************************************************************************
      ! Compute the derivatives of Psi parameters and store them in PsiDeriv%GaussPar
      ! ***************************************************************************************

      ALLOCATE( SbetaDotSinv(Psi%NrPrimGau,Psi%NrCfg), ProjGDeriv(Psi%NrPrimGau,Psi%NrPrimGau), &
                ProjHamilt(Psi%NrPrimGau,Psi%NrCfg), FirstMom(Psi%NrCfg,Psi%NrPrimGau) )
      ALLOCATE(  CMatrix(Psi%NrPrimGau,Psi%NrPrimGau), CMatrixInv(Psi%NrPrimGau,Psi%NrPrimGau), YVector(Psi%NrPrimGau) )

      CALL StartTimer(GaussDerivClock)
      CMatrixCondMin = 1.E+99

      DO GauSet = 1, Psi%NrGauSets     ! left loop over the GWP sets (1 for single set, NrStates for multi-set)

         ! Matrix with the first moments of the configurations products < g_alpha | x_i | g_beta >
         FirstMom(:,:) = CMPLX(0.0, 0.0)
         DO iCfgR = 1, Psi%NrCfg       ! loop over right configuration / diagonal configuration
            jPtr = (iCfgR-1)*Psi%GDim    ! pointer on M1 position for the combined iCfgR+iD index
            DO iD = 1, Psi%GDim          ! diagonal elements
               FirstMom(iCfgR,jPtr+iD) = PsiPrimMoment( 1, Psi,iCfgR,iCfgR,iD,GauSet,GauSet) * Psi%Overlap(iCfgR,iCfgR,GauSet,GauSet)
            END DO
            DO iCfgL = iCfgR+1, Psi%NrCfg     ! loop over left configuration
               iPtr = (iCfgL-1)*Psi%GDim      ! pointer on M1 position for the combined iCfgL+iD index
               DO iD = 1, Psi%GDim      ! off-diagonal elements
                  FirstMom(iCfgL,jPtr+iD) = PsiPrimMoment( 1, Psi,iCfgL,iCfgR,iD,GauSet,GauSet) * Psi%Overlap(iCfgL,iCfgR,GauSet,GauSet)
                  FirstMom(iCfgR,iPtr+iD) = CONJG( FirstMom(iCfgL,jPtr+iD)  )
               END DO
            END DO
         END DO

         ! Construct for current GWP set some preliminary matrix: S^{beta,0} * S^{-1} and  S^{beta,0} * S^{-1} * S^{0,alpha}
         CALL TheOneWithMatrixMultiplication(SbetaDotSinv, FirstMom(:,:), OverlapInv(:,:,GauSet), "C","N")
         CALL TheOneWithMatrixMultiplication(ProjGDeriv, SbetaDotSinv(:,:), FirstMom(:,:), "N","N")

         ! Initialize C Matrix and Y vector for the current Gaussian set
         CMatrix = CMPLX(0.0,0.0); YVector = CMPLX(0.0,0.0)

         ! Construct C Matrix
         DO iPrimR = 1, Psi%NrPrimGau           ! loop over the right gaussian derivative
            iDR = Prim2Dim(iPrimR, Psi%GDim); iCfgR = Prim2Cfg(iPrimR, Psi%GDim);
            DO iCfgL = 1, Psi%NrCfg             ! loop over the left gaussian configurations (separated from the derivative index for conveniency)
               IF ( Psi%WFType == vMCG_MULTISET ) THEN   ! DEFINE RHO ELEMENT: for multi-set it is the B^\star B of the selected GWP set
                  IF ( SimplifyCoefficients ) THEN
                     Rho = CMPLX(1.0,0.0)
                  ELSE IF ( .NOT. SimplifyCoefficients ) THEN
                     Rho = CONJG(Psi%Bvector(iCfgL, GauSet)) * Psi%Bvector(iCfgR, GauSet)
                  END IF
               ELSE                                      ! DEFINE RHO ELEMENT: for single-set it is defined as trace over the states
                  Rho = TheOneWithVectorDotVector(Psi%Bvector(iCfgL, :), Psi%Bvector(iCfgR, :))
               END IF
               DO iDL = 1,  Psi%GDim          ! loop over the dimensions for left derivative
                  iPrimL = Cfg2Prim(iCfgL, iDL, Psi%GDim)
                  IF ( iPrimL < iPrimR ) CYCLE    ! if on upper half triangle, it will be computed by hermitian symmetry
                  IF ( iDL == iDR ) THEN
                     CMatrix(iPrimL,iPrimR) = Rho*( Psi%Overlap(iCfgL, iCfgR, GauSet, GauSet) &
                        * PsiPrimMoment(2, Psi,iCfgL,iCfgR,iDL,GauSet,GauSet) - ProjGDeriv(iPrimL,iPrimR) )
                  ELSE
                     CMatrix(iPrimL,iPrimR) = Rho*( Psi%Overlap(iCfgL,iCfgR,GauSet,GauSet) * PsiPrimMoment(1,Psi,iCfgL,iCfgR,iDL,GauSet,GauSet) &
                        * PsiPrimMoment(1, Psi,iCfgL,iCfgR,iDR,GauSet,GauSet) - ProjGDeriv(iPrimL,iPrimR) )
                  END IF
                  IF ( iPrimL /= iPrimR ) THEN    ! when element is not diagonal, compute hermitial conjugated element
                     CMatrix( iPrimR, iPrimL) = CONJG( CMatrix( iPrimL, iPrimR) )
                  END IF
               END DO
            END DO
         END DO

         ! Handle separately the construction of the Y vector for single- and multi-set
         ! Start with multiset: each electronic state has its own Y array
         IF ( Psi%WFType == vMCG_MULTISET ) THEN

            DO iElR = 1, Psi%NrStates        ! loop over the right electronic states
               ! Compute the matrix with proj hamiltonian (S^(alpha,0) * S^-1 * H )
               CALL TheOneWithMatrixMultiplication(ProjHamilt, SbetaDotSinv(:,:), HMatrix(:,:,GauSet, iElR), "N","N")
               DO iCfgR = 1, Psi%NrCfg             ! loop over the gaussian configurations iCfgL, iCfgR
                  DO iCfgL = 1, Psi%NrCfg
                     ! Compute rho matrix for the current couple of electronic states and configurations
                     IF ( SimplifyCoefficients ) THEN
                        Rho = Psi%Bvector(iCfgR, iElR)
                     ELSE IF ( .NOT. SimplifyCoefficients ) THEN
                        Rho = CONJG(Psi%Bvector(iCfgL, GauSet)) * Psi%Bvector(iCfgR, iElR)
                     END IF
                     DO iDL = 1, Psi%GDim                ! iDL index loop over the derivatives
                        ! Compute the hamiltonian matrix between configuration and 1st derivative HAlpha
                        HAlpha = OperatorMatrixElement(Hamiltonian, Psi%GaussPar(:,CfgBeg(iCfgL,Psi%GDim):CfgEnd(iCfgL,Psi%GDim),GauSet), &
                                                Psi%GaussPar(:,CfgBeg(iCfgR,Psi%GDim):CfgEnd(iCfgR,Psi%GDim),iElR), GauSet, iElR, iDL, 1 )
                        ! Combine everything and accumulate the value of rho
                        YVector(Cfg2Prim(iCfgL,iDL,Psi%GDim)) = YVector(Cfg2Prim(iCfgL,iDL,Psi%GDim)) + &
                                                            Rho * ( HAlpha - ProjHamilt(Cfg2Prim(iCfgL,iDL,Psi%GDim),iCfgR) )
                     END DO
                  END DO
               END DO
            END DO

         ! Now single set: construct a single Y array by summing over electronic coordinates, weighted by the rho matrix
         ELSE IF ( Psi%WFType == vMCG_SINGLESET ) THEN

            DO iElL = 1, Psi%NrStates        ! loop over the left electronic states
               DO iElR = 1, Psi%NrStates        ! loop over the right electronic states
                  ! Compute the matrix with proj hamiltonian (S^(alpha,0) * S^-1 * H )
                  CALL TheOneWithMatrixMultiplication(ProjHamilt, SbetaDotSinv(:,:), HMatrix(:,:,iElL, iElR), "N","N")
                  DO iCfgR = 1, Psi%NrCfg             ! loop over the gaussian configurations iCfgL, iCfgR
                     DO iCfgL = 1, Psi%NrCfg
                        ! Compute rho matrix for the current couple of electronic states and configurations
                        Rho = CONJG(Psi%Bvector(iCfgL, iElL)) * Psi%Bvector(iCfgR, iElR)
                        DO iDL = 1, Psi%GDim                ! iDL index loop over the derivatives
                           ! Compute the hamiltonian matrix between configuration and 1st derivative HAlpha
                           HAlpha = OperatorMatrixElement(Hamiltonian, Psi%GaussPar(:,CfgBeg(iCfgL,Psi%GDim):CfgEnd(iCfgL,Psi%GDim),iElL), &
                                                   Psi%GaussPar(:,CfgBeg(iCfgR,Psi%GDim):CfgEnd(iCfgR,Psi%GDim),iElR), iElL, iElR, iDL, 1 )
                           ! Combine everything and accumulate the value of rho
                           YVector(Cfg2Prim(iCfgL,iDL,Psi%GDim)) = YVector(Cfg2Prim(iCfgL,iDL,Psi%GDim)) + &
                                                            Rho * ( HAlpha - ProjHamilt(Cfg2Prim(iCfgL,iDL,Psi%GDim),iCfgR) )
                        END DO
                     END DO
                  END DO
               END DO
            END DO

         END IF

         ! invert C-matrix
         CALL MatrixInversionDo( CMatrix(:,:), CMatrixInv(:,:), 2, InvCondNr, CMask(1:NInv) )
         IF ( InvCondNr < CMatrixCondMin ) CMatrixCondMin = InvCondNr

         ! set the Y vector of fixed gaussians
         DO iPrimL = NInv+1, Psi%NrPrimGau
            YVector(CMask(iPrimL)) = CMPLX(0.0,0.0)
         END DO

         ! Initialize derivatives array to zero
         PsiDeriv%GaussPar(:, :, GauSet) = CMPLX(0.0,0.0)

         ! build derivatives vector  ( use ZHEMV because only the upper part of the inverse C matrix is stored)
         NShort = Psi%NrPrimGau
         CALL ZHEMV("U", NShort, OverIFactor, CMatrixInv(:,:), NShort, YVector(:), 1, CMPLX(0.0,0.0), PsiDeriv%GaussPar(2, 1, GauSet), 3)
         IF ( SimplifyCoefficients ) THEN
            DO iPrimL = 1, NInv
               iCfgL = Prim2Cfg(CMask(iPrimL), Psi%GDim)
               IF ( ABS(Psi%Bvector(iCfgL, GauSet)) <= SimplifyThreshold ) THEN
                  PsiDeriv%GaussPar(2, CMask(iPrimL), GauSet) = CONJG(Psi%Bvector(iCfgL, GauSet))*PsiDeriv%GaussPar(2, CMask(iPrimL), GauSet) &
                       / ( ABS(Psi%Bvector(iCfgL, GauSet))**2 + SimplifyThreshold**2 )
               ELSE
                  PsiDeriv%GaussPar(2, CMask(iPrimL), GauSet) = PsiDeriv%GaussPar(2, CMask(iPrimL), GauSet) / Psi%Bvector(iCfgL, GauSet)
               END IF
            END DO
         END IF

         ! Compute derivative of Eta that keeps the gaussian normalization  (with phase of the gaussians = 0 )
         DO iPrimR = 1, Psi%NrPrimGau
!             PsiDeriv%GaussPar(3, iPrimR, GauSet) =  - Psi%qCoord(iPrimR,GauSet) * REAL(PsiDeriv%GaussPar(2,iPrimR,GauSet)) + &
!                     (0.25/REAL(Psi%GaussPar(1,iPrimR,GauSet)) - Psi%qCoord(iPrimR,GauSet)**2) * REAL(PsiDeriv%GaussPar(1,iPrimR,GauSet))
            PsiDeriv%GaussPar(3, iPrimR, GauSet) = - Psi%qCoord(iPrimR,GauSet) * (PsiDeriv%GaussPar(2, iPrimR, GauSet) )
         END DO

      END DO

      CALL StopTimer(GaussDerivClock)

      ! ***************************************************************************************
      ! Now compute the derivatives of Psi coefficients and store them in PsiDeriv%BVector
      ! ***************************************************************************************

      ALLOCATE( CoeffRHS(Psi%NrCfg) )

      CALL StartTimer(BVecDerivClock)

      DO iElL = 1, Psi%NrStates     ! left loop over electronic states

         ! Initialize vector to zero
         PsiDeriv%BVector(:,iElL) = CMPLX(0.0,0.0)

         ! Define the Gaussian set for the overlap and the tau matrix depending on the type of wavefunction
         IF ( Psi%WFType == vMCG_MULTISET ) THEN
            GauSet = iElL
         ELSE IF ( Psi%WFType == vMCG_SINGLESET ) THEN
            GauSet = 1
         END IF
         ! initialize the r.h.s. of the equation for the current value, temporary store for current value of iElL
         CoeffRHS(:) = CMPLX(0.0,0.0)
         DO iCfgL = 1, Psi%NrCfg
            DO iCfgR = 1, Psi%NrCfg

               ! Construct the iCfgL,iCfgR element of the Tau matrix for the GauSet gaussian set
               TauMatrix = CMPLX(0.0,0.0)
               DO iDR = 1,  Psi%GDim          ! loop over the dimensions for right derivative
                  iPrimR = Cfg2Prim( iCfgR, iDR, Psi%GDim )
                  ! This is the contribution coming from A, Xi and Eta derivative (with phase of the gaussians = 0 )
                  TauMatrix = TauMatrix +  Psi%Overlap(iCfgL, iCfgR, GauSet, GauSet) * ( PsiDeriv%GaussPar(3, iPrimR, GauSet) + &
                                 PsiPrimMoment( 1, Psi,iCfgL,iCfgR,iDR,GauSet,GauSet ) * PsiDeriv%GaussPar(2,iPrimR,GauSet)   + &
                                 PsiPrimMoment( 2, Psi,iCfgL,iCfgR,iDR,GauSet,GauSet ) * PsiDeriv%GaussPar(1,iPrimR,GauSet)       )
               END DO

               ! Sum over the right electronic state, 1/i*hbar (H - tau) B
               DO iElR = 1, Psi%NrStates
                  IF ( iElR == iElL ) THEN
                     CoeffRHS(iCfgL) = CoeffRHS(iCfgL) + ( OverIFactor * HMatrix(iCfgL,iCfgR,iElL,iElR) - TauMatrix ) * Psi%Bvector(iCfgR, iElR)
                  ELSE IF ( iElR /= iElL ) THEN
                     CoeffRHS(iCfgL) = CoeffRHS(iCfgL) + OverIFactor * HMatrix(iCfgL,iCfgR,iElL,iElR) * Psi%Bvector(iCfgR, iElR)
                  END IF
               END DO

            END DO
         END DO

         ! Now for fixed iElL the sum over iElR is over, and we can transform with the inverse overlap
         PsiDeriv%BVector(:,iElL) = TheOneWithMatrixVectorProduct( OverlapInv(:,:,GauSet), CoeffRHS )
      END DO

      CALL StopTimer(BVecDerivClock)

      ! If required, return the inverse condition number of the inverted matrices
      IF ( PRESENT( InvMatCondNr ) ) InvMatCondNr(:) = (/ OvMatrixCondMin , CMatrixCondMin /)

      DEALLOCATE( HMatrix, OverlapInv )
      DEALLOCATE( CoeffRHS, SbetaDotSinv, ProjGDeriv, ProjHamilt, FirstMom )
      DEALLOCATE( CMatrix, CMatrixInv, YVector )

   END SUBROUTINE ComputeDerivative_OldScheme

   SUBROUTINE ComputeDerivative_NovelScheme( Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)      :: Psi
      TYPE(Derivative), INTENT(INOUT)   :: PsiDeriv
      TYPE(OperatorData), INTENT(IN)   :: Hamiltonian
      REAL, INTENT(IN)                 :: ActualTime
      REAL, DIMENSION(2), OPTIONAL, INTENT(OUT)  :: InvMatCondNr

      COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE :: HMatrix
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: MMatrix, dHMatrix
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: TMatrix, TMatrixInv
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: STilde, STildeInv
      COMPLEX, DIMENSION(:), ALLOCATABLE   :: rhs, rhs2, temp, temp2
      
      INTEGER :: iEl
      INTEGER :: iCfg, iD, iPrm, iOcc
      INTEGER :: jCfg, jD, jPrm, jOcc

      COMPLEX :: OverIFactor
      REAL    :: CMatrixCondMin, OvMatrixCondMin
      COMPLEX :: HTildeValue, Rho
      
      INTEGER( SHORT_INTEGER_KIND )   :: NShort

      ! Compute 1/i (=-i) factor which is in the equations of motion
      OverIFactor = CMPLX(0.0,-1.0)

      CALL StartTimer(BVecDerivClock)

      ! Allocate memory to store the various integrals contained in the equations of motion
      ALLOCATE( HMatrix(1:Psi%NrCfg,1:Psi%NrCfg,1,1), MMatrix(1:NInv,1:Psi%NrCfg), dHMatrix(1:NInv,1:Psi%NrCfg) )
      ALLOCATE( TMatrix(1:NInv,1:NInv), TMatrixInv(1:NInv,1:NInv) )
      ALLOCATE( rhs(1:Psi%NrCfg), STilde(1:Psi%NrCfg,1:Psi%NrCfg), STildeInv(1:Psi%NrCfg,1:Psi%NrCfg) )
      ALLOCATE( rhs2(1:NInv), temp(1:NInv), temp2(1:NInv) )

      ! DEFINITIONS:

      ! GWPS are here rewritten in the form (definitions are the same, it is just rearrangment of 0th order terms)
      ! G = exp( a x**2 + xi * (x-x_0) + Norm(a,RE(xi))  )

      ! HMatrix = matrix element of the Hamiltonian with GWPs              < G_i | H | G_j >
      ! MMatrix = overlap between partial derivative and GWP               < partial_alpha G_i | G_j >
      ! dHMatrix = Hamiltonian matrix element between part deriv and GWP   < partial_alpha G_i | H | G_j >
      ! TMatrix = overlap between partial derivatives                      < partial_alpha G_i | partial_beta G_j >

      ! ***************************************************************************************
      ! Now compute the derivatives of Psi coefficients/parameters and store them in PsiDeriv
      ! ***************************************************************************************

      iEl = 1      ! single set calculation

      ! First compute some matrices which are used more than once:
      ! Construct HMatrix
      HMatrix = GetOperatorMatrix( Hamiltonian, Psi )

      ! Construct MMatrix and dHMatrix
      DO jCfg = 1, Psi%NrCfg       ! loop over right configuration / diagonal configuration
         DO iOcc = 1, NInv
            iPrm = CMask(iOcc); iCfg = Prim2Cfg( iPrm, Psi%GDim ); iD = Prim2Dim( iPrm, Psi%GDim )
            IF ( SimplifyCoefficients ) THEN
               Rho = CMPLX(1.0,0.0)
            ELSE IF ( .NOT. SimplifyCoefficients ) THEN
               Rho = CONJG(Psi%BVector(iCfg,iEl))
            END IF
            MMatrix(iOcc,jCfg) = Rho * ( PsiPrimMoment( 1, Psi,iCfg,jCfg,iD,iEl,iEl) - Psi%qCoord(iPrm,iEl) )             &
                                                * Psi%Overlap(iCfg,jCfg,iEl,iEl)
            dHMatrix(iOcc,jCfg)= Rho * ( OperatorMatrixElement(Hamiltonian,                                               &
               Psi%GaussPar(:,CfgBeg(iCfg,Psi%GDim):CfgEnd(iCfg,Psi%GDim),iEl), Psi%GaussPar(:,CfgBeg(jCfg,Psi%GDim):CfgEnd(jCfg,Psi%GDim),iEl),   &
               iEl, iEl, AugDim=iD, AugF=1 ) - HMatrix(iCfg,jCfg,1,1) * Psi%qCoord(iPrm,iEl) )
                  
         END DO
      END DO

      ! Construct TMatrix
      DO jOcc = 1, NInv
         jPrm = CMask(jOcc); jCfg = Prim2Cfg( jPrm, Psi%GDim ); jD = Prim2Dim( jPrm, Psi%GDim )
         DO iOcc = 1, jOcc
            iPrm = CMask(iOcc); iCfg = Prim2Cfg( iPrm, Psi%GDim ); iD = Prim2Dim( iPrm, Psi%GDim )

               IF ( SimplifyCoefficients ) THEN
                  Rho = CMPLX(1.0,0.0)
               ELSE IF ( .NOT. SimplifyCoefficients ) THEN
                  Rho = CONJG(Psi%BVector(iCfg,iEl)) * Psi%BVector(jCfg,iEl)
               END IF

               IF ( iD == jD ) THEN
                  TMatrix(iOcc,jOcc) = ( PsiPrimMoment( 2, Psi,iCfg,jCfg,iD,iEl,iEl) - &
      Psi%qCoord(iPrm,iEl)*PsiPrimMoment( 1, Psi,iCfg,jCfg,jD,iEl,iEl) - Psi%qCoord(jPrm,iEl)*PsiPrimMoment( 1, Psi,iCfg,jCfg,iD,iEl,iEl) + &
      Psi%qCoord(iPrm,iEl)*Psi%qCoord(jPrm,iEl) ) * Psi%Overlap(iCfg,jCfg,iEl,iEl) *  Rho
               ELSE
                  TMatrix(iOcc,jOcc) = ( PsiPrimMoment( 1, Psi,iCfg,jCfg,iD,iEl,iEl) - Psi%qCoord(iPrm,iEl)) * &
                  ( PsiPrimMoment( 1, Psi,iCfg,jCfg,jD,iEl,iEl) - Psi%qCoord(jPrm,iEl) ) * Psi%Overlap(iCfg,jCfg,iEl,iEl) *  Rho
               END IF

         END DO
      END DO

      ! Invert TMatrix and store the inverse in TMatrixInv (Only upper diagonal is computed by MatrixInversionDo!!!)
      CALL MatrixInversionDo( TMatrix, TMatrixInv, 1, OvMatrixCondMin )
      
      ! Construct RHS and STilde that define  the equation for the derivative of the coefficients

      ! initialize rhs to zero
      rhs(:) = CMPLX( 0.0, 0.0 )

      NShort = NInv
               
      DO jCfg = 1, Psi%NrCfg
      
         ! Compute T^-1 * dHMatrix
         CALL ZHEMV("U", NShort, CMPLX(1.0,0.0), TMatrixInv(:,:), NShort, dHMatrix(:,jCfg), 1, CMPLX(0.0,0.0), temp(:), 1)
         ! Compute T^-1 * MMatrix
         CALL ZHEMV("U", NShort, CMPLX(1.0,0.0), TMatrixInv(:,:), NShort, MMatrix(:,jCfg), 1, CMPLX(0.0,0.0), temp2(:), 1)
         
         DO iCfg = 1, Psi%NrCfg

            ! Compute HTilde
            HTildeValue = HMatrix(iCfg, jCfg, 1, 1)
            DO jOcc = 1, NInv
               HTildeValue = HTildeValue - CONJG(MMatrix(jOcc,iCfg)) * temp(jOcc)
            END DO
            
            ! Accumulate rhs with 1/ihbar * HTilde(iCfg,jCfg) * BVector(jCfg)
            rhs(iCfg) = rhs(iCfg) + OverIFactor * HTildeValue * Psi%BVector(jCfg,iEl)
            
            ! skip elements of the lower triangle
            IF ( iCfg > jCfg ) CYCLE
            
            ! compute STIlde
            STilde(iCfg,jCfg) = Psi%Overlap(iCfg,jCfg,iEl,iEl)
            DO jOcc = 1, NInv
               STilde(iCfg,jCfg) = STilde(iCfg,jCfg) - CONJG(MMatrix(jOcc,iCfg)) * temp2(jOcc)
            END DO
            
         END DO
      END DO

      ! Only upper diagonal is computed!
      CALL MatrixInversionDo( STilde, STildeInv, 2, CMatrixCondMin )
      
      NShort = Psi%NrCfg
      CALL ZHEMV("U", NShort, CMPLX(1.0,0.0), STildeInv(:,:), NShort, rhs(:), 1, CMPLX(0.0,0.0), PsiDeriv%BVector(:, iEl), 1)
       
      CALL StopTimer(BVecDerivClock)
      CALL StartTimer(GaussDerivClock)

      DO iOcc = 1, NInv
         rhs2(iOcc) = CMPLX( 0.0, 0.0 )
         DO jCfg = 1, Psi%NrCfg
            rhs2(iOcc) = rhs2(iOcc) + OverIFactor*dHMatrix(iOcc,jCfg)*Psi%BVector(jCfg,iEl) - MMatrix(iOcc,jCfg)*PsiDeriv%BVector(jCfg,iEl)
         END DO
      END DO

      NShort = NInv
      CALL ZHEMV("U", NShort, CMPLX(1.0,0.0), TMatrixInv(:,:), NShort, rhs2(:), 1, CMPLX(0.0,0.0), temp(:), 1)
      IF ( SimplifyCoefficients ) THEN
         DO iOcc = 1, NInv
            iPrm = CMask(iOcc); iCfg = Prim2Cfg(iPrm, Psi%GDim)
            IF ( ABS(Psi%Bvector(iCfg, iEl)) <= SimplifyThreshold ) THEN
               temp(iOcc) = CONJG(Psi%Bvector(iCfg, iEl))*temp(iOcc) / ( ABS(Psi%Bvector(iCfg, iEl))**2 + SimplifyThreshold**2 )
            ELSE
               temp(iOcc) = temp(iOcc) / Psi%Bvector(iCfg, iEl)
            END IF
         END DO
      END IF

      PsiDeriv%GaussPar(:, :, iEl) = CMPLX( 0.0, 0.0 )
      DO iOcc = 1, NInv
         iPrm = CMask(iOcc)
         PsiDeriv%GaussPar(2, iPrm, iEl) = temp(iOcc)
         PsiDeriv%GaussPar(3, iPrm, iEl) = - Psi%qCoord(iPrm,iEl) * (PsiDeriv%GaussPar(2, iPrm, iEl) )
      END DO

      CALL StopTimer(GaussDerivClock)
      
      ! If required, return the inverse condition number of the inverted matrices
      IF ( PRESENT( InvMatCondNr ) ) InvMatCondNr(:) = (/ OvMatrixCondMin , CMatrixCondMin /)

      ! Deallocate memory used to store the parts of the equations of motion
      DEALLOCATE( HMatrix, MMatrix, dHMatrix, TMatrix, TMatrixInv, rhs, STilde, STildeInv, rhs2, temp, temp2 )

   END SUBROUTINE ComputeDerivative_NovelScheme


   SUBROUTINE ComputeDerivative_SingleInversion( Psi, PsiDeriv, Hamiltonian, ActualTime, InvMatCondNr )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)             ::  Psi
      TYPE(Derivative), INTENT(INOUT)            ::  PsiDeriv
      TYPE(OperatorData), INTENT(IN)             ::  Hamiltonian
      REAL, INTENT(IN)                           ::  ActualTime
      REAL, DIMENSION(2), OPTIONAL, INTENT(OUT)  ::  InvMatCondNr

      COMPLEX, DIMENSION(:,:), ALLOCATABLE       ::  VMatrix, HMatrix, VMatrixInv
      COMPLEX, DIMENSION(:), ALLOCATABLE         ::  ZVector, Solution

      INTEGER :: iEl
      INTEGER :: iCfg, iD, iPrm, iOcc
      INTEGER :: jCfg, jD, jPrm, jOcc
      INTEGER( SHORT_INTEGER_KIND )   :: NShort  ! Integer index for LAPACK calls

      COMPLEX :: OverIFactor
      REAL    :: ConditionNr
      COMPLEX :: Rho, dHMatrix

      ! Compute 1/i (=-i) factor which is in the equations of motion
      OverIFactor = CMPLX(0.0,-1.0)

      CALL StartTimer(BVecDerivClock)

      ! Only single surface calculation
      iEl = 1

      ! Allocate memory to store the various integrals contained in the equations of motion
      ALLOCATE( VMatrix(1:Psi%NrCfg+NInv,1:Psi%NrCfg+NInv), ZVector(1:Psi%NrCfg+NInv), HMatrix(1:Psi%NrCfg,1:Psi%NrCfg) )
      ALLOCATE( VMatrixInv(1:Psi%NrCfg+NInv,1:Psi%NrCfg+NInv), Solution(1:Psi%NrCfg+NInv) )

      ! ***************************************************************************************
      ! Now compute the V matrix (only upper triangle is constructed)
      ! ***************************************************************************************

      ! - first block, overlap between Gaussian configurations
      DO jCfg = 1, Psi%NrCfg
         DO iCfg = 1, jCfg
            VMatrix(iCfg,jCfg) = Psi%Overlap(iCfg,jCfg,iEl,iEl)
         END DO
      END DO

      ! - second block, overlap between derivatives of Gaussian configurations
      DO jOcc = 1, NInv
         jPrm = CMask(jOcc); jCfg = Prim2Cfg( jPrm, Psi%GDim ); jD = Prim2Dim( jPrm, Psi%GDim )
         DO iOcc = 1, jOcc
            iPrm = CMask(iOcc); iCfg = Prim2Cfg( iPrm, Psi%GDim ); iD = Prim2Dim( iPrm, Psi%GDim )
            IF ( SimplifyCoefficients ) THEN
               Rho = CMPLX(1.0,0.0)
            ELSE IF ( .NOT. SimplifyCoefficients ) THEN
               Rho = CONJG(Psi%BVector(iCfg,iEl)) * Psi%BVector(jCfg,iEl)
            END IF
            IF ( iD == jD ) THEN
               VMatrix(Psi%NrCfg+iOcc,Psi%NrCfg+jOcc) = ( PsiPrimMoment( 2, Psi,iCfg,jCfg,iD,iEl,iEl) - &
      Psi%qCoord(iPrm,iEl)*PsiPrimMoment( 1, Psi,iCfg,jCfg,jD,iEl,iEl) - Psi%qCoord(jPrm,iEl)*PsiPrimMoment( 1, Psi,iCfg,jCfg,iD,iEl,iEl) + &
      Psi%qCoord(iPrm,iEl)*Psi%qCoord(jPrm,iEl) ) * Psi%Overlap(iCfg,jCfg,iEl,iEl) *  Rho
            ELSE
               VMatrix(Psi%NrCfg+iOcc,Psi%NrCfg+jOcc) = ( PsiPrimMoment( 1, Psi,iCfg,jCfg,iD,iEl,iEl) - Psi%qCoord(iPrm,iEl)) * &
                  ( PsiPrimMoment( 1, Psi,iCfg,jCfg,jD,iEl,iEl) - Psi%qCoord(jPrm,iEl) ) * Psi%Overlap(iCfg,jCfg,iEl,iEl) *  Rho
            END IF
         END DO
      END DO

      ! - off-diagonal block, overlap between GWPs and their derivatives
      DO iCfg = 1, Psi%NrCfg
         DO jOcc = 1, NInv
            jPrm = CMask(jOcc); jCfg = Prim2Cfg( jPrm, Psi%GDim ); jD = Prim2Dim( jPrm, Psi%GDim )
            IF ( SimplifyCoefficients ) THEN
               Rho = CMPLX(1.0,0.0)
            ELSE IF ( .NOT. SimplifyCoefficients ) THEN
               Rho = Psi%BVector(jCfg,iEl)
            END IF
            VMatrix(iCfg,Psi%NrCfg+jOcc) = Rho * (PsiPrimMoment(1, Psi,iCfg,jCfg,jD,iEl,iEl) - Psi%qCoord(jPrm,iEl)) * Psi%Overlap(iCfg,jCfg,iEl,iEl)
         END DO
      END DO

      ! ***************************************************************************************
      ! Now compute the Z vector
      ! ***************************************************************************************

      ZVector = CMPLX(0.0,0.0)

      ! First part of the vector is composed by the H matrix elements of the configurations SUM_j < G_i | H | G_j > A_j
      DO jCfg = 1, Psi%NrCfg
         ! Diagonal element
         HMatrix(jCfg,jCfg) = OperatorMatrixElement(Hamiltonian, Psi%GaussPar(:,CfgBeg(jCfg,Psi%GDim):CfgEnd(jCfg,Psi%GDim),iEl),    &
                  Psi%GaussPar(:,CfgBeg(jCfg,Psi%GDim):CfgEnd(jCfg,Psi%GDim),iEl), iEl, iEl )
         ZVector(jCfg) = ZVector(jCfg) + HMatrix(jCfg,jCfg) * Psi%BVector(jCfg,iEl)
         DO iCfg = 1, Psi%NrCfg
            IF ( iCfg >= jCfg ) CYCLE
            ! Off diagonal elements ( H(i,j)*A(j), H(j,i)*A(i), using Hermitian symmetry H(j,i) = H(i,j)^\star )
            HMatrix(iCfg,jCfg) = OperatorMatrixElement(Hamiltonian, Psi%GaussPar(:,CfgBeg(iCfg,Psi%GDim):CfgEnd(iCfg,Psi%GDim),iEl),    &
                Psi%GaussPar(:,CfgBeg(jCfg,Psi%GDim):CfgEnd(jCfg,Psi%GDim),iEl), iEl, iEl )
            HMatrix(jCfg,iCfg) = CONJG(HMatrix(iCfg,jCfg))
            ZVector(iCfg) = ZVector(iCfg) + HMatrix(iCfg,jCfg) * Psi%BVector(jCfg,iEl)
            ZVector(jCfg) = ZVector(jCfg) + CONJG(HMatrix(iCfg,jCfg)) * Psi%BVector(iCfg,iEl)
         END DO
      END DO

      ! Second part of the vector is composed by the H matrix elements of the derivatives SUM_j < \partial_alpha G_i | H | G_j > A_j
      DO jCfg = 1, Psi%NrCfg       ! loop over right configuration / diagonal configuration
         DO iOcc = 1, NInv
            iPrm = CMask(iOcc); iCfg = Prim2Cfg(iPrm, Psi%GDim); iD = Prim2Dim(iPrm, Psi%GDim)
            dHMatrix =  ( OperatorMatrixElement(Hamiltonian, Psi%GaussPar(:,CfgBeg(iCfg,Psi%GDim):CfgEnd(iCfg,Psi%GDim),iEl),            &
      Psi%GaussPar(:,CfgBeg(jCfg,Psi%GDim):CfgEnd(jCfg,Psi%GDim),iEl), iEl, iEl, AugDim=iD, AugF=1 ) - HMatrix(iCfg,jCfg) * Psi%qCoord(iPrm,iEl) )
            IF ( SimplifyCoefficients ) THEN
               Rho = Psi%BVector(jCfg,iEl)
            ELSE IF ( .NOT. SimplifyCoefficients ) THEN
               Rho = CONJG(Psi%BVector(iCfg,iEl))*Psi%BVector(jCfg,iEl)
            END IF
            ZVector(Psi%NrCfg+iOcc) = ZVector(Psi%NrCfg+iOcc) + Rho * dHMatrix
         END DO
      END DO

      ! ***************************************************************************************
      ! Now solve the equation
      ! ***************************************************************************************

      ! Invert TMatrix and store the inverse in TMatrixInv (Only upper diagonal is computed by MatrixInversionDo!!!)
      CALL MatrixInversionDo( VMatrix, VMatrixInv, 1, ConditionNr )

      ! Now multiply RHS by inverse matrix to get the solution of the linear equation
      NShort = Psi%NrCfg+NInv
      CALL ZHEMV("U", NShort, CMPLX(1.0,0.0), VMatrixInv(:,:), NShort, ZVector(:), 1, CMPLX(0.0,0.0), Solution(:), 1)

      ! Copy the derivatives of the coefficients to the PsiDeriv%BVector vector
      PsiDeriv%BVector(:,iEl) = OverIFactor * Solution(1:Psi%NrCfg)

      ! If coefficients have been simplified from the V^-1 * Z product, the solution has to be divided by the coefficients
      IF ( SimplifyCoefficients ) THEN
         DO iOcc = 1, NInv
            iPrm = CMask(iOcc); iCfg = Prim2Cfg(iPrm, Psi%GDim)
            IF ( ABS(Psi%Bvector(iCfg, iEl)) <= SimplifyThreshold ) THEN
              Solution(Psi%NrCfg+iOcc) = CONJG(Psi%Bvector(iCfg, iEl))*Solution(Psi%NrCfg+iOcc) / &
                                        ( ABS(Psi%Bvector(iCfg, iEl))**2 + SimplifyThreshold**2 )
            ELSE
               Solution(Psi%NrCfg+iOcc) = Solution(Psi%NrCfg+iOcc) / Psi%Bvector(iCfg, iEl)
            END IF
         END DO
      END IF

      ! Now copy the derivatives of the parameters to the PsiDeriv%GaussPar array
      PsiDeriv%GaussPar(:, :, iEl) = CMPLX( 0.0, 0.0 )
      DO iOcc = 1, NInv
         iPrm = CMask(iOcc)
         PsiDeriv%GaussPar(2, iPrm, iEl) = OverIFactor * Solution(Psi%NrCfg+iOcc)
         PsiDeriv%GaussPar(3, iPrm, iEl) = - Psi%qCoord(iPrm,iEl) * (PsiDeriv%GaussPar(2, iPrm, iEl) )
      END DO

      CALL StopTimer(BVecDerivClock)
      CALL StartTimer(GaussDerivClock)
      CALL StopTimer(GaussDerivClock)

      ! If required, return the inverse condition number of the inverted matrices
      IF ( PRESENT( InvMatCondNr ) ) InvMatCondNr(:) = (/ ConditionNr , 0.0 /)

      ! Deallocate memory used to store the parts of the equations of motion
      DEALLOCATE( VMatrix, ZVector, HMatrix, VMatrixInv, Solution )

   END SUBROUTINE ComputeDerivative_SingleInversion

   
!*******************************************************************************
!          FixGaussianConfigurationsWithSmallPopulations
!*******************************************************************************
!> Construct a mask that contains information on which gaussian configurations
!> are considered in a simplyfied way (e.g. kept fixed, or propagated with
!> classical ehrenfest forces).
!>
!> @param    _____          ________________________________________________
!*******************************************************************************
   SUBROUTINE FixGaussianConfigurationsWithSmallPopulations( Psi )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)                  ::  Psi

      REAL, DIMENSION(Psi%NrCfg) :: Populations
      INTEGER :: k,j, iConf

      ! Prepare memory to store information on which C matrix element to invert
      IF ( ALLOCATED(CMask) ) DEALLOCATE( CMask )
      ALLOCATE( CMask(Psi%NrPrimGau) )

      IF ( FreezeLowCoeffGaussians ) THEN

         ! Put at the beginning of the vector those elements which will be included in the inversion
         k = 0
         DO iConf = 1, Psi%NrCfg
            IF ( ABS(Psi%Bvector(iConf, 1)) >= CoeffMinThreshold ) THEN
               k = k+1
               CMask( (k-1)*Psi%GDim+1:k*Psi%GDim ) = (/ ( (iConf-1)*Psi%GDim+j, j = 1,Psi%GDim ) /)
            END IF
         END DO
         NInv = k*Psi%GDim

         ! the rest of the Cmatrix indices are put at the end of the vector, beyond the NInv index
         DO iConf = 1, Psi%NrCfg
            IF ( ABS(Psi%Bvector(iConf, 1)) < CoeffMinThreshold ) THEN
               k = k+1
               CMask( (k-1)*Psi%GDim+1:k*Psi%GDim ) = (/ ( (iConf-1)*Psi%GDim+j, j = 1,Psi%GDim ) /)
            END IF
         END DO

      ELSE IF ( FreezeLowPopGaussians ) THEN

         ! Compute Mulliken populations of the gaussian configurations
         Populations = WFGauPopulations( Psi, 1 )

         ! Put at the beginning of the vector those elements which will be included in the inversion
         k = 0
         DO iConf = 1, Psi%NrCfg
            IF ( ABS(Populations(iConf)) >= PopMinThreshold .AND. ABS(Populations(iConf)) <= PopMaxThreshold ) THEN
               k = k+1
               CMask( (k-1)*Psi%GDim+1:k*Psi%GDim ) = (/ ( (iConf-1)*Psi%GDim+j, j = 1,Psi%GDim ) /)
            END IF
         END DO
         NInv = k*Psi%GDim

         ! the rest of the Cmatrix indices are put at the end of the vector, beyond the NInv index
         DO iConf = 1, Psi%NrCfg
            IF ( ABS(Populations(iConf)) < PopMinThreshold .OR. ABS(Populations(iConf)) > PopMaxThreshold ) THEN
               k = k+1
               CMask( (k-1)*Psi%GDim+1:k*Psi%GDim ) = (/ ( (iConf-1)*Psi%GDim+j, j = 1,Psi%GDim ) /)
            END IF
         END DO

      ELSE

         NInv = Psi%NrCfg*Psi%GDim
         DO iConf = 1, Psi%NrCfg
            CMask( (iConf-1)*Psi%GDim+1:iConf*Psi%GDim ) = (/ ( (iConf-1)*Psi%GDim+j, j = 1,Psi%GDim ) /)
         END DO

      END IF

   END SUBROUTINE FixGaussianConfigurationsWithSmallPopulations


!*******************************************************************************
!          DisposeDerivative
!*******************************************************************************
!> Destructor of Derivative instance.
!>
!> @param    PsiDeriv   deallocate memory for the Derivative instance PsiDeriv
!*******************************************************************************
   SUBROUTINE DisposeDerivative( PsiDeriv )
      IMPLICIT NONE
      TYPE(Derivative), INTENT(INOUT) :: PsiDeriv

      ! Deallocate memory only when the derivative has already been set up
      IF ( .NOT. PsiDeriv%DerivativeIsSetup ) RETURN

      ! ALLOCATE MEMORY FOR DERIVATIVES
      SELECT CASE( PsiDeriv%WFType )
         CASE( vMCG_SINGLESET, vMCG_MULTISET)
            DEALLOCATE( PsiDeriv%BVector, PsiDeriv%GaussPar )
         CASE DEFAULT
            CALL AbortWithError( " DisposeDerivative: do not recognize input wavefunction type", ERR_OBJ_MISUSE )
      END SELECT

      ! Derivative is now well defined
      PsiDeriv%DerivativeIsSetup = .FALSE.

   END SUBROUTINE DisposeDerivative






!*******************************************************************************
!          GetOverlapMax
!*******************************************************************************
!> Compute the maximum absolute value of the entries of the overlap matrix.
!>
!> @param Psi      Input wavefunction
!> @returns        Maximum absolute value of the overlap
!*******************************************************************************
   REAL FUNCTION  GetOverlapMax( Psi )
      IMPLICIT  NONE
      TYPE(WaveFunct), INTENT(IN)  :: Psi
      INTEGER                      :: i,j

      IF ( Psi%NrCfg == 1 ) THEN
         GetOverlapMax = 0.0
         RETURN
      END IF

      GetOverlapMax = ABS(Psi%Overlap(2,1,1,1))
      DO j = 1, SIZE(Psi%Overlap,2)
         DO i = j+1, SIZE(Psi%Overlap,1)
            IF (ABS(Psi%Overlap(i,j,1,1)) > GetOverlapMax)  GetOverlapMax = ABS(Psi%Overlap(i,j,1,1))
         END DO
      END DO

   END FUNCTION  GetOverlapMax



!*******************************************************************************
!          WFNorm
!*******************************************************************************
!> Compute the norm of a given wavefunction.
!>
!> @param Psi      Input wavefunction
!> @returns        Norm of the input wavefunction
!*******************************************************************************
   REAL FUNCTION WFNorm( Psi )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(IN) :: Psi
      COMPLEX, DIMENSION(Psi%NrCfg) :: TmpVec
      INTEGER :: iEl

      ! Initialize norm
      WFNorm = 0.0
      ! Overlap is already available, and for each electronic state Norm is simply  Bvec^Dagger * SMatrix * Bvec
      IF ( Psi%WFType == vMCG_SINGLESET ) THEN
         DO iEl = 1, Psi%NrStates
            TmpVec = TheOneWithMatrixVectorProduct( Psi%Overlap(:,:,1,1), Psi%Bvector(:, iEl) )
            WFNorm = WFNorm + REAL(TheOneWithVectorDotVector( Psi%Bvector(:, iEl), TmpVec ))
         END DO
      ELSE IF ( Psi%WFType == vMCG_MULTISET ) THEN
         DO iEl = 1, Psi%NrStates
            TmpVec = TheOneWithMatrixVectorProduct( Psi%Overlap(:,:,iEl,iEl), Psi%Bvector(:, iEl) )
            WFNorm = WFNorm + REAL(TheOneWithVectorDotVector( Psi%Bvector(:, iEl), TmpVec ))
         END DO
      END IF
      WFNorm = SQRT( WFNorm )

   END FUNCTION WFNorm


!*******************************************************************************
!          WFEnergy
!*******************************************************************************
!> Compute the expectation value of a given operator for a given wavefunction.
!>
!> @param Psi      Input wavefunction
!> @returns        Norm of the input wavefunction
!*******************************************************************************
   COMPLEX FUNCTION WFExpectation( Psi, Op )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      TYPE(OperatorData), INTENT(IN)   :: Op
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrStates,Psi%NrStates)  ::  OpMatrix
      COMPLEX, DIMENSION(Psi%NrCfg) :: TmpVec
      INTEGER :: iEl, jEl

      ! Compute operator representation over the gaussian basis
      OpMatrix = GetOperatorMatrix( Op, Psi )

      ! Initialize expectation value
      WFExpectation = CMPLX( 0.0, 0.0 )

      ! Overlap is already available, and for each electronic state Norm is simply  Bvec^Dagger * SMatrix * Bvec
      DO jEl = 1, Psi%NrStates
         DO iEl = 1, Psi%NrStates
            TmpVec = TheOneWithMatrixVectorProduct( OpMatrix(:,:,iEl,jEl), Psi%Bvector(:, jEl) )
            WFExpectation = WFExpectation + TheOneWithVectorDotVector( Psi%Bvector(:, iEl), TmpVec)
         END DO
      END DO

   END FUNCTION WFExpectation


!*******************************************************************************
!          WFOverlap
!*******************************************************************************
!> Compute the overlap between two input wavefunctions.
!>
!> @param Psi      Input wavefunction
!> @returns        Norm of the input wavefunction
!*******************************************************************************
   COMPLEX FUNCTION WFOverlap( Psi1, Psi2 )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(IN) :: Psi1, Psi2
      COMPLEX :: Rho
      INTEGER :: Cfg2, g2Start, g2End, Cfg1, g1Start, g1End, iEl, GauSet

      ! check that the wavefunction are compatible and are defined on the same space
      CALL ERROR( Psi1%WFType /= Psi2%WFType,     " WFOverlap: Psi1 and Psi2 have different type ", ERR_SUBS_INPUT )
      CALL ERROR( Psi1%NrStates /= Psi2%NrStates, " WFOverlap: Psi1 and Psi2 have different number of el states ", ERR_SUBS_INPUT )
      CALL ERROR( Psi1%GDim /= Psi2%GDim, " WFOverlap: Psi1 and Psi2 have different number of dofs ", ERR_SUBS_INPUT )

      WFOverlap = CMPLX( 0.0, 0.0)

      DO GauSet = 1, Psi1%NrGauSets     ! left loop over the GWP sets (1 for single set, NrStates for multi-set)

         ! First compute matrix with the overlap of the gaussian functions
         DO Cfg2 = 1, Psi2%NrCfg
            ! Start and end of the configuration in the primitive gaussian list
            g2Start = (Cfg2-1)*Psi2%GDim+1
            g2End   = Cfg2*Psi2%GDim

            DO Cfg1 = 1, Psi1%NrCfg
               ! Start and end of the configuration in the primitive gaussian list
               g1Start = (Cfg1-1)*Psi1%GDim+1
               g1End   = Cfg1*Psi1%GDim

               IF ( Psi1%WFType == vMCG_SINGLESET ) THEN
                  ! Compute rho summed over the diagonal electronic states
                  Rho = CMPLX( 0.0, 0.0)
                  DO iEl = 1, Psi1%NrStates
                     Rho = Rho + CONJG(Psi1%Bvector(Cfg1, iEl)) * Psi2%Bvector(Cfg2, iEl)
                  END DO
               ELSE IF ( Psi1%WFType == vMCG_MULTISET ) THEN
                  ! State specific rho element
                  Rho = CONJG(Psi1%Bvector(Cfg1, GauSet)) * Psi2%Bvector(Cfg2, GauSet)
               END IF

               ! Compute overlap by summing rho * gaussian overlap
               WFOverlap = WFOverlap + Rho * GauConf_Overlap( Psi1%GaussPar(:, g1Start:g1End, GauSet), &
                                                              Psi2%GaussPar(:, g2Start:g2End, GauSet) )

            END DO ! loop over Cfg1
         END DO ! loop over Cfg2

      END DO

   END FUNCTION WFOverlap


!*******************************************************************************
!          WFDifference
!*******************************************************************************
!> Difference between two wavefunctions Psi1 and Psi2 computed according
!> to some defined metrics. The available metrics are :
!> 1) distance induced by usual scalar product of wavefunction in Hilbert space
!> 2) distance induced by l_1 norm, wavefunction seen as sequence of complex numbers
!> 3) distance induced by l_2 norm, wavefunction seen as sequence of complex numbers
!> 4) distance induced by l_inf norm, wavefunction seen as sequence of complex numbers
!>
!> @param Psi1      First wavefunction
!> @param Psi2      Second wavefunction
!> @param DiffType  Definition of distance adopted
!> @returns         Difference between the two wavefunctions
!*******************************************************************************
   REAL FUNCTION WFDifference( Psi1, Psi2, DiffType )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(IN) :: Psi1
      TYPE(WaveFunct), INTENT(IN) :: Psi2
      INTEGER, INTENT(IN)       :: DiffType
      INTEGER :: nAconfig, nEl, j
      REAL :: tmpDiff, weight

      ! Initialize difference
      WFDifference = 0.0

      ! When l_1, l_2 or l_inf distances are chosen, the wavefunction need to have the very same number of elements
      IF ( DiffType == WFDISTANCE_L1_NORM .OR. DiffType == WFDISTANCE_L2_NORM .OR. DiffType == WFDISTANCE_LINF_NORM ) &
         CALL ERROR( Psi1%WFType   /=  Psi2%WFType    .OR. &
                     Psi1%NrStates /=  Psi2%NrStates  .OR. &
                     Psi1%NrCfg    /=  Psi2%NrCfg     .OR. &
                     Psi1%GDim   /=  Psi2%GDim,  " WFDifference: incompatible Psi1 and Psi2", ERR_OBJ_MISUSE )

      IF ( DiffType == WFDISTANCE_H_SCALARPROD ) THEN

         ! Distance is Norm(Psi1) + Norm(Psi2) - 2 Re Ovlp(Psi1, Psi2)
         WFDifference = SQRT( ABS(WFNorm( Psi1 )**2 + WFNorm( Psi2 )**2 - 2.0 * REAL(WFOverlap( Psi1, Psi2 ))) )

      ELSE   ! In this case we need to cycle over all the numbers of the wavefunction

         ! Loop over all the complex numbers in the Bvector array
         DO nEl = 1, Psi1%NrStates
            DO nAconfig = 1,Psi1%NrCfg
               tmpDiff = ABS( Psi2%Bvector(nAconfig, nEl) - Psi1%Bvector(nAconfig, nEl) )
               SELECT CASE( DiffType )
                  CASE( WFDISTANCE_L1_NORM )
                     WFDifference = WFDifference  + tmpDiff
                  CASE( WFDISTANCE_L2_NORM )
                     WFDifference = WFDifference  + tmpDiff**2
                  CASE( WFDISTANCE_LINF_NORM )
                     IF ( tmpDiff > WFDifference ) WFDifference = tmpDiff
                  CASE DEFAULT
                     CALL AbortWithError( "WFDifference: undefined value of DiffType", ERR_OBJ_MISUSE  )
               END SELECT
            END DO
         END DO

         ! Loop over all the complex numbers in the GaussPar array
         DO nEl = 1, Psi1%NrGauSets
            DO nAconfig = 1, Psi1%NrPrimGau
               DO j = 1, 3
                  weight = 0.5*ABS(Psi2%Bvector(nAconfig, nEl) + Psi1%Bvector(nAconfig, nEl))
                  tmpDiff = weight * ABS( Psi2%GaussPar(j, nAconfig, nEl) - Psi1%GaussPar(j, nAconfig, nEl) )
                  SELECT CASE( DiffType )
                     CASE( WFDISTANCE_L1_NORM )
                        WFDifference = WFDifference  + tmpDiff
                     CASE( WFDISTANCE_L2_NORM )
                        WFDifference = WFDifference  + tmpDiff**2
                     CASE( WFDISTANCE_LINF_NORM )
                        IF ( tmpDiff > WFDifference ) WFDifference = tmpDiff
                     CASE DEFAULT
                        CALL AbortWithError( "WFDifference: undefined value of DiffType", ERR_OBJ_MISUSE  )
                  END SELECT
               END DO
            END DO
         END DO

         ! if l2 norm, a square root is needed
         IF ( DiffType == WFDISTANCE_L2_NORM ) WFDifference = SQRT(WFDifference)

      END IF

   END FUNCTION WFDifference


!*******************************************************************************
!          WFGauCenters
!*******************************************************************************
!> Get the p's and q's for the center of each gaussian configuration along
!> a certain degree of freedom for a given input wavefunction Psi.
!>
!> @param Psi       Input wavefunction
!> @param iCoord    Defines which of coordinate of the centers is returned by the function
!> @param iEl       Optional parameter to define the gaussian set in a multi-set WF (by default, iEl=1)
!> @returns         Real array with the list of coordinate values
!*******************************************************************************
   FUNCTION  WFGauCenters( Psi, iCoord, iEl )   RESULT(PQCenters)
      IMPLICIT  NONE
      TYPE(WaveFunct), INTENT(IN)        :: Psi
      INTEGER, OPTIONAL, INTENT(IN)      :: iCoord
      INTEGER, OPTIONAL, INTENT(IN)      :: iEl
      REAL, DIMENSION(2*Psi%NrCfg)       :: PQCenters
      INTEGER   :: iCfg, iPrm

      DO iCfg = 1, Psi%NrCfg
         iPrm = Cfg2Prim(iCfg, iCoord, Psi%GDim)
         IF ( PRESENT(iEl) ) THEN      
            PQCenters(2*iCfg-1) = Psi%qCoord( iPrm, iEl )
            PQCenters(2*iCfg) = AIMAG( Psi%GaussPar(2, iPrm, iEl) ) + 2.0 * AIMAG(Psi%GaussPar(1, iPrm, iEl)) * Psi%qCoord( iPrm, iEl )
         ELSE
            PQCenters(2*(iCfg-1)+1) = Psi%qCoord( iPrm, 1 )
            PQCenters(2*iCfg) = AIMAG( Psi%GaussPar(2, iPrm, 1) ) + 2.0 * AIMAG(Psi%GaussPar(1, iPrm, 1)) * Psi%qCoord( iPrm, 1 )
         END IF
      END DO

   END FUNCTION  WFGauCenters

!*******************************************************************************
!          WFGauWidth
!*******************************************************************************
!> Get the deltaQ for each gaussian configuration along
!> a certain degree of freedom for a given input wavefunction Psi.
!>
!> @param Psi       Input wavefunction
!> @param iCoord    Defines for which of coordinate the width is returned by the function
!> @param iEl       Optional parameter to define the gaussian set in a multi-set WF (by default, iEl=1)
!> @returns         Real array with the list of coordinate values
!*******************************************************************************
   FUNCTION  WFGauWidth( Psi, iCoord, iEl )   RESULT(DeltaQ)
      IMPLICIT  NONE
      TYPE(WaveFunct), INTENT(IN)        :: Psi
      INTEGER, OPTIONAL, INTENT(IN)      :: iCoord
      INTEGER, OPTIONAL, INTENT(IN)      :: iEl
      REAL, DIMENSION(Psi%NrCfg)         :: DeltaQ
      INTEGER   :: iCfg

      DO iCfg = 1, Psi%NrCfg
         IF ( PRESENT(iEl) ) THEN
            DeltaQ(iCfg) = SQRT( ( -1.0 / 4.0 / REAL(Psi%GaussPar(1, Cfg2Prim(iCfg, iCoord, Psi%GDim), iEl))) )
         ELSE
            DeltaQ(iCfg) = SQRT( ( -1.0 / 4.0 / REAL(Psi%GaussPar(1, Cfg2Prim(iCfg, iCoord, Psi%GDim), 1)) ) )
         END IF
      END DO

   END FUNCTION  WFGauWidth


!*******************************************************************************
!          WFGauPopulations
!*******************************************************************************
!> Get the Mulliken populations of the gaussian basis functions
!> for a given input wavefunction Psi.
!>
!> @param Psi       Input wavefunction
!> @param Entropy   Optional output value of the entropy corresponding to the populations
!> @param iEl       Optional parameter to define the gaussian set in a multi-set WF (by default, iEl=1)
!> @returns         Real array with the values of the populations
!*******************************************************************************
   FUNCTION  WFGauPopulations( Psi, iEl, Entropy )   RESULT( Populations )
      IMPLICIT  NONE
      TYPE(WaveFunct), INTENT(IN)        :: Psi
      REAL, DIMENSION(Psi%NrCfg)         :: Populations
      INTEGER, OPTIONAL, INTENT(IN)      :: iEl
      REAL, OPTIONAL, INTENT(OUT)        :: Entropy
      INTEGER   :: iCfg, jCfg
      COMPLEX   :: SumPop

      DO iCfg = 1, Psi%NrCfg
         SumPop = CMPLX( 0.0, 0.0 )
         DO jCfg = 1,Psi%NrCfg
            IF ( PRESENT(iEl) ) THEN
               SumPop = SumPop + CONJG(Psi%BVector(iCfg,iEl))*Psi%Overlap(iCfg, jCfg, iEl, iEl)*Psi%BVector(jCfg,iEl)
            ELSE
               SumPop = SumPop + CONJG(Psi%BVector(iCfg,1))*Psi%Overlap(iCfg, jCfg, 1, 1)*Psi%BVector(jCfg,1)
            END IF
         END DO
         Populations(iCfg) = REAL(SumPop)
      END DO

      IF ( PRESENT(Entropy) ) THEN
         Entropy = 0.0
         IF ( Psi%NrCfg > 1 ) THEN
            DO jCfg = 1,Psi%NrCfg
               IF ( Populations(jCfg) > 0.0 ) THEN
                  Entropy = Entropy - Populations(jCfg) * LOG(Populations(jCfg))
               END IF
            END DO
            Entropy = Entropy / LOG(REAL(Psi%NrCfg))
         END IF
      END IF
      
   END FUNCTION  WFGauPopulations

   
!===========================================================================================================
!                             OPERATIONS DEFINED WITH DERIVED DATA TYPES
!===========================================================================================================


!*******************************************************************************
!          CopyPsi
!*******************************************************************************
!> WaveFunct copy: if the target WaveFunct is already allocated,
!> copy the values of the arrays. Otherwise create a new Psi object
!> allocated and initialized with same data.
!> Note that at the moment when the target WaveFunct has been already
!> set up, there is no check on the consistency of the dimensions.
!> In future this check can be included, anyway a misuse of this
!> subroutine should result only in a array dimension mismatch in
!> the copy process.
!>
!> @param       PsiSrc            WaveFunct source of the copy process
!> @returns     PsiTar            WaveFunct target of the copy process
!*******************************************************************************
   SUBROUTINE CopyPsi( PsiTar, PsiSrc )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(OUT)  :: PsiTar
      TYPE(WaveFunct), INTENT(IN)   :: PsiSrc

      ! In case memory for PsiSrc is not allocated, error (there is nothing to copy)
      CALL ERROR( .NOT. PsiSrc%WaveFunctIsSetup, " CopyPsi: source WaveFunct is not set ", ERR_OBJ_MISUSE )

      ! When the target WaveFunct is already allocated, assume that the dimensions are the same and just copy everything
      ! Otherwise, when the target WaveFunct is not allocated, allocate it prior to copy

      IF ( .NOT. PsiTar%WaveFunctIsSetup ) THEN

         ! Store the values of the wavefunction type and dimensional parameters
         PsiTar%WFType =    PsiSrc%WFType
         PsiTar%NrStates =  PsiSrc%NrStates
         PsiTar%NrCfg =     PsiSrc%NrCfg
         PsiTar%GDim =    PsiSrc%GDim
         PsiTar%NrGauSets = PsiSrc%NrGauSets
         PsiTar%NrPrimGau = PsiSrc%NrPrimGau

         ! Allocate memory according to the defined dimensions
         ALLOCATE( PsiTar%Bvector(PsiTar%NrCfg, PsiTar%NrStates) )
         ALLOCATE( PsiTar%GaussPar(3, PsiTar%NrPrimGau, PsiTar%NrGauSets) )
         ALLOCATE( PsiTar%Overlap(PsiTar%NrCfg, PsiTar%NrCfg, PsiTar%NrGauSets, PsiTar%NrGauSets) )
         ALLOCATE( PsiTar%qCoord(PsiTar%NrPrimGau, PsiTar%NrGauSets) )

         ! Wavefunction is now setup
         PsiTar%WaveFunctIsSetup = .TRUE.

! #if defined(LOG_FILE)
!          __OPEN_LOG_FILE;
!          __WRITE_TO_LOG "CopyPsi: an instance of WaveFunct has been allocated with: "
!          SELECT CASE( PsiTar%WFType )
!             CASE( vMCG_SINGLESET )
!                __WRITE_TO_LOG " Single-Set type wavefunction "
!             CASE( vMCG_MULTISET )
!                __WRITE_TO_LOG " Multi-Set type wavefunction "
!          END SELECT
!          __WRITE_TO_LOG "* NrStates = ", NumberToString(PsiTar%NrStates)
!          __WRITE_TO_LOG "* NrCfg = ", NumberToString(PsiTar%NrCfg)
!          __WRITE_TO_LOG "* GauDim = ", NumberToString(PsiTar%GDim)
!          __WRITE_TO_LOG "* Bvector dimension ", NumberToString(PsiTar%NrCfg*PsiTar%NrStates)
!          __WRITE_TO_LOG "* GaussPar dimension ", NumberToString(3*PsiTar%NrPrimGau*PsiTar%NrStates)
!          __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
! #endif
      END IF

      ! Copy the values of the arrays of PsiSrc to PsiTar
      PsiTar%Bvector  = PsiSrc%Bvector
      PsiTar%GaussPar = PsiSrc%GaussPar
      PsiTar%Overlap  = PsiSrc%Overlap
      PsiTar%qCoord   = PsiSrc%qCoord

   END SUBROUTINE CopyPsi


!*******************************************************************************
!          PsiPlusDeltaPsi
!*******************************************************************************
!> Update an (alredy initialized) instance of WaveFunct with the gaussian
!> parameters and the coefficients coming from another instance of WaveFunct.
!> The gaussian parameters and coefficients are updated by
!> adding the values of the Derivative input argument.
!>
!> @param   PsiTar           WaveFunct target of the copy process
!> @param   PsiSrc           WaveFunct source of the copy process
!> @param   DeltaPsi         Changes of the wavefunction
!> @param   Coeffs           Corresponding weights
!*******************************************************************************
   SUBROUTINE PsiPlusDeltaPsi( PsiTar, PsiSrc, DeltaPsi, Coeffs )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)     :: PsiTar
      TYPE(WaveFunct), INTENT(IN)        :: PsiSrc
      TYPE(Derivative), DIMENSION(:), INTENT(IN)  :: DeltaPsi
      REAL, DIMENSION(:), INTENT(IN)              :: Coeffs
      INTEGER :: i

      ! Also here the compatibility of PsiTar, PsiSrc and DeltaPsi is not checked...
      ! later maybe I will add the check, for the time being paying attentio to psi definitions

      ! Initialize PsiTar with PsiSrc
      PsiTar%Bvector   = PsiSrc%Bvector
      PsiTar%GaussPar  = PsiSrc%GaussPar

      ! Increment PsiTar with the DeltaPsi, weighted by the input coefficients
      DO i = 1, SIZE(DeltaPsi)
         PsiTar%Bvector   = PsiTar%Bvector  + Coeffs(i) * DeltaPsi(i)%Bvector
         PsiTar%GaussPar  = PsiTar%GaussPar + Coeffs(i) * DeltaPsi(i)%GaussPar
      END DO

      ! Update the components of the WaveFunct instance
      CALL  UpdateWaveFunctInstance( PsiTar  )

      ! This subroutine is meant to be used as a frequent call, since it is the
      ! way to update the value of a wavefunction. Further, no allocation/deallocation
      ! is done, and the only things that this piece of code contains are
      ! dimensionality checks and arithmetic operations. For both these reasons,
      ! we do not include in the log file any information coming from this subroutine.

   END SUBROUTINE PsiPlusDeltaPsi

!*******************************************************************************
!          CopyDerivative
!*******************************************************************************
!> Derivative copy: if the target Derivative is already allocated,
!> copy the values of the arrays. Otherwise create a new Psi object
!> allocated and initialized with same data.
!>
!> @param       DeltaPsiSrc            Derivative source of the copy process
!> @returns     DeltaPsiTar            Derivative target of the copy process
!*******************************************************************************
   SUBROUTINE CopyDerivative( DeltaPsiTar, DeltaPsiSrc )
      IMPLICIT NONE
      TYPE(Derivative), INTENT(OUT)  :: DeltaPsiTar
      TYPE(Derivative), INTENT(IN)   :: DeltaPsiSrc

      ! In case memory for PsiSrc is not allocated, error (there is nothing to copy)
      CALL ERROR( .NOT. DeltaPsiSrc%DerivativeIsSetup, " CopyDerivative: source Derivative is not set ", ERR_OBJ_MISUSE )

      ! When the target WaveFunct is already allocated, assume that the dimensions are the same and just copy everything
      ! Otherwise, when the target WaveFunct is not allocated, allocate it prior to copy

      IF ( .NOT. DeltaPsiTar%DerivativeIsSetup ) THEN

         DeltaPsiTar%WFType =  DeltaPsiSrc%WFType
         ! Allocate memory according to the defined dimensions
         ALLOCATE( DeltaPsiTar%Bvector(SIZE(DeltaPsiSrc%Bvector,1),SIZE(DeltaPsiSrc%Bvector,2)) )
         ALLOCATE( DeltaPsiTar%GaussPar(SIZE(DeltaPsiSrc%GaussPar,1),SIZE(DeltaPsiSrc%GaussPar,2),SIZE(DeltaPsiSrc%GaussPar,3)) )

         ! derivative is now setup
         DeltaPsiTar%DerivativeIsSetup = .TRUE.

      END IF

      ! Copy the values of the arrays of DeltaPsiSrc to DeltaPsiTar
      DeltaPsiTar%Bvector  = DeltaPsiSrc%Bvector
      DeltaPsiTar%GaussPar = DeltaPsiSrc%GaussPar

   END SUBROUTINE CopyDerivative


!*******************************************************************************
!          SumDerivatives
!*******************************************************************************
!> Sum two instances of Derivative.
!>
!> @param    DeltaPsi1, DeltaPsi2  input variations to sum.
!> @returns  DeltaPsiSum           on output sum of the variations
!*******************************************************************************
   TYPE(Derivative) FUNCTION SumDerivatives( DeltaPsi1, DeltaPsi2 )  RESULT( DeltaPsiSum )
      IMPLICIT NONE
      TYPE(Derivative), INTENT(IN)          :: DeltaPsi1, DeltaPsi2

      ! In case derivative has not been setup, allocate memory with the correct dimensionality of the wavefunction
      CALL ERROR( .NOT. DeltaPsi1%DerivativeIsSetup .OR. .NOT. DeltaPsi2%DerivativeIsSetup, &
            " SumDerivatives: using non initialized Derivative instances as summands", ERR_OBJ_MISUSE )

      ! Allocate memory for the dummy output
      CALL CopyDerivative(DeltaPsiSum, DeltaPsi1)

      ! Checks on the derivative dimensions are not done here, anyway if different array size are used output error is obvious
      DeltaPsiSum%BVector  = DeltaPsi1%BVector  + DeltaPsi2%BVector
      DeltaPsiSum%GaussPar = DeltaPsi1%GaussPar + DeltaPsi2%GaussPar

   END FUNCTION SumDerivatives


!*******************************************************************************
!          MultiplyDerivativeByScalar
!*******************************************************************************
!> Compute the product of a instances of Derivative by a scalar number.
!>
!> @param    alpha                 scalar number
!> @param    DeltaPsi              input variations to multiply by scalar
!> @returns  DeltaPsiSum           on output product of DeltaPsi by the scalar
!*******************************************************************************
   TYPE(Derivative) FUNCTION MultiplyDerivativeByScalar( alpha, DeltaPsi )  RESULT( DeltaPsiProd )
      IMPLICIT NONE
      REAL, INTENT(IN)                      :: alpha
      TYPE(Derivative), INTENT(IN)          :: DeltaPsi

      ! In case derivative has not been setup, allocate memory with the correct dimensionality of the wavefunction
      CALL ERROR( .NOT. DeltaPsi%DerivativeIsSetup , &
            " MultiplyDerivativeByScalar: using non initialized Derivative instances as input", ERR_OBJ_MISUSE )

      ! Allocate memory for the dummy output
      CALL CopyDerivative(DeltaPsiProd, DeltaPsi)

      ! Checks on the derivative dimensions are not done here, anyway if different array size are used output error is obvious
      DeltaPsiProd%BVector  = alpha * DeltaPsi%BVector
      DeltaPsiProd%GaussPar = alpha * DeltaPsi%GaussPar

   END FUNCTION MultiplyDerivativeByScalar


!===========================================================================================================
!                               PRIVATE SUBROUTINES AND FUNCTIONS
!===========================================================================================================


!*******************************************************************************
!          NormalizeGaussians
!*******************************************************************************
!> Given an input Psi with a defined set of gaussian parameters, redefine
!> the eta parameters to get normalized primitive gaussians for all degrees
!> of freedom and for all configurations
!>
!> @param Psi            WaveFunct object
!*******************************************************************************
   SUBROUTINE NormalizeGaussians(Psi)
      TYPE(WaveFunct), INTENT(INOUT)   :: Psi
      INTEGER :: iGau, iState, iStart, iEnd

      ! Check that the WaveFunct object has been constructed already
      CALL ERROR( .NOT. Psi%WaveFunctIsSetup, " NormalizeGaussians: WaveFunct instance is not setup ", ERR_OBJ_MISUSE )

      ! Loop over the number of sets of gaussians and over the primitive gaussians
      DO iState = 1, Psi%NrGauSets
         DO iGau = 1, Psi%NrCfg
            iStart = (iGau-1)*Psi%GDim+1
            iEnd   = iGau*Psi%GDim
            CALL GauConf_Normalize(Psi%GaussPar(:, iStart:iEnd, iState))
         END DO
      END DO

  END SUBROUTINE NormalizeGaussians


!*******************************************************************************
!          UpdateOverlapMatrix
!*******************************************************************************
!> Given an input Psi, compute and returns the Overlap matrix
!> between the gaussian configurations, i.e.
!>            S(i,j)  =  < g_i | g_j >
!>
!> @param    Psi    WaveFunct object
!> @returns  S      Computed overlap matrix between the gaussian configurations
!*******************************************************************************
   FUNCTION UpdateOverlapMatrix(Psi) RESULT( S )
      TYPE(WaveFunct), INTENT(IN)               :: Psi
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrGauSets,Psi%NrGauSets) :: S
      INTEGER :: iCfg, jCfg, iEl, jEl, jStart, jEnd, iStart, iEnd

      ! Loop over the set of gaussians corresponding to different electronic states (in case of single set, loop is trivial)
      DO jEl = 1, Psi%NrGauSets
         DO iEl = jEl, Psi%NrGauSets   ! loop only with iEl >= jEl, since off diagonal symmetric values are related by conjugation

            DO jCfg = 1, Psi%NrCfg
               ! Start and end of the configuration in the primitive gaussian list
               jStart = (jCfg-1)*Psi%GDim+1
               jEnd   = jCfg*Psi%GDim

               ! CONFIGURATION DIAGONAL ELEMENT OF THE OVERLAP
               IF ( jEl == iEl ) THEN      ! gaussian configurations are normalized, then
                  S(jCfg,jCfg,iEl,jEl) = CMPLX(1.0,0.0)
               ELSE                        ! same cfg is different on diff el states, need to explicitly compute overlap
                  S(jCfg,jCfg,iEl,jEl) = &
                      GauConf_Overlap( Psi%GaussPar(:, jStart:jEnd, iEl), Psi%GaussPar(:, jStart:jEnd, jEl) )
                  S(jCfg,jCfg,jEl,iEl) = CONJG(S(jCfg,jCfg,iEl,jEl))
               END IF

               DO iCfg = 1, Psi%NrCfg

                  IF (iEl == jEl .AND. iCfg < jCfg ) CYCLE  ! skip values for same electronic subsets that can be related by conjg

                  ! Start and end of the configuration in the primitive gaussian list
                  iStart = (iCfg-1)*Psi%GDim+1
                  iEnd   = iCfg*Psi%GDim

                  ! CONFIGURATION OFF-DIAGONAL ELEMENT OF THE OVERLAP
                  S(iCfg,jCfg,iEl,jEl) = &
                      GauConf_Overlap( Psi%GaussPar(:, iStart:iEnd, iEl), Psi%GaussPar(:, jStart:jEnd, jEl) )
                  S(jCfg,iCfg,jEl,iEl) = CONJG(S(iCfg,jCfg,iEl,jEl))

               END DO ! loop over iCfg
            END DO ! loop over jCfg

         END DO ! loop over iEl
      END DO ! loop over jEl

   END FUNCTION UpdateOverlapMatrix


!*******************************************************************************
!          PsiPrimMoment
!*******************************************************************************
!> This function is meant only to give a simpler interface to GauPrim_Moment
!> in the case of primitive gaussians which are part of the Psi object
!>
!> @param      Order        order of the gaussian moment
!> @param      Psi          WaveFunct object
!> @param      lCfg, rCfg   index of the left and right configurations
!> @param      iD           coordinate of the primitive moment
!> @param      iElL, iElR   electronic state index of the left and right config
!> @returns                 the complex value of the moment
!*******************************************************************************
   COMPLEX FUNCTION PsiPrimMoment( Order, Psi, lCfg, rCfg, iD, iElL, iElR )
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: Order
      TYPE(WaveFunct), INTENT(INOUT)  :: Psi
      INTEGER, INTENT(IN)             :: lCfg, rCfg
      INTEGER, INTENT(IN)             :: iD
      INTEGER, INTENT(IN)             :: iElL, iElR

      INTEGER :: lPrim, rPrim

      ! Compute the indices of the primitive gaussians of the iD-th coordinate for the two gaussian config.s
      lPrim = (lCfg-1)*Psi%GDim+iD; rPrim = (rCfg-1)*Psi%GDim+iD
      ! Use GauPrim_Moment to compute the primitive moment
      PsiPrimMoment = GauPrim_Moment(  Psi%GaussPar(:, lPrim, iElL), Psi%GaussPar(:, rPrim, iElR), Order)

   END FUNCTION PsiPrimMoment


!*******************************************************************************
!          GetOperatorMatrix
!*******************************************************************************
!> Given an input Operator of type OperatorData and an input Psi of type
!> WaveFunct, compute and returns the matrix containing the representation
!> of the operator over the configurations of the wavefunction.
!>  H(iCfg,jCfg,iEl,jEl)  =  < Psi_iCfg,iEl | Op_iEl,jEl | Psi_jCfg,jEl >
!>
!> @param    Operator       OperatorData object with the operator
!> @param    Psi            WaveFunct object with the wavefunction
!> @returns  H              the representation of the operator Op on the cfg.s
!*******************************************************************************
   FUNCTION GetOperatorMatrix( Operator, Psi ) RESULT( H )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(IN)   :: Operator
      TYPE(WaveFunct), INTENT(IN)      :: Psi
      COMPLEX, DIMENSION(Psi%NrCfg,Psi%NrCfg,Psi%NrStates,Psi%NrStates) :: H
      INTEGER     :: iGauL, iGauR, LStart, LEnd, RStart, REnd, elGauL, elGauR
      INTEGER     :: iElL, iElR

      ! Initialize matrix element
      H = CMPLX(0.0, 0.0)

      DO iElR = 1,Psi%NrStates        ! loops over the electronic states
         DO iElL = 1,Psi%NrStates

            ! Depending on the type of wavefunction (multiset vs singleset) choose gaussian set for computing matrix element
            IF ( Psi%WFType == vMCG_SINGLESET ) THEN
               elGauL = 1; elGauR = 1
            ELSE IF ( Psi%WFType == vMCG_MULTISET ) THEN
               elGauL = iElL; elGauR = iElR
            END IF

            ! when the electronic states are different the loop runs over all possible i,j values
            IF ( iElR /= iElL ) THEN

               DO iGauR = 1, Psi%NrCfg         ! loops over the gaussian configurations
                  DO iGauL = 1, Psi%NrCfg

                     ! Set the position of the primitive gaussian corresponding to the iGauR and iGauL configurations
                     LStart = (iGauL-1)*Psi%GDim+1; LEnd = iGauL*Psi%GDim
                     RStart = (iGauR-1)*Psi%GDim+1; REnd = iGauR*Psi%GDim

                     ! Compute the expectation value for the two gaussian configurations
                     H(iGauR,iGauL,iElL,iElR) = OperatorMatrixElement( Operator,  Psi%GaussPar(:, LStart:LEnd, elGauL), &
                                   Psi%GaussPar(:, RStart:REnd, elGauR), iElL, iElR )

                  END DO
               END DO

            ! otherwise, hermitian symmetry of the matrix can be exploited
            ! and only upper triangle is computed (lower triangle is complex conj)
            ELSE

               DO iGauR = 1, Psi%NrCfg

                  ! Set the position of the primitive gaussian corresponding to the iGauR configurations
                  RStart = (iGauR-1)*Psi%GDim+1; REnd = iGauR*Psi%GDim
                  ! Diagonal elements
                  H(iGauR,iGauR,iElL,iElR) = OperatorMatrixElement( Operator,  Psi%GaussPar(:, RStart:REnd, elGauL), &
                                   Psi%GaussPar(:, RStart:REnd, elGauR), iElL, iElR  )

                  DO iGauL = iGauR+1, Psi%NrCfg

                     ! Set the position of the primitive gaussian corresponding to the iGauL configurations
                     LStart = (iGauL-1)*Psi%GDim+1; LEnd = iGauL*Psi%GDim

                     H(iGauL,iGauR,iElL,iElR) = OperatorMatrixElement( Operator,  Psi%GaussPar(:, LStart:LEnd, elGauL),     &
                                   Psi%GaussPar(:, RStart:REnd, elGauR), iElL, iElR  )
                     H(iGauR,iGauL,iElL,iElR) = CONJG( H(iGauL,iGauR,iElL,iElR) )

                  END DO
               END DO

            END IF

         END DO
      END DO

   END FUNCTION GetOperatorMatrix


END MODULE PsiObject
