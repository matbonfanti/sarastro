!***************************************************************************************
!*                              MODULE SharedData
!***************************************************************************************
!
!>  \brief     Input data
!>  \details   This module include the data to be used by the main, 
!>             possibly set by reading a formatted input file
!>             it contains also the input subroutine and
!>             all the subroutines which may be used to test whether the
!>             data has been assigned a valid value
!
!***************************************************************************************
MODULE SharedData
#include "preprocessoptions.cpp"
   USE InputField
   IMPLICIT NONE

   PUBLIC          ! This module contain data that is supposed to be shared

!=============================================================================================================

      ! PARAMETERS



!=============================================================================================================

      ! MISCELLANEOUS AND GENERAL VARIABLES

      !  Variable to handle the command line
      INTEGER :: NArgs                   ! number of command line arguments
      LOGICAL :: Help = .FALSE.          ! if true, print help and stop

      !> Input file name, set from command line arguments
      CHARACTER(120) :: InputFileName

      ! variable definition for preprocessing directives
      __TIME_STRING__

      !> total number of print steps of the dynamics and corresponding integer loop index
      INTEGER :: NPrintSteps, iPrintStep

      !> Actual time of the dynamical propagation
      REAL :: ActualTime

      ! Expectation values and similar stuff computed during propagation
      !> Actual energy and energy at initial time
      REAL :: Energy, Energy0
      !> Actual norm and norm at initial time
      REAL :: Norm, Norm0
      !> Actual value of the autocorrelation function
      COMPLEX :: Autocorrelation

      !> Output unit where expectation values are written
      INTEGER :: ExpectUnit
      !> Output units where the centers of the gaussians are written
      INTEGER, DIMENSION(:), ALLOCATABLE :: GauCentersUnits, GauWidthUnits
      !> Output unit to write the gaussian populations over time
      INTEGER :: PopUnit

      !> Status of the integrator when exit from the DoProgation routine
      INTEGER :: IntegratorStat


!=============================================================================================================

      ! VARIABLES SET FROM INPUT

      !> Number of gaussian functions of the vMCG Ansatz
      INTEGER :: nGaussian
      !> Number of degrees of freedom of the space
      INTEGER :: gDim
      !> Number of electronic states
      INTEGER :: nStates
      !> Logical variable to define whether a multiset calculations is required
      LOGICAL :: MultiSet

      ! names of the files to set hamiltonian and initial conditions
      !> operator file (in the vMCG format, see subroutine in OperatorDefine.f90)
      CHARACTER(30) :: HamiltonianFile
      !> file with the initial coefficients of the vMCG wavefunction
      CHARACTER(30) :: BVectorFile
      !> file with the initial gaussian parameters of the vMCG wavefunction
      CHARACTER(30) :: GaussianParFile

      ! type of potential energy surface to adopt
      INTEGER :: HamiltonianOpDefine
      INTEGER, PARAMETER :: HAMILTONIAN_FILE   = 1,  & ! An analytic hamiltonian is defined from input file
                            HAMILTONIAN_DDYN   = 2     ! Potential is computed on-the-fly (M.Filatov's interface)

      ! potential approximation order, to define local harmonic (LHA) and local cubic (LCA) approximations
      INTEGER :: PotentialApproxOrder
      INTEGER, PARAMETER :: POTAPPROX_FULL      = 0,  & ! Do not apply approximations, use full potential
                            POTAPPROX_HARMONIC  = 2,  & ! Potential is computed up to 2nd order (LHA)
                            POTAPPROX_CUBIC     = 3     ! Potential is computed up to 3rd order (LCA)

      !> Max value of the integration step
      REAL :: MaxIntegrationStep
      !> Min value of the integration step
      REAL :: MinIntegrationStep

      !> the output is written every PrintStep time intervals
      REAL :: PrintStep
      !> total time of the wavefunction time propagation
      REAL :: TotalTime

      ! Integrator setup
      !> in integrating, use fixed steps regardless of the estimated error (in this case step = MaxIntegrationStep )
      LOGICAL :: UseFixedSteps
      !> error tolerance for step adaptation (ignored in case of UseFixedSteps)
      REAL :: ErrorTolerance = 1.E+99

      ! Matrix inverse regularization scheme
      !> type of regularization used for Cmatrix and OvMatrix inversion
      INTEGER :: CMatrixRegType, OvMatrixRegType
      !> relevant threshold of the type of regularization chosen
      REAL :: CMatrixRegValue, OvMatrixRegValue

      ! Freeze gaussians with low populations 
      !> logical variable to activate the freezing of the low-populated gaussians
      LOGICAL ::  FreezeLowPopGaussians
      !> thresholds that define which gaussian to freeze
      REAL    ::  PopMinThreshold, PopMaxThreshold
      
      !> Type of equations of motion used in the dynamics
      INTEGER :: EquationType
      !> Logical flag, when true coeff.s are simplified from the parameter equations
      LOGICAL :: SimplifyCoefficients
      
!=============================================================================================================
   CONTAINS
!=============================================================================================================

   SUBROUTINE ReadInputFile()
      IMPLICIT NONE
      TYPE(InputFile) :: InputData         ! Derived type to handle input data
      LOGICAL :: DirectDynamics

      ! READ AND PROCESS INPUT VARIABLES

      ! Open and read from input file the input parameters of the calculation
      CALL OpenFile( InputData, InputFileName )

      ! ********************************************************************
      !                   SYSTEM AND POTENTIAL SETUP
      ! ********************************************************************

      ! Read the number of gaussian funtions to use and their dimension
      CALL SetFieldFromInput( InputData, "n",    nGaussian )
      CALL SetFieldFromInput( InputData, "gdim", gDim      )

      ! Define the number of electronic states and the type of wavefunction
      nStates = 1   ! this line is temporary... later nStates will be read from input
      IF ( nStates == 1 ) THEN        ! in this case there is no need for multiset
         MultiSet = .TRUE.
      ELSE
         MultiSet = .TRUE.     ! als this line is temporary... later MultiSet will be defined from input
      END IF

      ! Decide whether potential is computed on-the-fly
!       CALL SetFieldFromInput( InputData, "ddyn", DirectDynamics, .FALSE. )
      DirectDynamics = .FALSE.
      IF ( DirectDynamics ) THEN
         HamiltonianOpDefine = HAMILTONIAN_DDYN
      ELSE IF ( .NOT. DirectDynamics ) THEN
         HamiltonianOpDefine = HAMILTONIAN_FILE
      END IF

      ! Read the name of the operator file
      IF ( HamiltonianOpDefine == HAMILTONIAN_FILE ) THEN
         CALL SetFieldFromInput( InputData, "hamiltonian", HamiltonianFile )
         HamiltonianFile = ADJUSTL(HamiltonianFile)
      END IF

!       ! Optionally switch on a local approximation of the potential (no approx is default)
!       CALL SetFieldFromInput( InputData, "pot_localexp_order", PotentialApproxOrder, POTAPPROX_FULL )
      PotentialApproxOrder = POTAPPROX_FULL


      ! ********************************************************************
      !                   INITIAL CONDITIONS
      ! ********************************************************************

      ! Read the name of the file with the initial coefficients of the vMCG wavefunction
      CALL SetFieldFromInput( InputData, "read_avec", BVectorFile )
      BVectorFile = ADJUSTL(BVectorFile)
      ! Read the name of the file with the initial gaussian parameters of the vMCG wavefunction
      CALL SetFieldFromInput( InputData, "read_gaussians", GaussianParFile )
      GaussianParFile = ADJUSTL(GaussianParFile)


      ! ********************************************************************
      !             PROPAGATION TIME AND INTEGRATION STEPS
      ! ********************************************************************

      ! Time steps and similia
      CALL SetFieldFromInput( InputData, "maxstep", MaxIntegrationStep )
      MaxIntegrationStep = MaxIntegrationStep * MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "minstep", MinIntegrationStep )
      MinIntegrationStep = MinIntegrationStep * MyConsts_fs2AU

      CALL SetFieldFromInput( InputData, "maxtime", TotalTime )
      TotalTime = TotalTime * MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "writetime", PrintStep )
      PrintStep = PrintStep * MyConsts_fs2AU
      ! Compute number of steps between each log printout
      NPrintSteps = INT( TotalTime/ PrintStep )


      ! ********************************************************************
      !                 INTEGRATOR AND INVERSION SETUP
      ! ********************************************************************

      ! choose whether you want to use fixed integration steps without adaptation to error estimate
      CALL SetFieldFromInput( InputData, "use_fixed_steps", UseFixedSteps, .FALSE. )
      ! define error tolerance for step adaptation (ignored in case of UseFixedSteps)
      IF ( .NOT. UseFixedSteps ) CALL SetFieldFromInput( InputData, "tolerance", ErrorTolerance )

      ! Matrix regularization scheme: type of regularization for CMatrix and OVMAtrix ( default is EIGENVALUES )
      CALL SetFieldFromInput( InputData, "cmatrix_reg_type", CMatrixRegType, 1 )
      CALL SetFieldFromInput( InputData, "ovmatrix_reg_type", OvMatrixRegType, 1 )
      ! Matrix regularization scheme: numerical threshold appropriate for the method chosen
      CALL SetFieldFromInput( InputData, "cmatrix_eps", CMatrixRegValue, 1.E-6 )
      CALL SetFieldFromInput( InputData, "ovmatrix_eps", OvMatrixRegValue, 1.E-6 )

      ! define whether to activate the freezing of the gaussian with low populations 
      CALL SetFieldFromInput( InputData, "freeze_lowpopgau", FreezeLowPopGaussians, .FALSE. )
      ! define the corresponding thresholds 
      IF ( FreezeLowPopGaussians )  THEN
         CALL SetFieldFromInput( InputData, "popmin_threshold", PopMinThreshold, 1.E-06 )
         CALL SetFieldFromInput( InputData, "popmax_threshold", PopMaxThreshold, 1.E+99 )
      END IF

      ! define the type of equations of motion used in the dynamics
      CALL SetFieldFromInput( InputData, "eom_type", EquationType, 1 )
      
      ! define if the coefficients are simplified from the parameter equations
      CALL SetFieldFromInput( InputData, "eom_simplifycoeff", SimplifyCoefficients, .FALSE. )
      
      ! Close the input file
      CALL CloseFile( InputData )

      ! CHECK INPUT VARIABLES

      ! Check the value of the variable definying the type of hamiltonian to be used
      CALL CheckHamiltonianOpDefine( HamiltonianOpDefine )
      ! Check the value of the variable definying the approximation order of the potential
      CALL CheckPotentialApproxOrder( PotentialApproxOrder )
      ! Check that a full potential is not required when using the direct dynamics setup
      CALL ERROR( PotentialApproxOrder == POTAPPROX_FULL .AND. HamiltonianOpDefine == HAMILTONIAN_DDYN , &
       " ReadInputFile: full potential not available for direct dynamics run ", ERR_WRONG_INP )


      ! PRINT TO OUTPUT THE DEFINED VALUES OF THE VARIABLES

      WRITE(*,200)  gDim, nStates
      IF (MultiSet) THEN
         WRITE(*,210) nGaussian, gDim*nGaussian
      ELSE IF (.NOT. MultiSet) THEN
         WRITE(*,211) nGaussian, gDim*nGaussian
      END IF
      WRITE(*,220)  TRIM(ADJUSTL(BVectorFile)), TRIM(ADJUSTL(GaussianParFile))

      IF ( HamiltonianOpDefine == HAMILTONIAN_FILE ) THEN
         IF ( PotentialApproxOrder == POTAPPROX_FULL ) THEN
            WRITE(*,230) TRIM(ADJUSTL(HamiltonianFile)), " No Local Approx on the Potential        "
         ELSE
            WRITE(*,230) TRIM(ADJUSTL(HamiltonianFile)), &
               " Local Potential Approx of Order       = ", PotentialApproxOrder
         END IF
      ELSE IF ( HamiltonianOpDefine == HAMILTONIAN_DDYN ) THEN
         WRITE(*,231) PotentialApproxOrder
      END IF

   200 FORMAT  (/,' *******************  System information  *******************',/, &
                  " Nuclear degrees of freedom            = ",             I20,/, &
                  " Number of Electronic States           = ",             I20,/  )
   210 FORMAT  (  ' *******************  Wavefunction Ansatz *******************',/, &
                  " Number of Configurations              = ",             I20,/, &
                  " Number of Primitive Gaussians         = ",             I20,/, &
                  " Using the Multi-Set Formalism           ",                /  )
   211 FORMAT  (  ' *******************  Wavefunction Ansatz *******************',/, &
                  " Number of Configurations              = ",             I20,/, &
                  " Number of Primitive Gaussians         = ",             I20,/, &
                  " Using the Single-Set Formalism          ",                /  )
   220 FORMAT  (  ' ******************* Initial Wavefunction *******************',/, &
                  " Reading b-Vector from File            = ",             A20,/, &
                  " Reading Gaussian Parameters from File = ",             A20,/  )
   230 FORMAT  (  ' *******************  Hamiltonian Setup   *******************',/, &
                  " Hamiltonian Operator Read from File   = ",             A20,/, &
                  A41                                        ,             I20,/  )
   231 FORMAT  (  ' *******************  Hamiltonian Setup   *******************',/, &
                  " Potential is Computed On-The-Fly        ",             20X,/, &
                  " Local Potential Approx of Order       = ",             I20,/  )

   END SUBROUTINE ReadInputFile

!=============================================================================================================

   !> Subroutine to check the availability of a given Hamiltonian Operator option
   SUBROUTINE CheckHamiltonianOpDefine( IntNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IntNr
      LOGICAL :: Check

      Check = ( IntNr /= HAMILTONIAN_FILE .AND. &
                IntNr /= HAMILTONIAN_DDYN )
      CALL ERROR( Check, " CheckHamiltonianOpDefine: Invalid HamiltonianOpDefine option ", ERR_WRONG_INP )
   END SUBROUTINE CheckHamiltonianOpDefine


   !> Subroutine to check the availability of a given potential approximation order option
   SUBROUTINE CheckPotentialApproxOrder( IntNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IntNr
      LOGICAL :: Check

      Check = ( IntNr /= POTAPPROX_FULL     .AND. &
                IntNr /= POTAPPROX_HARMONIC .AND. &
                IntNr /= POTAPPROX_CUBIC )
      CALL ERROR( Check, " CheckPotentialApproxOrder: Invalid PotentialApproxOrder option ", ERR_WRONG_INP )
   END SUBROUTINE CheckPotentialApproxOrder

!=============================================================================================================

END MODULE SharedData
