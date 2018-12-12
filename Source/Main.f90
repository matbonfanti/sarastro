!***************************************************************************************
!*                              PROGRAM sarastro
!***************************************************************************************
!>  \mainpage      Program sarastro
!>
!>  quantum dynamics with a gaussian-type Ansatz wavefunction             \n
!>
!>  \author        Matteo Bonfanti
!>  \version       VERSIONTAG
!>  \date          1 July 2017
!>
!***************************************************************************************
PROGRAM sarastro
#include "preprocessoptions.cpp"
   USE SharedData
   USE PsiObject
   USE EmptyOptimize
   USE OperatorDefine
   USE MatrixInversion
   USE RungeKutta45

   USE GauConf

   IMPLICIT NONE

   !*************************************************************
   !   OBJECTS USED IN THE CODE  (they are here, not in SharedData)
   !*************************************************************

   ! Initial Wavefunction and actual wavefunction
   TYPE(WaveFunct) :: Psi0, PsiT
   ! Object which stores the hamiltonian definition
   TYPE(OperatorData) :: Hamiltonian

   INTEGER :: i
   CHARACTER(50) :: String

   !*************************************************************
   !   INITIAL MESSAGES AND MISCELLANOUS STUFF
   !*************************************************************

   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                               _________          ')"
   PRINT "(       '                    ==============================',/)"
   PRINT "(       '                       Author: Matteo Bonfanti'      )"
   PRINT "(       '                       Release: ',A)", VERSIONTAG
   PRINT "(       '                       Compilation: ',A,1X,A,/)", __DATE__, __TIME__

   PRINT "(       '                       << ______________________      ')"
   PRINT "(       '                         _______________________      ')"
   PRINT "(       '                        _______________________... >> '/)"
   PRINT "(       '                      [__________________________]  ',2/)"

   ! Set timers
   CALL InternalTimersSetup( )
   ! Set initial time
   CALL StartTimer(TotalClock)

   ! Set the machine precision at runtime
   CALL Calculate_Constant_EPS()

#if defined(LOG_FILE)
   __INIT_LOG_FILE
#endif

   !*************************************************************
   !         COMMAND LINE ARGUMENT
   !*************************************************************

   ! Check and read from command line the input file name
   NArgs = COMMAND_ARGUMENT_COUNT()
   IF (NArgs<1) THEN
      Help = .TRUE.
   ELSE
      CALL GET_COMMAND_ARGUMENT( 1, InputFileName )
      IF ( trim(InputFileName) == "help" ) Help = .TRUE.
   ENDIF
   IF (Help) THEN ! Call help
      PRINT*, ' Launch this program as:'
      PRINT*, ' % vmcg "InputFileName" '
      STOP
   ENDIF

   !*************************************************************
   !                 INPUT FILE
   !*************************************************************

   ! First, read the input variables from input file
   ! only general variables are read here (potential setup is done later)
   ! input variables are defined and stored in the SharedData module
   CALL ReadInputFile()

   ! set up the HAMILTONIAN OPERATOR
   SELECT CASE( HamiltonianOpDefine )

      ! in case of a hamiltonian defined via input file, set up the module OperatorDefine.f90
      ! which reads the hamiltonian input file and prepare the necessary data
      CASE( HAMILTONIAN_FILE )
         ! Read hamiltonian operator from file
         CALL SetOperatorFromFile( Hamiltonian, gDim, nStates, HamiltonianFile )
         ! Error when trying to use LHA or LCA which are not yet implemented
         CALL ERROR( PotentialApproxOrder /= POTAPPROX_FULL, " Local V approximation not yet implemented " , ERR_WRONG_INP )

      ! for a direct dynamics run, ... (not implememented yet)
      CASE ( HAMILTONIAN_DDYN )
         ! temporarily give error when trying to run direct dynamics ... it will implemented later
         CALL AbortWithError( " Direct dynamics not yet implemented", ERR_WRONG_INP )
         ! In this case the Hamiltonian needs to be setup with the kinetic operator in normal modes
         CALL SetOperatorKinOnly( Hamiltonian, gDim )
         ! **********************************************************
         ! Setup direct dynamics
         ! CALL DIRECT_DYNAMICS_SETUP()
         ! **********************************************************

   END SELECT
   ! Print infos on the Hamiltonian setup
   CALL LogOperator( Hamiltonian, "hamiltonian.log" )

   ! ****************************************************************************
   !   Define initial wavefunction
   ! ****************************************************************************

   ! Define a PsiData instance to store the initial wavefunction
   IF ( MultiSet ) THEN
      CALL SetPsi( Psi0, nStates, nGaussian, gDim, vMCG_MULTISET )
   ELSE
      CALL SetPsi( Psi0, nStates, nGaussian, gDim, vMCG_SINGLESET )
   END IF
   ! Define the initial wavefunction reading coefficients and parameters from input
   CALL ReadPsiFromFile( Psi0, (/ BVectorFile, GaussianParFile /) )

   ! Normalize initial wavefunction
   CALL NormalizePsi( Psi0 )

   ! Write to output file the actual (checked and normalized)
   ! initial values of the coeff.s and param.s
   CALL WritePsiToFile( Psi0, (/ "psi_input_bvec.inp     ", "psi_input_gaussians.inp" /) )

   ! ****************************************************************************
   ! Initialize dynamical propagation: set propagator details and
   ! compute the data which is needed for first expectation computation
   ! ****************************************************************************

   ! Regularized inversion setup
   CALL MatrixInversionSetup( 3, (/ "overlap      ", "cmatrix      ", "projemptyovlp" /), &
         (/ OvMatrixRegType, CMatrixRegType, TIKHONOV /), (/ OvMatrixRegValue, CMatrixRegValue, 1.E-12 /) )

      ! Optimize initial unoccupied states
!    CALL OptimizeEmptyGaussians( Psi0, Hamiltonian )

   ! Write to output file the optimized gaussians
!    CALL WritePsiToFile( Psi0, (/ "psi_opt_bvec.inp     ", "psi_opt_gaussians.inp" /) )

   ! Integrator setup for Runge-Kutta 4(5)
   CALL SetupRungeKutta45( UseFixedSteps, MinStep=MinIntegrationStep, MaxStep=MaxIntegrationStep, InpTolerance=ErrorTolerance )
   
   ! setup the equations of motion and optionally freeze gaussians with small-large populations / small coefficients
   ! the condition on the coefficient has priority over the condition on populations
   IF ( FreezeLowCoeffGaussians ) THEN
      CALL EquationsOfMotionSetup( EquationType, SimplifyCoefficients, InpCoeffThreshold=CoeffMinThreshold )
   ELSE IF ( FreezeLowPopGaussians ) THEN
      CALL EquationsOfMotionSetup( EquationType, SimplifyCoefficients, (/ PopMinThreshold, PopMaxThreshold /) )
   ELSE
      CALL EquationsOfMotionSetup( EquationType, SimplifyCoefficients )
   END IF
         
   ! Copy the initial wavefunction to the actual one
   PsiT = Psi0

   ! Define initial value of time, energy, norm, autocorrelation
   ActualTime = 0.0
   Energy0 = REAL( WFExpectation( Psi0, Hamiltonian ) )
   Norm0 = WFNorm( Psi0 )
   Autocorrelation = WFOverlap( Psi0, Psi0 )
   ! and write them to standard output
   WRITE(*,"(/,A,A,/)") " ****************************************************  ",   &
                                &   "STARTING THE PROPAGATION  ******************************************************** "
   WRITE(*,600)  ActualTime/MyConsts_fs2AU, Norm0, REAL(Energy0)*MyConsts_Hartree2eV, 0.0, 0.0

   ! Open expectation file and write first few lines and the values at initial time
   ExpectUnit = LookForFreeUnit()
   OPEN(UNIT=ExpectUnit, FILE="expectations.dat")
   WRITE(ExpectUnit,*) "# Expectation values "
   WRITE(ExpectUnit,"('#',A8,100(A20))") " Time ", "Norm", "Energy", "RE(Auto)", "IMAG(Auto)", "ABS(Auto)"
   WRITE(ExpectUnit,700)  ActualTime/MyConsts_fs2AU, Norm0, REAL(Energy0)*MyConsts_Hartree2eV, REAL(Autocorrelation),  &
                        AIMAG(Autocorrelation), ABS(Autocorrelation)

   ! Open output units in which the centers and width of the gaussian configurations will be written, write initial lines
   ALLOCATE( GauCentersUnits(gDim), GauWidthUnits(gDim) )
   DO i = 1, gDim
      GauCentersUnits(i) = LookForFreeUnit()
      WRITE(String,'(A,I0.3,A)') "centers_pq_",i,".dat"
      OPEN(UNIT=GauCentersUnits(i), FILE=TRIM(ADJUSTL(String)))
      WRITE(GauCentersUnits(i),*) "# p,q of the gaussian centers along coord ", i
      WRITE(GauCentersUnits(i),"('#',A8,100(A20))") " Time ", "Coord", "Momentum"
      WRITE(GauCentersUnits(i),700)  ActualTime/MyConsts_fs2AU, WFGauCenters( Psi0, i )
      GauWidthUnits(i) = LookForFreeUnit()
      WRITE(String,'(A,I0.3,A)') "width_",i,".dat"
      OPEN(UNIT=GauWidthUnits(i), FILE=TRIM(ADJUSTL(String)))
      WRITE(GauWidthUnits(i),*) "# delta_q of the gaussian centers along coord ", i
      WRITE(GauWidthUnits(i),"('#',A8,100(A20))") " Time ", "DeltaQ"
      WRITE(GauWidthUnits(i),700)  ActualTime/MyConsts_fs2AU, WFGauWidth( Psi0, i )
   END DO

   PopUnit = LookForFreeUnit()
   OPEN(UNIT=PopUnit, FILE="populations.dat")
   WRITE(PopUnit,*) "# gaussian mulliken populations "
   WRITE(PopUnit,"('#',A8,100(A20))") " Time ", "Populations"
   WRITE(PopUnit,700)  ActualTime/MyConsts_fs2AU, WFGauPopulations(Psi0)   

   ! ****************************************************************************
   ! Propagation wavefunction for NPrintSteps steps of time PrintStep, and
   ! every PrintStep write to output and to file the relevant information
   ! ****************************************************************************

   ! loop over the time propagation steps which are done between output printing
   DO iPrintStep = 1, NPrintSteps

      ! Integrate from time ActualTime to ActualTime+PrintStep
      IntegratorStat = DoPropagationRungeKutta45( PsiT, Hamiltonian, ActualTime, PrintStep )

      ! When IntegratorStat is not zero, propagation has been ended for small time step
      IF  (IntegratorStat /= 0 ) THEN
         WRITE(*,"(/,A,/)") " Run termination due to small integration steps! "
         EXIT
      END IF

      ! Update expectation values
      Energy = REAL( WFExpectation( PsiT, Hamiltonian ) )
      Norm = WFNorm( PsiT )
      Autocorrelation = WFOverlap( Psi0, PsiT )

      ! Write expectation values to screen and to file
      WRITE(*,600)  ActualTime/MyConsts_fs2AU, Norm, REAL(Energy)*MyConsts_Hartree2eV, &
           (Norm-Norm0)*1000., REAL(Energy-Energy0)*1000.*MyConsts_Hartree2eV
      WRITE(ExpectUnit,700)  ActualTime/MyConsts_fs2AU, Norm, REAL(Energy)*MyConsts_Hartree2eV, REAL(Autocorrelation), &
                           AIMAG(Autocorrelation), ABS(Autocorrelation)
      DO i = 1, gDim
         WRITE(GauCentersUnits(i),700)  ActualTime/MyConsts_fs2AU, WFGauCenters( PsiT, i )
         WRITE(GauWidthUnits(i),700)  ActualTime/MyConsts_fs2AU, WFGauWidth( PsiT, i )
      END DO
      WRITE(PopUnit,700)  ActualTime/MyConsts_fs2AU, WFGauPopulations(PsiT )   

   END DO

   ! Close output unit
   CLOSE(UNIT=ExpectUnit)
   DO i = 1, gDim
      CLOSE(UNIT=GauCentersUnits(i) )
   END DO
   DEALLOCATE( GauCentersUnits )
   CLOSE(UNIT=PopUnit)

   ! Deallocate memory for matrix inversion module
   CALL MatrixInversionDispose(  )
   ! Close units and deallocate memory for RK45
   CALL DisposeRungeKutta45( ActualTime )
   ! Deallocate memory for the hamiltonian operator object
   CALL DisposeOperator( Hamiltonian )
   ! Deallocate memory for the wavefunction objects
   CALL DisposePsi( PsiT )
   CALL DisposePsi( Psi0 )

   ! Stop main timer and log all the timings
   CALL StopTimer(TotalClock)

   ! Print to standard output a table with the timings recorded
   CALL LogTimings

#if defined(LOG_FILE)
   __END_LOG_FILE
#endif

600 FORMAT(" Time = ",F9.3,"   Norm = ",F14.8,"   Energy = ",F14.8,"   DeltaNorm/1000 = ",F14.8,"   DeltaEnergy/1000 = ",F14.8 )
700 FORMAT(F9.3,100(F20.12))


END PROGRAM
