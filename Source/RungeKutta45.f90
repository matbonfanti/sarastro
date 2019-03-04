!***************************************************************************************
!*                              MODULE RungeKutta45
!***************************************************************************************
!
!>  \brief     Time propagation of the wavefunction with Runge Kutta 4(5).
!>  \details   This module  implements the integration of the equations \n
!>             of motion with the Runge Kutta method at the 4th order,  \n
!>             adaptive time step chosen by estimating the error with   \n
!>             difference between 4th and 5th order.                    \n
!>             Two sets of numerical coefficients are available: the    \n
!>             ones of the Runge-Kutta-Fehlberg Integrator and the ones \n
!>             of the Cash-Karp Integrator.                             \n
!>             The module needs to be setup at the beginning with a     \n
!>             call to SetupRungeKutta45( ... ) which defines the time  \n
!>             step conditions and the error tolerance.                 \n
!>             Then the propagation can be done with the subroutine     \n
!>             DoPropagationRungeKutta45( ... ), while at the end       \n
!>             the module execution can be concluded with a call to     \n
!>             the subroutine DisposeRungeKutta45( ) which              \n
!>             deallocate memory and close open i/o units.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             25 July 2017
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          ______________________________________
!
!***************************************************************************************
MODULE RungeKutta45
#include "preprocessoptions.cpp"
   USE PsiObject
   USE OperatorDefine
   USE AdaptGBasis
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupRungeKutta45             !< setup module prior to the RK45 integration
   PUBLIC :: DoPropagationRungeKutta45     !< do a series of RK45 propagation steps for a fixed total time
   PUBLIC :: DisposeRungeKutta45           !< deallocate memory and close open i/o units

   ! Format of numbers in integrator.log file
   CHARACTER(100), PARAMETER :: LogNumbersFormat  = '(I10,1X,I7,1X,A6,2X,E23.16,2X,E23.16,2X,E23.16)'
   CHARACTER(100), PARAMETER :: LogNumbersFormat2 = '(E23.16,2X,E23.16)'
   CHARACTER(100), PARAMETER :: AddNumbersFormat  = '(2X,E23.16)'

   LOGICAL, SAVE :: PropagationIsSetup = .FALSE.      !< Logical flag to check whether the module is ready or not

   ! Type of coefficients of the RK45 integrator
   INTEGER, PUBLIC, PARAMETER :: RK45_FEHLBERG  = 0,  &   !  Runge–Kutta–Fehlberg method
                                 RK45_CASH_KARP = 1       !  Cash–Karp method
   INTEGER, SAVE              :: RK45_Type      = RK45_CASH_KARP

   ! Logical flag to keep the integration step fixed and switch off its adaptation
   LOGICAL :: KeepStepFixed = .FALSE.

   ! Numerical coefficients of the Runge-Kutta 4-5 Method
   REAL, SAVE :: c20, c2(1)
   REAL, SAVE :: c30, c3(2)
   REAL, SAVE :: c40, c4(3)
   REAL, SAVE :: c50, c5(4)
   REAL, SAVE :: c60, c6(5)
   REAL, SAVE :: a(6), b(6)
   INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: aIndex, bIndex

   ! Absolute number of RK45 steps which are used in a complete run of the code
   INTEGER, SAVE :: TotalNrRK45
   ! Absolute number of RK45 steps which are rejected in a complete run of the code
   INTEGER, SAVE :: TotalRejectNrRK45

   ! Actual stepsize and stepsize limits
   REAL, SAVE :: tStep
   REAL, SAVE :: tStepMin
   REAL, SAVE :: tStepMax

   ! Error tolerance for adaptive timestep
   REAL, SAVE :: Tolerance

   ! unit for integration log file
   INTEGER, SAVE :: LogUnit
   ! unit for time step analysis
   INTEGER, SAVE :: TStepUnit


   CONTAINS

!===========================================================================================================
!                               PUBLIC SUBROUTINES AND FUNCTIONS
!===========================================================================================================

!*******************************************************************************
!          SetupRungeKutta45
!*******************************************************************************
!> Setup module prior to the RK45 integration .
!>
!> @param FixedTime      Logical flag to set fixed time integration
!> @param MinStep        Min value of the integration time step
!> @param MaxStep        Max value of the integration time step
!>                       (if FixedTime this is the actual time step)
!> @param InpTolerance   Error tolerance of RK45 for time step adaptation
!*******************************************************************************
   SUBROUTINE SetupRungeKutta45( FixedTime, MinStep, MaxStep, InpTolerance )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: FixedTime
      REAL, INTENT(IN) :: MinStep, MaxStep, InpTolerance

      ! In case propagation has been already setup, give error
      CALL ERROR( PropagationIsSetup, " SetupRungeKutta45: propagation has been already set up", ERR_MODULE_SETUP )

      ! Initialize the total number of RK45 steps
      TotalNrRK45 = 0; TotalRejectNrRK45 = 0

      ! DEFINE TIME STEPS AND INTEGRATION THRESHOLDS
      IF ( FixedTime ) THEN
         KeepStepFixed = .TRUE.
         tStepMin = MaxStep
         tStepMax = MaxStep
         Tolerance  = 1.E+99
      ELSE IF ( .NOT. FixedTime ) THEN
         KeepStepFixed = .FALSE.
         tStepMin = MinStep
         tStepMax = MaxStep
         Tolerance  = InpTolerance
      END IF
      ! At the beginning of the integration, choose maximum step
      tStep = tStepMax

      ! Define an IO unit for printig log info on integration steps
      LogUnit = LookForFreeUnit()
      ! Open unit and write header of the file
      OPEN( FILE="integrator.log", UNIT=LogUnit )
      SELECT CASE( RK45_Type )
         CASE ( RK45_FEHLBERG )
            WRITE(LogUnit,*) "# Runge-Kutta-Fehlberg Integrator, 4th Order (with 5th Order error estimate)"
         CASE( RK45_CASH_KARP )
            WRITE(LogUnit,*) "# Cash-Karp Integrator, 4th Order (with 5th Order error estimate)"
      END SELECT
      WRITE(LogUnit,'("#",A10,1X,A7,1X,6X,2X,A23,2X,A23,2X,A23)',ADVANCE="no") "step", "attempt", "stepsize", "error", "time"
      WRITE(LogUnit,'(2X,A23)',ADVANCE="no") "overlap_inv-condnr"
      WRITE(LogUnit,'(2X,A23)',ADVANCE="no") "cmatrix_inv-condnr"
      WRITE(LogUnit,'(2X,A23)',ADVANCE="no") "overlap_max"
      WRITE(LogUnit,'()')

      ! Define an IO unit for printig log info on accepted steps
      TStepUnit = LookForFreeUnit()
      ! Open unit and write header of the file
      OPEN( FILE="accepted_steps.log", UNIT=TStepUnit )
      WRITE(TStepUnit,'("#",A23,2X,A23)',ADVANCE="no") "stepsize", "error"
      WRITE(TStepUnit,'(2X,A23)',ADVANCE="no") "overlap_inv-condnr"
      WRITE(TStepUnit,'(2X,A23)',ADVANCE="no") "cmatrix_inv-condnr"
      WRITE(TStepUnit,'()')

      ! Define coefficients of RungeKutta45 integrator
      SELECT CASE( RK45_Type )

      CASE ( RK45_FEHLBERG )
         c20 = 1.0/4.0;          c2 = (/ 1.0/4.0 /)
         c30 = 3.0/8.0;          c3 = (/ 3./32., 9./32. /)
         c40 = 12./13.;          c4 = (/ 1932./2197., -7200./2197., 7296./2197. /)
         c50 = 1.0;              c5 = (/ 439./216., -8., 3680./513., -845./4104. /)
         c60 = 0.5;              c6 = (/ -8./27., 2., -3544./2565., 1859./4104., -11./40. /)
         a = (/ 25./216., 0.,   1408./2565.,   2197./4104.,  -1./5., 0.     /)                   ! FOURTH ORDER
         b = (/ 16./135., 0.,  6656./12825., 28561./56430., -9./50., 2./55. /)                   ! FIFTH ORDER
         ALLOCATE(aIndex(4)); aIndex = (/ 1, 3, 4, 5 /)
         ALLOCATE(bIndex(5)); bIndex = (/ 1, 3, 4, 5, 6 /)

      CASE( RK45_CASH_KARP )
         c20 = 1./5. ;   c2 = (/ 1./5. /)
         c30 = 3./10.;   c3 = (/ 3./40. ,   9./40. /)
         c40 = 3./5. ;   c4 = (/ 3./10. ,  -9./10.,    6./5. /)
         c50 = 1.0   ;   c5 = (/ -11./54.,    5./2., -70./27.,  35./27. /)
         c60 = 7./8. ;   c6 = (/ 1631./55296.,  175./512.,  575./13824.,  44275./110592.,    253./4096. /)
         a = (/ 2825./27648., 0., 18575./48384., 13525./55296., 277./14336.,      1./4. /)   ! FOURTH ORDER
         b = (/ 37./378.    , 0.,     250./621.,     125./594.,          0., 512./1771. /)   ! FIFTH ORDER
         ALLOCATE(aIndex(5)); aIndex = (/ 1, 3, 4, 5, 6 /)
         ALLOCATE(bIndex(4)); bIndex = (/ 1, 3, 4, 6 /)

      END SELECT

      ! Propagation is now well defined
      PropagationIsSetup = .TRUE.

   END SUBROUTINE SetupRungeKutta45


!*******************************************************************************
!          DoPropagationRungeKutta45
!*******************************************************************************
!> Perform a series of Runge-Kutta integration steps to propagate the
!> wavefunction up to a fixed total time. The integration is done with
!> a Runge-Kutta 4th-Order method, with time adaptation based on error
!> estimate which is calculated as the difference to 5th-Order.
!>
!> @param   Psi            Wavefunction to propagate
!> @param   Hamiltonian    Hamiltonian of the propagation
!> @param   ActualTime     Actual time of the WF, updated during integration
!> @param   PropTime       Total amount of time propagation of the subroutine execution
!> @returns ExitStatus     Exit status definying the success or failure of the integration
!*******************************************************************************
   INTEGER FUNCTION DoPropagationRungeKutta45( Psi, Hamiltonian, ActualTime, PropTime )  RESULT( ExitStatus )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)               :: Psi               ! wavefunction at the beginning of the step
      TYPE(OperatorData), INTENT(IN)               :: Hamiltonian       ! Hamiltonian of the propagation
      REAL, INTENT(INOUT)                          :: ActualTime        ! starting time of the step, on exit final time
      REAL, INTENT(IN)                             :: PropTime          ! propagation step

      ! Temporary wavefunction to compute intermediate approximation and then 4th and 5th order solutions
      TYPE(WaveFunct) :: Psi4th, Psi5th
      ! Temporary memory to store intermediate derivatives giving the RK integral
      TYPE(Derivative), DIMENSION(6) :: Deriv
      ! Estimated error of integration
      REAL :: ErrorEst
      ! Estimated inverse condition number of the matrices which are inverted during propagation
      REAL, DIMENSION(2) :: InverseCondition

      REAL :: EndTime         ! ending time of the propagation step
      REAL :: NewStepSize     ! updated value of propagation step
      INTEGER :: nAttempts, i
      REAL    :: MaxOverlap

      ! In case propagation has been not set up yet, give error
      CALL ERROR( .NOT. PropagationIsSetup, " DoOneStepRungeKutta45: propagation has not been set up yet", ERR_MODULE_SETUP )

      ! Set time of propagation step
      EndTime = ActualTime + PropTime

      ! Initialize exit status of the function ( 0 = execution is ok )
      ExitStatus = 0

      TimePropagation: DO WHILE ( ActualTime < EndTime ) ! propagate for predefined time

         ! call the subroutine for the adaptation of the basis set during the propagation
!          CALL AddOrRemoveGaussianFunctions( Psi, Hamiltonian, ActualTime, LogUnit )

         ! call the subroutine to fix approximate propagation of some gaussian configurations depending on populations
         CALL FixGaussianConfigurationsWithSmallPopulations( Psi )

         ! Initialize temporary wavefunction and allocate memory by copying the initial wavefunction
         Psi4th = Psi; Psi5th = Psi

         ! Set initial time
         CALL StartTimer(RungeKuttaClock)

         ! evaluate derivatives at t1
         CALL ComputeDerivative( Psi, Deriv(1), Hamiltonian, ActualTime, InverseCondition )
         ! evaluate maximum of the abs non diagonal elements of the ovelap matrix
         MaxOverlap = GetOverlapMax( Psi )

         nAttempts = 0     ! Initialize the number of integration attempts
         DO           ! loop until the error on integration is accettable
            nAttempts = nAttempts + 1      ! increment the number of attempts that has been done so far

            ! estimate wavefunction at t2 and then compute derivative at t2
            CALL PsiPlusDeltaPsi( Psi4th, Psi, Deriv(1:1), tStep*c2(1:1) )
            CALL ComputeDerivative( Psi4th, Deriv(2), Hamiltonian, ActualTime+c20*tStep )

            ! estimate wavefunction at t3 and then compute derivative at t3
            CALL PsiPlusDeltaPsi( Psi4th, Psi, Deriv(1:2), tStep*c3(1:2) )
            CALL ComputeDerivative( Psi4th, Deriv(3), Hamiltonian, ActualTime+c30*tStep )

            ! estimate wavefunction at t4 and then compute derivative at t4
            CALL PsiPlusDeltaPsi( Psi4th, Psi, Deriv(1:3), tStep*c4(1:3) )
            CALL ComputeDerivative( Psi4th, Deriv(4), Hamiltonian, ActualTime+c40*tStep )

            ! estimate wavefunction at t5 and then compute derivative at t5
            CALL PsiPlusDeltaPsi( Psi4th, Psi, Deriv(1:4), tStep*c5(1:4) )
            CALL ComputeDerivative( Psi4th, Deriv(5), Hamiltonian, ActualTime+c50*tStep )

            ! estimate wavefunction at t6 and then compute derivative at t6
            CALL PsiPlusDeltaPsi( Psi4th, Psi, Deriv(1:5), tStep*c6(1:5) )
            CALL ComputeDerivative( Psi4th, Deriv(6), Hamiltonian, ActualTime+c60*tStep )

            !############### FOURTH-ORDER ESTIMATE #################
            CALL PsiPlusDeltaPsi( Psi4th, Psi, Deriv(aIndex), tStep*a(aIndex) )

            !################ FIFTH-ORDER ESTIMATE #################
            CALL PsiPlusDeltaPsi( Psi5th, Psi, Deriv(bIndex), tStep*b(bIndex) )

            !################ ERROR ESTIMATE #######################
            ErrorEst = WFDifference( Psi4th, Psi5th, WFDISTANCE_L2_NORM )

            IF ( .NOT. KeepStepFixed ) THEN
               ! Compute the new updated value of the time step
               NewStepSize = AdaptStepSize( tStep, Tolerance, ErrorEst )
               ! Check whether the new steps is higher or lower than the thresholds and act accordingly
               IF ( NewStepSize <= tStepMin ) THEN
                  ExitStatus = 1           ! ExitStatus = 1 means that the propagation is terminated because of too low timestep
                  CALL DisposePsi( Psi4th );  CALL DisposePsi( Psi5th )          ! clean up and stop running clock
                  CALL StopTimer(RungeKuttaClock)
                  EXIT TimePropagation
               ELSE IF ( NewStepSize > tStepMax ) THEN
                  NewStepSize = tStepMax
               END IF
            END IF

            ! Accept step when the error is small or when a calculation with fixed timestep is performed
            IF ( ErrorEst < Tolerance .OR. KeepStepFixed ) THEN

               ! The steps can be accepted, increment the total number of time steps
               TotalNrRK45 = TotalNrRK45 + 1
               ! Time has changed from current time to current time + timestep
               ActualTime = ActualTime + tStep
               ! Store updated Psi integrated with 4th order RK
               Psi = Psi4th

               ! Print info on accepted step to log file
               WRITE(LogUnit,LogNumbersFormat,ADVANCE="no") TotalNrRK45, nAttempts, "accept", &
                                           tStep/MyConsts_fs2AU, ErrorEst, ActualTime/MyConsts_fs2AU
               WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(1)
               WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(2)
               WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") MaxOverlap
               WRITE(LogUnit,'()')

               ! Print info on accepted step to log file
               WRITE(TStepUnit,LogNumbersFormat2,ADVANCE="no") tStep/MyConsts_fs2AU, ErrorEst
               WRITE(TStepUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(1)
               WRITE(TStepUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(2)
               WRITE(TStepUnit,'()')

               ! Adapt the time step with respect to the remaining time left for the integration period
               ! and store the new time step in tStep for the following time step
               IF ( .NOT. KeepStepFixed )  tStep = ConfrontStepWithRemainingTime( NewStepSize, EndTime-ActualTime )

               EXIT     ! Exit from inner loop and go on with following RK45 step

            ELSE

               ! Increment the counter of the rejected time steps
               TotalRejectNrRK45 = TotalRejectNrRK45 + 1
               ! Print info on failed step attempt to log file
               WRITE(LogUnit,LogNumbersFormat,ADVANCE="no") TotalNrRK45, nAttempts, "reject", &
                                             tStep/MyConsts_fs2AU, ErrorEst, ActualTime/MyConsts_fs2AU
               WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(1)
               WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(2)
               WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") MaxOverlap
               WRITE(LogUnit,'()')

               ! Update the time step (there is no need to confront with remaining time, since we are trying again the same step
               IF ( .NOT. KeepStepFixed )  tStep = NewStepSize
               CYCLE                           ! Continue within the inner DO loop (RK45 attempts)

            ENDIF

         END DO     ! end loop over the attempts to perform a single RK45 step

         ! Set final time
         CALL StopTimer(RungeKuttaClock)

         ! deallocate memory for the temporary storage of the wavefunction and the derivatives
         CALL DisposePsi( Psi4th );  CALL DisposePsi( Psi5th )
         DO i = 1, 6
            CALL DisposeDerivative(Deriv(i))
         END DO

!          CALL FLUSH()

      END DO TimePropagation ! end loop over the RK45 steps that give the time propagation for the desired time


   END FUNCTION DoPropagationRungeKutta45


!*******************************************************************************
!          DisposeRungeKutta45
!*******************************************************************************
!> Deallocate memory and close open i/o units.
!*******************************************************************************
   SUBROUTINE DisposeRungeKutta45( ActualTime  )
      IMPLICIT NONE
      REAL, INTENT(IN)    :: ActualTime        ! final time at which the integrator is disposed

      ! DEALLOCATE MEMORY
      DEALLOCATE( aIndex )
      DEALLOCATE( bIndex )

      ! CLOSE LOG UNIT
      CLOSE( UNIT=LogUnit )
      CLOSE( UNIT=TStepUnit )

      WRITE(*,200) TotalRejectNrRK45+TotalNrRK45, TotalNrRK45, TotalRejectNrRK45, &
                   ActualTime/REAL(TotalNrRK45)/MyConsts_fs2AU, &
                   ActualTime/REAL(TotalRejectNrRK45+TotalNrRK45)/MyConsts_fs2AU

      ! Propagation is now undefined
      PropagationIsSetup = .FALSE.

   200 FORMAT  (/,' ****************** Integration statistics  *****************',/, &
                  " Tot nr of attempted time steps        = ",             I20   ,/, &
                  " ** accepted                           = ",             I20   ,/, &
                  " ** rejected                           = ",             I20   ,/, &
                  " Average accepted time step /fs        = ",             F20.8 ,/, &
                  " Average overall tme step /fs          = ",             F20.8 )

   END SUBROUTINE DisposeRungeKutta45


!===========================================================================================================
!                               PRIVATE SUBROUTINES AND FUNCTIONS
!===========================================================================================================


!*******************************************************************************
!          AdaptStepSize
!*******************************************************************************
!> Adapt time step of Runge-Kutta 45 integration according to
!> the value of a function Delta of the relative error Error/Tolerance.
!>
!> @param      StepSize       Actual step size of the RK45 integration|
!> @param      Tolerance      Error threshold of RK45 integration
!> @param      Error          Error of the last integration step
!> @returns    NewStepSize    Updated value of the step size
!*******************************************************************************
   FUNCTION AdaptStepSize( StepSize, Tolerance, Error ) RESULT( NewStepSize )
      IMPLICIT NONE
      REAL, INTENT(IN) :: StepSize
      REAL, INTENT(IN) :: Tolerance, Error
      REAL             :: NewStepSize
      REAL, PARAMETER  :: Beta = 0.84
      REAL             :: Delta

      ! Define factor Delta which drives the adaptation of the timestep
      IF ( Error < Tolerance .OR. RK45_Type == RK45_FEHLBERG ) THEN
         Delta = Beta * (Tolerance/MAX(Error,MyConsts_EPS))**0.25
      ELSE
         Delta = Beta * (Tolerance/MAX(Error,MyConsts_EPS))**0.20
      ENDIF

      ! Adapt stepsize
      IF ( Delta < 0.1 ) THEN
         NewStepSize = 0.1*ABS(StepSize)
      ELSE IF ( Delta > 4.0 ) THEN
         NewStepSize = 4.0*ABS(StepSize)
      ELSE
         NewStepSize = Delta*ABS(StepSize)
      ENDIF

      ! If needed, give back to the stepsize its sign
      IF ( StepSize < 0.0 ) NewStepSize = -ABS(NewStepSize)

  END FUNCTION AdaptStepSize


!*******************************************************************************
!          ConfrontStepWithRemainingTime
!*******************************************************************************
!> Delimits stepsize of actual step to match the end of the DoPropagation
!> subroutines, in such a way that the execution terminates at the desired
!> final time.
!>
!> @param      StepSize       Actual step size of the RK45 integration
!> @param      TimeLeft       Remaining time to the end of propagation step
!> @returns    NewStepSize    Updated value of the step size
!*******************************************************************************
   FUNCTION ConfrontStepWithRemainingTime( StepSize, TimeLeft )  RESULT( NewStepSize )
      IMPLICIT NONE
      REAL, INTENT(IN) :: StepSize
      REAL, INTENT(IN) :: TimeLeft
      REAL             :: NewStepSize
      REAL, PARAMETER  :: StepMargin = 0.05

      IF ( TimeLeft < 2.0*MyConsts_EPS ) THEN    ! remaining time is lower than numerical precision
         NewStepSize = StepSize
      ELSE IF ( StepSize <  TimeLeft ) THEN      ! remaining time is larger than the time step
         IF ( ( TimeLeft - StepSize ) <= ( StepMargin * StepSize) ) THEN
            NewStepSize = TimeLeft                           ! slightly increase timestep to get to final time
         ELSE IF ( ( TimeLeft - StepSize ) <= ( 0.3*StepSize ) ) THEN
            NewStepSize = TimeLeft * 0.5                     ! decrease timestep to balance the last two time steps
         ELSE
            NewStepSize = StepSize                           !  we are far from final time, keep time step as it is
         END IF
      ELSE              ! remaining time is smaller than the time step
         NewStepSize = TimeLeft                  !  decrease timestep to match final time
      END IF

  END FUNCTION ConfrontStepWithRemainingTime


END MODULE RungeKutta45
