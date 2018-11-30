!***************************************************************************************
!*                              MODULE Clock
!***************************************************************************************
!
!>  \brief     CPU and WALL timing of code execution
!>  \details   This module store the WALL and CPU timings and compute
!>             the time spent in a given section of the code.
!>             The "Timer" object defined by this module stores the wall time and
!>             cpu time spent between the two start and stop methods ich can be used
!>             in the code with starting and stopping methods, in analogy with
!>             a stopwatch.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             April 2016
!>
!***************************************************************************************
!
!>   \remark         The module is based on the intrinsic subroutines  \n
!>                   CPU_TIME and DATE_AND_TIME and assume a serial    \n
!>                   execution (parallel extension is needed)
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg -- : --
!
!***************************************************************************************

MODULE Clock
   USE MyError

   PRIVATE
   PUBLIC :: Timer
   PUBLIC :: InternalTimersSetup, LogTimings
   PUBLIC :: StartTimer, StopTimer
   PUBLIC :: GetTotalTime, GetSectionFrequency, GetAverageTime

!********************************************************************************************************

   ! The following data defines timers internally defined by the module when
   ! the subroutine InternalTimersSetup is called
   ! they can be used outside with the public integer parameters below

   ! Nr of timers used in the code
   INTEGER, PARAMETER :: NrTimers = 6
   ! Labels for the timers
   CHARACTER(14), DIMENSION(NrTimers)  ::  TimersLabels = (/ "main          ", &
                                                             "runge-kutta45 ", &
                                                             "derivatives   ", &
                                                             "der_gausspar  ", &
                                                             "der_bvector   ", &
                                                             "mat_inversion " /)
   ! Timers identifier
   INTEGER, PARAMETER, PUBLIC :: TotalClock = 1
   INTEGER, PARAMETER, PUBLIC :: RungeKuttaClock = 2
   INTEGER, PARAMETER, PUBLIC :: DerivativesClock = 3
   INTEGER, PARAMETER, PUBLIC :: GaussDerivClock = 4
   INTEGER, PARAMETER, PUBLIC :: BVecDerivClock = 5
   INTEGER, PARAMETER, PUBLIC :: MatrixInversionClock = 6

!********************************************************************************************************

   !> Data type to define the Timer class. Wall time and cpu time are in seconds
   TYPE Timer
      PRIVATE
      LOGICAL :: StopwatchRunning = .FALSE.       ! Logical flag: timer is started or stopped
      REAL    :: TotalCPUTime = 0.0               ! Real variable to store CPU time (seconds)
      REAL    :: TotalWALLTime = 0.0              ! Real variable to store WALL time (seconds)
      INTEGER :: NrStart = 0                      ! Nr. of times the timer is started
      REAL    :: StartCPUTime                     ! Temporary variable to store start for CPU time
      INTEGER, DIMENSION(8) :: StartWALLTime      ! Temporary variable to store start for WALL time
      CHARACTER(20) :: Label = "none"             ! string definying what time is stored by the timer object
   END TYPE Timer

   ! Timers defined internally in the module
   TYPE(Timer), DIMENSION(:), ALLOCATABLE :: InternalTimers

   ! Logical variable to check the setup status of the module
   ! Module is ready to use when the initialisation subroutine has been called
   LOGICAL :: ModuleIsSetup = .FALSE.

   INTERFACE StartTimer
      MODULE PROCEDURE StartTimer_INT, StartTimer_OBJ
   END INTERFACE

   INTERFACE StopTimer
      MODULE PROCEDURE StopTimer_INT, StopTimer_OBJ
   END INTERFACE

   
!********************************************************************************************************
   CONTAINS
!********************************************************************************************************


!*******************************************************************************
!          InternalTimersSetup
!*******************************************************************************
!> Setup all the necessary timers which are defined by the parameters
!> of the module NrTimers and TimersLabels. These timers are supposed to
!> be designated outside the module by integer parameters defined
!> accordingly with the indices of TimersLabels.
!>
!*******************************************************************************
   SUBROUTINE InternalTimersSetup(  )
      IMPLICIT NONE
      INTEGER :: i

      ! give error if module has been already set up
      CALL ERROR( ModuleIsSetup, " InternalTimersSetup: module is already initialized ", ERR_MODULE_SETUP )

      ! Allocate internal array
      ALLOCATE( InternalTimers(NrTimers)  )

      ! Store the names of the types of inverted matrices
      DO i = 1, NrTimers
         InternalTimers(i)%Label = TRIM(ADJUSTL( TimersLabels(i) ))
      END DO

      ! Now the module is ready to be used
      ModuleIsSetup = .TRUE.

   END SUBROUTINE InternalTimersSetup


!*******************************************************************************
!          StartTimer
!*******************************************************************************
!> Start timer at the current line of code. The routine can be called
!> directly with an input object, or with an integer which refers to
!> one of the allocated internal timers.
!>
!> @param TimerObject     Timer data type / Integer index of the internal array of timers
!*******************************************************************************
   SUBROUTINE StartTimer_OBJ( TimerObject )
      IMPLICIT NONE
      TYPE(Timer), INTENT(INOUT) :: TimerObject

      CALL ERROR( TimerObject%StopwatchRunning, " Clock.StartTimer: timer is already started" )

      ! Store initial wall and cpu timings at present
      CALL CPU_TIME( TimerObject%StartCPUTime )
      CALL DATE_AND_TIME( VALUES=TimerObject%StartWALLTime )

      ! Start timer
      TimerObject%StopwatchRunning = .TRUE.

   END SUBROUTINE StartTimer_OBJ

   SUBROUTINE StartTimer_INT( TimerObject )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: TimerObject

      CALL ERROR( .NOT. ModuleIsSetup, " StartTimer: module has not been initialized yet ", ERR_MODULE_SETUP )
      CALL ERROR( InternalTimers(TimerObject)%StopwatchRunning, " Clock.StartTimer: timer is already started" )

      ! Store initial wall and cpu timings at present
      CALL CPU_TIME( InternalTimers(TimerObject)%StartCPUTime )
      CALL DATE_AND_TIME( VALUES=InternalTimers(TimerObject)%StartWALLTime )

      ! Start timer
      InternalTimers(TimerObject)%StopwatchRunning = .TRUE.

   END SUBROUTINE StartTimer_INT


!*******************************************************************************
!          StopTimer
!*******************************************************************************
!> Stop timer at the current line of code. The routine can be called
!> directly with an input object, or with an integer which refers to
!> one of the allocated internal timers.
!>
!> @param TimerObject     Timer data type / Integer index of the internal array of timers
!*******************************************************************************
      SUBROUTINE StopTimer_OBJ( TimerObject )
         IMPLICIT NONE
         TYPE(Timer), INTENT(INOUT) :: TimerObject
         REAL    :: StopCPUTime
         INTEGER, DIMENSION(8) :: StopWALLTime

         CALL ERROR( .NOT. TimerObject%StopwatchRunning, " Clock.StopTimer: timer has not been started" )

         ! Get final wall and cpu timings at present
         CALL CPU_TIME( StopCPUTime )
         CALL DATE_AND_TIME( VALUES=StopWALLTime )

         ! Stop timer
         TimerObject%StopwatchRunning = .FALSE.

         ! Compute time difference and add to the accumulated timings
         TimerObject%TotalCPUTime = TimerObject%TotalCPUTime + StopCPUTime - TimerObject%StartCPUTime
         TimerObject%TotalWALLTime = TimerObject%TotalWALLTime + TimeDifference(TimerObject%StartWALLTime,StopWALLTime)

         ! Increment the number of section recorded
         TimerObject%NrStart = TimerObject%NrStart + 1

      END SUBROUTINE StopTimer_OBJ

      SUBROUTINE StopTimer_INT( TimerObject )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: TimerObject
         REAL    :: StopCPUTime
         INTEGER, DIMENSION(8) :: StopWALLTime

         CALL ERROR( .NOT. ModuleIsSetup, " StartTimer: module has not been initialized yet ", ERR_MODULE_SETUP )
         CALL ERROR( .NOT. InternalTimers(TimerObject)%StopwatchRunning, " Clock.StopTimer: timer has not been started" )

         ! Get final wall and cpu timings at present
         CALL CPU_TIME( StopCPUTime )
         CALL DATE_AND_TIME( VALUES=StopWALLTime )

         ! Stop timer
         InternalTimers(TimerObject)%StopwatchRunning = .FALSE.

         ! Compute time difference and add to the accumulated timings
         InternalTimers(TimerObject)%TotalCPUTime = InternalTimers(TimerObject)%TotalCPUTime + &
                                                    StopCPUTime - InternalTimers(TimerObject)%StartCPUTime
         InternalTimers(TimerObject)%TotalWALLTime = InternalTimers(TimerObject)%TotalWALLTime + &
                                                   TimeDifference(InternalTimers(TimerObject)%StartWALLTime,StopWALLTime)

         ! Increment the number of section recorded
         InternalTimers(TimerObject)%NrStart = InternalTimers(TimerObject)%NrStart + 1

      END SUBROUTINE StopTimer_INT

!*******************************************************************************
!          LogTimings()
!*******************************************************************************
!> Write to standard output a scheme of the time recorded by the internal
!> timers.
!>
!*******************************************************************************
   SUBROUTINE LogTimings(  )
      IMPLICIT NONE
      INTEGER :: i

      ! header of the table
      WRITE(*,500)

      ! Write one line per each of the internal timers
      DO i = 1, SIZE(InternalTimers)
         WRITE(*,501) TRIM(ADJUSTL(InternalTimers(i)%Label)), GetTotalTime(InternalTimers(i)), &
                      GetSectionFrequency(InternalTimers(i)), GetAverageTime(InternalTimers(i))
      END DO

      ! end of the table
      WRITE(*,502)

500 FORMAT( /,"==================================================================================================",/, &
              " Code section  | tot WallTime/ s | tot CPUTime/ s  | Nr calls | avg WallTime/ s | avg CPUTime/ s  ",/, &
              "--------------------------------------------------------------------------------------------------" )
501 FORMAT(    A14,         " | ", F15.2,     " | ", F15.2,     " | ",I8,  " | ", F15.2,     " | ", F15.2          )
502 FORMAT(   "==================================================================================================",/ )

   END SUBROUTINE LogTimings


!*******************************************************************************
!          GetTotalTime
!*******************************************************************************
!> Returns the CPU and WALL time currently accumulated in the input object.
!>
!> @param TimerObject     Timer data type
!*******************************************************************************
   FUNCTION GetTotalTime( TimerObject )
      IMPLICIT NONE
      TYPE(Timer), INTENT(INOUT) :: TimerObject
      REAL, DIMENSION(2) :: GetTotalTime

      CALL ERROR( TimerObject%NrStart == 0, " Clock.GetTotalTime: timer has never been used" )

      ! Return the accumulated wall and cpu time
      GetTotalTime(1) = TimerObject%TotalCPUTime
      GetTotalTime(2) = TimerObject%TotalWALLTime

   END FUNCTION GetTotalTime


!*******************************************************************************
!          GetSectionFrequency
!*******************************************************************************
!> Returns the number of times the timer of the input object has been started.
!>
!> @param TimerObject     Timer data type
!*******************************************************************************
   INTEGER FUNCTION GetSectionFrequency( TimerObject )
      IMPLICIT NONE
      TYPE(Timer), INTENT(INOUT) :: TimerObject

      CALL ERROR( TimerObject%NrStart == 0, " Clock.GetSectionFrequency: timer has never been used" )

      ! Return the accumulated wall and cpu time
      GetSectionFrequency = TimerObject%NrStart

   END FUNCTION GetSectionFrequency


!*******************************************************************************
!          GetAverageTime
!*******************************************************************************
!> Returns the CPU and WALL average time per number of timer start
!> from the accumulated total time of the timer.
!>
!> @param TimerObject     Timer data type
!*******************************************************************************
   FUNCTION GetAverageTime( TimerObject )
      IMPLICIT NONE
      TYPE(Timer), INTENT(INOUT) :: TimerObject
      REAL, DIMENSION(2) :: GetAverageTime

      CALL ERROR( TimerObject%NrStart == 0, " Clock.GetAverageTime: timer has never been used" )

      ! Return the accumulated wall and cpu time
      GetAverageTime(1) = TimerObject%TotalCPUTime / REAL(TimerObject%NrStart)
      GetAverageTime(2) = TimerObject%TotalWALLTime / REAL(TimerObject%NrStart)

   END FUNCTION GetAverageTime


!*******************************************************************************
!     TimeDifference
!*******************************************************************************
!> Compute difference between final time and initial time formatted
!> according to the DATE_AND_TIME intrinsic function.
!>
!> @param InitialTime     Int array (dim 8) with the initial time
!> @param FinalTime       Int array (dim 8) with the final time
!> @returns               Real variable with seconds between initial and final time
!*******************************************************************************
   REAL FUNCTION TimeDifference( InitialTime, FinalTime )
      INTEGER, DIMENSION(8), INTENT(IN)    :: InitialTime, FinalTime

      TimeDifference =  REAL(FinalTime(3)-InitialTime(3))*24.0*60.0*60.0 + &
                        REAL(FinalTime(5)-InitialTime(5))*60.0*60.0 + &
                        REAL(FinalTime(6)-InitialTime(6))*60.0 + &
                        REAL(FinalTime(7)-InitialTime(7)) + &
                        REAL(FinalTime(8)-InitialTime(8))/1000.

   END FUNCTION TimeDifference



END MODULE Clock

