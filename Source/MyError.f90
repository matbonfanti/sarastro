!***************************************************************************************
!*                              MyError
!***************************************************************************************
!
!>  \brief     Warning, error and stop subroutines
!>  \details   Subroutines which show warning/error messages and stop the
!>             execution of the codes. These instructions are collected here
!>             to homogenize the appearence of the messages and to
!>             easily control the stop commands (which may be different in
!>             case of a non-standard execution - e.g. OMP and MPI).
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             April 2016
!>
!***************************************************************************************
MODULE MyError

   ! OS exit status corresponding to the type of error
   INTEGER, PARAMETER :: ERR_GENERIC      = 1,   & ! unspecified error
                         ERR_FILE_MISSING = 10,  & ! file does not exist
                         ERR_IO_OPEN      = 11,  & ! error when opening i/o file
                         ERR_INP_READ     = 12,  & ! error when reading input file
                         ERR_WRONG_INP    = 13,  & ! inconsistent input/wrong input value
                         ERR_MEM_ALLOC    = 20,  & ! memory allocation error
                         ERR_OBJ_MISUSE   = 30,  & ! incorrect use of object
                         ERR_MODULE_SETUP = 31,  & ! incorrect module setup
                         ERR_SUBS_INPUT   = 32,  & ! incorrect input value in function/subroutine
                         ERR_INT_ARITHM   = 40,  & ! integer arithmetic error
                         ERR_BLAS_LAPACK  = 41     ! blas/lapack error

   CONTAINS

!***************************************************************************************
!* Displays a warning message to stdout...
!***************************************************************************************
   SUBROUTINE ShowWarning( Message )
      IMPLICIT NONE
      CHARACTER(*) :: Message

      PRINT *
      PRINT '(A,A)', 'WARNING: ', TRIM(Message)
      PRINT *

   END SUBROUTINE ShowWarning

!***************************************************************************************
!* Test the test given and if true the warning is displayed...
!***************************************************************************************
   SUBROUTINE WARN( Test, Message )
      IMPLICIT NONE
      LOGICAL       :: Test
      CHARACTER(*)  :: Message

      IF( Test ) CALL ShowWarning(Message)

   END SUBROUTINE WARN

!***************************************************************************************
!*  Stop the execution of the program, with given error message and optional error status
!***************************************************************************************
  SUBROUTINE AbortWithError( ErrMsg, InpOSExitStatus, ErrCode )
      CHARACTER(LEN=*), INTENT(IN)   :: ErrMsg
      INTEGER, INTENT(IN), OPTIONAL  :: InpOSExitStatus
      INTEGER, INTENT(IN), OPTIONAL  :: ErrCode
      INTEGER :: OSExitStatus

      IF (.NOT. PRESENT(InpOSExitStatus)) THEN
         OSExitStatus = ERR_GENERIC
      ELSE
         OSExitStatus = InpOSExitStatus
      END IF

      SELECT CASE( OSExitStatus )
         CASE ( ERR_GENERIC )
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_GENERIC

         CASE ( ERR_FILE_MISSING )
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_FILE_MISSING

         CASE ( ERR_IO_OPEN )
            IF ( PRESENT(ErrCode) ) THEN
               WRITE(*,"(/,1X,A,A,I4)") trim(ErrMsg), " IOSTAT: ",ErrCode
            ELSE
               WRITE(*,"(/,1X,A)") trim(ErrMsg)
            END IF
            STOP ERR_IO_OPEN

         CASE ( ERR_INP_READ )
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_INP_READ

         CASE ( ERR_WRONG_INP )
            IF ( PRESENT(ErrCode) ) THEN
               WRITE(*,"(/,1X,A,A,I4)") trim(ErrMsg), " line ",ErrCode
            ELSE
               WRITE(*,"(/,1X,A)") trim(ErrMsg)
            END IF
            STOP ERR_WRONG_INP

         CASE ( ERR_MEM_ALLOC )
            IF ( PRESENT(ErrCode) ) THEN
               WRITE(*,"(/,1X,A,A,I4)") trim(ErrMsg), " Errorcode: ",ErrCode
            ELSE
               WRITE(*,"(/,1X,A)") trim(ErrMsg)
            END IF
            STOP ERR_MEM_ALLOC

         CASE ( ERR_OBJ_MISUSE )
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_OBJ_MISUSE

         CASE ( ERR_MODULE_SETUP )
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_MODULE_SETUP

         CASE ( ERR_SUBS_INPUT )
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_SUBS_INPUT

         CASE ( ERR_INT_ARITHM )
            IF ( PRESENT(ErrCode) ) THEN
               WRITE(*,"(/,1X,A,A,I4)") trim(ErrMsg), ErrCode
            ELSE
               WRITE(*,"(/,1X,A)") trim(ErrMsg)
            END IF
            STOP ERR_INT_ARITHM

         CASE ( ERR_BLAS_LAPACK )
            IF ( PRESENT(ErrCode) ) THEN
               WRITE(*,"(/,1X,A,A,I4)") trim(ErrMsg), " INFO: ", ErrCode
            ELSE
               WRITE(*,"(/,1X,A)") trim(ErrMsg)
            END IF
            STOP ERR_BLAS_LAPACK

         CASE DEFAULT
            WRITE(*,"(/,1X,A)") trim(ErrMsg)
            STOP ERR_GENERIC

      END SELECT

  END SUBROUTINE AbortWithError

!***************************************************************************************
!* Test the test given and if true, error message is displayed and program aborted...
!***************************************************************************************
   SUBROUTINE ERROR( Test, ErrMsg, InpOSExitStatus, ErrCode )
      IMPLICIT NONE
      LOGICAL                       :: Test
      CHARACTER(LEN=*), INTENT(IN)  :: ErrMsg
      INTEGER, INTENT(IN), OPTIONAL :: InpOSExitStatus
      INTEGER, INTENT(IN), OPTIONAL :: ErrCode
      INTEGER :: OSExitStatus

      IF (.NOT. PRESENT(InpOSExitStatus)) THEN
         OSExitStatus = ERR_GENERIC
      ELSE
         OSExitStatus = InpOSExitStatus
      END IF

      IF( Test ) CALL AbortWithError( ErrMsg, OSExitStatus, ErrCode )

   END SUBROUTINE ERROR

!***************************************************************************************

END MODULE MyError
