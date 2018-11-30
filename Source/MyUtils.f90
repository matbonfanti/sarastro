!***************************************************************************************
!*                              MODULE MyUtils
!***************************************************************************************
!
!>  \brief     Generic utilities subroutines
!>  \details   This module defines commonly used subroutines for I/O, sorting,
!>             type conversion, etc etc ...
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             July 2017
!>
!***************************************************************************************
MODULE MyUtils
   USE MyError

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: LookForFreeUnit         ! returns the smallest i>=20 avaiable as I/O unit
   PUBLIC :: NumberToString          ! convert to string an integer number
   PUBLIC :: CountLinesInFile        ! counts the number of lines of a file
   PUBLIC :: RemoveDups              ! returns the non repeated elements of an array
   PUBLIC :: Sort                    ! sort an array in ascending order
   PUBLIC :: Order                   ! returns the ordinal numbers of a given input array

   !> Wrapper for the unique elements extraction subroutine
   INTERFACE RemoveDups
      MODULE PROCEDURE RemoveDups1D_r, RemoveDups2D_r, RemoveDups1D_i
   END INTERFACE

   !> Wrapper for the sort subroutines
   INTERFACE Sort
      MODULE PROCEDURE Sort_r, Sort_i
   END INTERFACE

   !> Wrapper for the subroutine to convert number to string
   INTERFACE NumberToString
      MODULE PROCEDURE NumberToString_nonpadded, NumberToString_padded
   END INTERFACE

   CONTAINS


!*******************************************************************************
!                           LookForFreeUnit
!*******************************************************************************
!>  Set the smallest integer equal or greater than 20 that is
!>  available as unit i/o.
!>
!> @returns           The integer number of the smallest free unit.
!*******************************************************************************
   INTEGER FUNCTION LookForFreeUnit()
   IMPLICIT NONE
      LOGICAL               :: File_Opened
      INTEGER, PARAMETER    :: UnitMax = 300

      DO LookForFreeUnit = 20, UnitMax   ! Look for a non opened I/O unit
         INQUIRE( UNIT=LookForFreeUnit, OPENED=File_Opened )
         IF (.NOT. File_Opened) EXIT
      END DO
      ! If an available unit has been found, use as function results
      CALL ERROR((LookForFreeUnit == UnitMax), &
                                      "PrintTools: No free I/O unit available")
   END FUNCTION LookForFreeUnit


!*******************************************************************************
!                           NumberToString
!*******************************************************************************
!>  Given an integer number, five back a character string of the numbers
!>  with the input integer npadded, optionally pad the number with
!>  zeros giving a string of npadded digits.
!>
!> @param    N        Input integer number.
!> @param    npadded  Optional number of digits of the number padded with zeros.
!> @returns           A string with the number.
!*******************************************************************************
   FUNCTION NumberToString_nonpadded( N )  RESULT( String )
   IMPLICIT NONE
      INTEGER, INTENT(IN)        :: N
      CHARACTER(12) :: String

      WRITE(String,*) N
      String = ADJUSTL(String)
   END FUNCTION NumberToString_nonpadded

   FUNCTION NumberToString_padded( N, npadded ) RESULT( String )
   IMPLICIT NONE
      INTEGER, INTENT(IN)   :: N
      INTEGER, INTENT(IN)   :: npadded
      CHARACTER(npadded)    :: String
      CHARACTER(10) :: Fmt

      Fmt = '(I0.'//TRIM(ADJUSTL(NumberToString_nonpadded(N)))//')'
      WRITE(String,Fmt) N
   END FUNCTION NumberToString_padded


!*******************************************************************************
!                           CountLinesInFile
!*******************************************************************************
!>  Count the number of lines in a file with a given filename.
!>
!> @param    FileName     Input string with the file name.
!> @returns               The integer number of rows of file FileName.
!*******************************************************************************
   INTEGER FUNCTION CountLinesInFile( FileName )
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: FileName
      INTEGER             :: Stat, UnitNr
      CHARACTER(1)        :: String

      UnitNr=LookForFreeUnit()
      OPEN( FILE=TRIM(ADJUSTL(FileName)), UNIT=UnitNr )
      DO CountLinesInFile = 0,100000
         READ( UnitNr, *, IOSTAT=Stat ) String
         IF (Stat /= 0) EXIT
      END DO
      CLOSE( UnitNr )

   END FUNCTION CountLinesInFile


!*******************************************************************************
!                              RemoveDups
!*******************************************************************************
!>  Give the non repeated elements of an array . The input array is
!>  destroyed, instead the sub gives back an array with the non repeated entries
!>  in the elements 1:NrOfNonRepeated
!>
!> @param   InputArray        In input array
!> @param   NrOfNonRepeated   The nr of unique elements
!*******************************************************************************

SUBROUTINE RemoveDups1D_r( InputArray, NrOfNonRepeated )
   IMPLICIT NONE
   REAL, DIMENSION(:), INTENT(INOUT)  ::   InputArray
   INTEGER, INTENT(OUT)               ::   NrOfNonRepeated

   REAL, DIMENSION(SIZE(InputArray))  ::   TmpArray
   INTEGER :: i, j

   NrOfNonRepeated = 1
   TmpArray(1) = InputArray(1)
   Outer: DO i = 2, Size(InputArray)
      DO j = 1, NrOfNonRepeated
         IF ( TmpArray(j) == InputArray(i) )   CYCLE Outer
      END DO
      ! No match is found, so new element is added to non repeated array
      NrOfNonRepeated = NrOfNonRepeated + 1
      TmpArray( NrOfNonRepeated ) = InputArray(i)
   END DO Outer

   ! Store in InputArray the non repeated elements
   CALL Sort( TmpArray(1:NrOfNonRepeated) )
   InputArray(1:NrOfNonRepeated) = TmpArray(1:NrOfNonRepeated)

END SUBROUTINE RemoveDups1D_r

SUBROUTINE RemoveDups1D_i( InputArray, NrOfNonRepeated )
   IMPLICIT NONE
   INTEGER, DIMENSION(:), INTENT(INOUT)  ::   InputArray
   INTEGER, INTENT(OUT)                  ::   NrOfNonRepeated

   INTEGER, DIMENSION(SIZE(InputArray))  ::   TmpArray
   INTEGER :: i, j

   NrOfNonRepeated = 1
   TmpArray(1) = InputArray(1)
   Outer: DO i = 2, Size(InputArray)
      DO j = 1, NrOfNonRepeated
         IF ( TmpArray(j) == InputArray(i) )   CYCLE Outer
      END DO
      ! No match is found, so new element is added to non repeated array
      NrOfNonRepeated = NrOfNonRepeated + 1
      TmpArray( NrOfNonRepeated ) = InputArray(i)
   END DO Outer

   ! Store in InputArray the non repeated elements
   CALL Sort( TmpArray(1:NrOfNonRepeated) )
   InputArray(1:NrOfNonRepeated) = TmpArray(1:NrOfNonRepeated)

END SUBROUTINE RemoveDups1D_i

SUBROUTINE RemoveDups2D_r( InputArray, NrOfNonRepeated )
   IMPLICIT NONE
   REAL, DIMENSION(:,:), INTENT(INOUT)  ::   InputArray
   INTEGER, INTENT(OUT)                 ::   NrOfNonRepeated

   REAL, DIMENSION(SIZE(InputArray,1),SIZE(InputArray,2))  ::   TmpArray
   INTEGER :: i, j

   NrOfNonRepeated = 1
   TmpArray(1,:) = InputArray(1,:)
   Outer: DO i = 2, Size(InputArray,1)
      DO j = 1, NrOfNonRepeated
         IF ( ALL( TmpArray(j,:) == InputArray(i,:) ) )   CYCLE Outer
      END DO
      ! No match is found, so new element is added to non repeated array
      NrOfNonRepeated = NrOfNonRepeated + 1
      TmpArray( NrOfNonRepeated,: ) = InputArray(i,:)
   END DO Outer

   ! Store in InputArray the non repeated elements
   InputArray(1:NrOfNonRepeated,:) = TmpArray(1:NrOfNonRepeated,:)

END SUBROUTINE RemoveDups2D_r


!*******************************************************************************
!                              Sort
!*******************************************************************************
!> This subroutine receives an array x() and sorts it into ascending order.
!>
!> @param     x        In input real/integer array to sort
!> @returns            The sorted array
!*******************************************************************************

   SUBROUTINE  Sort_r(x)
      IMPLICIT  NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: x
      INTEGER    :: i,  Location
      REAL       :: Temp
      LOGICAL, DIMENSION(:), ALLOCATABLE :: Mask

      ALLOCATE(Mask(SIZE(x)))
      Mask = .TRUE.
      DO i = 1, SIZE(x)-1                  ! except for the last
         Location = MINLOC( x, 1, Mask )    ! find min from this to last
         Temp = x(i)
         x(i) = x(Location)
         x(Location) = Temp
         Mask(i) = .FALSE.
      END DO
      DEALLOCATE(Mask)
   END SUBROUTINE  Sort_r

   SUBROUTINE  Sort_i(x)
      IMPLICIT  NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: x
      INTEGER    :: i,  Location
      INTEGER    :: Temp
      LOGICAL, DIMENSION(:), ALLOCATABLE :: Mask

      ALLOCATE(Mask(SIZE(x)))
      Mask = .TRUE.
      DO i = 1, SIZE(x)-1                  ! except for the last
         Location = MINLOC( x, 1, Mask )    ! find min from this to last
         Temp = x(i)
         x(i) = x(Location)
         x(Location) = Temp
         Mask(i) = .FALSE.
      END DO
      DEALLOCATE(Mask)
   END SUBROUTINE  Sort_i


!*******************************************************************************
!                              Order
!*******************************************************************************
!>  This subroutine receives an array x() and gives an integer array
!>  of the same dimension with the ordinal numer of the array element
!>
!> @param     x   In input real array
!> @returns       An integer array with the ordinal numer of the array elements
!*******************************************************************************

   FUNCTION Order(x)
      IMPLICIT  NONE
      REAL, DIMENSION(:), INTENT(INOUT)        :: x
      INTEGER, DIMENSION(size(x))              :: Order
      LOGICAL, DIMENSION(size(x)) :: Mask
      INTEGER    :: i,  Location

      Mask = .TRUE.
      DO i = 1, SIZE(x)                     ! except for the last
         Location = MINLOC( x, 1, Mask )    ! find min from this to last
         Order(Location) = i
         Mask(Location)  = .FALSE.
      END DO
   END FUNCTION  Order

END MODULE MyUtils

! ! --------------------------------------------------------------------
! ! INTEGER FUNCTION  FindMinimum():
! !    This function returns the location of the minimum in the section
! ! between Start and End.
! ! --------------------------------------------------------------------
!
!    INTEGER FUNCTION  FindMinimum(x, Start, End)
!       IMPLICIT  NONE
!       REAL, DIMENSION(1:), INTENT(IN) :: x
!       INTEGER, INTENT(IN)             :: Start, End
!       INTEGER                         :: Minimum
!       INTEGER                         :: Location
!       INTEGER                         :: i
!
!       Minimum  = x(Start)               ! assume the first is the min
!       Location = Start                  ! record its position
!       DO i = Start+1, End               ! start with next elements
!          IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
!             Minimum  = x(i)             !      Yes, a new minimum found
!             Location = i                !      record its position
!          END IF
!       END DO
!       FindMinimum = Location            ! return the position
!    END FUNCTION  FindMinimum

! ! --------------------------------------------------------------------
! ! SUBROUTINE  Swap():
! !    This subroutine swaps the values of its two formal arguments.
! ! --------------------------------------------------------------------
!
!    SUBROUTINE  Swap(a, b)
!       IMPLICIT  NONE
!       REAL, INTENT(INOUT) :: a, b
!
!       Temp = a
!       a    = b
!       b    = Temp
!    END SUBROUTINE  Swap


!********************************************* END OF FILE *******************************
