!****************************************************************************************************
!*                              MODULE OperatorDefine
!****************************************************************************************************
!
!>  \brief     Define and compute an operator.
!>  \details   This module set up the necessary data to define an operator (arbitrary defined      \n
!>             operator from a file or a special direct-dynamics suited ones) written as usual     \n
!>             as a sum over product. The format of this file has been inherited from the          \n
!>             vMCG code by Pierre Eisenbradt.                                                     \n
!>                                                                                                 \n
!>             The storage of the sum-over-product expansion follows the following conventions:    \n
!>                * Operator%N                  is the number of sum over product terms            \n
!>             then, for  i = 1, ... , Operator%N                                                  \n
!>                * Operator%Coeffs(i)    is the linear coefficient of the i-th product operator   \n
!>                * Operator%OpType(:,i)  is a triple of integer specifying the type of operator:  \n
!>                    |--->  Operator%OpType(1,i) is the general form of the operator              \n
!>                    |           (-)  0 ... potential using only predefined functions             \n
!>                    |           (-)  1 ... kinetic type 1 (cartesian)                            \n
!>                    |           (-)  2 ... kinetic type 2 (not implemented)                      \n
!>                    |           (-) -1 ... potential using user defined function (or mixed)      \n
!>                    |---> Operator%OpType(2,i) is the left-hand-side electronic state            \n
!>                    \---> Operator%OpType(3,i) is the right-hand-side electronic state           \n
!>                              if they are equal the Operator element acts only on one state      \n
!>                              if not equal, this gives a coupling between two states!            \n
!>                * Operator%OpPar(:,i)   is a gDim vector with specific param.s of the operator   \n
!>                          in case of Operator%OpType(1,i) == 0, the vector is the list of powers \n
!>                          and the product operator is a power product of the coordinates         \n
!>                          i.e.  q_1 ^ Operator%OpPar(1,i) x q_2 ^ Operator%OpPar(2,i) x ...      \n
!>                          in case of Operator%OpType(1,i) == 1, the non-zero element of the vect \n
!>                          indicates on which coordinate is the cartesian operatro acting.        \n
!>             The variable Operator%LocalApproxOrder controls the enforcing of the local          \n
!>             approximation of the potential defined in Operator. When it is 0, the full sum-over-\n
!>             product potential is computed, otherwise for 2 or 3, all the potetnial terms are    \n
!>             substituted with their taylor expansion up to the corresponding order.              \n
!
!***************************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             17 July 2017
!>
!****************************************************************************************************
!
!>  \par Updates
!>  \arg 26 October 2017  : added
!
!>  \todo          add definition of masses for kinetic energy terms
!>  \todo          add local expansion of the potential defined analytically
!>  \todo          add local expansion of the potential computed on-the-fly
!
!*****************************************************************************************************
MODULE OperatorDefine
#include "preprocessoptions.cpp"
   USE GauConf
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: OperatorData                                            ! New data type to store the Operator definitions
   PUBLIC :: SetOperatorFromFile, SetOperatorKinOnly                 ! Constructors of OperatorData, from file, with kin only
   PUBLIC :: SetDerivPotentialOperator, SetHessianPotentialOperator
   PUBLIC :: LogOperator                                             ! Print to file infos on the setup Operator
   PUBLIC :: DisposeOperator                                         ! Destructor of a OperatorData instance
   PUBLIC :: OperatorMatrixElement                                   ! Compute matrix element of an operator with gaussians functions
   PUBLIC :: EvaluatePotential, EvaluateGradient, EvaluateHessian    ! Evaluate V/gradient/hessian of Potential contained in given operator at given point

   !> Data type to store the information on the Operator operator.                     \n
   TYPE OperatorData
      PRIVATE
      INTEGER                              :: gDim        !< Dimensionality of the space on which the operator is defined
      INTEGER                              :: N           !< Number of terms of the sum-over-product expansion
      REAL, DIMENSION(:), ALLOCATABLE      :: Coeffs      !< Coefficients of the sum-over-product expansion
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: OpType      !< Definition of the operators, type of the functions
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: OpPar       !< Definition of the operators, parameters of the functions
      INTEGER  :: LocalApproxOrder                        !< Order of the local approx of the potential defined in the datatype
   END TYPE OperatorData

   ! Definitions of the various type of operators implemented
   INTEGER, PUBLIC, PARAMETER :: OPTYPE_POTENTIAL = 0,  &   !  potential term, monomial function of the coordinate
                                 OPTYPE_KINETIC   = 1,  &   !  kinetic term in cartesian coordinate
                                 OPTYPE_UNDEFINED = 99      !  undefined potential term


!===========================================================================================================
   CONTAINS
!===========================================================================================================


!*******************************************************************************
!          SetOperatorFromFile
!*******************************************************************************
!> OperatorData constructor method: set the instance Operator of object
!> OperatorData reading data from input file. In case of an already setup
!> Operator, data is overwritten.
!>
!> @param Operator                   OperatorData object
!> @param gDim             Number of dimensions of the system
!> @param OperFile        Filename of the file with the A vector coefficients
! ! ! !> @param InpLocalApproxOrder   Integer local approximation order
!*******************************************************************************
   SUBROUTINE SetOperatorFromFile( Operator, gDim, nStates, OperFile )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(INOUT) :: Operator
      INTEGER, INTENT(IN)            :: gDim, nStates
      CHARACTER(*), INTENT(IN)       :: OperFile
!       INTEGER, INTENT(IN)            :: InpLocalApproxOrder

      INTEGER, PARAMETER   :: MaxLineNumber = 10000
      CHARACTER, PARAMETER :: Sep = ','
      LOGICAL              :: Success
      INTEGER              :: OperInp, iStat
      INTEGER              :: TermNumber, PreviousTermNumber, iLine, iGDim
      CHARACTER(LEN=2000)  :: Line, Substring, ErrorString

      ! In case memory is already allocated, deallocate (printing warning to be safe)
      IF ( ALLOCATED(Operator%Coeffs) .OR. ALLOCATED(Operator%OpType) .OR. ALLOCATED( Operator%OpPar ) ) THEN
         CALL ShowWarning( " OperatorDefine.SetOperatorFromFile: disposing previosly allocated for target OperatorData object" )
         CALL DisposeOperator( Operator )
      END IF

!       ! Check the value (no approx, LHA or LCA - 0,2,3 respectively) and store the local approximation order for the potential
!       CALL ERROR(  InpLocalApproxOrder /= 0 .AND. InpLocalApproxOrder /= 2 .AND. InpLocalApproxOrder /= 3,  &
!          " OperatorDefine.SetOperatorFromFile: incorrect input order of local approximation ", ERR_OBJ_MISUSE )
!       Operator%LocalApproxOrder = InpLocalApproxOrder

      ! Store the number of dimensions
      Operator%gDim = gDim

      ! Check the existence of the Operator input file, if it is there open it
      INQUIRE( FILE=TRIM(ADJUSTL(OperFile)), EXIST=Success )
      CALL ERROR( .NOT. Success, 'Operator input file -'//TRIM(ADJUSTL(OperFile))//'- not found', ERR_FILE_MISSING)
      OperInp = LookForFreeUnit()
      OPEN( UNIT=OperInp, FILE=TRIM(ADJUSTL(OperFile)), STATUS='old', ACTION='READ' )

      ! Count the line of the input file starting with consecutive integer numbers
      ! lines repeating the previous number will be ignored
      PreviousTermNumber = 0
      DO iLine = 1, MaxLineNumber

            ! Read line and break cycle when no more lines are available
            READ(OperInp,'(A2000)', IOSTAT=iStat) Line
            IF (iStat /= 0) EXIT

            ! Get first comma-separated value in the line, if read failed, stop program execution with error
            CALL get_next_substr(Line, Substring, Sep, Success)
            WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
                                 '  --> '//TRIM(ADJUSTL(Line))//NEW_LINE(Line), "  Wrong Operator input file at"
            CALL ERROR( .NOT. Success, ErrorString, ERR_WRONG_INP, iLine )

            ! the first field of the line is the operator term number
            READ(Substring,*) TermNumber
            IF ( TermNumber == PreviousTermNumber ) THEN   ! skip the line when the number is equal to the one of the previous line
               CYCLE
            ELSE IF ( TermNumber == PreviousTermNumber + 1 ) THEN   ! the line is ok, it contains the subsequent integer number
               Operator%N = TermNumber
               PreviousTermNumber = TermNumber
            ELSE         ! any other number is wrong
               WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
                        '  --> '//TRIM(Substring)//", "//trim(Line)//NEW_LINE(Line), " Non consecutive ascending term numbers at"
               CALL AbortWithError( ErrorString, ERR_WRONG_INP, iLine )
            END IF

      END DO
      ! Prepare input file for subsequent input variable parsing
      REWIND(OperInp)

      ! Allocate arrays Operator%Coeffs, Operator%OpType and Operator%OpPar
      ALLOCATE( Operator%Coeffs(Operator%N)     )
      ALLOCATE( Operator%OpType(3,Operator%N)   )
      ALLOCATE( Operator%OpPar(Operator%gDim,Operator%N) )

      ! Initialize all the arrays
      Operator%Coeffs(:)   = 0.0
      Operator%OpType(1,:) = OPTYPE_UNDEFINED
      Operator%OpType(2,:) = 1
      Operator%OpType(3,:) = 1
      Operator%OpPar(:,:)  = 0

      ! =========================================================================================
      ! Now read the Operator from File!                                                     !
      ! file contains comma separated values, one term per each line defined by the sequence:   !
      !        i, Operator%Coeffs(i), Operator%OpType(1...3,i), Operator%OpPar(1...gDim,i)      !
      ! all input data is integer except for Operator%Coeffs(i)                                 !
      ! =========================================================================================

      ! Read file line by line, ignoring line which repeats the previous number in the first field
      PreviousTermNumber = 0
      DO iLine = 1, MaxLineNumber

         ! Read line and break cycle when no more lines are available
         READ(OperInp,'(A2000)', IOSTAT=iStat) Line
         IF (iStat /= 0) EXIT

         ! Get first comma-separated value in the line (this step has been already done, so no need to check the success)
         CALL get_next_substr(Line, Substring, Sep, Success)
         READ(Substring,*) TermNumber

         ! Skip the line if the number is equal to the previous one
         IF ( TermNumber == PreviousTermNumber ) CYCLE
         ! Because of the first check, we are sure that the actual TermNumber is correct, we can update PreviousTermNumber
         PreviousTermNumber = TermNumber

         !read the Operator%Coeffs(i)
         CALL get_next_substr(Line, Substring, Sep, Success)
         IF ( Success ) THEN
            READ(Substring,*) Operator%Coeffs(TermNumber)
         ELSE
            WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
               "  --> primitive op nr.: ", TermNumber, NEW_LINE(Line), "  Incapable of reading Operator%Coeffs at"
            CALL AbortWithError( ErrorString, ERR_WRONG_INP, iLine )
         END IF

         ! When the calculation is not single state, read the electronic states corresponding to the operator element
         ! otherwise Operator%OpType(2,i) and Operator%OpType(3,i) are left with the initialization value (equal to 1)
         IF ( nStates > 1 ) THEN

            CALL get_next_substr(Line, Substring, Sep, Success)
            IF ( Success ) THEN
               READ(Substring,*) Operator%OpType(2,TermNumber)
            ELSE
               WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
               "  --> primitive op nr.: ", TermNumber, NEW_LINE(Line), "  Incapable of reading Operator%OpType(2,i) at"
               CALL AbortWithError( ErrorString, ERR_WRONG_INP, iLine )
            END IF
            CALL get_next_substr(Line, Substring, Sep, Success)
            IF ( Success ) THEN
               READ(Substring,*) Operator%OpType(3,TermNumber)
            ELSE
               WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
               "  --> primitive op nr.: ", TermNumber, NEW_LINE(Line), "  Incapable of reading Operator%OpType(3,i) at"
               CALL AbortWithError( ErrorString, ERR_WRONG_INP, iLine )
            END IF

         END IF

         ! Read the type of operator and save it to Operator%OpType(1,i)
         CALL get_next_substr(Line, Substring, Sep, Success)
         IF ( Success ) THEN
            READ(Substring,*) Operator%OpType(1,TermNumber)
         ELSE
            WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
                "  --> primitive op nr.: ", TermNumber, NEW_LINE(Line), "  Incapable of reading Operator%OpType(1,i) at"
            CALL AbortWithError( ErrorString, ERR_WRONG_INP, iLine )
         END IF

         ! Cycle over the number of degree of freedom to get the components of the operator on the primitive single-dof ones
         ! the cycle stops at gdim-1, because the last one is read from the remaining string after the last comma

         ! The error message here is the same for all the primitive operators
         WRITE(ErrorString,*) 'Error reading Operator input file '//TRIM(ADJUSTL(OperFile))//NEW_LINE(Line), &
                "  --> primitive op nr.: ", TermNumber, NEW_LINE(Line), "  Not enough primitive-dimension operators given at"

         DO iGDim = 1, Operator%gDim-1
            ! Read the portion of the line up to the first comma and give error message in case of failure
            CALL get_next_substr(Line, Substring, Sep, Success)
            CALL ERROR( .NOT. Success, ErrorString, ERR_WRONG_INP, iLine )
            ! Save the value to op_primdim
            READ(Substring,*) Operator%OpPar(iGDim,TermNumber)
         END DO
         ! The last portion of the line contains the last primitive operator
         READ(line,*) Operator%OpPar(Operator%gDim,TermNumber)

      END DO

      ! Close input file, everything has been done
      CLOSE( UNIT=OperInp )

   END SUBROUTINE SetOperatorFromFile

   SUBROUTINE get_next_substr(string,substring,separator,success)
      ! finds the next occurrance of the separator in the string
      ! returns as substring everything left of the separator
      ! returns a string everything right of the separator
      ! returns the string as it was if the separator was not found or if
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: separator
      CHARACTER(LEN=*), INTENT(INOUT) :: string
      CHARACTER(LEN=*), INTENT(OUT) :: substring
      LOGICAL, INTENT(OUT) :: success
      INTEGER :: sindex

      success=.FALSE.

      sindex=INDEX(string,separator)
      IF(sindex.GT.1)THEN
         success=.TRUE.
         substring=trim(adjustl(string(1:sindex-1)))
         string=trim(adjustl(string(sindex+1:)))
      END IF
   END SUBROUTINE get_next_substr



!*******************************************************************************
!          SetOperatorKinOnly
!*******************************************************************************
!> OperatorData constructor method: set the instance Operator of object
!> OperatorData including a standard kinetic energy part and
!> no potential. Usefull in case of a direct dynamics run.
!> In case of an already setup Operator, data is overwritten.
!>
!> @param Operator                   OperatorData object
!*******************************************************************************
   SUBROUTINE SetOperatorKinOnly( Operator, gDim )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(INOUT) :: Operator
      INTEGER, INTENT(IN)            :: gDim
      INTEGER              :: iGDim

      ! In case memory is already allocated, deallocate (printing warning to be safe)
      IF ( ALLOCATED(Operator%Coeffs) .OR. ALLOCATED(Operator%OpType) .OR. ALLOCATED( Operator%OpPar ) ) THEN
         CALL ShowWarning( " OperatorDefine.SetOperatorFromFile: disposing previosly allocated for target OperatorData object" )
         CALL DisposeOperator( Operator )
      END IF

      ! When only kinetic terms are defined, LocalApproxOrder is irrelevant
      Operator%LocalApproxOrder = 0

      ! Store the number of dimensions
      Operator%gDim = gDim
      ! The number of term is equal to the number of dimensions (one kinetic term per dimension)
      Operator%N = Operator%gDim

      ! Allocate arrays Operator%Coeffs, Operator%OpType and Operator%OpPar
      ALLOCATE( Operator%Coeffs(Operator%N)     )
      ALLOCATE( Operator%OpType(3,Operator%N)   )
      ALLOCATE( Operator%OpPar(Operator%gDim,Operator%N) )

      ! =================================================
      !            Dyrect Dynamics Hamiltonian          !
      ! kinetic energy operator in cartesian coordinate !
      ! =================================================

      ! Set Operator%Coeffs it is the set of linear coefficients of the potential terms, in this case is everywhere equal to one
      Operator%Coeffs(:) = 1

      ! Set Operator%OpType first number is the type of operator (kinetic/potential), second and third the electronic states
      ! all gdim elements are kinetic operator
      Operator%OpType(1,:) = 1
      ! for the time being, direct dynamics is single state only
      Operator%OpType(2,:) = 1
      Operator%OpType(3,:) = 1

      ! Set op_primdim : indicates on which coordinate is the the operator acting
      ! Initialize everything to zero
      Operator%OpPar = 0
      ! Kinetic energy is separable, with the i-th term acting on the i-th coordinate
      DO iGDim = 1, Operator%gDim
         Operator%OpPar(iGDim,iGDim) = 1
      END DO

   END SUBROUTINE SetOperatorKinOnly



!*******************************************************************************
!          SetDerivPotentialOperator
!*******************************************************************************
!> Defines the potential operator *DerivOperator* starting from a generic
!> operator *Operator* and by taking the derivative along the *iCoord*-th
!> coordinate of all the potential terms. In this way, the subroutine
!> defines the "force operator" corresponding to the potential contained
!> in Operator, which can be used to compute the classical Ehrenfest forces.
!>
!> @param DerivOperator   output gives the operator with the derivative of the potential
!> @param Operator        input, takes the OperatorData object which is differentiated
!> @param iCoord          input, takes the coordinate along which differentiate Operator
!*******************************************************************************
   SUBROUTINE SetDerivPotentialOperator( DerivOperator, Operator, iCoord )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(INOUT) :: DerivOperator
      TYPE(OperatorData), INTENT(IN)    :: Operator
      INTEGER, INTENT(IN)               :: iCoord

      INTEGER, DIMENSION(Operator%N) :: IndexMask
      INTEGER :: i, OpN, n

      ! In case memory is already allocated for DerivOperator, deallocate (printing warning to be safe)
      IF ( ALLOCATED(DerivOperator%Coeffs) .OR. ALLOCATED(DerivOperator%OpType) .OR. ALLOCATED( DerivOperator%OpPar ) ) THEN
         CALL ShowWarning( " OperatorDefine.SetDerivPotentialOperator: disposing previosly allocated for target OperatorData object (DerivOperator)")
         CALL DisposeOperator( DerivOperator )
      END IF
      ! In case memory of Operator is not allocated, exit with error
      CALL ERROR( .NOT. ALLOCATED(Operator%Coeffs) .OR. .NOT. ALLOCATED(Operator%OpType) .OR. .NOT. ALLOCATED( Operator%OpPar ),  &
         &        " OperatorDefine.SetDerivPotentialOperator: source OperatorData object is not defined (Operator)", ERR_OBJ_MISUSE )

      ! The number of dimension of the derivative operator is the same
      DerivOperator%gDim = Operator%gDim

      ! Count the number of potential energy operators contained in Operator in which Q_iCoord has at least power 1
      OpN = 0
      DO i = 1, Operator%N
        IF ( Operator%OpType(1,i) ==  OPTYPE_KINETIC .OR. Operator%OpPar(iCoord,i) == 0 ) CYCLE
        OpN = OpN+1; IndexMask(OpN) = i
      END DO
      ! OpN is now the number of potential terms to include in the derivative
      DerivOperator%N = OpN

      ! Now we can allocate memory to store the potential derivative operator
      ALLOCATE( DerivOperator%Coeffs(OpN)     )
      ALLOCATE( DerivOperator%OpType(3,OpN)   )
      ALLOCATE( DerivOperator%OpPar(DerivOperator%gDim ,OpN) )

      ! Define DerivOperator%OpType copying from Operator%OpType
      DerivOperator%Coeffs( 1:OpN ) = Operator%Coeffs( IndexMask(1:OpN) )
      DerivOperator%OpType( 1:3, 1:OpN ) = Operator%OpType( 1:3, IndexMask(1:OpN) )
      DerivOperator%OpPar( :, 1:OpN ) = Operator%OpPar( :, IndexMask(1:OpN) )

      DO n = 1, OpN
         i = IndexMask(n)
         DerivOperator%Coeffs(n) = Operator%Coeffs(i) * REAL(Operator%OpPar(iCoord , i))
         DerivOperator%OpPar(iCoord , n) = Operator%OpPar(iCoord , i)-1
      END DO

   END SUBROUTINE SetDerivPotentialOperator


!*******************************************************************************
!          SetHessianPotentialOperator
!*******************************************************************************
!> Defines the potential operator *HessianOperator* starting from a generic
!> operator *Operator* and by taking the second derivative along the *iCoord*-th
!> and the *jCoord*-th coordinates of all the potential terms. In this way, the subroutine
!> defines the "hessian operator" corresponding to the potential contained
!> in Operator.
!>
!> @param HessianOperator   output gives the operator with the derivative of the potential
!> @param Operator          input, takes the OperatorData object which is differentiated
!> @param iCoord            input, takes the first coordinate along which differentiate Operator
!> @param jCoord            input, takes the second coordinate along which differentiate Operator
!*******************************************************************************
   SUBROUTINE SetHessianPotentialOperator( HessianOperator, Operator, iCoord, jCoord )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(INOUT) :: HessianOperator
      TYPE(OperatorData), INTENT(IN)    :: Operator
      INTEGER, INTENT(IN)               :: iCoord, jCoord

      INTEGER, DIMENSION(Operator%N) :: IndexMask
      INTEGER :: i, OpN, n

      ! In case memory is already allocated for HessianOperator, deallocate (printing warning to be safe)
      IF ( ALLOCATED(HessianOperator%Coeffs) .OR. ALLOCATED(HessianOperator%OpType) .OR. ALLOCATED( HessianOperator%OpPar ) ) THEN
         CALL ShowWarning( " OperatorDefine.SetHessianPotentialOperator: disposing previosly allocated for target OperatorData object (HessianOperator)")
         CALL DisposeOperator( HessianOperator )
      END IF
      ! In case memory of Operator is not allocated, exit with error
      CALL ERROR( .NOT. ALLOCATED(Operator%Coeffs) .OR. .NOT. ALLOCATED(Operator%OpType) .OR. .NOT. ALLOCATED( Operator%OpPar ),  &
         &        " OperatorDefine.SetHessianPotentialOperator: source OperatorData object is not defined (Operator)", ERR_OBJ_MISUSE )

      ! The number of dimension of the derivative operator is the same
      HessianOperator%gDim = Operator%gDim

      ! Count the number of potential energy operators contained in Operator which are not zero upon differentiation
      OpN = 0
      DO i = 1, Operator%N
         IF ( Operator%OpType(1,i) ==  OPTYPE_KINETIC ) CYCLE    ! the term is kinetic
         IF ( iCoord /= jCoord ) THEN                                                        ! when i /= j
            IF ( Operator%OpPar(iCoord,i) == 0 .OR. Operator%OpPar(jCoord,i) == 0 ) CYCLE    ! remove terms which are constant along i or along j
         ELSE                                                                                ! when i == j
            IF ( Operator%OpPar(iCoord,i) == 0 .OR. Operator%OpPar(iCoord,i) == 1 ) CYCLE    ! remove terms which are linear or constant along i=j
         END IF
         OpN = OpN+1; IndexMask(OpN) = i
      END DO
      ! OpN is now the number of potential terms to include in the derivative
      HessianOperator%N = OpN

      ! Now we can allocate memory to store the potential derivative operator
      ALLOCATE( HessianOperator%Coeffs(OpN)     )
      ALLOCATE( HessianOperator%OpType(3,OpN)   )
      ALLOCATE( HessianOperator%OpPar(HessianOperator%gDim ,OpN) )

      ! Define HessianOperator%OpType copying from Operator%OpType
      HessianOperator%Coeffs( 1:OpN ) = Operator%Coeffs( IndexMask(1:OpN) )
      HessianOperator%OpType( 1:3, 1:OpN ) = Operator%OpType( 1:3, IndexMask(1:OpN) )
      HessianOperator%OpPar( :, 1:OpN ) = Operator%OpPar( :, IndexMask(1:OpN) )

      DO n = 1, OpN
         i = IndexMask(n)
         IF ( iCoord /= jCoord ) THEN
            HessianOperator%Coeffs(n) = Operator%Coeffs(i) * REAL(Operator%OpPar(iCoord , i)) * REAL(Operator%OpPar(jCoord , i))
            HessianOperator%OpPar(iCoord , n) = Operator%OpPar(iCoord , i)-1
            HessianOperator%OpPar(jCoord , n) = Operator%OpPar(jCoord , i)-1
         ELSE
            HessianOperator%Coeffs(n) = Operator%Coeffs(i) * REAL(Operator%OpPar(iCoord , i)) * REAL(Operator%OpPar(iCoord , i)-1)
            HessianOperator%OpPar(iCoord , n) = Operator%OpPar(iCoord , i)-2
         END IF
      END DO

   END SUBROUTINE SetHessianPotentialOperator


!*******************************************************************************
!          LogOperator
!*******************************************************************************
!> write to output a log file with the detailed description of the
!> operator defined in the input OperatorData instance
!>
!> @param Operator            OperatorData istance to describe in the log file
!*******************************************************************************
   SUBROUTINE LogOperator( Operator, FileName )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(INOUT) :: Operator
      CHARACTER(*), INTENT(IN)          :: FileName
      INTEGER :: OperInp, iStat, TermNumber, iGDim, i

      ! In case memory is not allocated, data is not setup and log cannot be done! abort with error!
      CALL ERROR((.NOT. ALLOCATED(Operator%Coeffs)).OR.(.NOT. ALLOCATED(Operator%OpType)).OR.(.NOT. ALLOCATED( Operator%OpPar )),&
         " OperatorDefine.LogOperator: cannot print info on operator which is not defined", ERR_OBJ_MISUSE )

      ! Now write a log of the defined Operator to an output file

      OperInp = LookForFreeUnit()
      OPEN(UNIT=OperInp, FILE=TRIM(ADJUSTL(FileName)), STATUS='REPLACE', ACTION='WRITE', IOSTAT=iStat)
      CALL ERROR( iStat /= 0, "OperatorDefine.LogOperator: could not write file Operator logfile", ERR_IO_OPEN )

      ! number of system degrees of freedom
      WRITE(OperInp,'("gdim=(")',ADVANCE='no')
      WRITE(OperInp,'(i2,")")',ADVANCE='yes') Operator%gDim

      ! write the number of operator terms
      WRITE(OperInp,'("me=",I3)') Operator%N

      ! write the coefficients of the operator terms
      WRITE(OperInp,'("-------------------")')
      WRITE(OperInp,'("The a_{r}")')
      WRITE(OperInp,*)''
      WRITE(OperInp,'("  r   a_{r}")',ADVANCE='yes')
      WRITE(OperInp,'(" -----------------")',ADVANCE='yes')
      DO TermNumber = 1, Operator%N
         WRITE(OperInp,'(I3,1X,f8.5)',ADVANCE='yes') TermNumber, Operator%Coeffs(TermNumber)
      END DO
      WRITE(OperInp,*) ''
      WRITE(OperInp,'("-------------------")')

      ! write the information about each operator term
      WRITE(OperInp,'("Definition of the operator:")')
      WRITE(OperInp,*)''
      WRITE(OperInp,'("no :   r    a_{r}    s1    s2    type ")',ADVANCE='no')
      DO iGDim = 1, Operator%gDim
         WRITE(OperInp,'("|",i2)',ADVANCE='no') mod(iGDim,10)
      END DO
      WRITE(OperInp,*) ''
      DO i=1,38
         WRITE(OperInp,'("-")',ADVANCE='no')
      END DO
      WRITE(OperInp,'("-------------------")')

      DO TermNumber = 1, Operator%N
         WRITE(OperInp,'(i3,":",1X,i3,2x,f8.5,1X,i3,2x,1X,i3,2x,1X,i3,4x)',ADVANCE='no') &
              TermNumber, TermNumber, Operator%Coeffs(TermNumber), &
              Operator%OpType(2,TermNumber), Operator%OpType(3,TermNumber), Operator%OpType(1,TermNumber)
         DO iGDim = 1, Operator%gDim
            IF ( Operator%OpPar(iGDim,TermNumber) /= 0 ) THEN
               WRITE(OperInp,'("|",i2)',ADVANCE='no') Operator%OpPar(iGDim,TermNumber)
            ELSE
               WRITE(OperInp,'("|  ")',ADVANCE='no')
            END IF
         END DO  ! d=1,gdim
         WRITE (OperInp,*) ''
      END DO ! re
      WRITE (OperInp,*) ''

   END SUBROUTINE LogOperator



!*******************************************************************************
!          DisposeOperator
!*******************************************************************************
!> OperatorData destructor: deallocate memory of an instance of OperatorData
!>
!> @param Operator            OperatorData istance to destruct
!*******************************************************************************
   SUBROUTINE DisposeOperator( Operator )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(INOUT) :: Operator

      ! In case memory is already allocated, deallocate
      IF ( ALLOCATED(Operator%Coeffs) ) DEALLOCATE( Operator%Coeffs )
      IF ( ALLOCATED(Operator%OpType) ) DEALLOCATE( Operator%OpType )
      IF ( ALLOCATED(Operator%OpPar) )  DEALLOCATE( Operator%OpPar )

#if defined(LOG_FILE)
      __OPEN_LOG_FILE;
      __WRITE_TO_LOG "DisposeOperator: an instance of OperatorData has been desctructed"
      __WHITELINE_TO_LOG; __CLOSE_LOG_FILE
#endif
   END SUBROUTINE DisposeOperator



!*******************************************************************************
!          OperatorMatrixElement
!*******************************************************************************
!> Given an input Operator of type OperatorData and two gaussian configurations
!> with integer indeces to label their electronic state,
!> compute the matrix element of the operator over the given bra and ket.
!> Optionally, augment the operator multiplying by a monomial in one dimension.
!>
!> @param    Operator       OperatorData object with the operator
!> @param    GParLeft       Bra gaussian configuration
!> @param    GParRight      Ket gaussian configuration
!> @param    iElL, iElR     Labels for the electronic state of the left and right gaussian cfg
!> @param    AugDim, AugF   augment the integral with x_d^n with d = aug_dim and n = aug_f
!> @returns  H              the matrix element
!*******************************************************************************
   COMPLEX FUNCTION OperatorMatrixElement( Operator, GParLeft, GParRight, iElL, iElR, AugDim, AugF ) RESULT( H )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(IN)                                   :: Operator
      COMPLEX, DIMENSION(:,:), INTENT(IN)                              :: GParLeft
      COMPLEX, DIMENSION(:,:), INTENT(IN)                              :: GParRight
      INTEGER, INTENT(IN)                                              :: iElL, iElR
      INTEGER, OPTIONAL, INTENT(IN)                                    :: AugDim, AugF

      COMPLEX :: OpTermMoment, FirstMom, SecondMom, ZeroMom, ACoeff, Eta
      INTEGER :: iOp, iD
      INTEGER, DIMENSION(Operator%gDim)  :: OrderVec
      
      CALL ERROR( SIZE(GParLeft,1) /= 3 .OR. SIZE(GParLeft,2) /= Operator%gDim,   " OperatorMatrixElement: left gaussian cfg dimension mismatch" )
      CALL ERROR( SIZE(GParRight,1) /= 3 .OR. SIZE(GParRight,2) /= Operator%gDim, " OperatorMatrixElement: left gaussian cfg dimension mismatch" )

      ! Initialize matrix element
      H = CMPLX(0.0, 0.0)

      ! Loop over the elements of the hamiltonian operator
      DO iOp = 1, Operator%N

         ! Use the current element iOp only if it couples the given electronic states
         IF ((Operator%OpType(2,iOp) /= iElL) .OR. (Operator%OpType(3,iOp) /= iElR))  CYCLE

         ! select construct handles the different type of operators defined in Operator%OpType(1,iOp)
         SELECT CASE (Operator%OpType(1,iOp))

         ! ===============================================================================
            CASE( OPTYPE_KINETIC )
         ! ===============================================================================

               ! Cycle over the degrees of freedom
               DO iD = 1, Operator%gDim

                  ! skip those degrees of freedom for which OpPar == 0
                  IF ( Operator%OpPar(iD,iOp) == 0 ) CYCLE

                  IF ( PRESENT(AugDim) .AND. PRESENT(AugF) ) THEN

                     IF ( AugDim == iD ) THEN ! in this case the derivation coordinate moments are increase by AugF
                        ! Compute the three necessary moments which sum up to give 2nd derivative
                        ZeroMom   = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), AugF )
                        FirstMom  = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), AugF+1 )
                        SecondMom = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), AugF+2 )

                     ELSE  ! in this case the deriv coord and the additional factor give different prim mom which are multiplied
                        ! Compute the moments
                        ZeroMom   = GauPrim_Moment( GParLeft(:,AugDim), GParRight(:,AugDim), AugF )
                        FirstMom  = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), 1 ) * ZeroMom
                        SecondMom = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), 2 ) * ZeroMom
                     END IF

                  ELSE
                     ! Compute the three necessary moments which sum up to give 2nd derivative
                     ZeroMom   = CMPLX(1.0,0.0)
                     FirstMom  = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), 1 )
                     SecondMom = GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), 2 )
                  END IF

                  ! extract the gaussian parameters of the right gaussian
                  ACoeff = GParRight(1, iD)
                  Eta    = GParRight(2, iD)

                  ! compute the second derivative combining the computed moments
                  OpTermMoment = -0.5 * (4.0*ACoeff**2*SecondMom + 4.0*ACoeff*Eta*FirstMom + (2.0*ACoeff+Eta**2) * ZeroMom )
                  EXIT        ! Exit from loop, only one coordinate is expected in this kinetic energy term

               END DO

         ! ===============================================================================
            CASE( OPTYPE_POTENTIAL )
         ! ===============================================================================

               ! *****************************************************************************************
               ! TODO: in case of LOCAL POTENTIAL EXPANSION, HERE THERE IS A CONDITION TO
               ! SKIP THIS TERM
               ! *****************************************************************************************

               ! In this case, iOp-th entry of OpPar gives the array of the powers of this potential term
               OrderVec = Operator%OpPar(:,iOp)
               ! if needed, increment the potential term with the arguments AugDim, AugF
               IF ( PRESENT(AugDim) .AND. PRESENT(AugF) )  OrderVec(AugDim) = OrderVec(AugDim) + AugF

               ! Initialize moment value
               OpTermMoment = CMPLX(1.0, 0.0)

               ! Cycle over the coordinates and multiply the initialized value of OpTermMoment by the
               ! appropriate coordinate component
               DO iD = 1, Operator%gDim
                  ! When zero-th order moment is required, the value to multiply is 1.0 and we can cycle
                  IF ( OrderVec(iD) == 0 ) CYCLE
                  ! compute the primitive moment for the actual iPtr and jPtr gaussians
                  OpTermMoment = OpTermMoment * &
                     GauPrim_Moment( GParLeft(:,iD), GParRight(:,iD), OrderVec(iD) )
               END DO

            CASE( OPTYPE_UNDEFINED )
               CALL AbortWithError( "HamiltonianElement: Undefined operator with iOp = "//TRIM(NumberToString(iOp)), ERR_GENERIC )

         END SELECT

         ! Multiply by the linear coefficient of the operator term and increment sum over potential terms
         H = H + Operator%Coeffs(iOp) * OpTermMoment

      END DO

      ! Multiply by the overlap (since GauPrim_Moment returns the value divided by the overlap!)
      H = H * GauConf_Overlap( GParLeft, GParRight )

      ! **********************************************************************************************************
      ! TODO: in case of LOCAL POTENTIAL EXPANSION, HERE THERE IS THE COMPUTATION OF THE MOMENTS IN THE LHA/LCA
      ! **********************************************************************************************************

   END FUNCTION OperatorMatrixElement



!    !##################################################################################
!    !#                                                                                #
!    !#                          Subroutine FuncV_LHA                                  #
!    !#                                                                                #
!    !# Calculates the potential hamiltionian. At the moment only polynomial-terms are #
!    !# implemented. Here local approximations are used, up to the third order         #
!    !# derivatives.                                                                   #
!    !# Again an additional x_d^n term can be added by using aug_dim and aug_f         #
!    !#                                                                                #
!    !##################################################################################
!
!    COMPLEX (KIND=dop) FUNCTION FuncV_LHA(temp_Gaussians,re,l_l,l_r,RefVector,aug_dim,aug_f)
!       IMPLICIT NONE
!       COMPLEX (KIND=dop), DIMENSION(:,:), INTENT(IN) :: temp_Gaussians
!       INTEGER(KIND=ilong), INTENT(IN)                :: re
!       INTEGER(KIND=ilong), INTENT(IN)                :: l_l, l_r
!       REAL(KIND=dop), DIMENSION(gdim), INTENT(IN)    :: RefVector
!       INTEGER(KIND=ilong), OPTIONAL, INTENT(IN)      :: aug_dim,aug_f
!
!       INTEGER(KIND=ilong):: el_l, el_r
!       INTEGER(KIND=ilong):: d1, d2, d3, te
!       COMPLEX(KIND=dop) :: d1_firstmom, d2_firstmom, aug_mom, d1_secmom
!
!       INTEGER(KIND=ilong), DIMENSION(gdim)  :: PowerList
!
!       REAL(KIND=dop)                        :: Potential
!       REAL(KIND=dop), DIMENSION(gdim)       :: Gradient
!       REAL(KIND=dop), DIMENSION(gdim,gdim)  :: Hessian
!
!       REAL(KIND=dop)                        :: Mom0
!       REAL(KIND=dop), DIMENSION(gdim)       :: Mom1
!       REAL(KIND=dop), DIMENSION(gdim,gdim)  :: Mom2
!
!       ! Initialize  variable to zero
!       FuncV_LHA = cmplx(dZero,dZero,KIND=dop)
!
!       ! Define the electronic indices of the operator
!       IF (multiset) THEN
!          el_l = op_re(re,2)
!          el_r = op_re(re,3)
!       ELSE
!          el_l = iOne
!          el_r = iOne
!       END IF
!
!       ! Compute potential, gradient and hessian
!       IF (DDYN) THEN
!          ! In the case of direct dynamics, copy variables into local arrays
!          Potential = ddyn_energy(l_r)
!          Gradient  = ddyn_gradient(:,l_r)
!          Hessian   = ddyn_hessian(:,:,l_r)
!
!          WRITE(*,*) " DEFINITION OF REFERENCE VECTOR IS REQUIRED! NOT YET IMPLEMENTED"
!          STOP
!       ELSE
!          PowerList = op_primdim(:,re)
!          ! In the case of analytic potential, compute the terms
!          ! potential is simply evaluated as zero-th order derivative
!          Potential = EvaluatePotential( PowerList, RefVector )
!          Gradient  = EvaluateGradient( PowerList, RefVector )
!          Hessian   = EvaluateHessian( PowerList, RefVector )
!       END IF
!
!       ! Define factors for evaluating potential
!       Mom0      = Potential
!       DO d1 = 1, gdim
!          Mom0     = Mom0 - Gradient(d1)*RefVector(d1)
!          Mom1(d1) = Gradient(d1)
!          DO d2 = 1, gdim
!             Mom0        = Mom0 + 0.5 * Hessian(d2,d1)*RefVector(d1)*RefVector(d2)
!             Mom1(d1)    = Mom1(d1) - 0.5*(Hessian(d2,d1)+Hessian(d1,d2))*RefVector(d2)
!             Mom2(d1,d2) = 0.5 * Hessian(d1,d2)
!          END DO
!       END DO
!
!       IF (present(aug_dim)) THEN
!
!          ! compute overlap between the two gaussians
!          aug_mom = moment_f_1d_new(temp_gaussians,l_l,l_r,aug_dim,aug_f,el_l,el_r)
!
!          ! Zeroth order term
!          FuncV_LHA = cmplx(Mom0,dZero,KIND=dop) * aug_mom
!
!          ! loop over the coordinates
!          DO d1 = 1, gdim
!
!             ! Compute moments of d1
!             IF ( d1 == aug_dim ) THEN
!                d1_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d1,aug_f+iOne,el_l,el_r)
!                d1_secmom   = moment_f_1d_new(temp_gaussians,l_l,l_r,d1,aug_f+2_ilong,el_l,el_r)
!             ELSE
!                d1_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d1,iOne,el_l,el_r)*aug_mom
!                d1_secmom   = moment_f_1d_new(temp_gaussians,l_l,l_r,d1,2_ilong,el_l,el_r)*aug_mom
!             END IF
!
!             ! Increment FuncV_LHA with first order term
!             FuncV_LHA = FuncV_LHA + cmplx(Mom1(d1),dZero,KIND=dop) * d1_firstmom
!             ! Increment FuncV_LHA with second order diagonal terms
!             FuncV_LHA = FuncV_LHA + cmplx(Mom2(d1,d1),dZero,KIND=dop) * d1_secmom
!
!             ! loop over the coordinates again, for hessian matrix
!             DO d2 = 1, gdim
!
!                ! Diagonal elements have already been computed
!                IF (d1 == d2) CYCLE
!
!                ! Compute moments of d2 and d1
!                IF ( aug_dim == d2 ) THEN
!                   d2_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d2,aug_f+iOne,el_l,el_r)
!                ELSE
!                   d2_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d2,iOne,el_l,el_r)*aug_mom
!                END IF
!
!                ! Compute moments of d2
!                d2_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d2,iOne,el_l,el_r)
!                ! Increment FuncV_LHA with second order term
!                FuncV_LHA = FuncV_LHA + cmplx(Mom2(d2,d1),dZero,KIND=dop) * d1_firstmom * d2_firstmom
!
!             END DO
!          END DO
!
!       ELSE ! aug_dim not present -> normal calculation
!
!          ! Zeroth order term
!          FuncV_LHA = cmplx(Mom0,dZero,KIND=dop)
!
!          ! loop over the coordinates
!          DO d1 = 1, gdim
!
!             ! Compute moments of d1
!             d1_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d1,iOne,el_l,el_r)
!             d1_secmom   = moment_f_1d_new(temp_gaussians,l_l,l_r,d1,2_ilong,el_l,el_r)
!             ! Increment FuncV_LHA with first order term
!             FuncV_LHA = FuncV_LHA + cmplx(Mom1(d1),dZero,KIND=dop) * d1_firstmom
!             ! Increment FuncV_LHA with second order diagonal terms
!             FuncV_LHA = FuncV_LHA + cmplx(Mom2(d1,d1),dZero,KIND=dop) * d1_secmom
!
!             ! loop over the coordinates again, for hessian matrix
!             DO d2 = 1, gdim
!
!                IF ( d1 == d2 ) CYCLE
!
!                ! Compute moments of d2
!                d2_firstmom = moment_f_1d_new(temp_gaussians,l_l,l_r,d2,iOne,el_l,el_r)
!                ! Increment FuncV_LHA with second order term
!                FuncV_LHA = FuncV_LHA + cmplx(Mom2(d2,d1),dZero,KIND=dop) * d1_firstmom * d2_firstmom
!
!             END DO
!          END DO
!
!       END IF ! aug_dim present ... else ...
!
!    END FUNCTION FuncV_LHA
!
!
!===========================================================================================================


!*******************************************************************************
!          EvaluatePotential
!*******************************************************************************
!> Given input Hamiltonian and a point in coordinate space, compute the
!> potential functions contained in the Hamiltonian in the given coordinate
!> point. The subroutine can be used then to compute the potential in the
!> center of a gaussian wavepacket.
!>
!> @returns Potential     output gives real value of the potential
!> @param Hamiltonian     input, takes the OperatorData object which is used as potential
!> @param Coord           input, takes the coordinates of the vector
!> @param ElState         input, electronic state to consider for the potential
!*******************************************************************************

   REAL FUNCTION EvaluatePotential( Hamiltonian, Coord, ElState ) RESULT( Potential )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(IN)                   :: Hamiltonian
      REAL, DIMENSION(Hamiltonian%gDim), INTENT(IN)   :: Coord
      INTEGER, INTENT(IN), OPTIONAL                    :: ElState
      REAL :: PowProduct
      INTEGER :: iOp, nEl, jDim

      ! In case memory of Hamiltonian is not allocated, exit with error
      CALL ERROR( .NOT. ALLOCATED(Hamiltonian%Coeffs) .OR. .NOT. ALLOCATED(Hamiltonian%OpType) .OR. .NOT. ALLOCATED( Hamiltonian%OpPar ),  &
         &        " OperatorDefine.EvaluatePotential: source OperatorData object is not defined (Hamiltonian)", ERR_OBJ_MISUSE )

      ! Define the electronic index of the potential energy functions which are included in the potetnial evaluation (diagonal terms only)
      IF ( PRESENT(ElState) ) THEN
         nEl = ElState
      ELSE
         nEl = 1
      END IF

      ! Initialize potential
      Potential = 0.0

      ! loop over the operators which are included in Hamiltonian
      DO iOp = 1, Hamiltonian%N

         ! use only potential operators
         IF ( Hamiltonian%OpType(1,iOp) /= OPTYPE_POTENTIAL ) CYCLE
         ! use only potential terms which are diagonal on the given electronic state
         IF ( Hamiltonian%OpType(2,iOp) /= nEl .OR. Hamiltonian%OpType(3,iOp) /= nEl ) CYCLE

         ! initialize product of coordinate powers
         PowProduct = 1.0
         ! accumulate product
         DO jDim = 1, Hamiltonian%gDim
            PowProduct = PowProduct * Coord(jDim)**Hamiltonian%OpPar(jDim,iOp)
         END DO
         ! increment potential, sum over products
         Potential = Potential + Hamiltonian%Coeffs(iOp) * PowProduct

      END DO

   END FUNCTION EvaluatePotential


!===========================================================================================================


!*******************************************************************************
!          EvaluateGradient
!*******************************************************************************
!> Given input Hamiltonian and a point in coordinate space, compute the
!> first order derivatives of the potential functions contained in the Hamiltonian
!> in the given coordinate point. The subroutine can be used then to compute
!> the gradient of the potential in the center of a gaussian wavepacket.
!>
!> @returns Gradient      output gives array of reals with the gradient
!> @param Hamiltonian     input, takes the OperatorData object which is used as potential
!> @param Coord           input, takes the coordinates of the vector
!> @param ElState         input, electronic state to consider for the potential
!*******************************************************************************

   FUNCTION EvaluateGradient( Hamiltonian, Coord, ElState ) RESULT( Gradient )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(IN)                  :: Hamiltonian
      REAL, DIMENSION(Hamiltonian%gDim), INTENT(IN)   :: Coord
      INTEGER, INTENT(IN), OPTIONAL                   :: ElState
      REAL, DIMENSION(Hamiltonian%gDim)               :: Gradient
      REAL :: PowProduct
      INTEGER :: iOp, nEl, jDim, lDim

      ! In case memory of Hamiltonian is not allocated, exit with error
      CALL ERROR( .NOT. ALLOCATED(Hamiltonian%Coeffs) .OR. .NOT. ALLOCATED(Hamiltonian%OpType) .OR. .NOT. ALLOCATED( Hamiltonian%OpPar ),  &
         &        " OperatorDefine.EvaluatePotential: source OperatorData object is not defined (Hamiltonian)", ERR_OBJ_MISUSE )

      ! Define the electronic index of the potential energy functions which are included in the potetnial evaluation (diagonal terms only)
      IF ( PRESENT(ElState) ) THEN
         nEl = ElState
      ELSE
         nEl = 1
      END IF

      ! Initialize gradient
      Gradient = 0.0

      ! loop over the operators which are included in Hamiltonian
      DO iOp = 1, Hamiltonian%N

         ! use only potential operators
         IF ( Hamiltonian%OpType(1,iOp) /= OPTYPE_POTENTIAL ) CYCLE
         ! use only potential terms which are diagonal on the given electronic state
         IF ( Hamiltonian%OpType(2,iOp) /= nEl .OR. Hamiltonian%OpType(3,iOp) /= nEl ) CYCLE

         ! jDim is the element of the gradient, which is the coordinate along which the derivation is done
         DO jDim = 1, Hamiltonian%gDim

            ! when deriving a constant function along jDim, overall derivative is zero
            IF ( Hamiltonian%OpPar(jDim,iOp)  == 0 ) THEN
               CYCLE

            ELSE IF ( Hamiltonian%OpPar(jDim,iOp) > 0 ) THEN
               ! Initialize the derivative with the numerical factor ( = exponent of the term to derive )
               PowProduct = REAL( Hamiltonian%OpPar(jDim,iOp) )

               ! Evaluate the monomial with the right power reduced by one
               DO lDim = 1, Hamiltonian%gDim
                  IF ( lDim == jDim ) THEN
                     PowProduct = PowProduct * Coord(lDim)**( Hamiltonian%OpPar(lDim,iOp) - 1 )
                  ELSE
                     PowProduct = PowProduct * Coord(lDim)**Hamiltonian%OpPar(lDim,iOp)
                  END IF
               END DO

               ! Increment gradient
               Gradient(jDim) = Gradient(jDim) + Hamiltonian%Coeffs(iOp) * PowProduct
            END IF

         END DO
      END DO

   END FUNCTION EvaluateGradient


!===========================================================================================================


!*******************************************************************************
!          EvaluateHessian
!*******************************************************************************
!> Given input Hamiltonian and a point in coordinate space, compute the
!> second order derivatives of the potential functions contained in the Hamiltonian
!> in the given coordinate point. The subroutine can be used then to compute
!> the hessian matrix of the potential in the center of a gaussian wavepacket.
!>
!> @returns Hessian       output gives matrix of reals with the hessian
!> @param Hamiltonian     input, takes the OperatorData object which is used as potential
!> @param Coord           input, takes the coordinates of the vector
!> @param ElState         input, electronic state to consider for the potential
!*******************************************************************************

   FUNCTION EvaluateHessian( Hamiltonian, Coord, ElState ) RESULT( Hessian )
      IMPLICIT NONE
      TYPE(OperatorData), INTENT(IN)                      :: Hamiltonian
      REAL, DIMENSION(Hamiltonian%gDim), INTENT(IN)       :: Coord
      INTEGER, INTENT(IN), OPTIONAL                       :: ElState
      REAL, DIMENSION(Hamiltonian%gDim,Hamiltonian%gDim)  :: Hessian
      REAL :: PowProduct
      INTEGER :: iOp, nEl, jDim, lDim, nDim

      ! In case memory of Hamiltonian is not allocated, exit with error
      CALL ERROR( .NOT. ALLOCATED(Hamiltonian%Coeffs) .OR. .NOT. ALLOCATED(Hamiltonian%OpType) .OR. .NOT. ALLOCATED( Hamiltonian%OpPar ),  &
         &        " OperatorDefine.EvaluatePotential: source OperatorData object is not defined (Hamiltonian)", ERR_OBJ_MISUSE )

      ! Define the electronic index of the potential energy functions which are included in the potetnial evaluation (diagonal terms only)
      IF ( PRESENT(ElState) ) THEN
         nEl = ElState
      ELSE
         nEl = 1
      END IF

      ! Initialize hessian
      Hessian = 0.0

      ! loop over the operators which are included in Hamiltonian
      DO iOp = 1, Hamiltonian%N

         ! use only potential operators
         IF ( Hamiltonian%OpType(1,iOp) /= OPTYPE_POTENTIAL ) CYCLE
         ! use only potential terms which are diagonal on the given electronic state
         IF ( Hamiltonian%OpType(2,iOp) /= nEl .OR. Hamiltonian%OpType(3,iOp) /= nEl ) CYCLE

         ! first compute diagonal elements...

         ! jDim,jDim is the index of the diagonal element of hessian, which equals the coordinates to derive two times
         DO jDim = 1, Hamiltonian%gDim
            ! when 2nd-deriv of a constant or of a linear function along jDim, overall derivative is zero
            IF ( Hamiltonian%OpPar(jDim,iOp) == 0 .OR. Hamiltonian%OpPar(jDim,iOp) == 1 ) THEN
               CYCLE

            ELSE IF ( Hamiltonian%OpPar(jDim,iOp) > 1 ) THEN
               ! Initialize the derivative with the numerical factor ( = exponent of the term to derive )
               PowProduct = REAL( Hamiltonian%OpPar(jDim,iOp) ) * REAL( Hamiltonian%OpPar(jDim,iOp)-1 )
               ! Evaluate the monomial with the right power reduced by two
               DO lDim = 1, Hamiltonian%gDim
                  IF ( lDim == jDim ) THEN
                     PowProduct = PowProduct * Coord(lDim)**( Hamiltonian%OpPar(lDim,iOp) - 2 )
                  ELSE
                     PowProduct = PowProduct * Coord(lDim)**Hamiltonian%OpPar(lDim,iOp)
                  END IF
               END DO
               ! Increment gradient
               Hessian(jDim, jDim) = Hessian(jDim, jDim) + Hamiltonian%Coeffs(iOp) * PowProduct
            END IF
         END DO

         ! ... and the off-diagonal ones taking advantage of hessian symmetry

         ! jDim, nDim are the indices of the element of the hessian, which equals the two coordinates to derive
         DO nDim = 1, Hamiltonian%gDim
            DO jDim = nDim+1, Hamiltonian%gDim
               ! when deriving a constant function along jDim or nDim, overall derivative is zero
               IF ( Hamiltonian%OpPar(jDim,iOp) == 0 .OR. Hamiltonian%OpPar(nDim,iOp) == 0 ) THEN
                  CYCLE

               ELSE IF ( Hamiltonian%OpPar(jDim,iOp) > 0 .AND. Hamiltonian%OpPar(nDim,iOp) > 0 ) THEN
                  ! Initialize the derivative with the correct numerical factor ( = exponent of the term to derive )
                  PowProduct = REAL( Hamiltonian%OpPar(jDim,iOp) ) * REAL( Hamiltonian%OpPar(nDim,iOp) )
                  ! Evaluate the monomial with the right power reduced by one
                  DO lDim = 1, Hamiltonian%gDim
                     IF ( lDim == jDim .OR. lDim == nDim ) THEN
                        PowProduct = PowProduct * Coord(lDim)**( Hamiltonian%OpPar(lDim,iOp) - 1 )
                     ELSE
                        PowProduct = PowProduct * Coord(lDim)**Hamiltonian%OpPar(lDim,iOp)
                     END IF
                  END DO
                  ! Increment gradient
                  Hessian(jDim, nDim) = Hessian(jDim, nDim) + Hamiltonian%Coeffs(iOp) * PowProduct
               END IF
            END DO ! loop over nDim
         END DO ! loop over jDim

      END DO ! loop over iOp

      ! now copy lower triangle to the upper one
      DO nDim = 1, Hamiltonian%gDim
         DO jDim = nDim+1, Hamiltonian%gDim
            Hessian(nDim, jDim) = Hessian(jDim, nDim)
         END DO
      END DO

   END FUNCTION EvaluateHessian

!===========================================================================================================

END MODULE OperatorDefine
