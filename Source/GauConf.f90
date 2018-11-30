!***************************************************************************************
!*                              MODULE GauConf
!***************************************************************************************
!
!>  \brief     Miscellanous utilities for handling gaussian functions
!>  \details   This module contains general purpose subroutines which can      \n
!>             be used to perform common tasks regarding gaussian functions,   \n
!>             e.g. normalization, integrals of different forms ... etc. etc.  \n
!>             For gaussian function here we mean a multi-dimensional function \n
!>             expressed as a product of 1D gaussians, with arbitrary complex  \n
!>             coefficients. Adopting a MCTDH-like naming convention, We call  \n
!>             GAUSSIAN CONFIGURATION such multidimensional function, and      \n
!>             PRIMITIVE GAUSSIAN FUNCTION the 1D gaussian function.           \n
!>             The standard definition of these gaussian configuration in this \n
!>             module is a list of coefficients GAUSSIAN(3,NDim) where NDim    \n
!>             is the dimension of the space and 3 is the number of parameters \n
!>             to specify: a, xi, eta such that the primitive gaussian is      \n
!>             g(x) = a * x^2 + xi * x + eta, with a,xi,eta complex.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             18 July 2017
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          ______________________________________
!
!***************************************************************************************
MODULE GauConf
#include "preprocessoptions.cpp"
   IMPLICIT NONE

   !> Order of the higher gaussian moments which are stored in the temporary array
   INTEGER, PARAMETER :: HigherOrder = 5

   PRIVATE

   PUBLIC :: MomentStorageArea
   PUBLIC :: EmptyMomentStorageArea

   PUBLIC :: GauConf_Normalize
   PUBLIC :: GauConf_Overlap
   PUBLIC :: GauPrim_Moment
   PUBLIC :: GauPrim_Normalize

   ! Two versions of the subroutine to compute moments of product gaussian distribution
   ! the first uses a storage arrays to recover previously computed values
   ! the second one computes them from scratch
   INTERFACE GauPrim_Moment
      MODULE PROCEDURE GauPrim_Moment_Storage, GauPrim_Moment_NoStorage
   END INTERFACE

   !> temporary storage of the moments of the primitive gaussians ( more precisely the integrals <gi|xf^n|gj>/<gi|gj> ),
   !> when a new value is computed then is added to the storage and every time the same integral is required,
   !> the relevant subroutine returns the stored value. The logical array CheckStore tells when the value
   !> has been already computed. Moments up the order HigherOrder are stored, the other - when required - are always recomputed.
   TYPE MomentStorageArea
      PRIVATE
      COMPLEX, DIMENSION(HigherOrder) :: Moment           !< array to store computed gaussian moments
      LOGICAL, DIMENSION(HigherOrder) :: Check = .FALSE.  !< check if a gaussian moment has already been computed
   END TYPE MomentStorageArea

   CONTAINS


!*******************************************************************************
!          EmptyMomentStorageArea
!*******************************************************************************
!> Returns an initialized TYPE(MomentStorageArea) with Moment = CMPLX(0.0,0.0)
!> and Check = .FALSE., which can be used to initialize the derived type
!> outside the module.
!>
!> @return    EmptyMomentStorageArea   initialized MomentStorageArea data type
!*******************************************************************************
   TYPE(MomentStorageArea) FUNCTION EmptyMomentStorageArea( )

      EmptyMomentStorageArea%Moment = CMPLX(0.0,0.0)
      EmptyMomentStorageArea%Check  = .FALSE.

   END FUNCTION EmptyMomentStorageArea


!*******************************************************************************
!          GauConf_Normalize
!*******************************************************************************
!> Given an input gaussian configuration, change the eta parameters in order
!> to normalize the function.
!>
!> @param   GaussCfg           Gaussian configuration to normalize
!> @param   NDim               Number of dimensions of the gaussian cfg
!*******************************************************************************
   SUBROUTINE GauConf_Normalize(GaussCfg)
      COMPLEX, DIMENSION(:,:), INTENT(INOUT) :: GaussCfg
      INTEGER :: iDim
      REAL    :: Eta_Real, A_Real, Xi_Real, Eta_Imag

      ! Loop over the primitive gaussians
      DO iDim = 1, SIZE(GaussCfg,2)
         ! Extract the values necessary for normalization
         A_Real  = real(GaussCfg(1, iDim))
         Xi_Real = real(GaussCfg(2, iDim))
         ! Compute the real part from normalization formula
         Eta_Real = 1.0/4.0*( (Xi_Real**2)/A_Real + log(-2.0*A_Real/MyConsts_PI) )
         Eta_Imag = aimag(GaussCfg(3, iDim))
         ! Substitute the value of eta
         GaussCfg(3, iDim) = CMPLX(Eta_Real,Eta_Imag)
      END DO

  END SUBROUTINE GauConf_Normalize


!*******************************************************************************
!          GauConf_Overlap
!*******************************************************************************
!> Given two gaussian configurations, calculate the overlap integral.
!> In formulas:  Ov = < g_left | g_right >
!> This subroutine uses the formula:
!>          <gleft|gright> = exp(q+p**2/(4o)) * sqrt(pi/o)
!> with: o=-(aleft*+aright), p=(xileft*+xiright), q=etaleft*+etaright
!>
!> @param   GParLeft           Left gaussian configuration
!> @param   GParRight          Right gaussian configuration
!> @param   NDim               Number of dimensions of the gaussian cfg
!*******************************************************************************

   COMPLEX FUNCTION GauConf_Overlap( Left, Right )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN) :: Left, Right
      INTEGER :: i
      COMPLEX:: p, q
      REAL :: o

      ! Initialize overlap
      GauConf_Overlap = CMPLX( 1.0, 0.0 )

      ! Accumulate product of the 1D overlap between primitive gaussians
      DO i = 1, SIZE(Left,2)
         o = - REAL(Left(1,i)) - REAL(Right(1,i))
         p =  CONJG(Left(2,i)) + Right(2,i)
         q =  CONJG(Left(3,i)) + Right(3,i)
         GauConf_Overlap = GauConf_Overlap * EXP(p**2/(4.0*o)+q)*SQRT(MyConsts_PI/o)
      END DO

  END FUNCTION GauConf_Overlap


!*******************************************************************************
!          GauPrim_Normalize
!*******************************************************************************
!> Given an input gaussian primitive function, change the eta parameters in order
!> to normalize the function.
!>
!> @param   GaussPrim          Gaussian function to normalize
!> @param   NDim               Number of dimensions of the gaussian cfg
!*******************************************************************************
   SUBROUTINE GauPrim_Normalize(GaussPrim)
      COMPLEX, DIMENSION(3), INTENT(INOUT) :: GaussPrim
      REAL    :: Eta_Real, A_Real, Xi_Real, Eta_Imag

      ! Extract the values necessary for normalization
      A_Real  = real(GaussPrim(1))
      Xi_Real = real(GaussPrim(2))
      ! Compute the real part from normalization formula
      Eta_Real = 1.0/4.0*( (Xi_Real**2)/A_Real + log(-2.0*A_Real/MyConsts_PI) )
      Eta_Imag = aimag(GaussPrim(3))
      ! Substitute the value of eta
      GaussPrim(3) = CMPLX(Eta_Real,Eta_Imag)

  END SUBROUTINE GauPrim_Normalize

!*******************************************************************************
!          GauPrim_Moment_Storage
!*******************************************************************************
!> Given two gaussian primitive functions, returns a moment of the product
!> distribution divided by the overlap of the gaussians.
!> In formulas, this function returns:
!>    < g_left(x) | x^n | g_right(x) > / < g_left(x) | g_right(x) >
!> where n specifies the order of the moment.
!> This specific subroutine uses a storage area to recover previously computed
!> moments and add new computed values to the storage itself.
!> If the moment needs to be really computed, the other GauPrim_Moment_NoStorage
!> is called.
!>
!> @param    GParLeft        Left gaussian configuration
!> @param    GParRight       Right gaussian configuration
!> @param    Order           Order of the moment to compute
!> @param    Storage         TYPE(MomentStorageArea) with the storage area for the moments
!> @returns  the complex value of the moment
!*******************************************************************************

   COMPLEX FUNCTION GauPrim_Moment_Storage( GParLeft, GParRight, Storage, Order ) RESULT( PrimMom )
      IMPLICIT NONE
      COMPLEX, DIMENSION(3), INTENT(IN)      :: GParLeft
      COMPLEX, DIMENSION(3), INTENT(IN)      :: GParRight
      TYPE(MomentStorageArea), INTENT(INOUT) :: Storage
      INTEGER, INTENT(IN)                    :: Order

      ! Initialize moment value
      PrimMom = CMPLX(1.0, 0.0)

      ! When zero-th order moment is required, the value is 1.0 and we can stop subroutine
      IF ( Order == 0 ) RETURN

      ! Decide whether the moment need to be actually computed or not
      IF ( Order > HigherOrder .OR. (.NOT. Storage%Check(Order) )) THEN

         ! CALCULATE MOMENT
         PrimMom = GauPrim_Moment_NoStorage( GParLeft, GParRight, Order )

         ! Store its value in case the order is smaller than the maximum
         IF ( Order <= HigherOrder ) THEN
            Storage%Moment( Order ) = PrimMom
            Storage%Check( Order )  = .TRUE.
         END IF

      ELSE

         ! This specific moment has been already computed... retrieve its value from storage
         PrimMom = Storage%Moment( Order )

      END IF

   END FUNCTION GauPrim_Moment_Storage


!*******************************************************************************
!          GauPrim_Moment_NoStorage
!*******************************************************************************
!> Given two gaussian primitive functions, calculate a moment of the product
!> distribution divided by the overlap of the gaussians.
!> In formulas, this function returns:
!>    < g_left(x) | x^n | g_right(x) > / < g_left(x) | g_right(x) >
!> where n specifies the order of the moment.
!> This subroutine uses the formula:
!>                   <gleft|x**n|gright> / <gleft|gright> =
!> \sum_{j=0}^{m=floor(n/2)[(n_over_2j) (p/2)**(n-2j) o**(j-n) \prod_{k=1}^{j} (2j-2k+1)/2]
!> with: o=-(aleft*+aright), p=(xileft*+xiright), q=etaleft*+etaright
!>
!> @param    Left               Left gaussian configuration
!> @param    Right              Right gaussian configuration
!> @param    Order              Order of the moment to compute
!> @returns  the complex value of the moment
!*******************************************************************************
  COMPLEX FUNCTION GauPrim_Moment_NoStorage( Left, Right, Order ) RESULT( GetMoment )
      IMPLICIT NONE
      COMPLEX, DIMENSION(3), INTENT(IN) :: Left
      COMPLEX, DIMENSION(3), INTENT(IN) :: Right
      INTEGER, INTENT(IN)               :: Order

      REAL    :: lo
      COMPLEX :: lp
      INTEGER :: j

      ! When zero-th order moment is required, the value is 1.0 and we can stop subroutine
      IF ( Order == 0 ) THEN
         GetMoment = CMPLX(1.0, 0.0)
         RETURN
      END IF
      
      lo = -( REAL(Left(1))   + REAL(Right(1)) )
      lp =    CONJG( Left(2)) + Right(2)

      SELECT CASE (Order)
      CASE(0)
         GetMoment = CMPLX(1.0, 0.0)
      CASE(1)
         GetMoment = lp / (2.0*lo)
      CASE(2)
         GetMoment = (lp**2 + 2.0*lo)/(4.0*lo**2)
      CASE(3)
         GetMoment = (6.0*lp*lo + lp**3)/(8.0*lo**3)
      CASE(4)
         GetMoment = (lp**4 + 12.0* lp**2 *lo + 12.0 * lo**2)/(16.0*lo**4)
      CASE DEFAULT
         GetMoment = CMPLX(0.0,0.0)
         DO j = 0, FLOOR( REAL(Order) / 2.0 )
            GetMoment = GetMoment + FactCoeff( Order, j ) * lp**(Order-2*j) / lo**(Order-j)
!             part = REAL( Binomial(Order, 2*j) ) * (lp/2.0)**(Order-2*j) * lo**(j-Order)
!             DO k = 1, j
!                part = part*REAL((2.0*REAL(j)-2.0*REAL(k)+1.0))/2.0
!             END DO
!             GetMoment = GetMoment + part
         END DO
         GetMoment = GetMoment / 2.0**Order
      END SELECT

!       o = - CONJG(GParLeft(1)) - GParRight(1)
!       p =   CONJG(GParLeft(2)) + GParRight(2)
!       m = FLOOR( REAL(Order) / 2.0 )
! 
!       PrimMom = CMPLX( 0.0, 0.0 )
! 
!       DO j = 0, m
!          part = REAL( Binomial(Order, 2*j) ) * (p/2.0)**(Order-2*j) * o**(j-Order)
!          DO k = 1, j
!             part = part*REAL((2.0*REAL(j)-2.0*REAL(k)+1.0))/2.0
!          END DO
!          PrimMom = PrimMom + part
!       END DO

   END FUNCTION GauPrim_Moment_NoStorage


! !*******************************************************************************
! !          GauConf_Moment
! !*******************************************************************************
! !> Given two gaussian configurations, calculate a moment of the product
! !> distribution divided by the overlap of the configurations.
! !> In formulas, this function returns:
! !>    < g_left | x_i^n | g_right > / < g_left | g_right >
! !> where i specifies one of the coordinates and n the order of the moment.
! !>
! !> @param    GParLeft           Left gaussian configuration
! !> @param    GParRight          Right gaussian configuration
! !> @param    iCoord             Coordinates of the moment to compute
! !> @param    Order              Order of the moment to compute
! !> @param    NDim               Number of dimensions of the gaussian cfg
! !> @returns  the complex value of the moment
! !*******************************************************************************
!   COMPLEX FUNCTION GauConf_Moment( GParLeft, GParRight, iCoord, Order, NDim )
!       IMPLICIT NONE
!       COMPLEX, DIMENSION(3,NDim), INTENT(IN) :: GParLeft
!       COMPLEX, DIMENSION(3,NDim), INTENT(IN) :: GParRight
!       INTEGER, INTENT(IN) :: iCoord
!       INTEGER, INTENT(IN) :: Order
!       INTEGER, INTENT(IN) :: NDim
! 
!       INTEGER :: m, j, k
!       COMPLEX :: o, p, part
! 
!       ! uses formula
!       ! <gleft|x**n|gright> =<gleft|gright>*\sum_{j=0}^{m=floor(n/2)[(n_over_2j) (p/2)**(n-2j) o**(j-n) \prod_{k=1}^{j} (2j-2k+1)/2]
!       ! with: o=-(aleft*+aright), p=(xileft*+xiright), q=etaleft*+etaright and
!       ! <gleft|gright>=exp(q+p**2/(4o))*dsqrt(pi/o)
! 
!       o=-(CONJG(GParLeft(1,iCoord))+GParRight(1,iCoord))
!       p=CONJG(GParLeft(2,iCoord))+GParRight(2,iCoord)
!       m = FLOOR( REAL(Order) / 2.0 )
!       GauConf_Moment = CMPLX( 0.0, 0.0 )
!       DO j = 0, m
!          part = REAL( Binomial(Order, 2*j) ) * (p/2.0)**(Order-2*j) * o**(j-Order)
!          DO k = 1, j
!             part = part*REAL((2.0*REAL(j)-2.0*REAL(k)+1.0))/2.0
!          END DO
!          GauConf_Moment = GauConf_Moment + part
!       END DO
! 
!    END FUNCTION GauConf_Moment


!    INTEGER FUNCTION factorial(a)
!       !   calculates the factorial of an integer
!       !   input parameter a, returns a! or -1 if a<0
!       IMPLICIT NONE
!       INTEGER, INTENT(IN)  :: a
!       INTEGER :: i
!
!       IF (a < 0) THEN
!          factorial = -1
!       END IF
!       factorial = 1
!       DO i = 1, a
!          factorial = factorial*i
!       END DO
!    END FUNCTION factorial


   REAL FUNCTION FactCoeff( n, k )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n, k
      INTEGER :: i

      FactCoeff = 1
      IF ( k /= 0 ) THEN
         DO i = n, n-2*k+1, -1
            FactCoeff = FactCoeff*REAL(i)
         END DO
      END IF
      IF ( k >= 2 ) THEN
         DO i = 2, k
            FactCoeff = FactCoeff/REAL(i)
         END DO
      END IF
     
   END FUNCTION FactCoeff
   
   INTEGER FUNCTION Binomial( n, k )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n, k
      INTEGER :: ksmall, i

      IF ( k > n/2 .AND. k <= n ) THEN
         ksmall = n - k
      ELSE IF ( k <= n/2 .AND. k >= 0 ) THEN
         ksmall = k
      ELSE
         CALL AbortWithError( "Binomial: wrong integer k = ", ERR_INT_ARITHM, k )
      END IF

      IF ( ksmall == 0 ) THEN
         Binomial = 1
      ELSE IF ( ksmall == 1 ) THEN
         Binomial = n
      ELSE
         Binomial = 1
         DO i = n-ksmall+1, n
            Binomial = Binomial * i
         END DO
         DO i = 2, ksmall
            Binomial = Binomial / i
         END DO
      END IF

   END FUNCTION Binomial



END MODULE GauConf
