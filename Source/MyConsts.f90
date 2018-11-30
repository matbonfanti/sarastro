!***************************************************************************************
!*                              MODULE MyConsts
!***************************************************************************************
!
!>  \brief     Often used mathematical and physical constants
!>  \details   This class defines commonly used mathematical and physical constants,
!>             and the numerical precision of the machine, that can be set dynamically
!>             at runtime with the subroutine Calculate_Constant_EPS().  \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          2.0
!>  \date             June 2012
!>
!***************************************************************************************
!
!>  \pre              Call the subroutine Calculate_Constant_EPS to set the
!>                    machine precision  at runtime
!
!***************************************************************************************
!
!>   \remark         Sources of the constants: \n
!>                   1) The NIST Reference on Constants, Units and
!>                      Uncertainty ( http://physics.nist.gov/cuu/index.html ) \n
!>                   2) The NIST Atomic Weights and Isotopic
!>                      Compositions ( http://www.nist.gov/pml/data/comp.cfm )  \par
!>   \remark         The module is based on Mark Somers' constants.f90
!>   \remark         Generic utility subroutines have been moved out of this
!>                   module to another module MyUtils (12 July 2017)
!
!***************************************************************************************
MODULE MyConsts
   IMPLICIT NONE

   !> \name NUMERICAL PRECISION
   !> Numerical precision of the machine \n
   !> It can be set dynamically through the subroutine Calculate_Constant_EPS()
   !> @{
   REAL               :: MyConsts_EPS      = 1E-12
   !> @}

   !> \name MATHEMATICAL CONSTANTS
   !> Values of frequently used mathematical constants
   !> @{
   REAL, PARAMETER    :: MyConsts_SQRT_2    = SQRT( 2.0 )                                   !< Square Root of 2
   REAL, PARAMETER    :: MyConsts_SQRT_3    = SQRT( 3.0 )                                   !< Square Root of 3
   REAL, PARAMETER    :: MyConsts_INVSQRT_3 = 1./MyConsts_SQRT_3                            !< Inverse of the square root of 3
   REAL, PARAMETER    :: MyConsts_PI        = 3.1415926535897932384626433832795028841971    !< Greek Pi
   COMPLEX, PARAMETER :: MyConsts_I         = CMPLX( 0.0, 1.0 )                             !< Imaginary Unit
   !> @}

   !> \name PHYSICAL CONSTANTS
   !> Values of frequently used physical constants
   !> @{
   REAL, PARAMETER    :: MyConsts_kb = 8.617332478e-5              !<  Boltzmann's constant in eV/K
   REAL, PARAMETER    :: MyConsts_uma = 1.660538921e-27            !<  Atomic Mass Unit (AMU) in kg (from NIST reference)
   REAL, PARAMETER    :: MyConsts_mel = 9.10938291e-31             !<  Electron mass in kg (from NIST reference)
   REAL, PARAMETER    :: MyConsts_NAvo = 6.02214129E23             !<  Avogadro's number (from NIST reference)
   REAL, PARAMETER    :: MyConsts_ThermoCal = 4.184                !<  Thermochemical Calorie ( from NIST reference )
   REAL, PARAMETER    :: MyConsts_AUofTime = 2.418884326502e-2     !<  Atomic unit of time in fs ( from NIST reference )
   REAL, PARAMETER    :: MyConsts_SpeedOfLight = 299792458         !<  Speed of light in the vacuum, in m/s ( from NIST reference )
   !> @}

   !> \name CONVERSION FACTORS
   !> Conversions factor between common units
   !> @{
   REAL, PARAMETER    :: MyConsts_Hartree2eV     = 27.21138505                 !< Conversion factor from Hartree to eV (from NIST reference)
   REAL, PARAMETER    :: MyConsts_Hatree2joule   = 4.35974434e-18              !< Conversion factor from Hartree to Joule (from NIST reference)
   REAL, PARAMETER    :: MyConsts_Hatree2kjmol   = MyConsts_Hatree2joule*MyConsts_NAvo/1000.   !< Conversion factor from Hartree to kJ/mol
   REAL, PARAMETER    :: MyConsts_Hatree2kcalmol = MyConsts_Hatree2kjmol/ MyConsts_ThermoCal   !< Conversion factor from Hartree to kcal/mol
   REAL, PARAMETER    :: MyConsts_Deg2Rad        = MyConsts_PI / 180.0         !< Conversion factor from Decimal degrees to radiants
   REAL, PARAMETER    :: MyConsts_Bohr2Ang       = 0.52917721092               !< Conversion factor from Bohr to Ang (from NIST reference)
   REAL, PARAMETER    :: MyConsts_Uma2Au         = MyConsts_uma / MyConsts_mel !< Conversion factor from UMA to Atomic Units
   REAL, PARAMETER    :: MyConsts_fs2AU          = 1.0 / MyConsts_AUofTime     !< Conversion factor from femtosecond to Atomic Units
   REAL, PARAMETER    :: MyConsts_K2AU           = MyConsts_kb / MyConsts_Hartree2eV   !< Conversion factor from Kelvin to Atomic Units
   REAL, PARAMETER    :: MyConsts_cmmin1tofsmin1 = 2.0*MyConsts_PI*MyConsts_SpeedOfLight*1e-13  !< Conversion from cm-1 to fs-1 (2PI for angular freq conversion)
   REAL, PARAMETER    :: MyConsts_cmmin1toAU     = MyConsts_cmmin1tofsmin1/MyConsts_fs2AU !< Conversion from cm-1 to Atomic Units (freq=energy)
   !> @}

   !> \name ELEMENTS MASSES
   !>  Masses of the elements (both pure and weighted for the isotopic composition) \n
   !>  Values have been taken from NIST:\n http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
   !> @{
   REAL, PARAMETER    :: MyConsts_mH  =  1.00782503207 * MyConsts_Uma2Au    !< Hydrogen atomic mass (from NIST)
   REAL, PARAMETER    :: MyConsts_mD  =  2.0141017778 * MyConsts_Uma2Au     !< Deuterium atomic mass (from NIST)
   REAL, PARAMETER    :: MyConsts_mCMix =  12.01078 * MyConsts_Uma2Au       !< Isotopic mix Carbon atomic mass (from NIST)
   REAL, PARAMETER    :: MyConsts_mCuMix = 63.5463 * MyConsts_Uma2Au        !< Isotopic mix Cupper atomic mass (from NIST)
   !> @}

   !> \name OTHER
   !> Other useful parameters, like special characters or numeric standard types \n
   !> @{
   CHARACTER(len=1), PARAMETER  :: NewLine = achar(10)
   ! The following kind parameters define the standard numeric kind, REGARDLESS of compiler options
   ! and coherently with the kind definition of the compiler
   ! hence SINGLE_PRECISION_KIND will be the kind of a true single precision real data for the compiler used
   ! and analogously for the other parameters
   INTEGER, PARAMETER           :: SINGLE_PRECISION_KIND = SELECTED_REAL_KIND(5,20)
   INTEGER, PARAMETER           :: DOUBLE_PRECISION_KIND = SELECTED_REAL_KIND(10,40)
   INTEGER, PARAMETER           :: LONG_INTEGER_KIND     = SELECTED_INT_KIND( 16 )
   INTEGER, PARAMETER           :: SHORT_INTEGER_KIND    = SELECTED_INT_KIND( 8 )
   !> @}


   CONTAINS

!*******************************************************************************
! Calculate_Constant_EPS
!*******************************************************************************
!>  Calculates the machines precision at runtime and stores
!>  it into the Constant_EPS variable so it can be used on any system
!>  or byte sizes of real types...
!*******************************************************************************
   SUBROUTINE Calculate_Constant_EPS( )
   IMPLICIT NONE
   REAL A

      MyConsts_EPS = 1.0
      DO
         A = 1.00 - MyConsts_EPS
         IF ( A >= 1.00 ) EXIT
         MyConsts_EPS = MyConsts_EPS / 2.00
      END DO

      MyConsts_EPS = ABS( MyConsts_EPS )
      IF( MyConsts_EPS <= 0.0000000 ) MyConsts_EPS = 1.0000E-12   ! added this to be sure ...

   END SUBROUTINE Calculate_Constant_EPS


END MODULE MyConsts

!********************************************* END OF FILE *******************************