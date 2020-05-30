
# ******************** MAKEFILE ****************

#----------------------------------------------------------------------------
#                         USER DEFINABLE OPTIONS
#----------------------------------------------------------------------------

# Executable name
EXENAME = vmcg

# Print relevant output to log file
LOGFILE = yes
LOGNAME = vmcg.log

# Compiler ( gfortran, ifort )
FC = gfortran

# Debugging options ( yes or no )
DEBUG = no

# Optimization level
OPTLEVEL = 3

# OpenMP libraries
OPENMP = no

# linking LAPACK and BLAS
LAPACK = yes

# Compile with standard real 8 (see details about the flags for each compiler...)
REAL8 = yes

#----------------------------------------------------------------------------
#                      LAPACK AND BLAS FINE DETAILS
#----------------------------------------------------------------------------

# Intel compiler version ( used only if FC=ifort and LAPACK=yes )
# 2016-SEQ, 2016-MULTI      -   2016 version, sequential / multithreaded
# 2013-SEQ, 2013-MULTI      -   2013 version, sequential / multithreaded
# 11-SEQ,   11-MULTI        -   11.x version, sequential / multithreaded
# 11-IA32                   -   11.x ia32 arch version, sequential
# 10-SEQ,   10-MULTI        -   10.x version, sequential / multithreaded
INTELVERS = 2016-SEQ

# gfortran lapack libraries
# GNU      - system default libraries
# ATLAS    - atlas libraries
GLAPACK = GNU

#----------------------------------------------------------------------------
#                             STRIP ALL SPACES
#----------------------------------------------------------------------------

# Strip leading and trailing spaces from all variables.
FC := $(strip ${FC})
DEBUG := $(strip ${DEBUG})
OPTLEVEL := $(strip ${OPTLEVEL})
LAPACK := $(strip ${LAPACK})
INTELVERS := $(strip ${INTELVERS})
GLAPACK := $(strip ${GLAPACK})
LOGFILE := $(strip ${LOGFILE})
LOGNAME := $(strip ${LOGNAME})
OPENMP := $(strip ${OPENMP})


#----------------------------------------------------------------------------
#                       Compiler specific statements.
#----------------------------------------------------------------------------

ifeq (${FC},gfortran)

   # Optimization flag
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flag
   DEBUGFLG =  -g -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow -Wall 

   # LAPACK AND BLAS flags
   ifeq (${GLAPACK},GNU)
      # GNU Lapack and Blas flags
      LAPACKFLG = -llapack -lblas
   endif
   ifeq (${GLAPACK},ATLAS)
      # ATLAS Lapack and Blas flags
      LAPACKFLG =  -L/usr/lib64/atlas/ -llapack -lf77blas -lcblas -latlas
   endif

   # OPENMP flags
   OPENMPFLG = -fopenmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
	DATAFLG = -fdefault-real-8
   endif

   # Flag to specify the position of mod files
   MODULEFLG = -fintrinsic-modules-path

   # Other miscellaneous compiling options
   MISCFLG = -ffree-line-length-none
endif

ifeq (${FC},ifort)

   # Optimization flags
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flags
   DEBUGFLG  =  -g -traceback -fpe-all=0 -debug all -check all

   # MKL flags
   ifeq (${INTELVERS},2016-SEQ)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -mkl=sequential
   endif
   ifeq (${INTELVERS},2016-MULTI)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -qopenmp -mkl=parallel
   endif
   ifeq (${INTELVERS},2013-SEQ)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -mkl=sequential
   endif
   ifeq (${INTELVERS},2013-MULTI)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -openmp -mkl=parallel
   endif
   ifeq (${INTELVERS},11-SEQ)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
      LAPACKCOMPILE = -I$(MKLROOT)/include
   endif
   ifeq (${INTELVERS},11-IA32)
      LAPACKFLG = -L${MKLROOT}/lib/ia32 -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -lm
      LAPACKCOMPILE = -I${MKLROOT}/include
   endif
   ifeq (${INTELVERS},11-MULTI)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
      LAPACKCOMPILE = -openmp -I$(MKLROOT)/include
   endif
   ifeq (${INTELVERS},10-SEQ)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
      LAPACKCOMPILE = -I$(MKLROOT)/include
   endif
   ifeq (${INTELVERS},10-MULTI)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
      LAPACKCOMPILE = -openmp -I$(MKLROOT)/include
   endif

   # OPENMP flags
   OPENMPFLG = -qopenmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
	DATAFLG = -r8 -i4
   endif

   # Flag to specify the position of mod files
   MODULEFLG = -I

   # Other miscellaneous compiling options
   MISCFLG =
endif

#----------------------------------------------------------------------------
#              Setup preprocessing options
#----------------------------------------------------------------------------

PPDEFINE =

# Define macro for writing current version (git tag) in the code
PPDEFINE += -DVERSIONTAG=\"$(shell git describe --tags)\"

# Write output to logfile
ifeq (${LOGFILE}, yes)
   PPDEFINE += -DLOG_FILE=\"${LOGNAME}\"
endif

# Preprocess with lapack calls
ifeq (${LAPACK}, yes)
   PPDEFINE += -DWITH_LAPACK
endif

# Preprocess with OPENMP
ifeq (${OPENMP}, yes)
   PPDEFINE += -DWITH_OPENMP
endif


#----------------------------------------------------------------------------
#              Setup linking and compilation flags
#----------------------------------------------------------------------------

# initialize flags
COMPILEFLG =
LINKFLG    =
LIBFLG     =
INCLUDEFLG =

# if debugging set the appropriate flags
ifeq (${DEBUG}, yes)
   COMPILEFLG += ${DEBUGFLG}
endif

# Set flags for defining standard variable kinds
COMPILEFLG += ${DATAFLG}

# If lapack, add the linking options
ifeq (${LAPACK}, yes)
   LIBFLG += ${LAPACKFLG}
   LINKFLG += ${LAPACKCOMPILE}
endif

# If OPENMP, add the linking options
ifeq (${OPENMP}, yes)
   LIBFLG += ${OPENMPFLG}
   COMPILEFLG += ${OPENMPFLG}
endif

COMPILEFLG += ${MISCFLG}

#----------------------------------------------------------------------------
#              Determine the optimization level to be used.
#----------------------------------------------------------------------------

# if debugging override input optimization level
ifeq (${DEBUG}, yes)
   OPTLEVEL = 0
endif

# Define optimization level
OPTFLAGS       = ${O0FLAGS}
ifeq (${OPTLEVEL},1)
  OPTFLAGS       = ${O1FLAGS}
endif
ifeq (${OPTLEVEL},2)
  OPTFLAGS       = ${O2FLAGS}
endif
ifeq (${OPTLEVEL},3)
  OPTFLAGS       = ${O3FLAGS}
endif

COMPILEFLG += ${OPTFLAGS}


#----------------------------------------------------------------------------
#                      List of directories
#----------------------------------------------------------------------------

SRCDIR  = Source
OBJDIR  = Objects
TESTDIR = Tests/Source
EXEDIR  = Executables

#----------------------------------------------------------------------------
#                      List of object files
#----------------------------------------------------------------------------

# Define list of object from the list of all f90 files in the directory
OBJSWITHMAIN =$(patsubst Source/%,Objects/%,$(patsubst %.f90,%.o,$(wildcard ${SRCDIR}/*.f90)))
OBJS =$(patsubst %/Main.o,,${OBJSWITHMAIN})


#----------------------------------------------------------------------------
#       Construct the compile, link, preprocess variables.
#----------------------------------------------------------------------------

# Compile command: ${COMPILE} <source>
COMPILE                 = ${FC} ${COMPILEFLG} ${MODULEFLG} ${OBJDIR} -c

# Link command: ${LINK} <exe name> <objects> <libflags>
LINK                    = ${FC} ${LINKFLG} ${MODULEFLG} ${OBJDIR} -o

# Preprocess commands, add to compilation
PREPROCESS              = -cpp ${PPDEFINE}

# Build static library : ${AR} <libraryname>.a <objects>
AR 			= ar cr

#----------------------------------------------------------------------------
#                       START OF MAKE RULES
#----------------------------------------------------------------------------

# Link objects to the produce the executable file
${EXENAME} : ${SRCDIR}/Main.f90 ${OBJS}
	${COMPILE} ${PREPROCESS} ${SRCDIR}/Main.f90
	${LINK} ${EXEDIR}/$@ Main.o $(OBJS) ${LIBFLG}
	rm Main.o

# Make target to build all the object files and assemble them
all : ${OBJS}

# Make a target object file by preprocessing and compiling the fortran code
${OBJDIR}/%.o : ${SRCDIR}/%.f90
	${COMPILE} ${PREPROCESS} ${SRCDIR}/$*.f90
	cp -p $*.o $(shell echo $* | tr A-Z a-z).mod ${OBJDIR}
	rm $*.o $(shell echo $* | tr A-Z a-z).mod

# Make target to build required directories
directories :
	mkdir -p ${OBJDIR} ${EXEDIR}

# Make documentation with doxygen
doc :
	doxygen Documentation/Doxyfile

# Remove compiled objects and related stuff
clean :
	rm -fr ${OBJDIR}/*

# Clean documentation
clean-doc :
	rm -fr Documentation/html
	rm -fr Documentation/latex

# --------------------------------------------------------------------------------------------
# ---------------------------------- START WITH DEPENDENCIES NOW -----------------------------
# --------------------------------------------------------------------------------------------

# Very basic files, which everything depends on
COMMONDEP = ${OBJDIR}/MyError.o  ${OBJDIR}/MyConsts.o ${OBJDIR}/MyUtils.o ${OBJDIR}/MyLinearAlgebra.o \
            ${OBJDIR}/Clock.o ${SRCDIR}/preprocessoptions.cpp Makefile

# Runge-Kutta 4(5) integrator
${OBJDIR}/RungeKutta45.o : ${SRCDIR}/RungeKutta45.f90 ${OBJDIR}/PsiObject.o ${OBJDIR}/SharedData.o \
                           ${OBJDIR}/OperatorDefine.o ${OBJDIR}/AdaptGBasis.o ${COMMONDEP}

# Input and shared data
${OBJDIR}/SharedData.o : ${SRCDIR}/SharedData.f90 ${OBJDIR}/InputField.o ${COMMONDEP}

# AdaptGBasis: additional methods for the Wavefunction object (basis set on-the-fly adaptation)
${OBJDIR}/AdaptGBasis.o : ${SRCDIR}/AdaptGBasis.f90 ${OBJDIR}/PsiObject.o ${OBJDIR}/GauConf.o ${OBJDIR}/MatrixInversion.o \
                          ${OBJDIR}/OperatorDefine.o ${OBJDIR}/SharedData.o ${COMMONDEP}

# Wavefunction object
${OBJDIR}/PsiObject.o : ${SRCDIR}/PsiObject.f90 ${OBJDIR}/GauConf.o ${OBJDIR}/MatrixInversion.o ${OBJDIR}/OperatorDefine.o \
                        ${OBJDIR}/SharedData.o ${COMMONDEP}

# Definition of operators
${OBJDIR}/OperatorDefine.o : ${SRCDIR}/OperatorDefine.f90 ${COMMONDEP}

# Matrix inversion with regularization
${OBJDIR}/MatrixInversion.o : ${SRCDIR}/MatrixInversion.f90 ${OBJDIR}/SharedData.o ${COMMONDEP}

# General purpose operations with gaussian functions
${OBJDIR}/GauConf.o : ${SRCDIR}/GauConf.f90 ${COMMONDEP}

# Machinery to read from input file
${OBJDIR}/InputField.o : ${SRCDIR}/InputField.f90 ${COMMONDEP}

# Machinery to deal with input/output units
${OBJDIR}/UnitConversion.o : ${SRCDIR}/UnitConversion.f90 ${COMMONDEP}

# Linear algebra subroutines
${OBJDIR}/MyLinearAlgebra.o : ${SRCDIR}/MyLinearAlgebra.f90 ${OBJDIR}/MyError.o ${OBJDIR}/MyConsts.o \
                              ${OBJDIR}/NRUtility.o Makefile

# Utility subroutines and functions of NR
${OBJDIR}/NRUtility.o : ${SRCDIR}/NRUtility.f90 Makefile

# General purpose subroutines
${OBJDIR}/MyUtils.o : ${SRCDIR}/MyUtils.f90 ${OBJDIR}/MyError.o Makefile

# Store timings of different sections of the code
${OBJDIR}/Clock.o : ${SRCDIR}/Clock.f90 ${OBJDIR}/MyError.o Makefile

# Values of physical and mathematical constants
${OBJDIR}/MyConsts.o : ${SRCDIR}/MyConsts.f90 Makefile

# Error and warning procedures
${OBJDIR}/MyError.o : ${SRCDIR}/MyError.f90 Makefile

