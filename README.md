# sarastro

*sarastro* is a FOTRAN 90/95 code to compute the solution
of the Time-Dependent Schr&ouml;dinger Equation (TDSE) with the
variational multiconfiguration Gaussian wavepacket (vMCG),
descring the wavefunction as a combination of multi-dimensional
Gaussian Wavepackets (GWPs).
The present implementation allows to integrate the TDSE for
simple analytical sum-over-products potential with 
polynomial terms.  This code is meant as a testing ground 
for alternative integration schemes for the dynamical vMCG 
equations.

## License

Copyright (c) 2019, AK Burghardt, Goethe-Universität Frankfurt am Main

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Installation

*sarastro* has been successfully compiled on Linux systems
with both the GNU and Intel FORTRAN compilers.
Here instructions are reported to install *sarastro*
with the GNU compiler gfortran and with standard Linux OS 
libraries.

### Requisites (for the *sarastro* main program):
* *gfortran* compiler
* [blas](http://www.netlib.org/blas/) and 
[lapack](http://www.netlib.org/lapack/) numerical libraries

the compiler and the libraries are usually available
as standard packages in most Linux distributions.
E.g., in Ubuntu 20.04 LTS, they can be installed with

    sudo apt install gfortran liblapack3 libblas3

### Requisites (for the analysis tools):
* Python 2.7 or 3.x
* *scipy* (~= 1.3.3)
* *numpy* (~= 1.17.4)
* *matplotlib* (~= 3.1.2)

### Quick installation

Choose the directory where you want to install *sarastro*.
In the following, this directory will be indicated as `$INSTALLDIR`.

Using the shell, move to the installation dir

    cd $INSTALLDIR
    
Clone the git repository of the code:

    git clone https://github.com/matbonfanti/sarastro.git
    
Now, you have to build the code:

    cd sarastro
    make
    
Upon sucessfull compilation, the `vmcg` executable will be
created in the `Executables` directory.

If you want to invoke the program and analysis tools
commands without path, you can add the directories
`$INSTALLDIR/sarastro/Executables` and `$INSTALLDIR/sarastro/Tools`
to the appropriate shell configuration file.
E.g. for `bash` you can add the following lines to 
the `$HOME/.bashrc` file:

    SARASTROPATH=$INSTALLDIR/sarastro
    export PATH=$PATH:$SARASTROPATH/Executables:$SARASTROPATH/Tools
    
### Different installation options

The procedure described above will install *sarastro*
with the GNU compiler and the standard OS blas and lapack 
libraries.

Different options can be specified by editing the file `Makefile`.
In particular, one can use the Intel compiler by changing
the value of the variable `FC`:

    FC = ifort
    
In this case, one should set the appropriate compiling and
linking commands for the MKL libraries, that will substitute 
the standard OS blas and lapack options.
The MKL options can be controlled by using the variable
 `INTELVERS`, and are set at line 107 of the Makefile.
 
## Testing

Two simple tests are available in the `Examples` directory.
The two subdirectories `harmonic1D` and `anharmonic1D` contain
 sample input files to solve the vMCG equations for a 
 simple 1D oscillator with an harmonic and anharmonic potential.
  
Move to the `harmonic1D` directory:

    cd $INSTALLDIR/sarastro/Examples/harmonic1D
    
Execute the `vmcg` code with the following command

    $INSTALLDIR/sarastro/Executables/vmcg vmcg.inp > vmcg.log

A file `vmcg.log` should be created, containing a list 
of information regarding the integration time, up to the 
final time of 500.0 femtoseconds.
A number of output files is created after a 
successfull run. One of the most important output files
is `expectations.dat`, that contains the norm, 
the energy and the autocorrelation function of the 
wavefunction over time. 

Now repeat the procedure for running the calculation 
with the input files contained in the  `anharmonic1D`
directory:

    cd $INSTALLDIR/sarastro/Examples/anharmonic1D
    $INSTALLDIR/sarastro/Executables/vmcg vmcg.inp > vmcg.log

This calculation will last a bit longer, as the
anharmonicity of the potential implies that 
a larger number of basis function is required to 
have a good description of the wavefunction dynamics 
(as a rule of thumb, the calculation ran for ~30 s 
on a machine with Intel Core i7-5500U CPU @ 2.40GHz).
 
Now we can visualize the difference between the dynamics
in the harmonic and anharmonic potential, using the
`AutoCorrelationFourierTransform.py` tool.
This script plots the autocorrelation function and its
Fourier transform on screen.

    cd $INSTALLDIR/sarastro/Examples
    $INSTALLDIR/sarastro/Tools/AutoCorrelationFourierTransform.py harmonic1D anharmonic1D 

The bottom panel of the graph represents the 
**autocorrelation of the wavefunction over time**.
The recurrences of this function describe the 
oscillatory behaviour of the system, and is a simple
sinusoidal function for the harmonic potential while
is more structured in the case of anharmonicities.
The decay of the function over a long ~ps time scale
is due to a exponential damping that is applied to 
the data to introduce a broadening in the Fourier 
transform.

The top panel is the Fourier transform of the function,
plotted as a function of frequency. It represents the
**energy spectrum of the system** and it is discretized
according the quantization of the 1D oscillators.

## Running *sarastro* 

To run, *sarastro* requires a main input file 
to specify the setup of the calculation and 
a few additional files that specify the 
 Hamiltonian operator and the initial conditions
of the basis functions.

The name of the **main input file** is given as 
first command line argument when running the 
`vmcg` executable. This input file is made of a
series of directives of the type

    keyword : value
    
The keywords are listed in the following table:

| keyword             |                                            |
|-------------------- |--------------------------------------------|
| `n`                 | number of Gaussian basis functions         |    
| `gdim`              | number of Gaussian dimensions              |    
| `hamiltonian`       | name of the Hamiltonian file               |    
| `read_avec`         | name of the initial A vector file          |    
| `read_gaussians`    | name of the initial GWP parameter file     |  
| `maxtime`           | maximum time of propagation (fs)           |    
| `maxstep`           | maximum integration time step (fs)         |    
| `minstep`           | minimum integration time step (fs)         |    
| `writetime`         | time step for writing output (fs)          |    
| `use_fixed_steps`   | use fixed time step                        |    
| `tolerance`         | error tolerance for adaptive timestep      |    
| `cmatrix_reg_type`  | type of regularization of the C matrix     |    
| `ovmatrix_reg_type` | type of regularization of the overlap      |    
| `cmatrix_eps`       | regularization eps for the C matrix        |    
| `ovmatrix_eps`      | regularization eps for the overlap         |    
| `freeze_lowpopgau`  | freeze GWP with too small/large population |    
| `popmin_threshold`  | threshold for small populations            |    
| `popmax_threshold`  | threshold for large populations            |    
| `eom_type`          | type of the vMCG equations of motion       |    
| `eom_simplifycoeff` | simplify A coeffs in the eq. of motion     |    

The names of the file for the Hamiltonian and initial GWPs 
are specified by the keywords `hamiltonian`, `read_avec` and `read_gaussians`
in the main input files. The format of these files has been inherited by 
the code *Sorbet* by Pierre Eisenbrandt. 

The **hamiltonian file** contains a sequence of lines, each representing
 an additive term of the Hamiltonian. Each line contains few comma-separated
elements:

    <i>, <coeff>, <type>, <order1>, <order2>, ...
    
Where:

* `<i>` is an ordinal number labelling the line (starting from 1)
* `<coeff>` is a floating point number that defines the coeffient 
  of the Hamiltonian term
* `<type>` is the type of Hamiltonian term, 1 for a kinetic energy term
  and 0 for a potential term
* `<order1>`, <order2>, ... should be a sequence of integer numbers,
  one for each degree of freedom of the system. The integers should all
  be 1 for a kinetic energy term, while in the case of a potential term
  they should specify the polynomial order in each of the variable. 
  should be equal to the polynomial order 

The **initial A vector file** should contain a sequence of lines, 
one for each GWP in the basis set of the system (`n`). Each line 
specify the initial amplitudes of the expansion of the wavefunction
over the GWP basis set. 
Each line contains three elements:

    <i>, <A coefficient real>, <A coefficient imag>
    
Where:

* `<i>` is an ordinal number labelling the line (starting from 1)
* `<A coefficient real>` is a floating point number with the 
  real part of the amplitude
* `<A coefficient imag>` is a floating point number with the 
  imaginary part of the amplitude
  
The **initial GWP parameter file** should contain `n` sequences
of `gdim` lines. Each section of `gdim` lines specifies the
GWP parameters for a single Gaussian functions. Then this 
is reiterated to describe the whole set of GWP functions.
Each line specify a set of one-dimensional GWP parameters

    <i>, <a real>, <a imag>, <xi real>, <xi imag>, <eta real>, <eta imag>
    
Where:

* `<i>` is an ordinal number labelling the line (starting from 1)
* `<a real>` and `<a imag>` are the real and imaginary part of 
  the *a* parameter, the coefficient of the quadratic term in 
  the Gaussian exponent
* `<xi real>` and `<xi imag>` are the real and imaginary part of 
  the *xi* parameter, the coefficient of the linear term in 
  the Gaussian exponent
* `<eta real>` and `<eta imag>` are the real and imaginary part of 
  the *eta* parameter, the coefficient of the zer-th order term in 
  the Gaussian exponent
  
