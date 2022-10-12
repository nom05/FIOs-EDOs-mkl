# FIOs-EDOs v2022.1
A Program to Calculate Electron Transport and Electric Response Properties from Electron Deformation Orbitals

## Authors
  - Nicolás Otero Martínez  -         nom05 (at) uvigo.es - University of Vigo
  - Marcos Mandado Alonso   -       mandado (at) uvigo.es - University of Vigo

## Collaborators
  - Nicolás Ramos Berdullas - nicolas.ramos (at) uvigo.es - University of Vigo
  - Sara Gil Guerrero       -           sgg (at) uvigo.es - University of Vigo

## How to use it
Legend for opts.:
     [CM] = common; [ED] = EDOs-only; [CO] = EDOs+Cond.; [FI] = FIOs-only; **\*\* = Mandatory**

Available options:

  -h  --help                        Print help.                                         [CM]
  
  -a  --calc-dmo                    alpha DMO will be computed instead of being read.   [FI]
  
  -C  --conduc  <arg>               Computation of EDOs + conductance.                  [CO]\*\*
  
  -c  --calc  <arg>                 Choose components to be calculated (FIOs).          [FI]
  
  -D  --dynamic                     Calculation of Dynamic alpha FIOs.                  [FI]
  
  -d  --debug                       Print DEBUG messages.                               [CM]
  
  -E  --edo  <arg>                  Computation of EDOs. Include 2nd fchk file.         [ED]\*\*
  
  -e  --estimate                    Estimate memory only. Stop before computation.      [CM]
  
  -F  --fchk  <arg>                 fchk file to be employed w/ or w/o extension.       [CM]\*\*
  
  -f  --field  <arg>                Set field by hand.                                  [FI]
  
  -G  --gau  <arg>                  G09 log file to be employed w/ or w/o extension.    [FI]\*\*
  
  -g  --mkgnuplot                   Make Gnuplot scripts to represent props vs # MO.    [FI]
  
  -i  --print-mat                   Print conductance matrices.                         [CO]
  
  -m  --mem-info                    Print memory info for each array (de)allocation.    [CM]
  
  -N  --nmkl  <arg>                 Number of threads for Intel(R) MKL libraries.       [CM]
  
  -n  --no-color                    Do not use colors on screen.                        [CM]
  
  -O  --print-occ                   Print def occup./MO contrib. to EDOs/FIOs,respectiv.[CM]
  
  -o  --output  <arg>               Change name of the output file (w/o extension)      [CM]
  
  -P  --proc  <arg>                 Number of threads for OpenMP.                       [CM]
  
  -p  --print-mat                   Print large matrices.                               [CM]
  
  -r  --read-mulm                   Read multipole matrices from G09 log file.          [FI]
  
  -S  --skip-ffch                   Skip FIOs/EDOs fchk files writing                   [CM]
  
  -s  --stealth                     Print on screen errors and warnings only.           [CM]
  
  -T  --trunc-dmo  <arg>            Truncate DMO with a set of MOs.                     [CM]
  
  -t  --test                        Perform some tests with DM derivatives and dip. mat.[FI]
  
  -W  --write-pop                   Write deformation orbital occupations and stop.     [FI]
  
  -w  --overwrite                   Overwrite output files.                             [CM]
  
Conductance input format:
  
   - "f=file"       : 2nd fchk file.
   - "d=real_number": electrode distance in A.
   - "e1=1-3,5"     : set of atoms for electrode 1.
   - "e2=12-20,22"  : set of atoms for electrode 2.
   - "print"        : print conductance matrices.
   - Arguments separated by colon(:) and all the text in quotes.

Computational details:
   * FIOs: fchk and G09 log files are mandatory.
   * EDOs: both fchk files are mandatory.
   * "#p"        -> to obtain correct calculation type (mandatory).
   * iop33(3=1)  -> to read multipole matrices from G09 log file (optional).
   * iop33(10=2) -> mandatory writing of density matrices derivatives.
   * field=read  -> mandatory field to perform numerical derivative of alpha.
   - We recommend: "integral=ultrafinegrid" for DFT calculations.
                   "scf=(conver=11,xqc)" for HF/DFT calculations
## GPLv3 license
See corresponding file called [`LICENSE`](LICENSE) for more details.

## Compile the code
We will include an extended manual in the following versions with more detailed information.

Dependencies: Intel(R) MKL libraries
We recommend to use Intel(R) MKL compiler (ifort) to use MKL libraries transparently.

After cloning the code of the repository:

```bash
 $ mkdir build && cd build/
 $ cmake ..
 $ make
```
  
To force CMake to compile the code with your favourite compiler, prepend the cmake line with e.g. FC=gfortran:

```bash
 $ FC=gfortran cmake ..
```
  
or use the explicit option:

```bash
 $ cmake -DCMAKE\_Fortran\_COMPILER=gfortran ..
```
