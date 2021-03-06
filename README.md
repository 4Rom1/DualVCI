# Dual Vibration Configuration Interaction (DVCI).

##### Author: 
Romain Garnier <rom1{dot}garnier{at}yahoo{dot}fr>.

##### Licensing provisions: 
GPL v3 or later. See <https://www.gnu.org/licenses/>.
    
##### Updated sources can be found at
https://github.com/4Rom1/DualVCI 

## Function :
Compute targeted vibrational spectrum of a wide range of molecules with a very good precision.
 The resources to be allocated directly depend on the number of nucleus 
 and the number of terms of the potential energy surface that should be provided. 
 The algorithm saves memory and computational time with a new factorization described in
  
  **[Dual vibration configuration interaction (DVCI). 
  An efficient factorisation of molecular Hamiltonian for
  high performance infrared spectrum computation,
  Computer Physics Communications, 2018](https://doi.org/10.1016/j.cpc.2018.07.008)** 
  
More precisely, this factorization allows fast computation of Hamiltonian matrix elements *on the fly*
and can be applied to any polynomial operator including derivative terms (e.g Coriolis terms) such as dipole moment, 
an example of implementation is given within the module Transitions [(see under).](#Transitions)
 
------------------------------------------------------------------------

## Languages:
 The code is written in C/C++ and uses ARPACK[6] and LAPACK[7]
 Fortran libraries (provided with the source files).

## Required material:
 This version of the program runs on a Linux machine.
 It should be compiled using gcc compiler,
 but it is possible to change it (e.g with intel compiler).
 It has successfully been compiled by using 
 gfortran and g++ commands of the gcc version 5.4.0.   

------------------------------------------------------------------------   
## Installation:

 To install it, unzip the DualVCI.zip, open a terminal and tape   
 **$ ./INSTALL**

 For a second installation on the same or an other machine, you need to go into the sub directory 
 ARPACK, open the file ARmake.inc and be sure that home directory writes "home = MyDir".
 Note that "home = MyDir" should be already written for the first installation, 
 and then overwritten by your current directory after you launch the command ./INSTALL. 

 If the installation is not running then make the installation file executable by typing   
 **$ sudo chmod 777 INSTALL**

 During the installation a link between the command DVCI and 
 the executable is created in ~/.bashrc.

 So if you open ~/.bashrc you'll see the line   
 alias DVCI='MyDir/DualVCI/source/./DVCI'   
 where MyDir is the directory you opened DualVCI.zip.  

----------------------------------------------------------------------------------  
## Execution:

 The program runs on a terminal with the command   
 **$ DVCI INPUT.in**  
 Or if for several applications you need to include the path of the executable   
 **$ MyDir/DualVCI/source/DVCI INPUT.in**  
 where INPUT.in contains the different keywords discussed in the user manual (DVCI_Manual.pdf).  
 If DVCI command is not recognized on your terminal, 
 then the bashrc file needs to be sourced with command   
 **$ source ~/.bashrc**

 The directory where is executed DVCI should also contain the potential energy file
 as explained in the user manual.  

 The informations of the frequencies, residues, assignment and CPU time will
 be displayed in output. Additionally you'll see 2 other files with respective extension 
 *-FinalBasis.in*, *-Vectors.bin* (if Printout=1).  
 The first one print all configurations of final basis set and the second one 
  the vectors components. 

## Transitions
 The module **Transitions** allows to get the transition moments <Psi_0 |Op|Psi_n> where
  Op is an operator having the same format as the PES and Phi_n are wave functions of
  the targets previously computed with DVCI (Psi_0 is the ground state).
 To run Transitions:  
 **$ Transitions INPUT-Op.in**  
  Where Input-Op.in should be the same file as Input.in except that PESName needs to be changed
  with the name of the file containing the operator Op 
  (one can also take the same PES name to do some verification as in the unitary tests). 
 If PrintOut=2 the matrix of the last iteration will be saved in binary format with extension
  *-Matrix.bin* to re-run the last diagonalization in the final basis set.
 To do so, you just need to use exactly the same input as DVCI and hit   
 **$ FinalVCI INPUT.in**

---------------------------------------------------------------------------------------    
## Test cases:

 If you go to the directory Molecules, you'll see subdirectories with different names of molecule.
 Each one of these directories contains necessary files to run directly. 
 You also will find a txt file named Inputs.txt consisting of the inputs
 used for the results shown in the manual.

 For N2H2 and C2H4 The potential energy files are extracted for PyPES librairy [1].

 [1] M. Sibaev, D. L. Crittenden, The PyPES library of high quality semi-global potential energy
 surfaces, Journal of Computational Chemistry 36 (29) (2015) 2200–2207.


 For C3H3NO, potential energy surface is generated with methods described in 
 [2] Sparta et al. JCTC, 6(10):3162-3175, 2010, 
 [3] Sparta et al. (J. Phys. Chem. A (2009), 113, 8712-8723), 
 [4] Sparta et al. in Theor. Chem. Acc. (2009), 123, 413.

 The force constants, equilibrium geometry and normal coordinates were taken from 

 [5] Madsen et al. "Efficient algorithms for solving the non-linear vibrational coupled-cluster
 equations using full and decomposed tensors." 
 The Journal of Chemical Physics 146(13) (2017) 134110.

 A converter script written in python is provided with the original files.
 To run it tape   
 **$./Convert.py oxazole_pes.mop 18**

 Unit tests are now available in directory Molecules with the 2 shell scripts 
 Test.sh and TestMemCheck.sh.
 The first script is a series of different combination of parameters,
 as the second one additionally containing memory checks with Valgrind.
 
---------------------------------------------------------------------------------------

## Source files

 The source files (+ corresponding headers) are the following:  
  
 DVCI.cc        -> Main program.  
 Shared.cc      -> Zeta coefficients[8] and initialization subroutines.  
 Assemblage.cc  -> Subroutines for assemblage of Hamiltonian matrix.  
 Generator.cc   -> Subroutines generating excitations and computing residual vectors on the fly.    
 Graph.cc       -> Sorting routines, initial subspace construction linear search and neighbors generator.    
 ReadData.cc    -> Reading input and potential energy files.  
 Basic.h        -> Basic subroutines.  
 Transitions.cc -> Transitions <Phi_0 |Op|Phi_n> where Phi_n are targets previously computed with DVCI   
 FinalVCI.cc    -> Re-run last iteration of DVCI when PrintOut=2   
  
## Additional sources
 
For the eigensolver, in file Solver2.cc the program uses
 DSAUPD, DSEUPD subroutines of the software ARPACK[6]  
  
 [6] R. B. Lehoucq, D. C. Sorensen, C. Yang, ARPACK Users ’ Guide : Solution of Large Scale
 Eigenvalue Problems with Implicitly Restarted Arnoldi Methods, Communication 6 (1998)
 147. doi:10.1137/1.9780898719628.  
  
 Additionally to the original source provided, a description of DSAUPD can be found [here](http://www.caam.rice.edu/software/ARPACK/UG/node136.html)   
  
--------------------------------------------------------------------  
  
In routine MomentInertie of file Shared.cc,  
 the computation of principal axes of inertia momentum requires to diagonalize a matrix
 with DSYEV subroutine of LAPACK librairy[7].  
   
 [7] E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra, J. Du Croz,
 A. Greenbaum, S. Hammarling, A. McKenney, D. Sorensen, LAPACK Users’ Guide,
 3rd Edition, Society for Industrial and Applied Mathematics, Philadelphia, PA, 1999.  
  
 Additionally to the original source provided, a description of DSYEV can be found [here](http://www.netlib.org/lapack/explore-3.1.1-html/dsyev.f.html)  

In routine Zeta,the rotational coefficients are computed according to the method of Meal and Polo [8]  

[8]Vibration—Rotation Interaction in Polyatomic Molecules. II. The Determination of Coriolis Coupling Coefficients
Janet Hawkins Meal and S. R. Polo, The Journal of Chemical Physics 1956 24:6, 1126-1133   

-------------------------------------------------------------------

The program uses the print banner routine taken from  
 www.cise.ufl.edu/class/cop4600/minix/src/commands/simple/banner.c  
 (By B.Wallis, 4 July 1988)
  
------------------------------------------------
  
In the file Graph.cc, Routines QuickSortMode, QuickSortModeShift are made
 from the Quicksort algorithm [9].  
  
 [9] C. A. R. Hoare. 1961. Algorithm 64: Quicksort. Commun. ACM 4, 7 (July 1961),   
 321-. DOI=http://dx.doi.org/10.1145/366622.366644   
  
 In the file Graph.cc, routines QsearchMode, QsearchModeShift are made 
 from half interval search algorithm also known as binary search [10].  
  
 [10] Willams, Jr., Louis F. (1975). A modification to the half-interval search 
 (binary search) method. Proceedings of the 14th ACM Southeast
 Conference. pp. 95–101. doi:10.1145/503561.503582.

------------------------------------------------

In Graph.cc, the routines nested_loop and nested_loop1 are depth variable nested loops sharing
 the same structure as the one presented in STACKOVERFLOW [11]  

 [11] How to set a variable for nested loops?   
 Posted on stackoverflow (by M Oehm) at the following address  
 http://stackoverflow.com/questions/25785863/how-to-set-a-variable-for-nested-loops  

