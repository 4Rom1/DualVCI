/*
    Dual Vibration Configuration Interaction (DVCI).
    A novel factorisation of molecular Hamiltonian for
    high performance infrared spectrum computation. 
    Copyright (C) 2018  Romain Garnier

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Updated sources can be found at 
    https://github.com/4Rom1/DualVCI 
*/
/*--------Additional sources-------------------

This part of the program also uses reverse communication interface subroutines 
DSAUPD, DSEUPD of the software ARPACK[1]

[1] R. B. Lehoucq, D. C. Sorensen, C. Yang, ARPACK Users ’ Guide : Solution of Large Scale
Eigenvalue Problems with Implicitly Restarted Arnoldi Methods, Communication 6 (1998)
147. doi:10.1137/1.9780898719628.
*/
#ifndef Solver2_H
#define Solver2_H

#include "Shared.h"

#define Tol0 1e-20
#define Sep 1e4
#define MaxIterDav 1000000

extern "C"
{
//Arpack subroutines,
//for more details see http://www.caam.rice.edu/software/ARPACK/UG/node136.html.
 void dsaupd_(int *ido, char *bmat, int *n, const char *which, int *nev, double *tol, double *resid,\
             int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl,\
             int *lworkl, int *info );
//
 void dseupd_(int *rvec, char *howmny, int *select, double *d, double* z, int *ldz, double *sigma,\
        char *bmat, int *n, const char *which, int *nev, double *tol, double *resid, int *ncv,\
        double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl,\
        int *lworkl, int *ierr);
//
}
//
int ConvertCharToFortran(char *StringInC, char* StringInF);
//Change the termination of string in c ('\0') by a blank for fortran equivalence
//
void SparseSymCSC(double *X,double *Y,CSC IJ, double *Val,SizeArray Size);
//Sparse symmetric matrix vector product: Y=H*X using CSC format.
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//Val[nn]=NNZ entries of matrix H.
//
uint32_t ShiftDiagCSC(SizeArray Size, CSC IJ, double *ValAct, double Shift);
//Shift the diagonal elements of active matrix Hb.
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//Digonal elements of ValAct[nn] are such that IJ.I[nn]=jj for column jj in [0,Size.DimAct[.
//
int SolverCSCSym(SizeArray *Size, int Iteration, int NEV, int NCV, double Shift, double tol,\
double *EigVal, double *EigVec, double *workd, double *workl,\
 CSC IJ, double *ValAct);
//Compute the eigenvalues and eigenvectors of matrix Hb stored in CSC format (IJ,ValAct),
//and store eigenvectors using column major order. 
//
#endif
