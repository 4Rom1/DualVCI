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


--------Additional sources-------------------

This part of the program also uses reverse communication interfaces subroutines 
DSAUPD, DSEUPD of the software ARPACK[1]

[1] R. B. Lehoucq, D. C. Sorensen, C. Yang, ARPACK Users ’ Guide : Solution of Large Scale
Eigenvalue Problems with Implicitly Restarted Arnoldi Methods, Communication 6 (1998)
147. doi:10.1137/1.9780898719628.

Additionally to the original source provided, a description of DSAUPD can be found here 
http://www.caam.rice.edu/software/ARPACK/UG/node136.html.
*/


#include "Solver2.h"
#include "Basic.h"
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cassert>
#include <new>          

//Call Arpack with c/fortran interface.

int SolverCSCSym(SizeArray *Size, int Iteration, int NEV, int NCV, double Shift, double tol,\
double *EigVal, double *EigVec, double *workd, double *workl, \
 CSC IJ, double *ValAct)
{
//Compute the eigenvalues and eigenvectors of matrix Hb stored in CSC format (IJ,ValAct),
//and store eigenvectors using column major order. 
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//Size[Iteration+1].NNZAct is the number of NNZ of active matrix Hb at current iteration.
//Size[Iteration+1].DimAct is the size of active space B at current iteration.
//ValAct is an array of size Size[Iteration+1].NNZAct containing NNZ of active matrix Hb 
//The matrix will be shifted as Hb=Hb-Shift*Id
//workd, workl, tol, NEV, NCV refer to the same variables defined in ARPACK
//see the header file solver2.h or http://www.caam.rice.edu/software/ARPACK/UG/node136.html
//for more details.
//
//-----------------------------------------------------------------------
//
//  %-----------------------------%
//  | Define leading dimensions   |
//  | for all arrays.             |
//  | MAXN:   Maximum dimension   |
//  |         of the A allowed.   |
//  | MAXNEV: Maximum NEV allowed |
//  | MAXNCV: Maximum NCV allowed |
//  %-----------------------------%
//  Here NEV is equal to MAXNEV
//   and NCV is equal to MAXNCV
//
      int             Ndim;
//
//  %--------------%
//  | Local Arrays |
//  %--------------%
//
       Ndim=Size[Iteration+1].DimAct;
//       
       double *resid=NULL;//Can be used as a workspace
       resid = new double[Ndim];
       InitTabInt(resid,Ndim)
//
      int *select=NULL;
      select = new int [NCV];
//
      int iparam[12]={0};
      int ipntr[12]={0};  
//  %---------------%
//  | Local Scalars |
//  %---------------%
      char bmat='I';
      char HOWMNY='A';
      char whichc[3]="LM";
      char whichf[3]={' '};
//      
      ConvertCharToFortran(whichc, whichf);
//
      int          ido=0, lworkl=0, info=0, ierr=0,\
                   nconv=0, maxitr=0, mode=0, ishfts=0, Iterate=0;
      int          rvec=0;
//     
      double       sigma=0;    
//      
      double        ShiftLoc=Shift;
//     
//  %-----------------------%
//  | Executable Statements |
//  %-----------------------%
//
//  %----------------------------------------------------%
//  | The number NX is the number of interior points     |
//  | in the discretization of the 2-dimensional         |
//  | Laplacian on the unit square with zero Dirichlet   |
//  | boundary condition.  The number N(=NX*NX) is the   |
//  | dimension of the matrix.  A standard eigenvalue    |
//  | problem is solved (BMAT = 'I'). NEV is the number  |
//  | of eigenvalues to be approximated.  The user can   |
//  | modify NEV, NCV, WHICH to solve problems of        |
//  | different sizes, and to get different parts of the |
//  | spectrum.  However, The following conditions must  |
//  | be satisfied:                                      |
//  |                   N <= MAXN,                       | 
//  |                 NEV <= MAXNEV,                     |
//  |             NEV + 1 <= NCV <= MAXNCV               | 
//  %----------------------------------------------------% 
//
//
//  %---------------------------------------------------%
//  | This program uses exact shifts with respect to    |
//  | the current Hessenberg matrix (IPARAM(1) = 1).    |
//  | IPARAM(3) specifies the maximum number of Arnoldi |
//  | iterations allowed.  Mode 1 of DSAUPD is used     |
//  | (IPARAM(7) = 1).  All these options may be        |
//  | changed by the user. For details, see the         |
//  | documentation in DSAUPD.                          |
//  %---------------------------------------------------%
//  
      ishfts = 1;
      maxitr = 3000;
      mode   = 1;
//   
      iparam[1] = ishfts; 
      iparam[3] = maxitr;
      iparam[7] = mode; 
//
// 
 if (!ShiftDiagCSC(Size[Iteration+1], IJ, ValAct, ShiftLoc))
{
     printf("\n!!!!!!!!IndexDiag mal assignee !!!!!!!!!! \n \n \n");
     return 0;
} 
//
//  %--------------------------------------------------%
//  | The work array WORKL is used in DSAUPD as        |
//  | workspace.  Its dimension LWORKL is set as       |
//  | illustrated below.  The parameter TOL determines |
//  | the stopping criterion.  If TOL<=0, machine      |
//  | precision is used.  The variable IDO is used for |
//  | reverse communication and is initially set to 0. |
//  | Setting INFO=0 indicates that a random vector is |
//  | generated in DSAUPD to start the Arnoldi         |
//  | iteration.                                       |
//  %--------------------------------------------------%
//
      lworkl = NCV*(NCV+8);
      info = 0;
      ido = 0;
      Iterate=1;
//  %-------------------------------------------%
//  | M A I N   L O O P (Reverse communication) |
//  %-------------------------------------------%   
       do{
//
//     %---------------------------------------------%
//     | Repeatedly call the routine DSAUPD and take | 
//     | actions indicated by parameter IDO until    |
//     | either convergence is indicated or maxitr   |
//     | has been exceeded.                          |
//     %---------------------------------------------%
//
  /*Fortran Call:

         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )*/

//Shifted pointer to 1 to have correspondance between c and fortran call.
// 
     dsaupd_( &ido, &bmat, &Ndim, whichf, &NEV, &tol, resid,\
             &NCV, &EigVec[1], &Ndim, &iparam[1], &ipntr[1], &workd[1], &workl[1],\
             &lworkl, &info);
//
         if ( (ido == -1) || (ido == 1) ) {
//
//        %--------------------------------------%
//        | Perform matrix vector multiplication |
//        |              y <--- OP*x             |
//        | The user should supply his/her own   |
//        | matrix vector multiplication routine |
//        | here that takes workd(ipntr(1)) as   |
//        | the input, and return the result to  |
//        | workd(ipntr(2)).                     |
//        %--------------------------------------%
//
     SparseSymCSC(&workd[ipntr[1]],&workd[ipntr[2]],IJ,ValAct,Size[Iteration+1]);
//
//        %-----------------------------------------%
//        | L O O P   B A C K to call DSAUPD again. |
//        %-----------------------------------------%

          } 
         else{Iterate=0;}
// 
   }while(Iterate);
//  %----------------------------------------%
//  | Either we have convergence or there is |
//  | an error.                              |
//  %----------------------------------------%
//
      if ( info < 0 )
      {
//     %--------------------------%
//     | Error message. Check the |
//     | documentation in DSAUPD. |
//     %--------------------------%
         printf("\n \n \n Error with _saupd, info = %i \n \n ",info);
       }
      else {
//
//     %-------------------------------------------%
//     | No fatal errors occurred.                 |
//     | Post-Process using DSEUPD.                |
//     |                                           |
//     | Computed eigenvalues may be extracted.    |  
//     |                                           |
//     | Eigenvectors may also be computed now if  |
//     | desired.  (indicated by rvec = .true.)    | 
//     %-------------------------------------------%
//        
         rvec = 1;
//       HOWMNY='A';
         sigma=0;
//Fortran call:
//
/*         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &        bmat, n, which, nev, tol, resid, NCV, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )*/
//              
//            
               dseupd_( &rvec, &HOWMNY, select, &EigVal[0], &EigVec[0], &Ndim, &sigma,\
                   &bmat, &Ndim, whichf, &NEV, &tol, &resid[0], &NCV, &EigVec[1], &Ndim,\
                   &iparam[1], &ipntr[1], &workd[1], &workl[1], &lworkl, &ierr);
//
//     %----------------------------------------------%
//     | Eigenvalues are returned in the first column |
//     | of the two dimensional array D and the       |
//     | corresponding eigenvectors are returned in   |
//     | the first NEV columns of the two dimensional |
//     | array V if requested.  Otherwise, an         |
//     | orthogonal basis for the invariant subspace  |
//     | corresponding to the eigenvalues in D is     |
//     | returned in V.                               |
//     %----------------------------------------------%
//
         if ( ierr != 0) {
//
//         %------------------------------------%
//         | Error condition:                   |
//         | Check the documentation of DSEUPD. |
//         %------------------------------------%
//
             printf("\n \n \n Error with _seupd, info =  %i \n \n ",ierr);
                         }
//Assign eigenvalues and eigenvectors 
         else 
            {
            nconv =  iparam[5];
            for(int jj=0; jj<nconv;jj++)
              {
            EigVal[jj]+=ShiftLoc;
              }
//Give back the Shift for all the diagonal entries
            ShiftDiagCSC(Size[Iteration+1], IJ, ValAct, -ShiftLoc);           
           }
         }
//
      delete [] resid;
      delete [] select;
      return nconv;
}
//
void SparseSymCSC(double *X,double *Y,CSC IJ, double *Val,SizeArray Size)
{
//Sparse symmetric matrix vector product: Y=H*X using CSC format.
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//Val[nn]=NNZ entries of matrix H.
for (uint32_t nn=0;nn<Size.DimAct;nn++){Y[nn]=0;}
for (uint32_t jj=0;jj<Size.DimAct;jj++)
{
for (uint64_t nn=IJ.NJ[jj];nn<IJ.NJ[jj+1];nn++)
 {
  Y[IJ.I[nn]]+=X[jj]*Val[nn];
  if(IJ.I[nn]-jj)
   {
  Y[jj]+=X[IJ.I[nn]]*Val[nn];
   }
  }
 }
}


uint32_t ShiftDiagCSC(SizeArray Size, CSC IJ, double *ValAct, double Shift)
{
//Shift the diagonal elements of active matrix Hb.
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//Digonal elements of ValAct[nn] are such that IJ.I[nn]=jj for column jj in [0,Size.DimAct[.
uint32_t Cmpt1=0;
for (uint32_t jj=0;jj<Size.DimAct;jj++)
{
for (uint64_t nn=IJ.NJ[jj];nn<IJ.NJ[jj+1];nn++)
 {
  if(!(IJ.I[nn]-jj))
   {
  ValAct[nn]-=Shift;
  Cmpt1++;
   }
  }
}
if (Cmpt1==Size.DimAct){return Cmpt1;}
else
 {
printf("Cmpt1 : %u, Size[Iteration+1].DimAct : %u \n",Cmpt1,Size.DimAct);
return 0;
 }
}
//
int ConvertCharToFortran(char *StringInC, char* StringInF) 
{
//Change the termination of string in c ('\0') by a blank for fortran equivalence
  int mon_i=0;
  while(StringInC[mon_i]!='\0')
    { 
     StringInF[mon_i]=StringInC[mon_i];
     mon_i++;
   }
   StringInF[mon_i]=' ';
   return mon_i;
}
