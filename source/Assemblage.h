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

#ifndef Assemblage_H
#define Assemblage_H
#include "Shared.h"
#include <pthread.h>

#define NBuffer 1000
#define MaxThread 100

//Q matrix
#define MatrixQ(Q,Dim)\
for (int i = 0 ; i < Dim ; i++){\
for (int j = 0 ; j < Dim ; j++){Q[i][j]=0;}}\
for (int i=0;i<Dim-1;i++){\
    Q[i][i+1]=sqrt((double)(i+1)/(double)(2));}\
for (int i=1;i<Dim;i++){\
    Q[i][i-1]=sqrt((double)(i)/(double)(2));}

//D2QMatrix
#define MatrixD2Q(Q,Dim)\
for (int i = 0 ; i < Dim ; i++){\
for (int j = i+1 ; j < Dim ; j++){Q[i][j]=0; Q[j][i]=0; }}\
for (int l=0; l<Dim; l++){ Q[l][l]=-((double)(l)+(double)(0.5)); }\
for (int i=2;i<Dim;i++){\
    Q[i-2][i]=sqrt((double)((i)*(i-1)))/(double)(2);}\
for (int i=0;i<Dim-2;i++){\
    Q[i+2][i]=sqrt((double)((i+2)*(i+1)))/(double)(2);}\

//D1QMatrix
#define MatrixD1Q(Q,Dim)\
for (int i = 0 ; i < Dim ; i++){\
for (int j = 0 ; j < Dim ; j++){Q[i][j]=0;}}\
for (int i=0;i<Dim-1;i++){\
    Q[i+1][i]=-sqrt((double)(i+1)/(double)(2));}\
for (int i=1;i<Dim;i++){\
    Q[i-1][i]=sqrt((double)(i)/(double)(2));}

//QD1QMatrix
#define MatrixQD1Q(Q,Dim)\
for (int i = 0 ; i < Dim ; i++){\
for (int j = 0 ; j < Dim ; j++){Q[i][j]=0;}}\
for (int i=2;i<Dim;i++){\
    Q[i-2][i]=sqrt((double)((i)*(i-1)))/(double)(2);}\
for (int i=0;i<Dim-2;i++){\
    Q[i+2][i]=-sqrt((double)((i+2)*(i+1)))/(double)(2);}\
for (int i=0;i<Dim;i++){\
Q[i][i]=-0.5;}

//D1QQMatrix
#define MatrixD1QQ(Q,Dim)\
for (int i = 0 ; i < Dim ; i++){\
for (int j = 0 ; j < Dim ; j++){Q[i][j]=0;}}\
for (int i=2;i<Dim;i++){\
    Q[i-2][i]=sqrt((double)((i)*(i-1)))/(double)(2);}\
for (int i=0;i<Dim-2;i++){\
    Q[i+2][i]=-sqrt((double)((i+2)*(i+1)))/(double)(2);}\
for (int i=0;i<Dim;i++){\
Q[i][i]=0.5;}

//Allocation +Fillin Matrix elements
#define MatricesElem(QQ,DegPol,dimax,NMatrices)\
MatrixElem *QQ=NULL;\
QQ=new MatrixElem[NMatrices];\
for (int i=0;i<=DegPol;i++){\
QQ[i].Coeff = NULL;\
QQ[i].Coeff = new double* [dimax+DegPol];\
for (int kk = 0; kk < dimax+DegPol; kk++){\
QQ[i].Coeff[kk] = NULL;\
QQ[i].Coeff[kk]=new double[dimax+DegPol];}\
MatrixQPower(QQ[i].Coeff, dimax+DegPol, i);}
//Allocation only
#define MatricesElemAllocate(QQ,dimax,NMatrices)\
MatrixElem *QQ=NULL;\
QQ=new MatrixElem[NMatrices];\
for (int i=0;i<NMatrices;i++){\
QQ[i].Coeff = NULL;\
QQ[i].Coeff = new double* [dimax];\
for (int kk = 0; kk < dimax; kk++){\
QQ[i].Coeff[kk] =  NULL;\
QQ[i].Coeff[kk] = new double[dimax];}}
//Free MatricesElements
#define FreeMatricesElem(QQ,NMatrices,DimMat)\
for (int kk = 0; kk < NMatrices; kk++){\
for (int i = 0; i < DimMat; i++)\
{delete [] QQ[kk].Coeff[i];}\
delete [] QQ[kk].Coeff;}\
delete [] QQ;
//
 struct MatrixElem
 {
double **Coeff;
 };
//
void MatrixQPower(double **APow, int Dim, int Pow);
//Compute the square matrix Q^(Pow) (Dimension Dim x Dim) and return it in APow
//The matrix Q beeing <Psi_n(Q),Q,Psi_m(Q)>, (n,m)in [0,Dim[
//
double MatrixElement (int NMode, LocalFF LFF, int DegrePolP1,int NPES,\
 uint8_t *ModeLin, uint8_t *ModeCol, MatrixElem *QQ, KForce KFC, double **ZetaXYZ, double *Omega);
/*Compute <\Phi_Lin | H |\Phi_Col> for Lin = ModeLin and Col = ModeCol 
(Lin and Col are multi-indexes of length NMode).
Corresponding force constants are the one associated with index of excitation |Lin-Col| and are located in LFF.*/
//
double MatrixEvalHarm (int NMode,int DegrePolP1,int NCPol, int NPES, int *DegreCoupl,\
 ConfigId ModeLin, ConfigId ModeCol,  KForce KFC, MatrixElem *QQ, double **ZetaXYZ, double *Omega,\
 unsigned int NXDualHPlus,ConfigId *DualHPos, LocalFF *LFF, unsigned int *Permuter, int LMC, int IncK2);
/*Compute matrix element between ModeLin(=multi-index Lin) and ModeCol(multi-index Col).
Here matrix element <Phi_Lin | H |Phi_Col> involves the force constants included in
{LFF(Lin-Col)}={K_[c1, ... , cn], |Lin_i-Col_i|=ci-2ti}, 
First is computed the differences Lin-Col, check if the set LFF non void and add the harmonic energy E_Col^0 when Col=Lin. 
Lin-Col corresponds to an integer xx assigned to an excitation DualHPos[xx].
IncK2:Say if the the harmonic energy should be added (IncK2=0) or not (IncK2=1)
*/ 
//
void FlyRezMVP2(ConfigId *ModeRez, ConfigId *ModeAct,ConfigId *DualHPos, MatrixElem *QQ,\
     unsigned int *Permuter, unsigned int *TabNull, KForce KFC,SizeArray *Size, int NMode,\
     CSC IJRez,unsigned int NXDualHPlus,  int DegrePol,\
     int NPES, int Iteration, int NScreen, int *TabScreen, SizeArray SizeMax,\
     float *RezVect,LocalFF *LFF,double *EigVec, double **ZetaXYZ, double *Omega);
//
/*Complete the matrix vector product Hsb*X for elements (Col,Lin) in (B-A)x(H*(B-A)-B)
 thanks to the graph of residual matrix stored in IJRez. */
unsigned int AssembleHarmCSC (int NMode,int DegrePol, int NCPol, int Iteration, int NPES,\
 ConfigId *ModeAct,SizeArray *Size,KForce KFC,double *ValAct, CSC IJAct, MatrixElem *QQ, int *DegreCoupl,\
 double **ZetaXYZ, double *Omega, double ThrMat, uint64_t NNZActMax,ConfigId *DualHPos, LocalFF *LFF, \
 unsigned int *Permuter, unsigned int NXDualHPlus, int IncK2);
//Assemble the graph of active sparse matrix in CSC format in IJAct, and compute NNZ values in ValAct.
//Return the new number of non zero elements
//IncK2:Say if the the harmonic energy should be added (IncK2=0) or not (IncK2=1)
//
double CorElem(uint8_t *ModeLin, uint8_t *ModeCol, int ni, int nj, int nk, int nl, MatrixElem *QQ, int DegrePolP1,\
double *Omega, int NMode);
//Evaluate corriolis integral terms between ModeLin(=m) and ModeCol(=n):
//<Phi_m|(\sqrt{\nu_j/\nu_i}Qi*dQj- \sqrt{\nu_i/\nu_j} QjdQi)(\sqrt{\nu_l/\nu_k}Qk*dQl- \sqrt{\nu_k/\nu_l} Ql*dQk)|Phi_n>
//m_i=n_i +- 1 or 2
//
#endif
