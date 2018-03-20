/*
    Dual Vibration Configuration Interaction (DVCI).
    A novel factorisation of molecular Hamiltonian for
    high performance infrared spectrum computation.
    Copyright (C) 2017  Romain Garnier

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
//
//
#include "Assemblage.h"
#include "Shared.h"
#include "Graph.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Basic.h"
#include <cassert>
//
void MatrixQPower(double **APow, int Dim, int Pow)
{
//Compute the square matrix Q^(Pow) (Dimension Dim x Dim) and return it in APow
//The matrix Q beeing <Psi_n(Q),Q,Psi_m(Q)>, (n,m)in [0,Dim[
double **Q=NULL;
double **B=NULL;
//Set APow to zero
InitMat(APow, Dim, Dim)
//Allocate temporary matrix B 
allocate_MatDb(B, Dim, Dim)
//APow^(0)=Identity matrix 
for (int ii=0; ii<Dim; ii++)
APow[ii][ii]=1;
allocate_MatDb(Q, Dim, Dim)
MatrixQ(Q,Dim)
/*Loop begin at power 1 
1-Copy matrix APow into B
2-Set APow to zero
3-Multiply Q by B.
*/
for (int pp=1;pp<=Pow; pp++)
{
CopyMat(B, APow, Dim, Dim)
for (int ii=0 ; ii<Dim ; ii++)
  {
for (int jj=0 ; jj<Dim ; jj++)
   {
APow[ii][jj]=0;
for(int ll=0;ll<Dim;ll++)
     {
APow[ii][jj]=APow[ii][jj]+(B[ii][ll] * Q[ll][jj]);
     }
   }
  }
 }
FreeMat(B,Dim)
FreeMat(Q,Dim)
}
//
double CorElem(uint8_t *ModeLin, uint8_t *ModeCol, int ni, int nj, int nk, int nl, MatrixElem *QQ, int DegrePolP1,\
double *Omega, int NMode)
{
//Evaluate corriolis integral terms between ModeLin(=m) and ModeCol(=n):
//<Phi_m|(\sqrt{\nu_j/\nu_i}Qi*dQj- \sqrt{\nu_i/\nu_j} QjdQi)(\sqrt{\nu_l/\nu_k}Qk*dQl- \sqrt{\nu_k/\nu_l} Ql*dQk)|Phi_n>
//m_i=n_i +- 1 or 2
//Omega : Harmonic force constants in Hartree : used as conversion factor i.e \sqrt{\nu_j/\nu_i},\sqrt{\nu_l/\nu_k}
//in coriolis terms.
/*QQ[dd].Coeff[ii][jj] are matrix elements of operator Q^dd for dd in [0,DegrePol],
QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of operator D2Q (second order derivative)
QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order derivative)
QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q 
QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ*/

int IG=22;
int ID=33;
//
int **LocalMonm=NULL;
allocate_MatInt(LocalMonm, 4, NMode);
InitMat(LocalMonm,4,NMode)
//
//First term Qi DQj Qk DQl
// X Omega for derivative and /Omega for power
LocalMonm[0][ni]+=IG;
LocalMonm[0][nj]-=IG;
LocalMonm[0][nk]+=ID;
LocalMonm[0][nl]-=ID;
//
//Second term Qi DQj Ql DQk
LocalMonm[1][ni]+=IG;
LocalMonm[1][nj]-=IG;
LocalMonm[1][nl]+=ID;
LocalMonm[1][nk]-=ID;
//
//Third term Qj DQi Qk DQl
LocalMonm[2][nj]+=IG;
LocalMonm[2][ni]-=IG;
LocalMonm[2][nk]+=ID;
LocalMonm[2][nl]-=ID;
//
//Forth term Qj DQi Ql DQk
LocalMonm[3][nj]+=IG;
LocalMonm[3][ni]-=IG;
LocalMonm[3][nl]+=ID;
LocalMonm[3][nk]-=ID;
//
double Elem[4]={1,1,1,1};
//InitTabVal(Elem,4,1)
int PowD;
//
 for (int kk=0; kk<4; kk++)
    {
     for (int mm=0; mm<NMode; mm++)   
        {
          if(LocalMonm[kk][mm]==-55)
          {
          PowD=DegrePolP1;
          Elem[kk]*=pow(Omega[mm],2);
          Elem[kk]*=QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]];
          }
          else if(LocalMonm[kk][mm]==55)
          {
          PowD=2;
          Elem[kk] /=pow(Omega[mm],2);
          Elem[kk]*=QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]];
          }
          else if(LocalMonm[kk][mm]==22 || LocalMonm[kk][mm]==33 )
          {
          PowD=1;
          Elem[kk] /=Omega[mm];
          Elem[kk]*=QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]];
          }
          else if(LocalMonm[kk][mm]==-22 || LocalMonm[kk][mm]==-33 )
          {
          PowD=DegrePolP1+1;
          Elem[kk]*=Omega[mm];
          Elem[kk]*=QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]];
          }
          else if(LocalMonm[kk][mm]==-11)
          {
          PowD=DegrePolP1+2;
          Elem[kk]*=QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]];
          }
          else if(LocalMonm[kk][mm]==11)
          {
          PowD=DegrePolP1+3;
          Elem[kk]*=QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]];
          } 
        }
     }
FreeMat(LocalMonm, 4)
double Tot=Elem[0]-Elem[1]-Elem[2]+Elem[3];
return Tot;
//
}

double MatrixEvalHarm (int NMode,int DegrePolP1,int NCPol, int NPES, int *DegreCoupl,\
 ConfigId ModeLin, ConfigId ModeCol,  KForce KFC, MatrixElem *QQ, double **ZetaXYZ, double *Omega,\
 uint32_t NXDualHPlus,ConfigId *DualHPos, LocalFF *LFF, uint32_t *Permuter, int LMC)
{
/*Compute matrix element between ModeLin(=multi-index Lin) and ModeCol(multi-index Col).
Here matrix element <Phi_Lin | H |Phi_Col> involves the force constants included in
{LFF(Lin-Col)}={K_[c1, ... , cn], |Lin_i-Col_i|=ci-2ti}, 
First is computed the differences Lin-Col, check if the set LFF non void and add the harmonic energy E_Col^0 when Col=Lin. 
Lin-Col corresponds to an integer xx assigned to an excitation DualHPos[xx].
*/ 
//Force constants:
//KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
//Local force field:
//LFF[xx] : local force field associated with positive excitation DualHPos[xx].
//LFF[xx].Idx[ii] < 0 are indexes of non coupled force constants:
//KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
//0 <= LFF[xx].Idx[ii] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF[xx].Idx[ii]].
//LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
//with a unique corresponding (ni,nj,nk,nl)  
//LFF[xx].Idx[ii] are defined for ii in [0,LFF[xx].Num[ 
//Permuter : give a sorted order (according to memcmp) of elements of DualHPos.
//ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}
//Omega : Harmonic force constants in Hartree, used as conversion factor to compute rotational matrix elements
//NCPol : Maximal number of couplings
//DegreCoupl[xx] : degree of PES for coupling xx-1 (start at one)
/*QQ[dd].Coeff[ii][jj] are matrix elements of operator Q^ii for ii in [0,DegrePol],
QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of operator D2Q (second order derivative)
QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order derivative)
QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q 
QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ*/
//LMC=Lin - Col : if zero then the surface IndexEx corresponds to the zero excitation
//located in position 0. 
//           
   double ValMat=0; 
   uint8_t *Tester; //Variable multi index during the loop
   allocate_TabUint8(Tester,NMode)
   int DistanceModal;
   int IndexEx;
//Sum over all the monomial coefficients  prod_mm H_{mm,kk} : outside the loop mm.
ValMat=0;
//
if(!LMC)
{
  IndexEx=0;
  ValMat=MatrixElement (NMode, LFF[IndexEx],DegrePolP1,\
       NPES, ModeLin.Degrees, ModeCol.Degrees, QQ,KFC, ZetaXYZ, Omega);
  ValMat+=GetEHarmonic0(ModeLin.Degrees, NMode, KFC);
}
else
{
//
for(int ii=0;ii<NMode;ii++)
 {
if(ModeLin.Degrees[ii]>ModeCol.Degrees[ii])
  {Tester[ii]=ModeLin.Degrees[ii]-ModeCol.Degrees[ii];}
else
  {Tester[ii]=ModeCol.Degrees[ii]-ModeLin.Degrees[ii];}
 }
//
int Cpld=CmptNNULL(Tester, NMode);
//
if ( Cpld<=NCPol ) 
 {
//
DistanceModal=sum_array<uint8_t>(Tester,NMode);
//         
if(DistanceModal<=DegreCoupl[Cpld-1])//Only for Harmonic basis set : Start with 0=1 modal!
  {
 IndexEx=QsearchMode(Tester, DualHPos,Permuter,NXDualHPlus, NMode);
 if(IndexEx>0)//Zero treated separately
   {
  ValMat=MatrixElement (NMode, LFF[IndexEx],DegrePolP1,\
       NPES, ModeLin.Degrees, ModeCol.Degrees, QQ,KFC, ZetaXYZ, Omega);
   }
  }
 } 
}
//
delete [] Tester; 
return ValMat;
}
//
uint32_t AssembleHarmCSC (int NMode,int DegrePol, int NCPol, int Iteration, int NPES,\
 ConfigId *ModeAct,SizeArray *Size,KForce KFC,double *ValAct, CSC IJAct, MatrixElem *QQ, int *DegreCoupl,\
 double **ZetaXYZ, double *Omega, double ThrMat, uint64_t NNZActMax,ConfigId *DualHPos, LocalFF *LFF, \
 uint32_t *Permuter, uint32_t NXDualHPlus)
{
//Assemble the graph of sparse matrix Hb in CSC format in IJAct, and compute NNZ values in ValAct.
//Return the new number of non zero elements
//Force constants:
//KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
//Local force field:
//LFF[xx] : local force field associated with positive excitation DualHPos[xx].
//LFF[xx].Idx[ii] < 0 are indexes of non coupled force constants:
//KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
//0 <= LFF[xx].Idx[ii] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF[xx].Idx[ii]].
//LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
//with a unique corresponding (ni,nj,nk,nl)  
//LFF[xx].Idx[ii] are defined for ii in [0,LFF[xx].Num[  
//Permuter : give a sorted order (according to memcmp) of elements of DualHPos.
//ModeAct[Lin].Degrees multi-index array of length NMode corresponding to one element Lin in active space B.
//ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}
//Omega : Harmonic force constants in Hartree : used as conversion factor to compute \sqrt{\frac{\nu_j}{\nu_i}},\sqrt{\frac{\nu_p}{\nu_m}}
//in coriolis terms
//NCPol : Maximal number of couplings
//DegreCoupl[xx] : degree of PES for coupling xx-1 (start at one)
/*QQ[dd].Coeff[ii][jj] are matrix elements of operator Q^ii for ii in [0,DegrePol],
QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of operator D2Q (second order derivative)
QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order derivative)
QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q 
QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ
Size[Iteration].NNZAct: NNZ of active matrix Hb already incremented in a previous iteration;
Size[Iteration+1].NNZAct: Size of active space B beeing modified by the new NNZ values;
//LMC=Lin - Col : if zero then the surface IndexEx corresponds to the zero excitation. 
//located in position 0. 
*/
Size[Iteration+1].NNZAct=Size[Iteration].NNZAct;
//
   int DegrePolP1=DegrePol+1;        
   double ValMat;
   uint32_t Lin,Col;   
   int LMC;  
   int AssignFirstCol;
//
for (Col=Size[Iteration].DimAct ; Col < Size[Iteration+1].DimAct ; Col++)
{
AssignFirstCol=1;
for (Lin=0 ; Lin <= Col ; Lin++)
 {
 LMC=Col-Lin;
 ValMat=MatrixEvalHarm(NMode, DegrePolP1, NCPol, NPES, DegreCoupl,\
  ModeAct[Col],ModeAct[Lin], KFC, QQ, ZetaXYZ, Omega,NXDualHPlus, DualHPos, LFF, Permuter, LMC);
//
if( (SIGN<double>(ValMat)*ValMat>ThrMat)  )
    {
    if(Size[Iteration+1].NNZAct < NNZActMax)     
       {
         ValAct[Size[Iteration+1].NNZAct]=ValMat;
//
         IJAct.I[Size[Iteration+1].NNZAct]=Lin;
          if(AssignFirstCol)
           {
          IJAct.NJ[Col]=Size[Iteration+1].NNZAct;
          AssignFirstCol=0;
           }
//Increase the dimension of the NNZ elements 
         Size[Iteration+1].NNZAct++;
       }
    else
       {
       printf("***NNZ(Hb) overflow: To poursue, increase Memory or KNNZ***\n\n");
       return 0;
       }
     }
   } 
  }
//
IJAct.NJ[Size[Iteration+1].DimAct]=Size[Iteration+1].NNZAct;
//
return Size[Iteration+1].DimAct;
}
//
double MatrixElement (int NMode, LocalFF LFF, int DegrePolP1, int NPES,\
 uint8_t *ModeLin, uint8_t *ModeCol, MatrixElem *QQ, KForce KFC, double **ZetaXYZ, double *Omega)
{
/*
Compute <\Phi_Lin | H |\Phi_Col> for Lin = ModeLin and Col = ModeCol 
(Lin and Col are multi-indexes of length NMode).
Corresponding force constants are the one associated with index of excitation |Lin-Col| and are located in LFF.
Force constants:
KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
Local force field:
LFF[xx] : local force field associated with positive excitation DualHPos[xx].
LFF[xx].Idx[ii] < 0 are indexes of non coupled force constants:
KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
0 <= LFF[xx].Idx[ii] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF[xx].Idx[ii]].
LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
with a unique corresponding (ni,nj,nk,nl)  
LFF[xx].Idx[ii] are defined for ii in [0,LFF[xx].Num[ 
ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}
Omega : Harmonic force constants in Hartree : used as conversion factor to compute \sqrt{\frac{\nu_j}{\nu_i}},\sqrt{\frac{\nu_p}{\nu_m}}
in coriolis terms.
QQ[dd].Coeff[ii][jj] are matrix elements of operator Q^ii for ii in [0,DegrePol],
QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of operator D2Q (second order derivative)
QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order derivative)
QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q 
QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ
*/
  double Valk;
  double ValMat=0; 
  int PowD;      
  int dd,mm,kk;
  int ni,nj,nk,nl;
//  
for (int ii=0; ii<LFF.Num; ii++)
   {     
      if(LFF.Idx[ii]<0)
      {
       dd=-LFF.Idx[ii]/NMode;
       mm=-LFF.Idx[ii]-dd*NMode;
       PowD=dd+1;
       ValMat+=(QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]])*KFC.KijNCpld[dd][mm];//Should be a common factor
      }          
      else if(LFF.Idx[ii]<NPES && LFF.Idx[ii]>=0)
      {
       kk=LFF.Idx[ii];
       Valk=1;
       for (mm=0; mm<NMode; mm++)
        {   
      PowD=KFC.Monm[kk][mm];                   
      Valk*=(QQ[PowD].Coeff[ModeLin[mm]][ModeCol[mm]]);//Should be a common factor
        }
        ValMat+=KFC.KijCpld[kk]*Valk;
      }
       else if (LFF.Idx[ii]>=NPES)
       {
       GetIJKL(LFF.Idx[ii]-NPES,NMode, &ni, &nj, &nk, &nl);
       if(ni>nj || nk>nl)
        {
       printf( "ni %i, nj %i, nk %i, nl %i\n",ni,nj,nk,nl);
       printf( "Problem with Coriolis indexes in MatrixElement\n");
       return 0; 
        }
       ValMat-=CorElem(ModeLin, ModeCol, ni, nj, nk, nl, QQ, DegrePolP1,Omega,  NMode)*\
       ZetaXYZ[ni + ((nj-1)*nj)/2][nk + ((nl-1)*nl)/2];
       }
   }
     return ValMat;
}
//
void FlyRezMVP2(ConfigId *ModeRez, ConfigId *ModeAct,ConfigId *DualHPos, MatrixElem *QQ,\
     uint32_t *Permuter, uint32_t *TabNull, KForce KFC,SizeArray *Size, int NMode,\
     CSC IJRez,uint32_t NXDualHPlus, int DegrePol,\
     int NPES, int Iteration, int NScreen, int *TabScreen, SizeArray SizeMax,\
     float *RezVect,LocalFF *LFF,double *EigVec, double **ZetaXYZ, double *Omega)
{
/*Complete the matrix vector product Hsb*X for elements (Col,Lin) in (B-A)x(H*(B-A)-B)
 thanks to the graph of residual matrix stored in IJRez. 
The compressed sparse column format is used:
IJRez.NJ[jj] = number of the first nnz elmement of column jj
IJRez.I[nn] = line number of nnz element number nn 
Permuter : permutation index array of elements of DualHPos (memcmp order)
Force constants:
KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
Local force field:
LFF[xx] : local force field associated with positive excitation DualHPos[xx].
LFF[xx].Idx[ii] < 0 are indexes of non coupled force constants:
KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
0 <= LFF[xx].Idx[ii] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF[xx].Idx[ii]].
LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
with a unique corresponding (ni,nj,nk,nl)  
LFF[xx].Idx[ii] are defined for ii in [0,LFF[xx].Num[ 
ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}.
Omega : Harmonic force constants in Hartree : used as conversion factor i.e \sqrt{\frac{\nu_j}{\nu_i}},\sqrt{\frac{\nu_p}{\nu_m}}
in coriolis terms.
If TabNull[LinRez]=0 then element LinRez is no longer part of residual space.
Matrix elements:
QQ[dd].Coeff[ii][jj] are matrix elements of operator Q^ii for ii in [0,DegrePol],
QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of operator D2Q (second order derivative)
QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order derivative)
QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q 
QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ
Size[Iteration].DimAct: Size of active space B at previous iteration;
Size[Iteration+1].DimAct: Size of active space B at current iteration;
Size[Iteration+1].NNZRez: NNZ of residual matrix beeing updated for TabNull[LinRez]=0,
indicating that LinRez is not a member of the residual space anymore.
*/
int DegrePolP1=DegrePol+1;
uint32_t LinRez=0;
uint8_t *Tester=NULL; //Variable multi index during the loop     
double MonomProd;
long IndexEx;
//
uint32_t vv, Col;
//
allocate_TabUint8(Tester,NMode)
//
for (Col=0 ; Col < Size[Iteration].DimAct ; Col++)
{ 
 for(vv=IJRez.NJ[Col];vv<IJRez.NJ[Col+1];vv++)
 {
  LinRez=IJRez.I[vv];
  if(TabNull[LinRez])
  {
 for(int ii=0;ii<NMode;ii++)
   {
 if(ModeAct[Col].Degrees[ii]>=ModeRez[LinRez].Degrees[ii])
    {Tester[ii]=ModeAct[Col].Degrees[ii]-ModeRez[LinRez].Degrees[ii];}
 else
    {Tester[ii]=ModeRez[LinRez].Degrees[ii]-ModeAct[Col].Degrees[ii];}
   }
  IndexEx=QsearchMode(Tester, DualHPos,Permuter,NXDualHPlus, NMode);
//
  MonomProd=MatrixElement (NMode, LFF[IndexEx] ,DegrePolP1,\
       NPES, ModeAct[Col].Degrees, ModeRez[LinRez].Degrees, QQ,KFC, ZetaXYZ, Omega);
//
        for (int ll=0; ll < NScreen; ll++)
        {  
        RezVect[LinRez+SizeMax.DimRez*ll]+=MonomProd*\
        EigVec[TabScreen[ll]*Size[Iteration+1].DimAct+Col];
        }
   } 
  }
 }
   delete [] Tester;
//
   uint32_t *NJTmp=new uint32_t[Size[Iteration+1].DimAct];
   uint32_t NNNZ=0;
   int AssignFirstCol;
// 
for (Col=0; Col < Size[Iteration+1].DimAct; Col++)
{
AssignFirstCol=1;
 for(vv=IJRez.NJ[Col];vv<IJRez.NJ[Col+1];vv++)
 {
  LinRez=IJRez.I[vv];
  if(TabNull[LinRez])
  {
  IJRez.I[NNNZ]=LinRez;
  if(AssignFirstCol)//First nnz of column number Col has to be assigned
    {
    NJTmp[Col]=NNNZ;
    AssignFirstCol=0;
    }
  NNNZ++;//First should be zero because NNNZ has not been incremented
  }//But not necessary Col not necessary
 }
// 
  if(AssignFirstCol)//Fillin of holes
  {
   NJTmp[Col]=NNNZ;
  }  
}   
//
for (Col=0 ; Col < Size[Iteration+1].DimAct ; Col++)
 {
 IJRez.NJ[Col]=NJTmp[Col];
 }
//
IJRez.NJ[Size[Iteration+1].DimAct]=NNNZ;
Size[Iteration+1].NNZRez=NNNZ;
//
delete [] NJTmp;//Works if whole between columns
}
