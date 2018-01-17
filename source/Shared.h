/*
    Dual Vibration Configuration Interaction (DVCI).
    A novel factorisation of molecular Hamiltonian for
    high performance infra-red spectrum computation.
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

In routine MomentInertie, the computation of principal axes of inertia momentum requires to diagonalize a matrix
with DSYEV subroutine of LAPACK librairy[1].

[1] E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra, J. Du Croz,
A. Greenbaum, S. Hammarling, A. McKenney, D. Sorensen, LAPACK Users’ Guide,
3rd Edition, Society for Industrial and Applied Mathematics, Philadelphia, PA, 1999.

Additionally to the source provided, a description of DSYEV can be found here 
http://www.netlib.org/lapack/explore-3.1.1-html/dsyev.f.html 

-------------------------------------------------

In routine Zeta

The rotational coefficients are computed according to the method of Meal and Polo[2]

[2] Vibration—Rotation Interaction in Polyatomic Molecules. II. The Determination of Coriolis Coupling Coefficients
Janet Hawkins Meal and S. R. Polo, The Journal of Chemical Physics 1956 24:6, 1126-1133 */
//
#ifndef Shared_H
#define Shared_H
#include "stdint.h"
#include <cstdio>
#include <cmath>

#define HA_TO_CM 219474.6435
#define AMU_TO_ME 1.660538921E-27/9.10938291E-31
#define ME_TO_AMU 9.10938291E-31/1.660538921E-27
#define BOHR_TO_ANGST 0.52917721092
#define PI 3.141592653589793
#define Plank_JS 6.62607004081E-34
#define BOHR_TO_NM 0.05291772109217
#define BOHR_TO_M 0.05291772109217E-9
#define ME_TO_KG 9.1093829140E-31
#define CLight_MS 299792458
#define MEBOHR2_TO_CM 6.62607004081E9/(8*PI*PI*5.291772109217*5.291772109217*2.99792458*9.1093829140)
//#define MEBOHR2_TO_CM 5807.049070575313504001047865451953183834319123869964251972560085559754
#define AMU_TO_KG 1.66053904020E-27
#define AMUANG2_TO_CM 6626.07004/(8*PI*PI*1.66053904020*2.99792458)
#define IterMax 200
#define MaxDegrePI 40
#define MaxCpld 20
#define KNNZRez 0.9
#define MaxTarget 200
#define MaxNormal 100
#define MaxRef 2000

#define MaxChar 10000
#define MinChar 10

//Allocation routines for array and structures
#define allocate_MatDb(arr, n, m)\
arr=NULL;\
arr = new double*[n];\
for (int j = 0; j < n; j++){\
arr[j]=NULL;\
arr[j] = new double[m];}
#define allocate_MatInt(arr, n, m)\
arr=NULL;\
arr = new int32_t* [n];\
for (int i = 0; i < n; i++){\
arr[i]=NULL;\
arr[i] = new int32_t[m];}
#define allocate_TabDb(arr, n)\
arr=NULL;\
arr = new double[n];
#define allocate_TabInt(arr, n)\
arr = new int32_t[n];
#define allocate_TabSize(Size1,n)\
SizeArray * Size1;\
Size1 = new SizeArray [n];\
for (int k = 0; k < n; k++)\
{Size1[k].DimAct=0;Size1[k].DimRez=0;\
Size1[k].NNZAct=0;Size1[k].NNZRez=0;}
#define AllocateCSC(IJ,n,m)\
CSC IJ;\
InitCSC(IJ,n,m);
#define allocate_MatUInt(arr, n, m)\
uint32_t **arr = NULL;\
arr= new uint32_t* [n];\
for (int i = 0; i < n; i++){\
arr[i]=NULL;\
arr[i] = new uint32_t[m];}
#define allocate_TabUint8(arr, n)\
arr=NULL;\
arr = new uint8_t[n];
#define allocate_TabUint16(arr, n)\
arr=NULL;\
arr = new uint16_t[n];
#define allocate_MatUint8(arr, n, m)\
arr=NULL;\
arr = new uint8_t* [n];\
for (int i = 0; i < n; ++i){\
arr[i]=NULL;\
arr[i] = new uint8_t[m];}
#define allocU32MatU8(arr, n, m)\
arr=NULL;\
arr = new uint8_t* [n];\
for (uint32_t i = 0; i < n; ++i){\
arr[i]=NULL;\
arr[i] = new uint8_t[m];}
#define allocate_TabUInt(arr, n)\
arr=NULL;\
arr = new uint32_t[n];
//Delete routines
#define FreeMat(arr,n)\
for (int i = 0; i < n; i++){\
delete [] arr[i];}\
delete [] arr;
#define FreeUMat(arr,n)\
for (uint32_t i = 0; i < n; i++){\
delete [] arr[i];}\
delete [] arr;
#define FreeTab(arr) delete [] arr;
#define FreeCSC(IJ)\
delete [] IJ.I; delete [] IJ.NJ;
#define FreeLFF(LFF,NLFF)\
for (uint32_t i = 0 ; i < NLFF ; i++)\
{DeleteLFF(LFF[i]);}\
delete [] LFF;
//Display routines
#define afficheMatInt(Mat, lin, col)\
for (int i = 0 ; i < lin ; i++){\
for (int j = 0 ; j < col ; j++){\
printf("%d  ", Mat[i][j]);} printf(" \n");}
#define afficheMatDbl(Mat, lin, col)\
for (int i = 0 ; i < lin ; i++){\
for (int j = 0 ; j < col ; j++){\
printf("%f  ", Mat[i][j]);}printf(" \n");}
#define afficheTabDb(Tab, lin)\
for (int i = 0 ; i < lin ; i++)\
{printf(" %f \n", Tab[i]);}        
#define afficheTabInt(Tab, lin)\
for (int i = 0 ; i < lin ; i++)\
{printf(" %d , ", Tab[i]);}\
printf("\n"); 
#define afficheTabUI(Tab, lin)\
for (uint32_t i = 0 ; i < lin ; i++)\
{printf(" %u , ", Tab[i]);}\
printf("\n");
//Display Size
#define DisplaySize(Size1)\
printf("DimAct %u, DimRez %u, NNZAct %u, NNZRez %u \n",\
Size1.DimAct,Size1.DimRez-Size1.DimAct,Size1.NNZAct,Size1.NNZRez);
//Check size
#define TestSize(Size1,SizeActMax,SizeRezMax,NNZActMax,NNZRezMax)\
if (Size1.DimAct >= SizeActMax)\
{printf("DimAct trop grande %u", Size1.DimAct);\
return 0;}\
else if (Size1.DimRez >= SizeRezMax)\
{printf("DimRez trop grande %u", Size1.DimRez);\
return 0;}\
else if (Size1.NNZAct >= NNZActMax)\
{printf("NNZAct trop grande %u", Size1.NNZAct);\
return 0;} \
else if (Size1.NNZRez >= NNZRezMax)\
{printf("NNZRez trop grande %u", Size1.NNZRez);\
return 0;}    
//Initialisation array 
#define InitTab(Tab, lin)\
for (uint32_t i = 0 ; i < lin ; i++){Tab[i]=0;}
#define InitTabInt(Tab, lin)\
for (int i = 0 ; i < lin ; i++){Tab[i]=0;}
#define InitTabOne(Tab, lin)\
for (uint32_t i = 0 ; i < lin ; i++){Tab[i]=1;}
#define InitTab1(Tab, lin)\
for (int i = 0 ; i < lin ; i++){Tab[i]=1;}
#define InitTabVal(Tab, lin, Val)\
for (int i = 0 ; i < lin ; i++){Tab[i]=Val;}
//Copy
#define CopyTab(Tab1,Tab2, lin)\
for (int i = 0 ; i < lin ; i++){Tab1[i]=Tab2[i];}
#define CopyTabUint(Tab1,Tab2, lin)\
for (uint32_t i = 0 ; i < lin ; i++){Tab1[i]=Tab2[i];}
#define InitMat(Mat, lin, col)\
for (int i = 0 ; i < lin ; i++){\
for (int j = 0 ; j < col ; j++){Mat[i][j]=0;}}
#define InitMatUI(Mat, lin, col)\
for (uint32_t i = 0 ; i < lin ; i++){\
for (uint32_t j = 0 ; j < col ; j++){Mat[i][j]=0;}}
#define InitMatId(Mat, lin, col)\
for (int i = 0 ; i < lin ; i++){\
for (int j = 0 ; j < col ; j++){\
if(i!=j){Mat[i][j]=0;}\
else{Mat[i][j]=1;}}}
#define CopyMat(Mat1, Mat2, lin, col)\
for (int i = 0 ; i < lin ; i++){\
for (int j = 0 ; j < col ; j++){Mat1[i][j]=Mat2[i][j];}}
//Sum of arrays
#define SumTab(Tab1,Tab2,lin)\
for (int i = 0 ; i < lin ; i++)\
{Tab1[i]=Tab1[i]+Tab2[i];}
#define SumTabSquare(Tab1,Tab2,lin)\
for (int i = 0 ; i < lin ; i++)\
{Tab1[i]=Tab1[i]+pow(Tab2[i],2);}
#define SumTabAbs(Tab1,Tab2,lin)\
for (uint32_t i = 0 ; i < lin ; i++)\
{Tab1[i]=Tab1[i]+SIGN<double>(Tab2[i])*Tab2[i];}
#define SumKTabAbs(Tab1,Tab2,lin,Kte)\
for (uint32_t i = 0 ; i < lin ; i++)\
{Tab1[i]=Tab1[i]+SIGN<float>(Tab2[i])*Tab2[i]*Kte;}
//
//Maximal number of assignments to show
#define MaxAssign 100
//
/*Measure time in seconds */
#define TIME_DIFFS(t1, t2) \
        t2.tv_sec - t1.tv_sec 
//
//Measure CPU time in seconds
#define TIME_CPUS(t1, t2) \
        ((double) (t2 - t1)) / CLOCKS_PER_SEC; 
//
#define FinalMessaj()\
printf("\nDVCI Copyright (C) 2018 Romain Garnier.\n");\
printf("This program comes with ABSOLUTELY NO WARRANTY;\n\
This is free software, and you are welcome\n\
to redistribute it under certain conditions.\n\
See the GNU General Public License for more details.\n\n");
struct KForce
{
//Force constants:
//KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
double **KijNCpld;
double *KijCpld;
int32_t **Monm;
};
//
struct SizeArray
{
uint32_t DimAct;
uint64_t DimRez;
uint64_t NNZAct;
uint64_t NNZRez;
};
//
//typedef struct ConfigId ConfigId ;
 struct ConfigId
 {
uint8_t *Degrees;
 };
//
 struct LocalFF
 {
//LFF[xx] : local force field associated with positive excitation DualHPos[xx].
//LFF[xx].Idx[ii] < 0 are indexes of non coupled force constants:
//KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
//0 <= LFF[xx].Idx[ii] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF[xx].Idx[ii]].
//LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
//with a unique corresponding (ni,nj,nk,nl)  
//LFF[xx].Idx[ii] are defined for ii in [0,LFF[xx].Num[ 
int32_t Num;
int32_t *Idx;
 };
//
//Compressed sparse column format
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
struct CSC
{
uint32_t *I;
uint64_t *NJ;
};
//
extern "C"
{
void dsyev_(char* JOBZ, char* UPLO, int *N,double *A, int *LDA, double *W, double *WORK, int *LWORK,int *INFO );
}
//
void InitMode(ConfigId & M1,int Nmode);
void InitLFF(LocalFF& LFF,int NFF);
double factorial(uint32_t n);
//
double GetEHarmonic0(uint8_t* Harm,int NMode, KForce KFC);
////Return Harmonic energy sum_NM (d_i+0.5)*\nu_i
//
void InitKforce(KForce KFC, double **KijNCpld, double *KijCpld,\
 int **Monm, int DegrePol, int NMode, int NPES);
void DeleteKForce(KForce& KFC, int DegrePol, int NPES);
//
void CrossProd(double *Vec1, double *Vec2, double *Vec3);
//
void MassCenter(double *Mas, double *CoorCart, int NA, double *MassC); 
//Compute the center of Mass. 
//X=(\sum_i m_ix_i)/(\sum_i m_i)
//Y=(\sum_i m_iy_i)/(\sum_i m_i)
//Z=(\sum_i m_iz_i)/(\sum_i m_i)
//
//x_1 =CoorCart[0], y_1 =CoorCart[1], z_1 =CoorCart[2], 
//x_2 =CoorCart[3], y_2 =CoorCart[4], z_2 =CoorCart[5],
//....
//x_n =CoorCart[3*(n-1)], y_n =CoorCart[3*(n-1)+1], z_n =CoorCart[3*(n-1)+2],
//
double MomentInertie(double *Mas, double *CoorCart, int NA, double *MassC, double *PrincAx, double *PrincEig);
//Compute the principal moments of inertia and corresponding axes.
//i.e eigenvalue and eigenvectors of inertia momentum matrix. 
//Shift coordinates with center of masse MassC.
//x_i=X_i-X_c
//x_i=Y_i-Y_c
//x_i=Z_i-Z_c
//Diagonal entries of inertia matrix
//Ixx=\sum_i m_i(x_i^2+z_i^2)
//Iyy=\sum_i m_i(x_i^2+z_i^2)
//Izz=\sum_i m_i(x_i^2+y_i^2)
//Non diagonal entries
//Ixy=Iyx=-\sum_i m_i x_iy_i
//Ixz=Izx=-\sum_i m_i x_iz_i
//Izy=Izx=-\sum_i m_i z_iy_i
//
int Zeta(double **CoorMode, int NMode, int NA, double *Mas, double *PrincAx, double *PrincEig, double **ZetaXYZ, int DoRot);
//Routines that compute the Zeta coefficients according to the method of Meal and Polo
//PrincAxes are eigenvectors of inertia momentum matrix
//PrincEig are eigenvalues of inertia momentum matrix
//Normal coordinate vectors are 
//CoorMode[ii][3*(jj)], CoorMode[ii][3*(jj)+1], CoorMode[ii][3*(jj)+2], 0 <= jj< NA, 0<= ii <NMode
//
//DoRot specifies non mass weighted normal coordinates for odd value and mass weighted for even value 
//
int DetectMaxCoordinate(double *Vec, int* PositionTarget, int NTargetStates, double ThrCoor);
//Return 1 if a coordinate listed in PositionTarget is bigger than ThrCoor in absolute value.
//Return 0 if no such coordinate have been detected when NTargetStates>0
//Return 1 when NTargetStates=0
//
float EffectiveMemory2(char PidName[MaxChar]);
//Execute a terminal command ps aux | grep PidName
//and convert the output in megabytes
//
uint32_t GetSizeActMax2(int MaxEv,int MAXNCV, int NMode, int MaxQLevel, int DegrePol, int NPES, int Verbose,\
 int NXDualHPlus, double Mem, double KNNZ, double KNREZ, double KNZREZ, int SumTarget, uint32_t NIdx32, uint32_t NXDualHTrunc);
//Return the maximal size of the active space SizeActMax for a given memory Mem(in megabytes),
/*
NXDualHTruncPos-1: Total number of raising excitations in H* (after truncation through ThrKX)
NXDualHTrunc : Total number of excitations in H* (after truncation through ThrKX)
SizeRezMax=SizeActMax*(NXDualHTruncPos-1)*KNREZ
NNZActMax=SizeActMax*(KNNZ)*(NXDualHTrunc)
NNZRezMax=SizeActMax*(KNZREZ)*(NXDualHTrunc) 
*/
//
void swap (uint32_t* a,uint32_t* b );
//
float Norm2F(float *X, uint32_t Dim);
//
double Norm2D(double *X, uint32_t Dim);
//
void ConvertTimeS(uint32_t seconds);
//
double GetClosest(double TabRef[MaxRef], double Val);
//Get the closest entry to Val from the values in Tabref
//
double GetFreq0(uint8_t* Harm,int NMode, KForce KFC);
//Return Harmonic frequency sum_NM (d_i)*\nu_i
//
uint32_t MaxContribAll(float *Vec, uint32_t* IndexAssign, int *TabScreen, uint32_t *TabNull, int DoGraph,\
const uint64_t DimRez,const int NAdd, int MaxAdd, const float Trh, int NNotConv, int NScreenTot, uint32_t SizeRezMax);
//Return the most contributive indexes of residual vectors in IndexAssign.
//TabScreen contains the indexes of screened eigenvectors
//TabNull[nn] is equal to zero when residual contribution nn has been chosen for enrichment in previous iteration (else 1)
//Vec[nn*SizeRezMax] contains the eigenvector number nn.
//Nadd = number of elments to add at next iteration per non converged eigenpairs NNotConv.
//MaxAdd = Limit of number of elments to add at next iteration.
//Thr = Threshold above which a coordinate will be selected for enrichment
//
void GetIJKL(int32_t Rest, int NMode, int *ni, int *nj, int *nk, int *nl);
//Get rotational indexes ni, nj, nk, nl
//For a given number Rest=nl*NMode^3+nk*NMode^2+nj*NMode+ni
//
void ContribRot(uint8_t *Tester, int ni, int nj, int nk, int nl, int NMode);
//Increment by one each normal coordinate involved in rotational coefficients.
//computation and return an equivalent excitation array in Tester.
//Tester is the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors).
//
int32_t AssignRot(uint8_t *Tester, int NMode, int ni, int nj, int nk, int nl);
//Store the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors) in Tester  
//and return the equivalent number nl*NMode^3+nk*NMode^2+nj*NMode+ni
//
int CmptNNULL(uint8_t *Tab, int NMode);
//
double ScaleVect(float *Vec, int Iteration , SizeArray *Size, uint32_t SizeRezMax, int NScreenTot);
//
void InitCSC(CSC& IJ, uint32_t NNZ, uint32_t SizeActMax);
//
void DeleteLFF(LocalFF &LFF);
//
double ScualProd(double *X, double *Y, uint32_t Dim);
//
#endif
