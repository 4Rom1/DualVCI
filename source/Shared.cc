/* 
    Dual Vibration Configuration Interaction (DVCI) 
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
    along with this program. If not, see <https://www.gnu.org/licenses/>.

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

The rotational coefficients are computed according to the method of Meal and Polo [2]

[2]Vibration—Rotation Interaction in Polyatomic Molecules. II. The Determination of Coriolis Coupling Coefficients
Janet Hawkins Meal and S. R. Polo, The Journal of Chemical Physics 1956 24:6, 1126-1133 

*/
//
#include "Shared.h"
#include "Basic.h"
#include <string.h>
#include <sys/time.h>
//
void InitMode(ConfigId& M1,int Nmode)
{
   	int mm;
   	M1.Degrees=NULL;
   	M1.Degrees=new uint8_t [Nmode];
   	for (mm = 0 ; mm < Nmode; mm++)
    	{
     	M1.Degrees[mm] = 0;
    	} 
}
void InitLFF(LocalFF& LFF,int NFF)
{
   	int mm;
   	LFF.Num=0;
        LFF.Idx=NULL;
   	LFF.Idx=new int32_t [NFF];
   	for (mm = 0 ; mm < NFF; mm++)
    	{
     	LFF.Idx[mm] = 0;
    	}         
}
//
double GetEHarmonic0(uint8_t* Harm,int NMode, KForce KFC)
{//Return Harmonic energy sum_NM (d_i+0.5)*\nu_i
double Val=0;
for (int mm=0;mm<NMode;mm++)
 {
Val+=(Harm[mm]+0.5)*KFC.KijNCpld[1][mm];
 }
return Val;
}
//
double GetFreq0(uint8_t* Harm,int NMode, KForce KFC)
{//Return Harmonic frequency sum_NM (d_i)*\nu_i
double Val=0;
for (int mm=0;mm<NMode;mm++)
 {
Val+=(Harm[mm])*KFC.KijNCpld[1][mm];
 }
return Val;
}
//
 double factorial(uint32_t n)
 {
 	double val = 1; 	

      for (uint32_t ii = 1; ii <= n; ii++)
 		{val *= ii;}
 	
       return val;
 }
//
void InitKforce(KForce KFC, double **KijNCpld, double *KijCpld,\
int **Monm, int DegrePol, int NMode, int NPES)
{
   CopyMat(KFC.KijNCpld,KijNCpld, DegrePol, NMode)
   CopyMat(KFC.Monm,Monm,NPES, NMode)
   CopyTab(KFC.KijCpld,KijCpld,NPES)
}
//
void DeleteKForce(KForce& KFC, int DegrePol, int NPES)
{  
FreeMat(KFC.KijNCpld,DegrePol)
FreeMat(KFC.Monm,NPES)
FreeTab(KFC.KijCpld)
}
//
void InitCSC(CSC& IJ, uint32_t NNZ, uint32_t SizeActMax)
{
   uint32_t mm;
   IJ.I=NULL;
   IJ.NJ=NULL;
   IJ.I=new uint32_t [NNZ];
   IJ.NJ=new uint64_t [SizeActMax+1];
//
   for (mm = 0 ; mm <SizeActMax+1; mm++)
    {
     IJ.NJ[mm] = 0;
    }

   for (mm = 0 ; mm < NNZ; mm++)
    {
     IJ.I[mm] = 0;    
    } 
}
//
int DetectMaxCoordinate(double *Vec, int* PositionTarget, int NTargetStates, double ThrCoor)
{
//Return 1 if a coordinate listed in PositionTarget is bigger than ThrCoor in absolute value.
//Return 0 if no such coordinate have been detected when NTargetStates>0
//Return 1 when NTargetStates=0
for (int tt=0;tt<NTargetStates;tt++)
 {
  if(fabs(Vec[PositionTarget[tt]])>ThrCoor)
  {return 1;}
 }
  if(NTargetStates >0)
  {return 0;}
  else{return 1;}
}
//
uint32_t MaxContribAll(float *Vec, uint32_t* IndexAssign, int *TabScreen, uint32_t* TabNull, int DoGraph, \
const uint64_t DimRez,const int NAdd, int MaxAdd, const float Thr, int NNotConv, int NScreenTot,  uint32_t SizeRezMax)
{
//Return the most contributive indexes of residual vectors in IndexAssign.
//TabScreen contains the indexes of screened eigenvectors
//TabNull[nn] is equal to zero when residual contribution nn has been chosen for enrichment in previous iteration (else 1)
//Vec[nn*SizeRezMax] contains the eigenvector number nn.
//Nadd = number of elments to add at next iteration per non converged eigenpairs NNotConv.
//MaxAdd = Limit of number of elments to add at next iteration.
//Thr = Threshold above which a coordinate will be selected for enrichment
   int NScreenNC=0;
   int NNAdd=Min<int>(MaxAdd,NAdd);
   NScreenNC=0; 
   int ll=0;  
   int cmpt=0;
//   
   while (cmpt<NNAdd && NScreenNC < MaxAdd)
    {    
    ll=0;
    while (ll<NNotConv && NScreenNC < MaxAdd)
     {
    IndexAssign[NScreenNC]=FindMaximumAbs<float>(&Vec[TabScreen[ll]*SizeRezMax], DimRez);
// 
   for (int tt=0;tt<NScreenTot;tt++)
      {
     Vec[IndexAssign[NScreenNC]+tt*SizeRezMax]=0;
      }
     if(DoGraph)
      {TabNull[IndexAssign[NScreenNC]]=0;}
     NScreenNC++;
     ll++;
     }
    cmpt++;
    }

float Thrtmp=0;
uint32_t Indextmp=0;
cmpt=NScreenNC;
//uint32_t TmpScreen;
ll=0;
int go=1; 
   while (ll<NNotConv && cmpt < MaxAdd)
   {
    go=1;
    while(go)
    {
     Indextmp=FindMaximumAbs<float>(&Vec[TabScreen[ll]*SizeRezMax], DimRez);
     Thrtmp=Vec[TabScreen[ll]*SizeRezMax+Indextmp]*\
     SIGN<float>(Vec[TabScreen[ll]*SizeRezMax+Indextmp]);
    if(Thrtmp>Thr && cmpt < MaxAdd)
    {
    go=1;
    IndexAssign[cmpt]=Indextmp;
    for (int mm=0;mm<NScreenTot;mm++)
     {
    Vec[Indextmp+mm*SizeRezMax]=0;
     }
    if(DoGraph)
      {TabNull[Indextmp]=0;}
      cmpt++;
    }
    else{go=0;}
   }
    ll++;
  }
// 
   for (ll=0;ll<NScreenTot;ll++)
    {
     for (uint32_t nn=0; nn<DimRez; nn++)
     {
    Vec[nn+ll*SizeRezMax]=0;
     }
    }
return cmpt;
//
}
//
void CrossProd(double *Vec1, double *Vec2, double *Vec3){
Vec3[0]=Vec1[1]*Vec2[2]-Vec2[1]*Vec1[2];
Vec3[1]=Vec1[2]*Vec2[0]-Vec2[2]*Vec1[0];
Vec3[2]=Vec1[0]*Vec2[1]-Vec2[0]*Vec1[1];
}
//
void MassCenter(double *Mas, double *CoorCart, int NA, double *MassC)
{
//Compute the center of Mass. 
//X=(\sum_i m_ix_i)/(\sum_i m_i)
//Y=(\sum_i m_iy_i)/(\sum_i m_i)
//Z=(\sum_i m_iz_i)/(\sum_i m_i)
//
//x_1 =CoorCart[0], y_1 =CoorCart[1], z_1 =CoorCart[2], 
//x_2 =CoorCart[3], y_2 =CoorCart[4], z_2 =CoorCart[5],
//....
//x_n =CoorCart[3*(n-1)], y_n =CoorCart[3*(n-1)+1], z_n =CoorCart[3*(n-1)+2],
double MasTot=0;
for (int aa=0;aa<NA;aa++)
 {
MasTot+=Mas[aa];
 }
//Intialisation to zero
InitTab(MassC, 3)
for (int jj=0;jj<3;jj++)
{
for (int ii=0;ii<NA;ii++)
 {
MassC[jj]+=Mas[ii]*CoorCart[3*ii+jj];
 }
}
//
for (int xyz=0;xyz<3;xyz++)
 {
MassC[xyz]/=MasTot;
 }
}
//
double MomentInertie(double *Mas, double *CoorCart, int NA, double *MassC, double *PrincAx, double *PrincEig)
{
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
double Ixx=0;
double Iyy=0;
double Izz=0;
double Iyz=0;
double Ixz=0;
double Ixy=0;
//X_1 =CoorCart[0], Y_1 =CoorCart[1], Z_1 =CoorCart[2], 
//X_2 =CoorCart[3], Y_2 =CoorCart[4], Z_2 =CoorCart[5],
//....
//X_n =CoorCart[3*(n-1)], Y_n =CoorCart[3*(n-1)+1], Z_n =CoorCart[3*(n-1)+2],
for (int ii=0;ii<NA;ii++)
{
Ixx+=Mas[ii]*( pow(CoorCart[3*(ii)+1]-MassC[1],2)+ pow(CoorCart[3*(ii)+2]-MassC[2],2)  );
Iyy+=Mas[ii]*( pow(CoorCart[3*(ii)]-MassC[0],2)+ pow(CoorCart[3*(ii)+2]-MassC[2],2)  );
Izz+=Mas[ii]*( pow(CoorCart[3*(ii)]-MassC[0],2)+ pow(CoorCart[3*(ii)+1]-MassC[1],2)  );
//
Iyz-=Mas[ii]*( (CoorCart[3*(ii)+1]-MassC[1])*(CoorCart[3*(ii)+2]-MassC[2])    );
Ixz-=Mas[ii]*( (CoorCart[3*(ii)]-MassC[0])*(CoorCart[3*(ii)+2]-MassC[2])    );
Ixy-=Mas[ii]*( (CoorCart[3*(ii)+1]-MassC[1])*(CoorCart[3*(ii)]-MassC[0])    );
}
// T[jj][ii] = (ii*3+jj)=Iij=Iji 
PrincAx[0]=Ixx;
//
PrincAx[4]=Iyy;
//
PrincAx[8]=Izz;
//
PrincAx[1]=Ixy;
//
PrincAx[3]=Ixy;
//
PrincAx[2]=Ixz;
//
PrincAx[6]=Ixz;
//
PrincAx[5]=Iyz;
//
PrincAx[7]=Iyz;
//
double WORK[9];
int INFO;
int Dim=3;
int LWORK=9;
char JOBZ='V';
char UPLO='U';
//Get eigenpairs of moment of inertia matrix
dsyev_(&JOBZ,&UPLO, &Dim,&PrincAx[0], &Dim, &PrincEig[0], WORK, &LWORK, &INFO );
//Tensor constant to shift the fundamental value by U=-1/8\sum _{\alpha =1}^{3}\mu_{\alpha \alpha }
 return -(1.0/8.0)*(1/PrincEig[0]+1/PrincEig[1]+1/PrincEig[2])*HA_TO_CM;
}
//
int Zeta(double **CoorMode, int NMode, int NA, double *Mas, double *PrincAx, double *PrincEig, double **ZetaXYZ, int DoRot)
{
//Routines that compute the Zeta coefficients according to the method of Meal and Polo
//PrincAxes are eigenvectors of inertia momentum matrix
//PrincEig are eigenvalues of inertia momentum matrix
//Normal coordinate vectors are 
//CoorMode[ii][3*(jj)], CoorMode[ii][3*(jj)+1], CoorMode[ii][3*(jj)+2], 0 <= jj< NA, 0<= ii <NMode
//
//DoRot specifies non mass weighted normal coordinates for odd value and mass weighted for even value 
double Ui[3]={0};
double Vj[3]={0};
double W[3]={0};
//Rotational constants
double RotC[3]={0};
double Zeta1[3]={0};
double Zeta2[3]={0};
double ZetaT1[3]={0};
double ZetaT2[3]={0};
//-1/2 in front of the sum
RotC[0]=1/(2*PrincEig[0]);
RotC[1]=1/(2*PrincEig[1]);
RotC[2]=1/(2*PrincEig[2]);
//
double *ZetaX=new double[NMode*NMode];
double *ZetaY=new double[NMode*NMode];
double *ZetaZ=new double[NMode*NMode];
//
InitTabInt(ZetaX, NMode*NMode)
InitTabInt(ZetaY, NMode*NMode)
InitTabInt(ZetaZ, NMode*NMode)
//
double LocalNormii=0;
double ScalTest=0;
//
 for (int ii=0;ii<NMode;ii++)
  {   
        LocalNormii=0;
        for (int aa=0;aa<NA;aa++)
	 {
          if(DoRot%2) //2 possible type of normal coordinates : mass weighted or not.
          {
          CoorMode[ii][3*(aa)]*=sqrt(Mas[aa]);
          CoorMode[ii][3*(aa)+1]*=sqrt(Mas[aa]);
          CoorMode[ii][3*(aa)+2]*=sqrt(Mas[aa]);
          }  
	 LocalNormii+=pow(CoorMode[ii][3*(aa)],2);
	 LocalNormii+=pow(CoorMode[ii][3*(aa)+1],2);
	 LocalNormii+=pow(CoorMode[ii][3*(aa)+2],2);
        }
//
   if(fabs(LocalNormii-1)>1e-5)
   {
    printf("Normal coordinates eigenvectors are not normal : change DoRot\n\n");
    return 0;
   }
 }
 for (int ii=0;ii<NMode;ii++)
  {   
  for (int jj=ii+1;jj<NMode;jj++)
   {      
        ScalTest=0;
        for (int aa=0;aa<NA;aa++)
	 {      
	 ScalTest+=CoorMode[ii][3*(aa)]*CoorMode[jj][3*(aa)]+\
	 CoorMode[ii][3*(aa)+1]*CoorMode[jj][3*(aa)+1]+\
	 CoorMode[ii][3*(aa)+2]*CoorMode[jj][3*(aa)+2];        
         } 
   if(fabs(ScalTest)>1e-5)
   {
    printf("Normal coordinates eigenvectors are not orthogonals.\n\n");
    return 0;
   }
  }
 } 
//
for (int ii=0;ii<NMode;ii++)
 {
  for (int jj=ii+1;jj<NMode;jj++)
   {
     for (int aa=0;aa<NA;aa++)
      {
//Mass weighted displacements vectors. Careful the transformation is orthogonal
//Then the weighted vectors should have total norm equal to one.
//Zeta matrix definite with mass weighted coordinates.
//We have the GF transformation without mass weight.
//Then we have to scale the GF matrix by [sqrt(M1)*I3,...,sqrt(MN)*I3]
//
	Ui[0]=(CoorMode[ii][3*(aa)]);
	Ui[1]=(CoorMode[ii][3*(aa)+1]);
	Ui[2]=(CoorMode[ii][3*(aa)+2]);
//
	Vj[0]=(CoorMode[jj][3*(aa)]);
	Vj[1]=(CoorMode[jj][3*(aa)+1]);
	Vj[2]=(CoorMode[jj][3*(aa)+2]);
/*
(\zeta_{ii,jj}^{x}, \zeta_{ii,jj}^{y}, \zeta_{ii,jj}^{z})=(UixVj), 
(\zeta_{jj,ii}^{x}, \zeta_{jj,ii}^{y}, \zeta_{jj,ii}^{z})=(UjxVi), 
*/
//Vectorial product stored in W
	CrossProd(Ui, Vj, W);
//
	ZetaX[ii*NMode+jj]+=W[0];
//
	ZetaY[ii*NMode+jj]+=W[1];
//
	ZetaZ[ii*NMode+jj]+=W[2];
//Anti symmetry
	ZetaX[jj*NMode+ii]-=W[0];
//
	ZetaY[jj*NMode+ii]-=W[1];
//
	ZetaZ[jj*NMode+ii]-=W[2];
	}
     }
  }
//Init Matrix Zeta
InitMat(ZetaXYZ, (NMode*(NMode+1))/2, (NMode*(NMode+1))/2 )
//if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
 for (int ii=0;ii<NMode;ii++)
 { 
  for (int jj=ii+1;jj<NMode;jj++)
   {
    for (int kk=0;kk<NMode;kk++)
     {
      for (int ll=kk+1;ll<NMode;ll++)
       {
Zeta1[0]=ZetaX[ii*NMode+jj];
Zeta2[0]=ZetaX[kk*NMode+ll];
//
Zeta1[1]=ZetaY[ii*NMode+jj];
Zeta2[1]=ZetaY[kk*NMode+ll];
//
Zeta1[2]=ZetaZ[ii*NMode+jj];
Zeta2[2]=ZetaZ[kk*NMode+ll];
//
InitTab(ZetaT1,3)
InitTab(ZetaT2,3)
//
//ZetaT1=P^T Zeta1
      for (int ss=0;ss<3;ss++)
	{
	for (int mm=0;mm<3;mm++)
	 {
	ZetaT1[ss]+=PrincAx[ss*3+mm]*Zeta1[mm];
	ZetaT2[ss]+=PrincAx[ss*3+mm]*Zeta2[mm];
	 }
	}
// D^-1 P^T Zeta1
	for (int ss=0;ss<3;ss++)
	  {
	ZetaT1[ss]*=RotC[ss];
	  }
//Zeta2^T P D^-1 (P^T Zeta1)
//=zeta_ij^x zeta_kl^x RotC_x + zeta_ij^y zeta_kl^y RotC_y + zeta_ij^z zeta_kl^z RotC_z
	ZetaXYZ[ii + ((jj-1)*jj)/2][kk + ((ll-1)*ll)/2]=ScualProd(ZetaT2,ZetaT1, 3)*HA_TO_CM;
//ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}
	}
      }
    }
 }
//
delete [] ZetaY;
delete [] ZetaX;
delete [] ZetaZ;
//
return 1;
}
//
uint32_t GetSizeActMax2(int MaxEv,int MAXNCV, int NMode, int MaxQLevel, int DegrePol, int NPES, int Verbose,\
 int NXDualHPlus, double Mem, double KNNZ, double KNREZ, double KNZREZ, int SumTarget, uint32_t TotalLFF, uint32_t NXDualHTrunc)
{
//Return the maximal size of the active space B SizeActMax for a given memory Mem(in megabytes),
/*
NXDualHTruncPos-1: Total number of raising excitations in H* (after truncation through ThrKX)
NXDualHTrunc : Total number of excitations in H* (after truncation through ThrKX)
SizeRezMax=SizeActMax*(NXDualHTruncPos-1)*KNREZ
NNZActMax=SizeActMax*(KNNZ)*(NXDualHTrunc)
NNZRezMax=SizeActMax*(KNZREZ)*(NXDualHTrunc) 
*/
//TotalFF :Total number of force constants per excitation in dual operator
//MaxNCV : Maximal number of Lanczos basis set vectors
uint32_t UI=0;
double dble=0;
uint8_t U8=0;
int32_t I32=0;
float fl=0; 
int Targets=SumTarget;
double Diff;
if(SumTarget<1){Targets=MaxEv;}
double Memoir=(double)(Mem*(pow(1000,2)));
int NCV;
if(!MAXNCV)
{
 NCV=2*MaxEv;
}
else
{
 NCV=MAXNCV;//Maximal number of Lanczos basis set vectors
}
/*
if(Verbose)
{
printf("\n\n Size of int %lu, double %lu, float %lu, uint32_t %lu:  \n\n",\
sizeof(Intt),sizeof(dble),sizeof(fl),sizeof(UI));
}*/
if(Verbose)
{
printf("\n******Memory (bytes):%e*******\n",Memoir);
}
Diff=(sizeof(I32)*(TotalLFF+NXDualHPlus)\
+sizeof(U8)*(NXDualHPlus*NMode)+2*sizeof(UI)*NXDualHPlus); //LocalFF + Permuter
Memoir-=Diff;
if(Verbose)
{
printf("\n**Memory - Local Force Field=%e, substracted: %e**\n",Memoir,Diff);
}
Diff=sizeof(dble)*(NPES);
Memoir-=Diff; //PES
if(Verbose)
{
printf("\n**Memory - KForce:%e, substracted: %e**\n",Memoir,Diff);
}
Diff=sizeof(I32)*(NPES)*NMode;
Memoir-=Diff; //PES
if(Verbose)
{
printf("\n**Memory - monomial indexes=%e, substracted: %e**\n",Memoir,Diff);
}
Diff=sizeof(dble)*(DegrePol+1+NMode+4)*((MaxQLevel+DegrePol)*(MaxQLevel+DegrePol));//QQ
Memoir-=Diff;//QQ
if(Verbose)
{
printf("\n**Memory - 1d matrix elements=%e, substracted: %e**\n",Memoir,Diff);
}
Diff=sizeof(dble)*((NMode*(NMode+1))*(NMode*(NMode+1))/4);//Zeta
Memoir-=Diff;//Zeta
if(Verbose)
{
printf("\n**Memory - zeta coefficients=%e, substracted: %e**\n",Memoir,Diff);
}
Diff=2*sizeof(dble)*((MaxEv));//EigVal+EigvalOld
Memoir-=Diff;//EigVal
if(Verbose)
{
printf("\n**Memory - eigenvalues=%e, substracted: %e**\n",Memoir,Diff);
}
Diff=sizeof(dble)*((NCV)*(NCV+8));//EigVec
Memoir-=Diff;//EigVec
if(Verbose)
{
printf("\n**Memory - eigenvectors=%e, substracted: %e**\n",Memoir,Diff);
}
if(Memoir<0){return 0;}
double C=floor(Memoir);
//double A=(sizeof(UI)+sizeof(dble)/2)*KNNZ;
double B=ceil((sizeof(UI)+sizeof(dble))*KNNZ*(double)NXDualHTrunc)+\
ceil(sizeof(UI)*KNZREZ*(double)NXDualHTrunc)+\
ceil(sizeof(fl)*(double)(NXDualHPlus*Targets)*KNREZ)+\
ceil(sizeof(dble)*(6+NCV))+\
ceil(sizeof(UI)*(3*(double)NXDualHPlus*KNREZ+6))+\
ceil(sizeof(U8)*(NMode+(double)(NXDualHPlus*NMode)*KNREZ));
if(KNZREZ<=0)
{B-=ceil(sizeof(UI)*((double)NXDualHPlus*KNREZ));}
if(Verbose)
{
printf("\n**SizeActMax=C/B with**\n");
printf("**C=%e, B=%e, C/B=%e**\n",C,B,C/B);
}
//ValAct + IJAct = S2*(sizeof(UI)+sizeof(dble)/2)*KNNZ)
//+S*(Tester out + RezVec)
//+S*(WorkD+Eigenvectors)
//+S*(Max PermutAct + Nulify+Permuter+Corres, IJRez)
//+S*(Active space, Residual space)
//IndexDiag? : 4*NXDualHPlus=5*NXDualHPlus
//CSC U64*2=4UI+3*Permut
double Sol=floor((double)C/(double)B);
//
return (uint32_t)Sol;
}
//
float EffectiveMemory2(char PidName[MaxChar])
{
//Execute a terminal command ps aux | grep PidName
//and convert the output in megabytes
  FILE *FileCommand=NULL;
  char path[1035]={0};
//
  char Exec[MaxChar]={0};
  sprintf(Exec,"ps aux | grep %s | grep -v grep | awk 'BEGIN { sum=0 } {sum=sum+$6;} END\
 { print sum / 1024 }'",PidName);
  float Memoire=0;
  /*Execute the command specified by the string Exec.*/
  FileCommand = popen(Exec, "r");
  if (FileCommand == NULL) {
    printf("Failed to execute command\n" );
    return 0;
   }
  /* Read one output line at a time. */
  while (fgets(path, sizeof(path)-1, FileCommand) != NULL) 
  {  
    sscanf(path, "%f", &Memoire);
  }
//
  pclose(FileCommand);
  return Memoire; 
}
//
void swap (uint32_t* a,uint32_t* b )
{
    uint32_t t = *a;
    *a = *b;
    *b = t;
}
//
double Norm2D(double *X, uint32_t Dim)
{
double nrm=0;
for (uint32_t nn=0;nn<Dim;nn++)
nrm=nrm+pow(X[nn],2);
//
return sqrt(nrm);
}
//
float Norm2F(float *X, uint32_t Dim)
{
float nrm=0;
for (uint32_t nn=0;nn<Dim;nn++)
nrm=nrm+pow(X[nn],2);
//
return sqrtf(nrm);
}
//
void ConvertTimeS(uint32_t seconds)
{
int days = seconds / (60 * 60 * 24);
seconds -= days * (60 * 60 * 24);
int hours = seconds / (60 * 60);
seconds -= hours * (60 * 60);
int minutes = seconds / 60;
seconds -= minutes * (60);
printf("%i days %i hours %i minutes %i seconds\n",days,hours,minutes,seconds);
}
//
double GetClosest(double TabRef[MaxRef], double Val)
{
//Get the closest entry to Val from the values in Tabref
double Eps=4000;
double ValOut=TabRef[0];
for(int ii=0; ii<MaxRef; ii++)
 {
 if(fabs(Val-TabRef[ii])<Eps)
  {
  ValOut=TabRef[ii];
  Eps=fabs(Val-TabRef[ii]);
  }
 }
return ValOut;
}
//
void GetIJKL(int32_t Rest, int NMode, int *ni, int *nj, int *nk, int *nl)
{
//Get rotational indexes ni, nj, nk, nl
//For a given number Rest=nl*NMode^3+nk*NMode^2+nj*NMode+ni
*nl=Rest/Powo<int32_t>(NMode,3);
Rest-=*nl*Powo<int32_t>(NMode,3); 
*nk=Rest/Powo<int32_t>(NMode,2);
Rest-=*nk*Powo<int32_t>(NMode,2); 
*nj=Rest/NMode;
Rest-=*nj*NMode; 
*ni=Rest;
}
//
void ContribRot(uint8_t *Tester, int ni, int nj, int nk, int nl, int NMode)
{
//Increment by one each normal coordinate involved in rotational coefficients.
//computation and return an equivalent excitation array in Tester.
//Tester is the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors).
InitTabInt(Tester, NMode)
Tester[ni]++;
Tester[nj]++;
Tester[nk]++;
Tester[nl]++;
}
//
int32_t AssignRot(uint8_t *Tester, int NMode, int ni, int nj, int nk, int nl)
{
//Store the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors) in Tester  
//and return the equivalent number nl*NMode^3+nk*NMode^2+nj*NMode+ni
InitTabInt(Tester, NMode)
Tester[ni]++;
Tester[nj]++;
Tester[nk]++;
Tester[nl]++;
return ni+NMode*nj+Powo<int32_t>(NMode,2)*nk+Powo<int32_t>(NMode,3)*nl;
}
//
int CmptNNULL(uint8_t *Tab, int NMode) {
   int mon_i=0;
   int NKey=0;   
   for (mon_i=0; mon_i<NMode;mon_i++)
    { 
     if(Tab[mon_i]!=0)
     {NKey++;}
    }
    return NKey;
}
//
void DeleteLFF(LocalFF &LFF)
{   
  delete[] LFF.Idx; 
}
//
double ScualProd(double *X, double *Y, uint32_t Dim)
{
double Scual=0;
for (uint32_t nn=0;nn<Dim;nn++)
{Scual+=X[nn]*Y[nn];}
return Scual;
}
