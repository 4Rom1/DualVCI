/*
    Dual Vibration Configuration Interaction (DVCI).
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "Shared.h"
#include "Graph.h"
#include "ReadData.h"
#include "Assemblage.h"
#include "Basic.h"
#include "Solver2.h"
#include "Generator.h"
#include <cassert>
//#include <omp.h> might be used in parallel version
#include <limits.h>
#include <vector>
#include <float.h>
#include <algorithm> 
//
#define IncK2 1
//IncK2:Say if the the quadratic term should be included or not in the local force field,
//also say if the harmonic energy should be added to diagonal elements (IncK2=0) or not (IncK2=1)
int main(int argc, char *argv[])
{
//Files to read    
    FILE *FileKey = NULL;
    FILE *FilePES = NULL;
//
    long LinesKijCpld=0;//Number of coupled force constants
    uint8_t TargetState[MaxTarget][MaxNormal]={0};
//To show time spend    
    struct timeval BeginAll, EndAll;
    clock_t CPUBeginAll, CPUEndAll, CPUTime;
    uint32_t time_spent=0;
//
PrintBanner ("--dvci--");
//
printf("Dual vibration configuration interaction (DVCI)\n\
A novel factorisation of molecular Hamiltonian for\n\
high performance infrared spectrum computation,\n\
Romain Garnier, submitted in 2018.\n");
//
     if(argc != 2)
     { 
      printf("\n***Usage: ../../source/Transitions Key.in***\n");
      printf("One argument requested : Key file missing.  \n");
      printf("Dipole surface file should also be in the same directory. \n\n");
      FinalMessaj()
      return 0;
     } 
//
     FileKey=fopen(argv[1],"r");//Key words to read
//        
	if(FileKey==NULL)
	{
	printf("\n***Wrong or inexistant value of Key file***\n\n");
        FinalMessaj()
	return 0;
	} 
//    
int NMode=0;//Number of normal coordinates
int DegrePol=0;//Maximal degree in PES : not given as input will be set after reading of PES 
int DoRot=1;//If > 0 then Coriolis terms are computed
// with non mass weighted normal coordinates for odd and mass weighted for even 
int PESType=1;//if non 0 derivatives in a.u else force constants in cm-1
int NAtome=0;//Number of atomes
int MAXNCV=0;//Number of ritz vectors for the eigensolver
int MaxQLevel=15;//Maximal quantum level in each direction 
//Be carefull when increasing potentials might be very poluted around 20 
double MinFreq=-100;// [MinFreq,MaxFreq]= frequency window target
double MaxFreq=4000;// 
double GroundState=-1;//Groundstate should be given by the user.
//except when it is actually computed then it has negative value.
double Kappa=1.2;//Multiplicatif factor for intial subspace construction 
double Tol=1e-8;//Threshold for the eigensolver (same as tol in DSAUPD) 
uint32_t SizeActMax=0;//Maximal size for active space B
int NAdd=10;//Number of added elements at each iteration  
int MaxAdd=1000;//Maximal number of added elements at each iteration  
double EtaComp=1.5;//New basis set will be selected from the maximal component
// of error vector divided by EtaComp*NNotConv, where NNotConv is the number of
// non converged eigenvectors at a given iteration
char OutName[MaxChar]="NoName";
double ThrMat=DBL_EPSILON; //Threshold for matrix elements
double ThrCoor=0.6; //Threshold for component selection of targeted vectors
int AddTarget=2; //Number of additional residual vectors to correct the potential increasing of targets
int DegreEx[MaxCpld]={0};//Degree of excitation in each direction for the generator used here only to show values 
char PESName[MaxChar]={0}; //Name of file of PES to be read.
char RefName[MaxChar]={0};//Reference txt file containing frequencies in increasing order
int MaxEv=30;//Maximal number of wanted eigenvalues for the Arnoldi eigensolver
int DeltaNev=1000; //Adjustment number for MaxEV that will be equal to Min(MaxEV,  MaxScreen+DeltaNev) 
//where PosEig is the eigenvalue position of maximal target
int NCPol=0;//Maximal couplings in polynomial: not given as input will be set after reading of PES. 
double Memory=0;//Total allocated Memory in Megabytes
int PrintOut=0; //If >0 print final basis set and eigenvectors for post treatment 
double EpsRez=6e-3; //Convergence criteria for global residue
int Verbose=0;//If > 0 printout additional informations between iterations
float ThrKX=1; //Threshold for contributions of sum force constants in dual operator
double ThrPES=DBL_EPSILON; //Threshold for derivative or force constants of the PES
double Freq0Max=30000;//Frequency wall above groundState
double KNNZ=0.03;// Shrinking factor for the NNZ of active matrix Hb
double KNREZ=0.2;//Shrinking factor for residual space size
double KNZREZ=0.3;//Shrinking factor for residual matrix NNZ
int EvalDeltaE=0;//If positif add correction vector and correction energy to targets
/*
NXDualHTruncPos : Total number of raising excitations for the dual of H (after truncation through ThrKX)
NXDualHTrunc : Total number of excitations for the dual of H (after truncation through ThrKX)
SizeRezMax=SizeActMax*(NXDualHTruncPos-1)*KNREZ
NNZActMax=SizeActMax*(KNNZ)*(NXDualHTrunc)
NNZRezMax=SizeActMax*(KNZREZ)*(NXDualHTrunc) 
*/
int DegreCoupl[MaxCpld]={0}; //Maxdegree of PES for each coupling: starts at 0 for coupling 1
int DoGraph=1; // If not equal to zero, it will store the graph of residual matrix in CSC format
//to gain computational time and then significantly increases memory requirement
//
printf("\n*****Reading Datas*****\n");
//
int NTargetStates=GetKeyWords(FileKey, &NMode, &DoRot, &MaxEv, &DeltaNev, &AddTarget,\
    &PESType, TargetState, &DoGraph, &MAXNCV,\
    &MaxQLevel, &NAdd, &MaxAdd, &PrintOut, &EvalDeltaE, &EpsRez, &KNREZ, &KNNZ, &KNZREZ, &EtaComp,\
     &Tol, &Kappa, &ThrMat, &ThrPES, &ThrCoor, &GroundState, &MinFreq, &MaxFreq, &Freq0Max, &ThrKX,\
    &Memory, &Verbose, OutName, PESName, RefName);
//
if(NMode==0)
{
printf("\n***NMode must be specified***\n\n");
FinalMessaj()
return 0;
}
//
NAtome=(NMode+6)/3;
//
printf("\n Output file name:%s \n",OutName);
//
printf("\n Number of targeted states %i \n",NTargetStates);
printf("\n PES file name %s \n",PESName);
//
      printf("\n Number of normal coordinates NMode = %d \n",NMode);
      printf("\n PESType %d \n", PESType);
      MaxQLevel++; //Because < in Matrices
      printf("\n Maximal quantum level for each direction MaxQLevel = %d \n",MaxQLevel-1);
      if(ThrPES >0 && ThrMat>0)
      {
      //Cannot be under error machine.
      ThrPES=Max<double>(DBL_EPSILON,ThrPES);
      ThrMat=Max<double>(DBL_EPSILON,ThrMat);
      printf("\n Threshold for PES ThrPES = %e,\n Threshold for Hamiltonian matrix elements ThrMat = %e \n", ThrPES, ThrMat);}
      else{printf("\n Threshold for PES and for Hamiltonian matrix elements have to be > 0 \n");FinalMessaj() return 0;
      }
//
      if(NTargetStates)
      {
       if((ThrCoor >= 1)  || (ThrCoor<0)) 
       {printf("\n ***Threshold for coordinate must be lower than one and positive*** \n");FinalMessaj() return 0;}
       else 
       {printf("\n Threshold for coordinate selection ThrCoor = %e \n", ThrCoor);}
      }

//
       printf("\n EpsRez = %.4f -> maximal allowed value of relative residues |R|/E \n",EpsRez);
//
       printf("\n EtComp: component selected bigger than |R|/(%.1f*NNotConv) \n",EtaComp);
//
       printf("\n NAdd:%d, NAdd*(Iteration+1) components added\n \
per non converged eigenpairs at each iteration, bounded by %d \n",NAdd,MaxAdd);
//
       printf("\nMaximal allowed harmonic frequency in B+Bs Freq0Max = %.1f cm-1\n", Freq0Max);
//
       if(!DoGraph)
       {
       printf("\nDoGraph Null means no storage of residual graph: more cpu time, less memory\n");
       KNZREZ=0;
       }
       else
       {
       printf("\nDoGraph non null means storage of residual graph: more memory, less cpu time\n");
       }
//
       printf("\nShrinking factor for residual space KNREZ = %.1f,\n\
 and shrinking factor for the maximal number of indexes in the residual graph KNZREZ =  %.1f\n", KNREZ, KNZREZ);
//
if(EtaComp>100)
{
printf("\n*** EtaComp too big then lowered to 100 ***\n\n");
EtaComp=100;
}
//
if(KNREZ<=0 || KNREZ >1)
{
printf("\n***KNREZ must be in ]0,1]***\n\n");
FinalMessaj()
return 0;
}
//
if(KNNZ<=0 || KNNZ >1)
{
printf("\n***KNNZ must be in ]0,1]***\n\n");
FinalMessaj()
return 0;
}
//
if(MaxQLevel>=MaxDegrePI)
{
printf("\nMaxQLevel too big, only handle until %i \n\n",MaxDegrePI);
FinalMessaj()
return 0;
}
//
if( DoRot < 0)
{
printf("\nDoRot must be positif odd or even \n\n");
FinalMessaj()
return 0;
}
//
if(NTargetStates>=MaxTarget)
{
printf("\nMaximal number of targets limited to %i \n\n",MaxTarget);
FinalMessaj()
return 0;
}
//
if(NMode>=MaxNormal)
{
printf("\nMaximal number of normal coordinates limited to %i \n\n",MaxNormal);
FinalMessaj()
return 0;
}
//
if(!PrintOut)
{
printf("\n***PrintOut must be equal to 1 and basis set saved from previous DVCI run*** \n\n");
FinalMessaj()
return 0;
}
      if(Memory<=0)
       {
           printf("****Need to indicate a positive value of Memory in megabytes****\n\n");
           FinalMessaj()
           return 0;
       }
//  
       printf("\n Allocated memory %.1f Mo \n",Memory);

      if(GroundState > 0 )
      {
       printf("\n Groundstate provided %f \n", GroundState);   
      }
//
FilePES=fopen(PESName,"r");
//
if(FilePES==NULL)
{
printf("***Wrong or inexistant value of PESName***\n");
FinalMessaj()
return 0;
}
//
if(NTargetStates)
{printf("\n Target states readed : \n");}
//
for (int tt=0;tt<NTargetStates;tt++)
{
AfficheNu(TargetState[tt], NMode);
if(tt<NTargetStates-1)
{printf(" , ");}
if(!((tt + 1) % 3)){printf("\n");}
}
//
printf("\n");
//
//    
      if(PESType)//Derivatives in a.u
      {LinesKijCpld=GetNumKijCpld(FilePES, NMode, ThrPES, &NCPol, &DegrePol, DegreCoupl);}
      else//Force constants in cm-1
      {LinesKijCpld=GetNumKijCpld0(FilePES, NMode, ThrPES, &NCPol, &DegrePol, DegreCoupl);}
//
      printf("\nMax couplings NCPOL %i and max degree (sum of monomial indexes) DegrePol %i\n",NCPol,DegrePol);
//
if(NCPol>=MaxCpld)
{
printf("\nNCPol too big, only handle until %i \n\n",MaxCpld);
FinalMessaj()
return 0;
}
//
      if(LinesKijCpld<0){FinalMessaj() return 0;}
//
if(pow(NMode,4)+LinesKijCpld>INT32_MAX)
{
printf("\nNumber of NNUL PES too big use int64_t to store them instead of int32_t\n\n");
FinalMessaj()
return 0;
} 
    const int32_t NPES(LinesKijCpld);
//Force constant variables
    int32_t** Monm=NULL;
    double** KijNCpld=NULL;
    double* KijCpld=NULL;
//
allocate_MatInt(Monm, NPES, NMode)
allocate_MatDb(KijNCpld, DegrePol, DegrePol*NMode)
allocate_TabDb(KijCpld, NPES)
//
//Init Arrays
InitMat(Monm, NPES, NMode)
InitMat(KijNCpld, DegrePol, DegrePol*NMode)
InitTabInt(KijCpld, NPES)      
//
double **CoorMode=NULL;
allocate_MatDb(CoorMode, NMode, 3*NAtome)
double *CoorCart=NULL;
CoorCart=new double[3*NAtome];
double *Mas=NULL;
Mas=new double[3*NAtome];
//center of masse, principal axes and principal moments  
double MassC[3]={0};
double PrincAx[9]={0};
double PrincEig[3]={0};
//Get force constants of the PES and check orthogonality of normal modes
 int DoRotCheck=Min<int>(GetDataPES(FilePES, KijCpld,KijNCpld, Monm, CoorCart,Mas, CoorMode,\
 NMode, NAtome, PESType, Verbose, ThrPES),DoRot);
//
 if(!DoRotCheck){DoRot=0;}
//
 if(DoRot<0){FinalMessaj() return 0;}
//
  if(DoRot)
   {printf("\n *** Coordinates, masses and equilibrium geometry readed *** \n");}
//
//DoRot must be set up to zero anyway, because Coriolis not part of operator
  DoRot=0;
  double *Omega=new double[NMode];//Omega is used as converter from Hartree to cm-1
//
uint8_t *Pid=NULL;//Degrees in each direction
Pid=new uint8_t[NMode];
//
  printf("\n-Here the harmonic frequencies \n");
for (int mm=0; mm<  NMode ; mm ++) 
  {   
    if(PESType)
    {
     Omega[mm]=pow(KijNCpld[1][mm],1.0/4.0);
     printf(" \\nu_{%d} : %.1f cm-1, ",mm+1,pow(Omega[mm],2)*HA_TO_CM);
     Pid[mm]=Min<uint32_t>((uint32_t)ceil(Freq0Max/(pow(Omega[mm],2)*HA_TO_CM)),MaxQLevel);
    }
    else
    {
     Omega[mm]=2*KijNCpld[1][mm];
     printf(" \\nu_{%d} : %.1f cm-1, ",mm+1,Omega[mm]);
     Pid[mm]=Min<uint32_t>((uint32_t)ceil(Freq0Max/Omega[mm]),MaxQLevel);  
     Omega[mm]=sqrt(Omega[mm]/HA_TO_CM);//Inverse conversion for Coriolis terms
    }
    if( ( (mm+1)%3==0) && (mm>0))
    {printf("\\\\\n");}
  }
  printf("\n");
//
printf("Quantum levels by coordinate\n");
printf("[");
for (int mm=0; mm<  NMode-1; mm ++) 
  { 
       printf("%u, ",Pid[mm]-1);  
  }
printf("%u]\n",Pid[NMode-1]-1);
//
 double KFMax=AdjustCoeff(KijNCpld, KijCpld,DegrePol, NPES, PESType, Monm, NMode);
//   
  printf("\nMaximal Force constant %f (exclude Harmonic)\n",KFMax); 
//   
      if((ThrKX > 0) && (ThrKX < KFMax) )
      {printf("\n Threshold for acceptable excitations in H* : ThrKX = %e ,\n", ThrKX);}
      else
      {printf("\n ***Threshold for H* must be lower than KFMax and strictly positif*** \n"); FinalMessaj() return 0;}
//Force constants:
//KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES. 
   KForce KFC;  
   allocate_MatDb(KFC.KijNCpld,  DegrePol, NMode)
   allocate_TabDb(KFC.KijCpld, NPES)
   allocate_MatInt(KFC.Monm, NPES, NMode)
//
   InitKforce(KFC, KijNCpld, KijCpld,\
   Monm,  DegrePol, NMode, NPES);//Gather Force constant data's in KFC structure. 
//
       FreeMat(Monm,NPES)  
       FreeTab(KijCpld)
       FreeMat(KijNCpld,DegrePol) 
//
if(PrintOut >2)
{
printf("\n***Compute size of pruned basis set \n");
int NMillion=1;
uint8_t MultiDegrees[MaxNormal]={0};
uint8_t MultiSurf[MaxNormal]={0};
uint8_t MaxDegrees[MaxNormal]={0};
uint8_t MaxDegSurf[MaxNormal]={0};
//uint8_t Tester[MaxNormal]={0};
uint64_t SizePrun=0;
double SubFreq=0;
double ffreq=0;
//
double AdjustFreq; 
int ss;
//
for (int NEccit=1;NEccit<=NMode;NEccit++)
 {
 for (int dd=0; dd< NEccit; dd++)
   {
  MaxDegSurf[dd]=1;
   }
//
 for (int dd=NEccit; dd< NMode; dd++)
   {
  MaxDegSurf[dd]=0;
   }
//
//
//Combination i<j<k<l...
//Surface i<j<l<m... until NEccit ={ MultiSurf[ss], 0<=ss<NEccit }
//Given by the indexes xx of MaxDegSurf[xx]!=0
//MultiSurf[ss]<NMode contains the indexes of surface for 0<=ss<NEccit
//MaxDegrees[ss] will be the maximal degree of modal excitation on coordinate MultiSurf[ss]
printf("\n------------%d-Modes couplings surface-------------\n",NEccit);
do{
  SubFreq=0;
//  InitTabInt(Tester, NMode)  
    for (int mm=0; mm< NMode; mm++)
    {
   if(MaxDegSurf[mm])
     {
     SubFreq+=KFC.KijNCpld[1][mm];
     }
    }
  ss=0;  
  AdjustFreq=Freq0Max-SubFreq;
  if(AdjustFreq>=0)
  {    
  for (int mm=0; mm< NMode; mm++)
    {
   if(MaxDegSurf[mm])
     {
//Freq0Max-SubFreq+KFC.KijNCpld[1][mm] explains the +2 yes but if 0.5 maximal 1 can be included
   MaxDegrees[ss]=Min<int>((int)ceil((AdjustFreq)/KFC.KijNCpld[1][mm])+1,MaxQLevel);
//printf("Pid[%d]:%d,",mm+1,MaxDegrees[ss]);  
   MultiDegrees[ss]=1;
   MultiSurf[ss]=mm;
   ss++;
     }
    }
//   
// printf("\n");
   do{    
//Sum_ee Multidegrees[ee] is all the possible sum of NEccit uplets lower than MaxDegrees+1.
    ffreq=0;
    for (int ee=0;ee<NEccit;ee++)
     {
    ffreq+=MultiDegrees[ee]*KFC.KijNCpld[1][MultiSurf[ee]];
     }    
   if(ffreq<Freq0Max)
      {   
      SizePrun++;
        if(!(SizePrun % (int)1e6))
        {
     printf("%i millions\n",NMillion);
     NMillion++;
        }
       }  
     }while(nested_loop1(MultiDegrees, MaxDegrees, NEccit));  
    }//Adjust freq need to be positif 
  }while(std::prev_permutation(MaxDegSurf, MaxDegSurf+NMode));
}
//
printf("\n\n***The size is %lu***\n\n",SizePrun+1);//Do not forget 0
FinalMessaj() return 0;
}


if(DoRot)
{
//Compute center of Mass
MassCenter(Mas, CoorCart, NAtome, MassC);
//Compute principal axes of moment of inertia: Axes and values
double Ull=MomentInertie(Mas, CoorCart, NAtome, MassC, PrincAx, PrincEig);//Watson Correction
printf("\n*** Watson term = -1/8*(Sum U_ll) %.4f ***\n\n",Ull); 
if(Verbose)
{
printf("\n*** Center of masses ***\n");
for (int jj=0;jj<3;jj++)
{
printf(" Coordinates %d , %e\n",jj+1,MassC[jj]);
}
//
printf("\nMoment of Inertia principal eigenvalues\n");
afficheTabDb(PrincEig,3)
printf("\nconverted from me.bhor^2 to amu.A^2\n");
    for (int ss=0;ss<3;ss++)
    {
    printf(" %e\n",PrincEig[ss]*pow(BOHR_TO_ANGST,2)*ME_TO_AMU);
    }
//
printf("\nRotational constant in cm-1\n");
  for (int ss=0;ss<3;ss++)
  {
  printf(" %f\n",MEBOHR2_TO_CM/PrincEig[ss]);
  }
 }
}
//Compute Zeta elements
// ZetaXYZ[i + (j-1)*j/2][k + (l-1)*l/2]=
//\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta} , i<j, k<l
double **ZetaXYZ=NULL;
allocate_MatDb(ZetaXYZ, ((NMode+1)*NMode)/2, ((NMode+1)*NMode)/2)
if(DoRot)
{
 if(Zeta(CoorMode, NMode, NAtome, Mas, PrincAx, PrincEig, ZetaXYZ, DoRot)<=0)
  {
  FinalMessaj();
  return 0;
  }
}
//Maximal number of neighboors and maximal size of the arrays
    uint32_t NGenPoz=NMode*DegreCoupl[0]+1;//One Modal excitation+Zero excitation
//
   printf("\nMaximal degree for 1-mode couplings = %d\n",DegreCoupl[0]); 
//        
   for (int nn=2;nn<=MaxCpld;nn++)
    {
    if(DoRot && nn<=4) //Rotational terms may increase the couplings, except for water...
     {
     if(NMode>3)
      {
      DegreCoupl[nn-1]=Max<int>(DegreCoupl[nn-1],4);
      NCPol=Max<int>(NCPol,4);
      }
     }
     if(DegreCoupl[nn-1]) 
     {
     NGenPoz+=GetSizeEccitPoz(NMode,nn,DegreCoupl[nn-1]);  
     printf("Maximal degree for %d-mode couplings = %d\n",nn,DegreCoupl[nn-1]); 
     }     
    }
   //Allocation for generator XPolSupPos[Idm(xx)+mm], (xx,mm) in [0,NGenPoz[ x [0,NMode[
    allocate_TabMode(XPolSupPos,NGenPoz) //Set of upper limit of positive excitations in H*
    uint32_t SizeGenOld=GeneratorTotPoz (NMode, XPolSupPos, DegreCoupl, NCPol);
      if(NGenPoz != SizeGenOld)
       { 
      printf("**Modal Excitations not equal between each other SizeGenOld %d, NGenPoz %u**\n\n",SizeGenOld,NGenPoz);
      FinalMessaj()
      return 0;
       }
    if(Verbose)
    {
    printf("\nUpper limit of raising excitations in H* NGenPoz=%u\n\n",NGenPoz);
    }     
    int *NLFFPerEx=NULL;
    NLFFPerEx=new int[NGenPoz]; 
    InitTab(NLFFPerEx, NGenPoz)  
    uint32_t NXDualHPlus=CptLFFPerEx (DegrePol, NMode, NPES, DoRot, NGenPoz,\
    KFC, XPolSupPos, NLFFPerEx, ZetaXYZ, ThrPES, IncK2);
//
if(!NLFFPerEx[0]){NXDualHPlus++;}//If zero excitation not included 
//
 for (uint32_t ll=0; ll<NGenPoz; ll++)
 { 
 if(NLFFPerEx[ll])
  {
//""DegreEx[CmptNNULL(XPolSupPos[ll].Degrees, NMode)]++; 
    DegreEx[CmptNNULL(&XPolSupPos[Idm(ll)],NMode)]++;    
  }
 } 
    if(Verbose)
    {
       for (int dd=1; dd < MaxCpld; dd++)
        {
        if(DegreEx[dd])
         {
        printf("Number of %d-mode excitations %d\n",dd,DegreEx[dd]);
         }
        }
     } 
//
 allocate_TabMode(DualHPos,NXDualHPlus)  
//DualHPos[Idm(xx)+mm], (xx,mm) in [0,NXDualHPlus[ x [0,NMode[.
//""    
 LocalFF LFF;
 LFF.Num=NULL;
 LFF.Num=new uint32_t[NXDualHPlus+1];
//
/*LFF[xx] : local force field associated with positive excitation DualHPos[xx].
LFF[xx].Idx[ii] < 0 are indexes of non coupled force constants:
KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
0 <= LFF[xx].Idx[ii] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF[xx].Idx[ii]].
LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
with a unique corresponding (ni,nj,nk,nl)  
LFF[xx].Idx[ii] are defined for ii in [0,LFF[xx].Num[

LFF.Idx[ii(xx)] are defined for ii [LFF.Num[xx],LFF.Num[xx+1][*/
//
   NXDualHPlus=0;
   uint32_t TotalLFF=0;//Total number of force constants per excitation in dual operator
//
//Guaranty zero excitation
//
   if(!NLFFPerEx[0])
   {
   printf("*******No zero excitation**********\n");
//""CopyTab(DualHPos[NXDualHPlus].Degrees,XPolSupPos[0].Degrees,NMode)
   memcpy(&DualHPos[Idm(NXDualHPlus)],&XPolSupPos[0],(size_t)NMode);
   LFF.Num[NXDualHPlus]=0;
   NXDualHPlus++;
   TotalLFF++;
   }   
//   
int MaxTerm=0;
//
 for (uint32_t ll=0; ll<NGenPoz; ll++)
 { 
 if(NLFFPerEx[ll])
  {
//""CopyTab(DualHPos[NXDualHPlus].Degrees,XPolSupPos[ll].Degrees,NMode) 
   memcpy(&DualHPos[Idm(NXDualHPlus)],&XPolSupPos[Idm(ll)],NMode);
//""InitLFF(LFF[NXDualHPlus],NLFFPerEx[ll]);
//LFF.Idx[ii(xx)] are defined for ii [LFF.Num[xx],LFF.Num[xx+1][ then
   LFF.Num[NXDualHPlus]=TotalLFF;
   NXDualHPlus++;
   TotalLFF+=NLFFPerEx[ll];
//
   if(NLFFPerEx[ll]>MaxTerm && ll)
    {
    MaxTerm=NLFFPerEx[ll];
    }
   }
  }

//""Add last term
   LFF.Num[NXDualHPlus]=TotalLFF;
   LFF.Idx=new int32_t[TotalLFF];
//

 int AvgTerm=TotalLFF/NXDualHPlus;
//Evaluated number of operation per matrix element
 uint32_t AvgOp=(uint32_t)log((double)(NXDualHPlus))*(double)NMode+AvgTerm*NMode;
 uint32_t MaxOp=(uint32_t)log((double)(NXDualHPlus))*(double)NMode+MaxTerm*NMode;
  if(Verbose)
  {
  printf("\nAverage number of terms in the local force fields, AvgTerm = %u anf total %u \n",AvgTerm,TotalLFF);
  printf("\nMaximal number of terms in the local force fields (exclude 0), MaxTerm = %u\n",MaxTerm);
  printf("\nNumber of raising excitations in H* NXDualHPlus = |LFF*| = %u\n",NXDualHPlus);
  printf("\nAverage number of operation per matrix element\n (log(|LFF*|)+AvgTerm)*NMode = %u\n",AvgOp);
  printf("\nMaximal number of operation per matrix element\n (log(|LFF*|)+MaxTerm)*NMode = %u\n",MaxOp);
  }
// 
//
 uint32_t *Corres=NULL;
 Corres=new uint32_t[NXDualHPlus];
//
//Permutation for generator
//
      uint32_t *PermutGen=NULL;
      PermutGen=new uint32_t [NXDualHPlus];
//
 uint32_t NXDualHTruncPos=AssignLFFPerEx (DegrePol, NMode, DoRot, NPES, NXDualHPlus,\
 KFC, DualHPos, LFF, ZetaXYZ, ThrKX, ThrPES, PermutGen, Corres, IncK2);
// 
 if(Verbose)
 {
 printf("\nTotal number of raising excitations in H*\n(after truncation through ThrKX) NXDualHTruncPos-1 = %u\n",NXDualHTruncPos-1);
 }
//
uint32_t NXDualHTrunc=GetNXDualHTrunc(DualHPos,NXDualHTruncPos,Corres, NMode);
 if(Verbose)
 {
  printf("\nMaximal number of excitations in H*\n(after truncation through ThrKX) NXDualHTrunc = %u \n",NXDualHTrunc);
 }
//
   delete [] NLFFPerEx;
//""
   delete [] XPolSupPos;
//FreeTabMode(XPolSupPos,NGenPoz)  
//  
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Here the construction of the matrices elements done outside the routine
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
*/ 
   MatricesElem(QQ,DegrePol,MaxQLevel,DegrePol+1+4)
/*
QQ[dd].Coeff[ii][jj] are matrix elements of operator Q^ii for ii in [0,DegrePol],
QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of operator D2Q (second order derivative)
QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order derivative)
QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q 
QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ
*/
//
//Derivative terms
     for (int mm=0;mm<4;mm++) 
     {
     allocate_MatDb(QQ[mm+DegrePol+1].Coeff, MaxQLevel+DegrePol, MaxQLevel+DegrePol)
     }
//
MatrixD2Q(QQ[DegrePol+1].Coeff,MaxQLevel+DegrePol)
MatrixD1Q(QQ[DegrePol+2].Coeff,MaxQLevel+DegrePol)
MatrixQD1Q(QQ[DegrePol+3].Coeff,MaxQLevel+DegrePol)
MatrixD1QQ(QQ[DegrePol+4].Coeff,MaxQLevel+DegrePol)
//
//Get The total time.
gettimeofday(&BeginAll,NULL);
CPUBeginAll = clock();
//
uint32_t FinalSize;//Final space size
//******************************************
//Correction of MaxFreq according to targets
//Protection in case of wrong data entries 
//******************************************    
     uint8_t *ZeroSt=new uint8_t[NMode];     
     InitTabInt(ZeroSt, NMode)     
      double FMax=0;     
      double MinMaxFreq=MaxVal<double>(KFC.KijNCpld[1],NMode)/2;
//     
      if ((MaxFreq< MinMaxFreq))
      {
       MaxFreq=MinMaxFreq;
       MinFreq=-100;
      }
      if ((MinFreq>0) && (GroundState <=0)) 
      {
       printf("\n***GroundState value must be specified for window target when MinFreq>0***\n\n");
       FinalMessaj()
       return 0;
      }      
//
     if(NTargetStates)
    {
     for (int nn=0;nn<NTargetStates;nn++)
      {
     if(FMax <= GetFreq0(TargetState[nn],NMode, KFC))
       FMax=GetFreq0(TargetState[nn],NMode, KFC);
      }
     if(FMax > MinMaxFreq)
      {MaxFreq=FMax;}
     else
      {
      MaxFreq=MinMaxFreq;
      MinFreq=-100;
      } //In case of ground state target
       //Should be larger for fundamentals frequencies.
       //To avoid overlap with the last eigenvalue of the solver.
     }
     //Kappa or MaxEv should be increased when close to infrared region
       if(Kappa<1.3 && MaxFreq > 3000)
       {Kappa=1.3;}      
// When the ground state is given <=0 it means that we actally want to compute it.  
     printf("\nMaximum targeted frequency %.1f, elongation factor Kappa %.1f \n", MaxFreq, Kappa);
//          
     int SumTarget=0;
//
     if(NTargetStates){SumTarget=NTargetStates+AddTarget;}//Add target to allocate proper memory
//
     SizeActMax=GetSizeActMax2(MaxEv, MAXNCV, NMode, MaxQLevel, DegrePol,Verbose,\
 NPES,NXDualHTruncPos-1,Memory,KNNZ,KNREZ,KNZREZ,SumTarget,TotalLFF, NXDualHTrunc);//Compute maximal size from parameters
//                
         printf("\nMaximal number of elements in B SizeActMax = %u \n",SizeActMax);
         if(SizeActMax<(uint32_t)10*NMode)
         {
          printf("Not enough allocated memory %u \n\n",SizeActMax);
          FinalMessaj()
          return 0;
         }         
     double TargetFreq=Kappa*MaxFreq; // Target frequency for initial subspace construction.
//
//     printf("\nMaximal energy for initial subspace construction\n %f cm-1 = kappa*(MaxFreq) \n", TargetFreq);        
//
       double TargetMax=TargetFreq;       
       if(TargetMax>Freq0Max)
       {printf("******Maximal energy too small*******\n\n"); FinalMessaj() return 0;}
//     
  if(NTargetStates)
     {
     if(Verbose)
      {
     MaxFreq=Freq0Max;
     printf("\nMaximum targeted frequency increased to %f\n\
to avoid null intersection with targets\n",MaxFreq);
      }
     }

FILE *FileRef;//To compare with reference frequencies in output 
FileRef=fopen(RefName,"r");
if(FileRef==NULL)
{
printf("\n FileRef not detected \n");
}
else{printf("\n FileRef detected \n");}
//
//Files to write Final basis set (FileBasis) and eigenvectors in binary form(FileVec)
//Take the previous name generated by DVCI
  char OutVec[MaxChar]={0};
//
  char OutBasis[MaxChar]={0};
//
  char BasisOut[30]="-FinalBasis.bin";//Final basis file name if PrintOut = 1
//
  char VectOut[30]="-Vectors.bin";//Binary file name with targeted eigenvectors if PrintOut = 1
//
  sprintf(OutVec, "%s", OutName);
  sprintf(OutBasis, "%s", OutName);
//
  strcat(OutVec, VectOut);
  strcat(OutBasis, BasisOut);
//
  FILE *FileVec=NULL;
  FILE *FileBasis=NULL;
//     
  FileVec=fopen(OutVec,"rb");
  FileBasis=fopen(OutBasis,"rb");   
//
  if (FileVec == NULL || FileBasis==NULL)
   {
   printf("Error opening files \n");
   FinalMessaj()
   return 0;
   }
//  

//Get The total time.
     gettimeofday(&BeginAll,NULL);
     CPUBeginAll = clock();
//
//
 int NScreen=0;//NUmber of screened states in output
//
uint8_t *ModeAct=NULL;
//
    FinalSize=GetConfsBin(FileBasis, ModeAct, NMode, &NScreen);
    if(!FinalSize){
                  printf("No basis set detected.\n");
                  printf("To use this programm first run DVCI to save final basis set.\n");
                  FinalMessaj()
                  return 0;
                  }
                  else
                  {
                  printf("Final basis set detected \n");
                  }   
//         
      printf("\nFinal basis set size %u , number of targets %i \n",FinalSize,NScreen);
//
if(GroundState > 0)
{
NScreen++;
printf("\nGroundState not into targets then NScreen incremented\n");
} 
     int *PositionTarget=new int[NScreen+1]; //For zero
//
      if(GroundState<=0)
        {
      PositionTarget[0]=0; //Always from initiale subspace construction
        }
  if(NTargetStates && Verbose)
   {printf("Position of targets in final space\n");}
//
 long detect=-1;
//
if(NTargetStates)
{
 for (int tt=0;tt<NTargetStates;tt++)
 {
 detect=LinearModeSearch(ModeAct, TargetState[tt], NMode*sizeof(uint8_t),  FinalSize);
 if(detect>=0)
      {
 PositionTarget[tt]=detect;
 // AfficheNu(ModeAct[detect].Degrees, NMode)
  if(Verbose)
       {
       AfficheNu(TargetState[tt], NMode);
       printf(": Pos %lu, ",detect);
          if(!((tt+1)%3))
          {printf("\n");}
       }
      }
 else{
     printf("Intersection nul between target");
     AfficheNu(TargetState[tt], NMode);
     printf(" and initial subspace\n\n");
     FinalMessaj()
     return 0;
      }
  }
printf("\n\n");
//Check zeros state
 detect=LinearU8Search(TargetState, ZeroSt, NMode*sizeof(uint8_t), NTargetStates);
 if(detect<0 && GroundState <= 0)
      {
     printf("\n**Ground state value must be provided when not included in targets **\n\n");
     FinalMessaj()
     return 0;
      }
 else if(detect>=0 && GroundState > 0)
      {
     printf("\n**Ground state cannot be among targets when GroundState > 0 **\n\n");
     FinalMessaj()
     return 0;
      }
//Check for doublons
for (int ii=0;ii<NTargetStates;ii++)
 {
for (int jj=0;jj<NTargetStates;jj++)
 {
 if(ii!=jj)
  {
   if(PositionTarget[ii]==PositionTarget[jj])
    {
     printf("Target %i and Target %i equal -> not possible all targets must be differents \n\n",ii,jj);
     FinalMessaj()
     return 0;  
    }  
   }
  }
 }
}
else
{
 for (int tt=0;tt<NScreen;tt++)
 {
 PositionTarget[tt]=tt;
 }
} 
      allocate_TabSize(Size,2);
      Size[1].DimAct=FinalSize;
//    Maximal NNZs  
      const uint64_t NNZActMax=(uint64_t)((double)(SizeActMax*NXDualHTrunc)*(KNNZ));
//    Allocate CSC pointers  
      AllocateCSC(IJAct,NNZActMax,FinalSize)//IJAct.NJ[jj],IJ.I[nn] in [0,FinalSize[ x [0,NNZActMax[
//    CSC values
      double *ValAct=NULL;
      ValAct=new double[NNZActMax];
//Initialize it to zero
      InitTab(ValAct,NNZActMax)  
//
/*Assemble the matrix of operator in final basis set*/
       uint32_t CheckOutputAssemble=AssembleHarmCSC (NMode,DegrePol, NCPol, 0, NPES,\
       ModeAct,Size,KFC,ValAct, IJAct, QQ, DegreCoupl, ZetaXYZ, Omega, ThrMat, NNZActMax,DualHPos, LFF,\
       PermutGen, NXDualHPlus, IncK2);
//     
        if(!CheckOutputAssemble)
        {
        printf("Faile assemble matrix operator\n");
        FinalMessaj()
        return 0;
        }
//
//Read Vectors from binary file and corresponding eigenvalue in position Size[1].DimAct*ii
      double *EigVec=NULL;
      EigVec=new double[(Size[1].DimAct+1)*NScreen];
//
      for (int ii=0;ii<NScreen; ii++)
       {
       fread (&EigVec[ii*(Size[1].DimAct+1)] , sizeof(double) , Size[1].DimAct+1 , FileVec );
       }
//
       double *Y=new double[Size[1].DimAct];
       char DipoleChar[MaxChar];       
//
       SubStringChar (PESName, DipoleChar, '.');
//  
//PrintBanner ("transition");
//EigVec[0] contains the ground state eigenvector if included into targets
//else take Harmonic ground state
//
//   
//
//ZPE
       double E0;
       if(GroundState <=0)
       {E0=EigVec[Size[1].DimAct];}
       else
       {
       E0=GroundState;
       printf("E0 calculated when not included into targets : %f\n",EigVec[Size[1].DimAct]); 
       int *PositionTTmp=new int[NScreen+1]; //For zero
       PositionTTmp[0]=0;
       for (int ii=0;ii<NScreen-1; ii++)
        {
       PositionTTmp[ii+1]=PositionTarget[ii];
        }
       for (int ii=0;ii<NScreen; ii++)
        {
       PositionTarget[ii]=PositionTTmp[ii];
        }
       delete [] PositionTTmp;
       }
//Compute Op*VecE0
       SparseSymCSC(EigVec,Y,IJAct, ValAct,Size[1]);
//
       for (int ii=0;ii<NScreen; ii++)
       {
       double Trans=ScualProd(&EigVec[ii*(Size[1].DimAct+1)],Y, Size[1].DimAct);
       printf("< Phi_0 | %s |  Phi_%d > : %f ,\n",DipoleChar,PositionTarget[ii],Trans);
       printf("Corresponding energy transition E_%d - E_0 : %f\n\n",PositionTarget[ii],\
       EigVec[(ii+1)*(Size[1].DimAct)+ii]-E0);
       }
//
     gettimeofday(&EndAll,NULL);
     time_spent = TIME_DIFFS(BeginAll,EndAll);
     CPUEndAll = clock();
     CPUTime=TIME_CPUS(CPUBeginAll, CPUEndAll); 

    if(Verbose)
    {printf("Total effective time spent: %u s and CPU wall time %ld \n \n",time_spent,CPUTime);}
//    
    printf("Total effective time spent: ");
    ConvertTimeS(time_spent);
//
    printf("Total CPU wall time: ");
    ConvertTimeS(CPUTime);
//
//
//Free
    fclose(FileBasis);
    fclose(FileVec);
//
    fclose(FileKey);
    fclose(FilePES);
//
//     FreeTabMode(ModeAct,FinalSize)
//  
//     FreeTabMode(DualHPos,NXDualHPlus) 
//
     delete [] ModeAct;
     delete [] DualHPos;
     FreeMat(ZetaXYZ, ((NMode+1)*NMode)/2)
//  
     FreeMat(KFC.KijNCpld,  DegrePol)
//
     FreeMat(KFC.Monm, NPES)
//
     FreeMat(CoorMode, NMode)
//
     delete [] KFC.KijCpld;  
//  
     FreeCSC(IJAct)
//""
     FreeLFF(LFF)
//
     delete [] Size;
     delete [] Omega;
     delete [] PermutGen;
     delete [] Mas;
     delete [] CoorCart;
     delete [] Corres;
     delete [] ZeroSt;
     delete [] EigVec;
     delete [] Pid;  
     delete [] ValAct;    
     delete [] Y; 
     delete [] PositionTarget; 
//
FreeMatricesElem(QQ,DegrePol+1+4,DegrePol+MaxQLevel)
// 
FinalMessaj()
//
 return 0;
//
}//End of main
