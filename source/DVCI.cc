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
#include <cstring>
#include "Shared.h"
#include <sys/time.h>
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
//
#define IncK2 0
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
    struct timeval begin, end, BeginAll, EndAll;
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
      printf("\n***Usage: DVCI Key.in***\n");
      printf("One argument requested : Key file missing.  \n");
      printf("PES file should also be in the same directory. \n\n");
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
/*
NXDualHTruncPos : Total number of raising excitations for the dual of H (after truncation through ThrKX)
NXDualHTrunc : Total number of excitations for the dual of H (after truncation through ThrKX)
SizeRezMax=SizeActMax*(NXDualHTruncPos-1)*KNREZ
NNZActMax=SizeActMax*(KNNZ)*(NXDualHTrunc)
NNZRezMax=SizeActMax*(KNZREZ)*(NXDualHTrunc) 
*/
int DegreCoupl[MaxCpld]={0}; //Maxdegree of PES for each coupling: starts at 0 for coupling 1
int DoGraph=1; // If not equal to zero, it will store the graph of residual matrix in CSC format
int EvalDeltaE=0;//If positif add correction vector and correction energy to targets
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
int NAddVar=NAdd;
int NAddIter=NAdd;
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
      printf("\n Threshold for PES ThrPES = %e,\n Threshold for Hamiltonian matrix elements ThrMat = %e \n", ThrPES, ThrMat);
      }
      else{printf("\n Threshold for PES and for Hamiltonian matrix elements have to be > 0 \n");FinalMessaj() return 0;}
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
       if(EvalDeltaE)
       {
       printf("\nEvalDeltaE non null correction of order 2 will be printed out at final stage\n");
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
FILE *FileRef;//To compare with reference frequencies in output 
FileRef=fopen(RefName,"r");
if(FileRef==NULL)
{
printf("\n FileRef not detected \n");
}
else{printf("\n FileRef detected \n");}
//
  char OutBasis[MaxChar]={0};
  char OutVec[MaxChar]={0};
//
  char OutMat[MaxChar]={0};
//
  char BasisOut[30]="-FinalBasis.bin";//Final basis file name if PrintOut = 1
//
  char VectOut[30]="-Vectors.bin";//Binary file name with targeted eigenvectors if PrintOut = 1
//
  char MatOut[30]="-Matrix.bin";//Binary file name with targeted eigenvectors if PrintOut > 1 
//
  sprintf(OutVec, "%s", OutName);
  sprintf(OutBasis, "%s", OutName);
  sprintf(OutMat, "%s", OutName);
//
  strcat(OutBasis, BasisOut);
  strcat(OutVec, VectOut);
  strcat(OutMat, MatOut);
//
//Files to write Final basis set (FileBasis) and eigenvectors in binary form(FileVec)
//      
      if(PESType)//Derivatives in a.u
      {LinesKijCpld=GetNumKijCpld(FilePES, NMode, ThrPES, &NCPol, &DegrePol, DegreCoupl);}
      else//Force constants in cm-1
      {LinesKijCpld=GetNumKijCpld0(FilePES, NMode, ThrPES, &NCPol, &DegrePol, DegreCoupl);}
//
      if(!LinesKijCpld || !DegrePol)
      {
      printf("\n*****Problem reading force field*****\n");
      printf("*****Check if key words written at beginning of lines*****\n");
      FinalMessaj() 
      return 0;
      }
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
printf("\n***Quantum levels by coordinate***\n");
printf("[");
for (int mm=0; mm<  NMode-1; mm ++) 
  { 
       printf("%u, ",Pid[mm]-1);  
  }
printf("%u]\n",Pid[NMode-1]-1);
//
printf("\n***Pruning condition***\n");
printf("Sum ");
double FreqMin=MinimumVal(KFC.KijNCpld[1],NMode);
for (int mm=0; mm<  NMode-1; mm ++) 
  { 
       printf("%ib_{%i} + ",(int)round(KFC.KijNCpld[1][mm]/FreqMin),mm+1);  
  }
printf("%ib_{%i} <= %i\n",(int)round(KFC.KijNCpld[1][NMode-1]/FreqMin),NMode,(int)round(Freq0Max/FreqMin));
//
//
double Ull=0;//Watson term
if(DoRot)
{
//Compute center of Mass
MassCenter(Mas, CoorCart, NAtome, MassC);
//Compute principal axes of moment of inertia: Axes and values
Ull=MomentInertie(Mas, CoorCart, NAtome, MassC, PrincAx, PrincEig);//Watson Correction
printf("\n*** Watson term = -1/8*(Sum mu_ll) %.4f ***\n\n",Ull); 
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
//""DualHPos[Idm(xx)+mm], (xx,mm) in [0,NXDualHPlus[ x [0,NMode[.
//    
 LocalFF LFF;
 LFF.Num=NULL;
 LFF.Num=new uint32_t[NXDualHPlus+1];
//
/*LFF[ii(xx)] : local force field associated with positive excitation &DualHPos[Idm(xx)].
LFF.Idx[ii(xx)] are defined for ii [LFF.Num[xx],LFF.Num[xx+1][
LFF.Idx[ii(xx)] < 0 are indexes of non coupled force constants:
KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
0 <= LFF.Idx[ii(xx)] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF.Idx[ii(xx)]].
LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
with a unique corresponding (ni,nj,nk,nl)  
LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
with a unique corresponding (ni,nj,nk,nl) 
*/
//
   NXDualHPlus=0;
   uint32_t TotalLFF=0;//Total number of force constants per excitation in dual operator
//
//Guaranty zero excitation
//
   if(!NLFFPerEx[0])
   {
   printf("*******No zero excitation**********\n");
//""
   memcpy(&DualHPos[Idm(NXDualHPlus)],&XPolSupPos[0],(size_t)NMode);
   LFF.Num[NXDualHPlus]=0;
   NXDualHPlus++;
   TotalLFF++;
   }   
//   
int MaxTerm=0;
 for (uint32_t ll=0; ll<NGenPoz; ll++)
 { 
 if(NLFFPerEx[ll])
  {
//""
   memcpy(&DualHPos[Idm(NXDualHPlus)],&XPolSupPos[Idm(ll)],(size_t)NMode);
//""
//LFF.Idx[ii(xx)] are defined for ii [LFF.Num[xx],LFF.Num[xx+1][ then
   LFF.Num[NXDualHPlus]=TotalLFF;
   NXDualHPlus++;
   TotalLFF+=NLFFPerEx[ll];
/*
   CopyTab(DualHPos[NXDualHPlus].Degrees,XPolSupPos[ll].Degrees,NMode) 
   InitLFF(LFF[NXDualHPlus],NLFFPerEx[ll]);
   NXDualHPlus++;
   TotalLFF+=NLFFPerEx[ll];
*/
   if(NLFFPerEx[ll]>MaxTerm && ll)
    {
    MaxTerm=NLFFPerEx[ll];
    }
   }
  }
//""/Add last term
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
uint32_t SizeInit;//Intial subspace size
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
//Initial subspace construction    
     allocate_TabMode(ModeAct,SizeActMax)
//ModeAct : Multi-dimensional array for active space B : ModeAct[Idm(nn)+mm], (nn,ii) in [0,SizeActMax[ x [0,NMode[.
//Changed to ModeAct[Idm(nn)+mm]
     double TargetFreq=Kappa*MaxFreq; // Target frequency for initial subspace construction.
//
     printf("\nMaximal energy for initial subspace construction\n %f cm-1 = kappa*(MaxFreq) \n", TargetFreq);        
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
//
printf("\n******Starting program******\n\n");
// 
     gettimeofday(&begin,NULL);
//
     SizeInit=InitB0(ModeAct, KFC, Pid, NMode, SizeActMax, TargetMax);
//             
     if(!SizeInit){printf("\n\n"); FinalMessaj() return 0;}
//     
     gettimeofday(&end,NULL);
     time_spent=TIME_DIFFS(begin, end); 
//Define the maximum size of active space B and residual space
     printf("Initial space size %u, time spent to build it %u s\n", SizeInit,time_spent);
// 
     int *PositionTarget=new int[NTargetStates+1]; //For zero
//
      if(GroundState<=0)
        {
      PositionTarget[0]=0; //Always from initiale subspace construction
        }
//
  if(NTargetStates && Verbose)
   {printf("Position of targets in initial space\n");}
//
 long detect=-1;
//
if(NTargetStates)
{
 for (int tt=0;tt<NTargetStates;tt++)
 {
 detect=LinearModeSearch(ModeAct, TargetState[tt], NMode*sizeof(uint8_t),  SizeInit);
 if(detect>=0)
      {
 PositionTarget[tt]=detect;
 //'' AfficheNu(&ModeAct[Idm(detect)], NMode)
  if(Verbose)
       {
       AfficheNu(TargetState[tt], NMode);
       printf(": Pos %lu, ",detect);
          if(!((tt+1)%3))
          {printf("\n");}
       }
      }
 else{
     printf("Intersection nul between target ");
     AfficheNu(TargetState[tt], NMode);
     printf(" and initial subspace\n\n");
     FinalMessaj()
     return 0;
      }
  }
printf("\n");
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
     printf("\n**Ground state cannot be among targets when GroundState > 0**\n\n");
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
//
    const uint64_t SizeRezMax=(uint64_t)((double)(SizeActMax*(NXDualHTruncPos-1))*KNREZ);
//Allocation of the dimensions of the different arrays and configurations
     allocate_TabSize(Size,IterMax);    
     Size[1].DimAct=SizeInit; 
//Definition of the maximal nnz elements of the active and residual matrices.    
     const uint64_t NNZActMax=(uint64_t)((double)(SizeActMax*NXDualHTrunc)*(KNNZ));
     const uint64_t NNZRezMax=(uint64_t)((double)(SizeActMax*NXDualHTrunc)*(KNZREZ));      
    printf("\nMaximal residual space size\n SizeRezMax (=SizeActMax*(NXDualHTruncPos-1)*KNREZ) : %lu \n",SizeRezMax);
//
    if(SizeRezMax>UINT32_MAX && DoGraph)
    {
    printf("\nMaximal residual space size too big for construction of the graph of\n");
    printf(" residual matrix, set DoGraph to 0 or decrease KNRez \n\n");
    FinalMessaj()
    return 0;
    }
//
    printf("\nMaximal NNZ for matrix Hb,\n NNZActMax (=SizeActMax*(KNNZ)*(NXDualHTrunc)) : %lu \n",NNZActMax);
    printf("Maximal indexes in graph of residual matrix Hsb,\n NNZRezMax (= SizeActMax*(KNZREZ)*(NXDualHTrunc)) : %lu \n\n",NNZRezMax);
//
//ModeRez : Multi-dimensional array for residual space : ModeRez[Idm(nn)+mm], (nn,mm) in [0,SizeRezMax[ x [0,NMode[. 
//Changed to ModeRez[Idm(nn)+mm]
    allocate_TabMode(ModeRez,SizeRezMax)   
//
// Allocation and intialisation of the array and configurations    
//Compressed sparse  column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
AllocateCSC(IJAct,NNZActMax,SizeActMax)//IJAct.NJ[jj],IJ.I[nn] in [0,SizeActMax[ x [0,NNZActMax[
AllocateCSC(IJRez,NNZRezMax,SizeActMax) //IJRez.NJ[jj],IJRez.I[nn] in [0,SizeRezMax[ x [0,NNZRezMax[
//
//Assign maximal sizearray (structure for dimensions)
SizeArray SizeMax;
SizeMax.DimRez=SizeRezMax;
SizeMax.DimAct=SizeActMax;
SizeMax.NNZAct=NNZActMax;
SizeMax.NNZRez=NNZRezMax;
//Permutation index for residual configurations
     uint32_t *PermutRez=NULL;
     PermutRez=new uint32_t [SizeRezMax]; 
//     
       uint32_t Alloc=1;
//       
       if(!DoGraph){Alloc=0;}
//If TabNull[Lin]=0 then element Lin is no longer part of residual space.
//Used only to update the graph of residual matrix
       uint32_t *TabNull=NULL;
       TabNull=new uint32_t [SizeRezMax*Alloc];
//
     for (uint64_t nn = 0; nn < SizeRezMax; nn++)
      {
      //PermutRez[nn]=nn:before but now avoid size overhead 
       PermutRez[nn]=nn-(nn/MaxQSortSeq)*MaxQSortSeq;    
      if(DoGraph)
       {TabNull[nn]=1;}
      }
//Permutation index for active configurations in B
     uint32_t *PermutAct=NULL;
     PermutAct=new uint32_t [SizeActMax]; 
//
     for (uint32_t nn = 0; nn < SizeActMax; nn++)
      {
      PermutAct[nn]=nn;
      }
//
 double Shift=0; //Shift for the matrix to compute the eigenvalues
//Assignation of the shift outside the loop
     for(int mm=0;mm<NMode;mm++)
     {
     Shift=Shift+(KFC.KijNCpld[1][mm])*(Pid[mm]-0.5);
     }
//Allocation of the non zero values of the active matrix Hb
 double *ValAct=NULL;
 ValAct=new double[NNZActMax];
//Initialize it to zero
 InitTab(ValAct,NNZActMax)  
//***********
//Eigensolver
//************
//Correction of MaxEv to maximal allowed eigenvalue position.
MaxEv=Min<int>(MaxEv,SizeInit-1);
int NEV=MaxEv;
//Allocate the eigenvalues
double *EigVal=NULL;
EigVal=new double[NEV];
//Allocate the relative residues
float *RezRel=NULL;
RezRel=new float[NEV];
double Ground=GroundState;//Ground state may be variable if we actually compute it
 int Iteration=0;
 uint32_t CheckOutputPolyX=0; //Check good execution of PolyX
 uint32_t CheckOutputAssemble=0; //Check good execution of assemblage 
// int NotConverj=1; //Turn to zero if all the eigenpairs have converged.
 int nconv=0; //Number of converged eigenvalues 
//Max number of ritz vectors
if(!MAXNCV)
{MAXNCV=(NEV)*2;} //Automatic choice of number of Lanczos basis vectors : 2*NEV.
//
int NCV=MAXNCV; //Maximal number of Lanczos basis vectors
//External allocation of arrays for eigensolver
      double *workd=NULL;
      workd=new double [3*SizeMax.DimAct+1];
      InitTab(workd,3*SizeMax.DimAct+1) 
//
      double *workl=NULL;
      workl = new double[(MAXNCV)*(MAXNCV+8)+1]; 
       InitTabInt(workl,(MAXNCV)*(MAXNCV+8)+1)  
      double *EigVec=NULL;
      EigVec=new double[SizeMax.DimAct*NCV+1];   
      InitTab(EigVec,SizeMax.DimAct*NCV+1) 
//
//Component jj of Eigenvector ii is accessible at EigVecii[ii*Size[Iteration+1].DimAct + jj]
    int NScreen=NEV; //Maximal number of screened eigenvalues 
    if(NTargetStates)
      {
       NScreen=Min<int>(NEV,NTargetStates+AddTarget);
       for (int tt=0;tt<NTargetStates;tt++)
        {
       EigVec[tt*SizeInit+PositionTarget[tt]]=1;
        }
      }   
//****************
//Assignment arrays
//****************  
    int NScreenTot=NScreen;//Maximal number of screened eigenvalues    
    int *TabScreen;
    TabScreen=new int[NScreen];//Indexes of screened eigenvalues 
//Initialise TabScreen
       for (int tt=0;tt<NScreen;tt++)
        {
        TabScreen[tt]=tt;
        }
int MaxScreen=0;//Higher eigenvalue position of the targets
//
uint32_t *AssignRez= new uint32_t[MaxAdd];//Store the indexes of residual vectors selected components 
//
double DeltaFreq=0;//Maximal difference between previous and current eigenvalues
//
double MaxRelRez=100;//Maximal relative residu
//
double *EigValOld=NULL;//Old eigenvalues to compare to new ones
EigValOld=new double[NEV];
  for (int nn = 0 ; nn < NEV; nn++)
   { EigValOld[nn]=0;
     for (int mm = 0 ; mm < NMode; mm++)
    {
//""EigValOld[nn]=EigValOld[nn]+(ModeAct[Idm(nn)+mm]+0.5)*KFC.KijNCpld[1][mm];
    EigValOld[nn]=EigValOld[nn]+(ModeAct[Idm(nn)+mm]+0.5)*KFC.KijNCpld[1][mm];
    } 
   }
   SortAsc2<double>(EigValOld, NEV);
//   
//Remind parameters before beginning of the loop
  if(Verbose)
   {
   printf("Tol %e and Shift for the eigensolver %.1f \n \n", Tol,Shift);
   printf("MaxEv = %i wanted eigenvalues, MAXNCV = %i Lanczos basis vectors,\n\
 Maximal variation of eigenvalues DeltNEV = %i, \nFrequency window : [%f , %f] .\n", MaxEv, MAXNCV, DeltaNev, MinFreq, MaxFreq);
   printf("SizeActMax %u,SizeRezMax %lu,NNZActMax %lu, NNZRezMax %lu \n",\
   SizeActMax,SizeRezMax,NNZActMax, NNZRezMax);
   }
//***************
//Residual vectors
//*************** 
    float *RezVect=NULL;   
    RezVect=new float[SizeMax.DimRez*NScreen];
    InitTab(RezVect,SizeMax.DimRez*NScreen)
//Component jj of RezVect ii is accessible at RezVect[ii*SizeMax.DimRez + jj]       
    size_t SizeBit=NMode*sizeof(uint8_t);//sizeof(uint8_t)=1!
    int ExitSuccess=1;
    int MaxAddVar=0;//MaxAddVar=Min<int>(MaxAdd,SizeActMax-Size[Iteration+1].DimAct);
    long DetectZero=0;//If 1 the ground state has been detected
    long PosZero=-1;//If > 0 give the position of the groundstate that might be non equal to zero  
//
 while(MaxRelRez > EpsRez) //Main loop
 { 
  printf("\n***** Iteration number : %i *****\n",Iteration);

   if(Verbose)
     {
    //Memory usage
    printf("\n---Effective memory usage %.1f Mo---\n", EffectiveMemory2(argv[0]));
     } 
       gettimeofday(&begin,NULL);
/*Assemble the matrix from which will be computed the eigenvalues*/
       CheckOutputAssemble=AssembleHarmCSC (NMode,DegrePol, NCPol, Iteration, NPES,\
       ModeAct,Size,KFC,ValAct, IJAct, QQ, DegreCoupl, ZetaXYZ, Omega, ThrMat, NNZActMax,DualHPos, LFF,\
       PermutGen, NXDualHPlus, IncK2);
//
        if(!CheckOutputAssemble)
        {
        printf("Current size |B|: %u \n",Size[Iteration].DimAct);
        MaxRelRez=0;Iteration--;
        ExitSuccess=0;
        break;
        }
//    
   gettimeofday(&end,NULL);
   time_spent=TIME_DIFFS(begin, end);
//
      if(Verbose)
           {
      printf("Time spent to assemble the active matrix Hb: %u  s\n",time_spent);
           }
//
         if(Verbose)
         {
         printf("Number of wanted eigenvalues : %i  \n",NEV);
         }
//
     gettimeofday(&begin,NULL);
//Set maximal number of Lanczos basis vectors for diagonalization in eigensolver     
   NCV=(int)Min<uint32_t>(Size[Iteration+1].DimAct,(uint32_t)MAXNCV);    
//    
//Compute the smallest eigenvalue of the reduced Hamiltonian matrix    
   nconv=SolverCSCSym(Size, Iteration, NEV, NCV, Shift, Tol, EigVal,\
   EigVec, workd, workl,IJAct, ValAct);
//
        if(!nconv)
        {
        printf("\n**Non convergence on eigensolver**\n");
        MaxRelRez=0;Iteration--;
        ExitSuccess=0;
        break;
        }
//
          DeltaFreq=0;   
// 
//Get the varation of eigenvalues between current and last iteration DeltaFreq=Max_l|E_l^(j)-E_l^(j-1)|
          for (int cc=0;cc<nconv;cc++)
          {
             if(SIGN<double>(EigVal[cc]-EigValOld[cc])*(EigVal[cc]-EigValOld[cc]) >DeltaFreq)
             {DeltaFreq=SIGN<double>(EigVal[cc]-EigValOld[cc])*(EigVal[cc]-EigValOld[cc]);}
             EigValOld[cc]=EigVal[cc];
          }
//    
   gettimeofday(&end,NULL);
   time_spent=TIME_DIFFS(begin, end);
//
    if(Verbose)
         {
   printf("Value of DeltaFreq=Max_l|E_l^(j)-E_l^(j-1)|: %.1f \n",DeltaFreq);
   printf("Time spent to solve eigenvalue problem: %u s\n",time_spent);
   printf("First eigenvalues\n");
   afficheTabDb(EigVal, Min<int>(nconv,10))
         }
//Next iteration, the active space B has got minimum Size[Iteration+1].DimAct elements 
//Size[Iteration+2].DimAct will be incremented for next iteration when adding new basis functions.
Size[Iteration+2].DimAct=Size[Iteration+1].DimAct;
//
    if(Verbose)
    {
    printf("\nNumber | Eigenvalue | Frequency   | Assignment ||\n");
    printf("       |            | En-E0 or E0 | (Component)  ||\n");
    printf("------------------------------------ \n");
    }
    NScreen=0;
    MaxScreen=0;
    DetectZero=0; //If >0 the groundstate has been detected
    detect=1; //If it remains positive one of the Target has been detected 
    PosZero=-1; //Position of groundstate if detected (might be not equal to zero) 
    for (int PosEig=0;PosEig<nconv;PosEig++)
     {
//Sparse matrix-vector product HR*X 
//X eigenvector PosEig stored in &EigVec[PosEig*Size[Iteration+1].DimAct], HR residual matrix
        if(GroundState<=0)
        {    
        DetectZero=DetectMaxCoordinate(&EigVec[PosEig*Size[Iteration+1].DimAct],PositionTarget, 1,ThrCoor);
        if(DetectZero)
         {
         PosZero=PosEig;         
         Ground=EigVal[PosEig];
         } 
        }
        if(NTargetStates && !DetectZero)
        {detect=DetectMaxCoordinate(&EigVec[PosEig*Size[Iteration+1].DimAct], PositionTarget, NTargetStates,ThrCoor);}
//Freq[PosEig]=EigVal[PosEig]-Ground should be in the interval of interest [MinFreq,MaxFreq+DeltaFreq]                 
           if( ((EigVal[PosEig]-Ground <= MaxFreq+DeltaFreq) && (EigVal[PosEig]-Ground >= MinFreq)) && (detect || DetectZero) )
              { 
               if(PosEig>=MaxScreen)
               {
               MaxScreen=PosEig;
               } 
               if(NTargetStates && (NScreen>=AddTarget+NTargetStates))
               {
               NScreen++;
               break;                   
               } 
               TabScreen[NScreen]=PosEig;
               NScreen++;    
              if(Verbose)
              {
//Print eigenvalue corresponding frequency, residues and ground state for PosEig=0
              if(PosEig==PosZero)
               {
              printf("   %d   |   %.2f   |  %.2f  | ",\
              PosEig,EigVal[PosEig],EigVal[PosEig]);
               }
              else
               {
              printf("   %d   |  %.2f  |  %.2f  | ",\
              PosEig,EigVal[PosEig],EigVal[PosEig]-Ground);
               }
//Assignment display        
            AfficheNuSup(ModeAct, &EigVec[PosEig*Size[Iteration+1].DimAct], NMode, Size[Iteration+1].DimAct, ThrCoor);
//           
             printf(" || \n \n");
             }              
          }
     }
//Update MaxEv to gain computational cost for eigensolver
NEV=Min<int>(MaxEv,MaxScreen+DeltaNev);   
//
      if(!NScreen){
                   printf("\nNot enough eigenvalues to reach targets -> increase MaxEV or DeltaNev,\n");
                   printf("or maximal components are too small -> decrease ThrCoor,\n");
                   printf("or initial subspace too small -> increase Kappa,\n");
                   printf("or could also be the potential badly approximated -> decrease MaxQLevel.\n\n");     
                   FinalMessaj()
                   ExitSuccess=0;              
                   break;                   
                   }
     if(NTargetStates && (NScreen>AddTarget+NTargetStates))
                   {
                   printf("\nToo many targets -> increase AddTarget or ThrCoor\n\n");
                   FinalMessaj()
                   ExitSuccess=0;
                   break;                   
                   }
     else if(NTargetStates && (NScreen<NTargetStates))
                   {
                   printf("\nNot enough targets screened -> decrease ThrCoor\n or increase MaxEv, DeltaNev and/or Kappa\n\n");
                   FinalMessaj()
                   ExitSuccess=0;
                   break;                   
                   }
    gettimeofday(&begin,NULL);
    if(DoGraph)
    {//The space H*(A) is fetched at each iteration where A is the set of added basis functions.
    CheckOutputPolyX=PolyX(ModeRez,ModeAct,DualHPos, QQ, PermutRez, PermutAct, Corres, KFC, Size,\
    Pid, Freq0Max, TargetMax, NMode, NXDualHTruncPos, NCPol, DegrePol, NPES,Iteration, NScreen,  TabScreen,\
    SizeMax, RezVect,LFF,EigVec,ZetaXYZ,Omega,IJRez);
    }
    else
    {//The whole space H*(B) is fetched at each iteration.
    CheckOutputPolyX=PolyXNoG(ModeRez,ModeAct,DualHPos, QQ, PermutRez, PermutAct, Corres, KFC, Size,\
    Pid, Freq0Max, TargetMax, NMode, NXDualHTruncPos, NCPol, DegrePol, NPES,Iteration, NScreen, TabScreen,\
    SizeMax, RezVect,LFF,EigVec,ZetaXYZ,Omega);
    }
    if(CheckOutputPolyX==0)
    {
     printf("\n|B|: %u , |Bs|: %lu, \n NNZ(Hb) : %lu , NNZ(Hsb) : %lu \n",\
     Size[Iteration+1].DimAct, Size[Iteration+1].DimRez, Size[Iteration+1].NNZAct, Size[Iteration+1].NNZRez);
     MaxRelRez=0;
     ExitSuccess=0;
     break;
    }
   gettimeofday(&end,NULL);
   time_spent=TIME_DIFFS(begin, end);
//
     printf("\n|B|: %u , |Bs|: %lu, \n NNZ(Hb) : %lu , NNZ(Hsb) : %lu \n",\
     Size[Iteration+1].DimAct,Size[Iteration+1].DimRez, Size[Iteration+1].NNZAct, Size[Iteration+1].NNZRez);
//     
     if(Verbose)
     {
     printf("Time spent bulding residual space: %u \n",time_spent);
     }
//
    if(Iteration && DoGraph)
      {
     gettimeofday(&begin,NULL);
//
//Matrix vector product of residual matrix on the fly completed thanks to graph stored in IJRez
     FlyRezMVP2(ModeRez, ModeAct,DualHPos, QQ,PermutGen,TabNull, KFC,Size,NMode,\
     IJRez,NXDualHPlus, DegrePol, NPES, Iteration, NScreen, TabScreen,SizeMax,\
     RezVect, LFF,EigVec, ZetaXYZ, Omega);
//
     gettimeofday(&end,NULL);
//
      time_spent=TIME_DIFFS(begin, end);
//
      if(Verbose)
       {
      printf("Time spent for residual MVP : %u \n",time_spent);
       }
      }    
// 
    MaxRelRez=0;
//
    printf("\nNumber | Frequency   | Relative | Assignement || \n");
    printf("       | En-E0 or E0 |  Residue | (Component) || \n");
    printf("------------------------------------ \n");

     for (int ll=0;ll<NScreen;ll++)
     { 
              RezRel[ll]=Norm2F(&RezVect[SizeRezMax*ll],Size[Iteration+1].DimRez)/EigVal[TabScreen[ll]];           
//
              if(TabScreen[ll]==PosZero)//Last position of zero readed
               {
              printf("   %d   |  %.4f  |  %.4f  | ",\
              TabScreen[ll],EigVal[TabScreen[ll]],RezRel[ll]);           
               }
              else
               {
              printf("   %d   |  %.4f  |  %.4f  | ",\
              TabScreen[ll],EigVal[TabScreen[ll]]-Ground,RezRel[ll]);
               }
//Assignment display        
             AfficheNuSup(ModeAct, &EigVec[TabScreen[ll]*Size[Iteration+1].DimAct], NMode, Size[Iteration+1].DimAct, ThrCoor);
             printf(" || \n \n");
//
        if(MaxRelRez<RezRel[ll])
        {
        MaxRelRez=RezRel[ll];
        }  
     }
//
     printf("Maximal relative residue %f \n", MaxRelRez); 
//
      if (MaxRelRez>EpsRez)
      {
        float MaxComp=0;//Sum of maximal component of maximal component of residual vectors
//
        int NNotConv=0;//number of non converged eigenvector at this iteration
//
        for (int ll=0;ll<NScreen;ll++)
         {   
          if(RezRel[ll]>EpsRez) 
           {
           MaxComp+=FindMaximumValPoz<float>(&RezVect[SizeRezMax*ll],Size[Iteration+1].DimRez);
           TabScreen[NNotConv]=ll;
           NNotConv++;
           }
         }
//The entries of the RezVec greater than to ||R||/(EtaComp*NNotConv) will be considered for enrichment
          MaxComp/=(EtaComp*NNotConv);
//
          if(Verbose)
          {printf("Maximal component of error vector divided by EtaComp*NNotConv %f \n", MaxComp);}
//
         NAddIter=NAdd*(Iteration+1);
//                  
         if(Size[Iteration+1].DimAct<SizeActMax)
          {
          MaxAddVar=Min<int>(MaxAdd,SizeActMax-Size[Iteration+1].DimAct);
          }
          else
          {
         MaxAddVar=0;
          }        
          if(MaxAddVar)
          {
   gettimeofday(&begin,NULL);
//                   
            NAddVar=MaxContribAll(RezVect,AssignRez,TabScreen,TabNull, DoGraph,\
            Size[Iteration+1].DimRez,NAddIter, MaxAddVar, MaxComp, NNotConv, NScreenTot, SizeRezMax);  
//    
   gettimeofday(&end,NULL);
   time_spent=TIME_DIFFS(begin, end);
//
           printf("Number of added basis functions |A| : %i, Time spent to get them : %u s\n", NAddVar,time_spent); 
          }
         else
          {
     printf("\n|B|: %u , |Bs|: %lu, \n NNZ(Hb) : %lu , NNZ(Hsb) : %lu \n",\
     Size[Iteration+1].DimAct, Size[Iteration+1].DimRez, Size[Iteration+1].NNZAct, Size[Iteration+1].NNZRez);
     printf("Maximum number of basis functions reached at previous iteration\n\
-> Algorithm stop at iteration %i \n",Iteration);
     MaxRelRez=0;
     break;
          }
         for (int aa=0;aa < NAddVar ;aa++)
          {//Add elements for next iteration and increment Size[Iteration+2].DimAct
           if (Size[Iteration+2].DimAct < SizeActMax)
             {
//""
          memcpy(&ModeAct[Idm(Size[Iteration+2].DimAct)], &ModeRez[Idm(AssignRez[aa])],SizeBit);         
          Size[Iteration+2].DimAct++;
             }
           }
        }    
 else
     {
     printf("\n|B|: %u , |Bs|: %lu, \n NNZ(Hb) : %lu , NNZ(Hsb) : %lu \n",\
     Size[Iteration+1].DimAct, Size[Iteration+1].DimRez, Size[Iteration+1].NNZAct, Size[Iteration+1].NNZRez);
     printf("No added elements at iteration %i\n-> Stopping criterion verified MaxRelRez = %.4f \n",Iteration,MaxRelRez);     
     MaxRelRez=0;
     break;
     }
//
     Iteration++;     
// 
             if (Iteration >= IterMax-1)
             {
             printf("Maximal number of iteration reached at iteration %i -> algorithm stop\n \n",Iteration);
             MaxRelRez=0;
             break;
             }
 }//End of iterations 
//
     gettimeofday(&EndAll,NULL);
     time_spent = TIME_DIFFS(BeginAll,EndAll);
     CPUEndAll = clock();
     CPUTime=TIME_CPUS(CPUBeginAll, CPUEndAll);    
//
    if(Verbose)
    {printf("Total effective time spent: %u s and CPU wall time %ld \n \n",time_spent,CPUTime);}
//    
    printf("Total effective time spent: ");
    ConvertTimeS(time_spent);
//
    printf("Total CPU wall time: ");
    ConvertTimeS(CPUTime);
//
//Delta E correction energies
double *VPTE=NULL;
//
if(ExitSuccess)
 {
 
 if(EvalDeltaE)
  {
     VPTE=new double[NScreen];
//
     printf("\n****Starting computation of correction energies****\n");
     gettimeofday(&begin,NULL);
//
     InitTabInt(VPTE,NScreen)
     CorrectionEnergy(RezVect,Size[Iteration+1].DimRez, EigVal, ModeRez, NScreen, NPES,\
     DegrePol,NMode, TabScreen, QQ, SizeMax, LFF, ZetaXYZ, Omega, VPTE,KFC);
//
     gettimeofday(&end,NULL);
//    
     time_spent=TIME_DIFFS(begin, end); 
//Define the maximum size of active space B and residual space
     printf("\n**Total time for energy correction evaluation %u s**\n",time_spent);
  }   

      double TabRef[MaxRef]={0}; //Array to store the reference values
          if(FileRef!=NULL)
          {
                printf("\n %d lines detected in Reference file\n",GetValRef(FileRef, TabRef));
          }
          if(EvalDeltaE)
          { 
                printf("\n Assignment  | Frequency | Relativ |  Correction | Error  || \n");
                printf(" (Component) | (number)  | Residue |   Energy   | Ref-Here || \n");
                printf("------------------------------------ \n");              
          }
          else
          { 
                printf("\n Assignment  | Frequency | Relativ | Error  || \n");
                printf(" (Component) | (number)  | Residue | Ref-Here || \n");
                printf("------------------------------------ \n");              
          }
     int ss=0;    
     int SaveZero=0;//Last position of zero readed      
     for (int PosEig=0;PosEig<nconv;PosEig++)
     {   
    DetectZero=0; //If >0 the groundstate has been detected
    detect=1; //If it remains positive one of the Target has been detected 
    PosZero=-1; //Position of groundstate if detected (might be not equal to zero)             
        if(GroundState<=0)
        {
        DetectZero=DetectMaxCoordinate(&EigVec[PosEig*Size[Iteration+1].DimAct], PositionTarget, 1,ThrCoor);
        if(DetectZero)
         {
         PosZero=PosEig; 
         SaveZero=PosEig;      
         Ground=EigVal[PosEig];
/*
         if(EvalDeltaE)
          {
            Ground=EigVal[PosEig]+VPTE[ss];
          }*/
         } 
        }
        if(NTargetStates && !DetectZero)
        {detect=DetectMaxCoordinate(&EigVec[PosEig*Size[Iteration+1].DimAct], PositionTarget, NTargetStates,ThrCoor);}    
          if( ((EigVal[PosEig]-Ground <= MaxFreq+DeltaFreq) && (EigVal[PosEig]-Ground >= MinFreq)) && (detect || DetectZero) )
          { 
           //Assignment display
           AfficheNuSupTex(ModeAct, &EigVec[PosEig*Size[Iteration+1].DimAct], NMode, Size[Iteration+1].DimAct, ThrCoor);
           //""                  
           if(PosEig==PosZero)
           {
           if(EvalDeltaE)
            {
           if(!DoRot)
             {
            printf(" | %.2f(%d) | %.4f | %.4f  | %.4f  |  // \n \n", \
            Ground,PosEig,RezRel[ss],VPTE[ss],GetClosest(TabRef,Ground)-(Ground));
             }
           else//Add Watson
             {
           printf(" | %.2f(%d) (mu : %.2f) | %.4f | %.4f | %.4f  // \n \n", \
           Ground,PosEig,Ull,RezRel[ss],VPTE[ss],GetClosest(TabRef,Ground)-(Ground));//Add watson
             }
            }
          else
           {
           if(!DoRot)
             {
            printf(" | %.2f(%d) | %.4f | %.4f // \n \n", \
            Ground,PosEig,RezRel[ss],GetClosest(TabRef,Ground)-(Ground));
             }
           else//Add Watson
             {
           printf(" | %.2f(%d) (mu : %.2f) | %.4f | %.4f  // \n \n", \
           Ground,PosEig,Ull,RezRel[ss],GetClosest(TabRef,Ground)-(Ground));//Add watson
             }
            }
//
           }
           else
           {
           if(EvalDeltaE)
            {
//           EigVal[PosEig]+=VPTE[ss];
           printf(" | %.2f(%d) | %.4f | %.4f | %.4f // \n \n", \
           EigVal[PosEig]-Ground,PosEig,RezRel[ss],VPTE[ss],GetClosest(TabRef,EigVal[PosEig]-Ground)-(EigVal[PosEig]-Ground));         
            }
           else
            {
           printf(" | %.2f(%d) | %.4f | %.4f // \n \n", \
           EigVal[PosEig]-Ground,PosEig,RezRel[ss],GetClosest(TabRef,EigVal[PosEig]-Ground)-(EigVal[PosEig]-Ground));
            }
           }         
           ss++;                                                                                                          
           }
          }
//
if(PrintOut==1)//Print final basis set in FileBasis
 {
PrintConfsBin(OutBasis,ModeAct, NMode, NScreen, Size[Iteration+1].DimAct);
//""
PrintVecBin(OutVec, EigVec, EigVal, NScreen,TabScreen, SaveZero, GroundState, Size[Iteration+1].DimAct);
 }
if(PrintOut>1)//Matrix in binary
 {
 PrintMatCSCBin(OutMat, Size, Iteration, NEV, NCV, SizeInit, Shift, Tol, Ull, IJAct, ValAct);
 PrintConfsBin(OutBasis,ModeAct, NMode, NScreen, Size[Iteration+1].DimAct);
 }
}//If exitsuccess
else if(!ExitSuccess && PrintOut)//""Indicate non success with null value into binarry files
{
FILE *FileMat=fopen(OutMat,"wb");
FILE *FileBasis=fopen(OutBasis,"wb");
//
fwrite (&ExitSuccess, sizeof(int) , 1 , FileMat); 
fwrite (&ExitSuccess, sizeof(int) , 1 , FileBasis);
// 
fclose(FileMat);
fclose(FileBasis);
}
//
    fclose(FileKey);
    fclose(FilePES);
//
    if(FileRef!=NULL){fclose(FileRef);}
//
    if(EvalDeltaE){delete [] VPTE;}
//
//""FreeTabMode(ModeAct,SizeActMax)

  delete [] ModeAct;
  delete [] ModeRez;
  delete [] DualHPos;

//  
//""FreeTabMode(DualHPos,NXDualHPlus) 
//
//""FreeTabMode(ModeRez,SizeRezMax)  
//
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
     FreeCSC(IJRez)
//""
     FreeLFF(LFF)
//
     delete [] Size;
     delete [] RezRel;
     delete [] Omega;
     delete [] PermutGen;
     delete [] Mas;
     delete [] CoorCart;
     delete [] Corres;
     delete [] ZeroSt;
     delete [] workd;
     delete [] workl;
     delete [] EigVec;
     delete [] EigVal;
     delete [] Pid;
     delete [] RezVect;
     delete [] TabScreen;
     delete [] PositionTarget;
     delete [] PermutRez;
     delete [] PermutAct;
     delete [] TabNull;
     delete [] EigValOld;   
     delete [] AssignRez;
     delete [] ValAct;      
//
    FreeMatricesElem(QQ,DegrePol+1+4,DegrePol+MaxQLevel)
// 
FinalMessaj()
//
 return 0;
}
