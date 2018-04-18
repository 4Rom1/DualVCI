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
#include "Shared.h"
#include <sys/time.h>
#include "Graph.h"
#include "ReadData.h"
#include "Assemblage.h"
#include "Basic.h"
#include "Solver2.h"
#include "Generator.h"
#include <cassert>
//
int main(int argc, char *argv[])
{
//Files to read    
    FILE *FileKey = NULL;
//
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
      printf("\n***Usage: ../../source/FinalVCI Key.in***\n");
      printf("One argument requested : Key file missing.  \n");
      printf("DVCI should have been previously run with PrintOut>1. \n\n");
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
int DoRot=1;//If > 0 then Coriolis terms are computed
// with non mass weighted normal coordinates for odd and mass weighted for even 
int PESType=1;//if non 0 derivatives in a.u else force constants in cm-1
int MAXNCV=0;//Number of ritz vectors for the eigensolver
int MaxQLevel=15;//Maximal quantum level in each direction 
//Be carefull when increasing potentials might be very poluted around 20 
double MinFreq=-100;// [MinFreq,MaxFreq]= frequency window target
double MaxFreq=4000;// 
double GroundState=-1;//Groundstate should be given by the user.
//except when it is actually computed then it has negative value.
double Kappa=1.2;//Multiplicatif factor for intial subspace construction 
double Tol=1e-8;//Threshold for the eigensolver (same as tol in DSAUPD) 
int NAdd=10;//Number of added elements at each iteration  
int MaxAdd=1000;//Maximal number of added elements at each iteration  
double EtaComp=1.5;//New basis set will be selected from the maximal component
// of error vector divided by EtaComp*NNotConv, where NNotConv is the number of
// non converged eigenvectors at a given iteration
char OutName[MaxChar]="NoName";
double ThrMat=1e-20; //Threshold for matrix elements
double ThrCoor=0.6; //Threshold for component selection of targeted vectors
int AddTarget=2; //Number of additional residual vectors to correct the potential increasing of targets
char PESName[MaxChar]={0}; //Name of file of PES to be read.
char RefName[MaxChar]={0};//Reference txt file containing frequencies in increasing order
int MaxEv=30;//Maximal number of wanted eigenvalues for the Arnoldi eigensolver
int DeltaNev=1000; //Adjustment number for MaxEV that will be equal to Min(MaxEV,  MaxScreen+DeltaNev) 
//where PosEig is the eigenvalue position of maximal target
double Memory=0;//Total allocated Memory in Megabytes
int PrintOut=0; //If >0 print final basis set and eigenvectors for post treatment 
double EpsRez=6e-3; //Convergence criteria for global residue
int Verbose=0;//If > 0 printout additional informations between iterations
float ThrKX=1; //Threshold for contributions of sum force constants in dual operator
double ThrPES=1e-20; //Threshold for derivative or force constants of the PES
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
int DoGraph=1; // If not equal to zero, it will store the graph of residual matrix in CSC format
int DoVPT=0;//If positif add correction vector and correction energy to targets
//to gain computational time and then significantly increases memory requirement
//
printf("\n*****Reading Datas*****\n");
//
int NTargetStates=GetKeyWords(FileKey, &NMode, &DoRot, &MaxEv, &DeltaNev, &AddTarget,\
    &PESType, TargetState, &DoGraph, &MAXNCV,\
    &MaxQLevel, &NAdd, &MaxAdd, &PrintOut, &DoVPT, &EpsRez, &KNREZ, &KNNZ, &KNZREZ, &EtaComp,\
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
      {printf("\n Threshold for PES ThrPES = %e,\n Threshold for Hamiltonian matrix elements ThrMat = %e \n", ThrPES, ThrMat);}
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
if(PrintOut <=1)
{
printf("\n***PrintOut must be > 1 and matrix saved from previous DVCI run***\n\n");
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
//Get The total time.
gettimeofday(&BeginAll,NULL);
CPUBeginAll = clock();


FILE *FileRef;//To compare with reference frequencies in output 
FileRef=fopen(RefName,"r");
if(FileRef==NULL)
{
printf("\n FileRef not detected \n");
}
else{printf("\n FileRef detected \n");}
//
  char OutMat[MaxChar]={0};
//
  char OutBasis[MaxChar]={0};
//
  char BasisOut[30]="-FinalBasis.bin";//Final basis file name if PrintOut = 2
//
  char MatOut[30]="-Matrix.bin";//Binary file name with matrix coefficients if PrintOut = 2
//
  sprintf(OutMat, "%s", OutName);
  sprintf(OutBasis, "%s", OutName);
//
  strcat(OutMat, MatOut);
  strcat(OutBasis, BasisOut);
//
  FILE *FileMat=NULL;
  FILE *FileBasis=NULL;
//     
  FileMat=fopen(OutMat,"rb");
  FileBasis=fopen(OutBasis,"rb");   
//
  if (FileMat == NULL || FileBasis==NULL)
   {
   printf("Error opening files \n");
   FinalMessaj()
   return 0;
   }
//
uint64_t NNZAct; 
double Shift, tol, Ull;
CSC IJAct;
int NEV,NCV;
double *ValAct;
uint32_t DimAct,SizeInit;
//
 if(!ReadMatCSCBin(FileMat, &NEV, &NCV, &SizeInit, &DimAct, &NNZAct, &Shift, &tol, &Ull, IJAct, ValAct))
  {
  FinalMessaj()
  return 0;
  }
//
//
int NScreen=0;//NUmber of screened states in output
//
uint8_t *ModeAct=NULL;
//ConfigId *ModeAct=NULL;
////Configuration space
 if(!GetConfsBin(FileBasis, ModeAct, NMode, &NScreen))
  {
  FinalMessaj()
  return 0;
  }
//
      allocate_TabSize(Size,2);
      Size[1].DimAct=DimAct;
      Size[1].NNZAct=NNZAct;
//
//External allocation of arrays for eigensolver
      double *workd=NULL;
      workd=new double [3*DimAct+1];
      InitTab(workd,3*DimAct+1) 
//
      double *workl=NULL;
      workl = new double[(NCV)*(NCV+8)+1]; 
//      InitTabInt(workl,(NCV)*(NCV+8)+1)  
      double *EigVec=NULL;
      EigVec=new double[DimAct*NCV+1];   
//      InitTab(EigVec,DimAct*NCV+1) 
      double *EigVal=NULL;
      EigVal=new double[NEV];   
//
//Detect position Tragets for printing out
     uint8_t *ZeroSt=new uint8_t[NMode];     
     InitTabInt(ZeroSt, NMode)    
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
   printf("\n---Effective memory usage %.1f Mo---\n", EffectiveMemory2(argv[0]));
//
   int nconv=SolverCSCSym(Size, 0, NEV, NCV, Shift, Tol, EigVal,\
   EigVec, workd, workl,IJAct, ValAct);
//
        if(!nconv)
        {
        printf("\n**Non convergence on eigensolver**\n");        
        }
//
 if(nconv)
   {
    double TabRef[MaxRef]={0}; //Array to store the reference values
    int ss=0;     
    double Ground=GroundState; 
    int DetectZero=0; //If >0 the groundstate has been detected
    int detect=1; //If it remains positive one of the Target has been detected 
    long PosZero=-1; //Position of groundstate if detected (might be not equal to zero) 
    if(FileRef!=NULL)
     {
                printf("\n %d lines detected in Reference file\n",GetValRef(FileRef, TabRef));
                printf("\n Assignment  | Frequency | Error  || \n");
                printf(" (Component) | (number)    | Here-Ref || \n");
                printf("------------------------------------ \n");              

     for (int PosEig=0;PosEig<nconv;PosEig++)
     {   
    DetectZero=0; //If >0 the groundstate has been detected
    detect=1; //If it remains positive one of the Target has been detected 
    PosZero=-1; //Position of groundstate if detected (might be not equal to zero)             
        if(GroundState<=0)
        {
        DetectZero=DetectMaxCoordinate(&EigVec[PosEig*Size[1].DimAct], PositionTarget, 1,ThrCoor);
        if(DetectZero)
         {
         PosZero=PosEig;     
         Ground=EigVal[PosEig];
         } 
        }
        if(NTargetStates && !DetectZero)
        {detect=DetectMaxCoordinate(&EigVec[PosEig*Size[1].DimAct], PositionTarget, NTargetStates,ThrCoor);}    
          if( ((EigVal[PosEig]-Ground <= MaxFreq) && (EigVal[PosEig]-Ground >= MinFreq)) && (detect || DetectZero) )
          { 
           //Assignment display
           AfficheNuSupTex(ModeAct, &EigVec[PosEig*Size[1].DimAct], NMode, Size[1].DimAct, ThrCoor);
//                  
           if(PosEig==PosZero)
           {
           if(!DoRot)
            {
            printf(" | %.2f(%d) | %.4f // \n \n", \
            Ground,PosEig,GetClosest(TabRef,Ground)-(Ground));
            }
           else//Add Watson
            {
           printf(" | %.2f(%d) | %.4f  | %.4f (Watson) // \n \n", \
           Ground,PosEig,GetClosest(TabRef,Ground)-(Ground),Ull);//Add watson
            }
           }
           else
           {
           printf(" | %.2f(%d) | %.4f // \n \n", \
           EigVal[PosEig]-Ground,PosEig,GetClosest(TabRef,EigVal[PosEig]-Ground)-(EigVal[PosEig]-Ground));
           }         
           ss++;                                                                                                          
           }
        }
      } 
   else
     {
                printf("\n Assignment  | Frequency  || \n");
                printf(" (Component) | (number)     || \n");
                printf("------------------------------------ \n");              

     for (int PosEig=0;PosEig<nconv;PosEig++)
     {   
    DetectZero=0; //If >0 the groundstate has been detected
    detect=1; //If it remains positive one of the Target has been detected 
    PosZero=-1; //Position of groundstate if detected (might be not equal to zero)              
        if(GroundState<=0)
        {
        DetectZero=DetectMaxCoordinate(&EigVec[PosEig*Size[1].DimAct], PositionTarget, 1,ThrCoor);
        if(DetectZero)
         {
         PosZero=PosEig;     
         Ground=EigVal[PosEig];
         } 
        }
        if(NTargetStates && !DetectZero)
        {detect=DetectMaxCoordinate(&EigVec[PosEig*Size[1].DimAct], PositionTarget, NTargetStates,ThrCoor);}    
          if( ((EigVal[PosEig]-Ground <= MaxFreq) && (EigVal[PosEig]-Ground >= MinFreq)) && (detect || DetectZero) )
          { 
           //Assignment display
           AfficheNuSupTex(ModeAct, &EigVec[PosEig*Size[1].DimAct], NMode, Size[1].DimAct, ThrCoor);
//                  
           if(PosEig==PosZero)
           {
           printf(" | %.2f(%d)  // \n \n", \
           Ground,PosEig);
           }
           else
           {
           printf(" | %.2f(%d) // \n \n", \
           EigVal[PosEig]-Ground,PosEig);
           }         
           ss++;                                                                                                          
           }
         }
      }  
    }
//
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
    if(FileRef!=NULL){fclose(FileRef);}
//
fclose(FileMat);
//
     FreeCSC(IJAct)
//
//'' FreeTabMode(ModeAct,DimAct)
// 
     delete [] ModeAct;
     delete [] Size;
     delete [] workd;
     delete [] workl;
     delete [] EigVec;
     delete [] EigVal;
     delete [] PositionTarget;
     delete [] ValAct;  
     delete [] ZeroSt;
}
