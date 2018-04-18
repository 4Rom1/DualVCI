/*
    Dual Vibration Configuration Interaction (DVCI).
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


--------Additional sources-------------------

Routines QuickSortMode, QuickSortModeShift are made from the Quicksort algorithm [1]

[1] C. A. R. Hoare. 1961. Algorithm 64: Quicksort. Commun. ACM 4, 7 (July 1961), 321-. DOI=http://dx.doi.org/10.1145/366622.366644 

Routines QsearchMode, QsearchModeShift are made from the half interval search algorithm also known as binary search [2]

+

[2] Willams, Jr., Louis F. (1975). A modification to the half-interval search (binary search) method. Proceedings of the 14th ACM Southeast Conference. pp. 95–101. doi:10.1145/503561.503582.


----------------------------------------------

The routines nested_loop and nested_loop1 are depth variable nested loops sharing
 the same structure as the one presented in STACKOVERFLOW [3]

[3] How to set a variable for nested loops? 
Posted on stackoverflow (by M Oehm) at the following address
http://stackoverflow.com/questions/25785863/how-to-set-a-variable-for-nested-loops
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "Graph.h"
#include "Shared.h"
#include "Basic.h"
#include <sys/time.h>
//
int nested_loop(uint8_t *ix, uint8_t *im, int depth)
{
//Simulation of a nested loop of a given depth to generate multi-indexes stored
//in ix with ix[Cpt] bounded by im[Cpt].
//ix must be initialized by 0.
//
    int Cpt = 0;
//
    do {
        ix[Cpt]++;
        if(ix[Cpt] < im[Cpt]) {return 1;}
        ix[Cpt] = 0;
        Cpt++;
       } while (Cpt < depth);
//
    return 0;
}
//
int nested_loop1(uint8_t *ix, uint8_t *im, int depth)
{
//Same as nested_loop but ix is initialized to one instead
    int Cpt = 0;
    do {
        ix[Cpt]++;
        if(ix[Cpt] < im[Cpt]) {return 1;}
        ix[Cpt] = 1;
        Cpt++;
    } while (Cpt < depth);
//
    return 0;
}
//
//
long LinearModeSearch(uint8_t *TabMode, uint8_t *ModeT, size_t NMode,  uint32_t SizeTab)
{
//Look if the multi-indexes ModeT is located in the array of multi-indexes TabMode.
//Size of TabMode is SizeTab*SizeBit.
//Don't forget the case where the size of array is zero.
//size_t SizeBit=NMode*sizeof(uint8_t)=NMode;
   uint32_t Cmpt=0;
 //memcmp return 0 if the configurations are identical else it returns a non nul value.
//''
 while(memcmp(&TabMode[Idm(Cmpt)], ModeT, NMode))
    {   
//"" 
//Need to add NMode to parameters
     Cmpt++;
     if(Cmpt >= SizeTab)
      {return -1;} 
    }
    return Cmpt;
//
}
//
  long QsearchMode(uint8_t *Key, uint8_t *base, uint32_t *Permuter, uint32_t size, int NMode)
  {
//Binary search of configurations "Key" in basis "base" 
//already sorted in the memcmp order given by the indexes of Permuter array
          uint32_t debut = 0, fin = size-1;
          uint32_t milieu=0;
          int Cmp=0;  
//""
          int Cmp0=memcmp(Key, &base[Idm(Permuter[debut])], NMode);
          int Cmp1=memcmp(Key, &base[Idm(Permuter[fin])], NMode);
//         
          if(Cmp0<0 || Cmp1 >0)
             {return -1;}   
          else if (!Cmp0)
             {return Permuter[debut];}
          else if (!Cmp1)
             {return Permuter[fin];}   
//
          while (debut < fin) 
           {
                  milieu = debut + (fin - debut) / 2;  
                    Cmp = memcmp(Key, &base[Idm(Permuter[milieu])], NMode);
                  if (Cmp < 0)
                       {
                          fin = milieu;
                       }
                  else if (Cmp > 0)
                      {
                          debut = milieu + 1;
                      }
                  else
                      {
                          return Permuter[milieu];
                      } 
         }
//  
          return -1;
  }
//
 long QsearchModeShift(uint8_t *Key, uint8_t *base, uint32_t *Permuter, uint32_t size, int NMode, uint64_t Shift)
  {
//Same as previous subroutine : Binary search of configurations "Key" in basis "base" 
//already sorted in the memcmp order given by the indexes of Permuter array.
//But indexes are shifted by the number Shift in order to avoid troubles when they oversize MAXUINT32
//Because Permuter is an array of UINT32_t
          uint32_t debut = 0, fin = size-1;
          uint32_t milieu=0;
          int Cmp=0;     
//""
          int Cmp0=memcmp(Key, base+Idm(Permuter[debut]+Shift), NMode);
          int Cmp1=memcmp(Key, base+Idm(Permuter[fin]+Shift), NMode);
//
//    
          if(Cmp0<0 || Cmp1 >0)
             {return -1;}   
          else if (!Cmp0)
             {return Permuter[debut]+Shift;}
          else if (!Cmp1)
             {return Permuter[fin]+Shift;}   
//
          while (debut < fin) 
           {
                  milieu = debut + (fin - debut) / 2;  
//""
                  Cmp = memcmp(Key, base+Idm(Permuter[milieu]+Shift), NMode);
                  if (Cmp < 0)
                       {
                          fin = milieu;
                       }
                  else if (Cmp > 0)
                      {
                          debut = milieu + 1;
                      }
                  else
                      {
                          return Permuter[milieu]+Shift;
                      } 
         } 
          return -1;
  }
//
int32_t partition (uint8_t *arr, uint32_t Permuter[], int ll, int rr, int NMode)
{
/*
Permuter[rr] is the pivot element.
All the smaller and bigger elements than the pivot are respectively 
placed to the left and the right of the latter.
Comparison of elements is done with memcmp c subroutine.
 */
    int32_t ii = (ll - 1);
 
    for (int32_t jj = ll; jj <= rr- 1; jj++)
    {
//""    
        if(memcmp(&arr[Idm(Permuter[jj])],&arr[Idm(Permuter[rr])], NMode)<=0) 
        {
            ii++;
            swap (&Permuter[ii], &Permuter[jj]);
        }
    }
    swap (&Permuter[ii + 1], &Permuter[rr]);
    return (ii + 1);
}
//
int32_t PartShift (uint8_t *arr, uint32_t Permuter[], int ll, int rr, int NMode, uint64_t Shift)
{
/*
Permuter[rr]+Shift is the pivot element.
All the smaller and bigger elements than the pivot are respectively 
placed to the left and the right of the latter.
There is an additional shift to avoid memory overflow of the residual indexes than can be bigger than UINT32MAX.
Comparison of elements is done with memcmp c subroutine.
*/
    int32_t ii = (ll - 1);
 
    for (int32_t jj = ll; jj <= rr- 1; jj++)
    {
//      
        if(memcmp(&arr[Idm(Permuter[jj]+Shift)],&arr[Idm(Permuter[rr]+Shift)], NMode)<=0) 
        {
//""
            ii++;
            swap (&Permuter[ii], &Permuter[jj]);
        }
    }
    swap (&Permuter[ii + 1], &Permuter[rr]);
    return (ii + 1);
}
//
void QuickSortMode (uint8_t *MIndex, int lo,int hi, unsigned int *Permuter, int NMode)
{
//Sort multi-indexes in MIndex according to memcmp order and return ascending order indexes in Permuter
  if (lo < hi)
    {
        // pt is partitioning index of MIndex[Permuter[pt]] 
        int pt = partition(MIndex, Permuter, lo, hi, NMode);

        QuickSortMode(MIndex,lo, pt - 1, Permuter, NMode);  // Before pt
        QuickSortMode(MIndex, pt + 1, hi, Permuter, NMode); // After pt
    }
}
//
void QuickSortModeShift (uint8_t *MIndex, int lo,int hi, uint32_t *Permuter, int NMode, uint64_t Shift)
{
//Sort multi-indexes in &MIndex[+Shift] according to memcmp order and return ascending order indexes in Permuter
  if (lo < hi)
    {
// pt is the partitioning index of MIndex[Permuter[pt]+Shift] 
        int pt = PartShift(MIndex, Permuter, lo, hi, NMode, Shift);
//
        QuickSortModeShift(MIndex,lo, pt - 1, Permuter, NMode, Shift);  // Before pt
        QuickSortModeShift(MIndex, pt + 1, hi, Permuter, NMode, Shift); // After pt
    }
}
//
uint32_t InitB0(uint8_t  *ModeAct, KForce KFC,\
   uint8_t *Pid, int NMode, uint32_t SizeActMax, double TargetMax)
{
//Compute the initial subspace En<TargetMax, where En harmonic energy
//and n are configuration arrays stored in ModeAct.
//Use a greedy algorithm that sweeps all the region until TargetMax.
//Here is used only one modal nested excitations on the ground state (0,...,0).
//Pid[mm] is the maximal quantum level + 1 in each direction mm 
//Works because of order of energies.
//Force constants:
//KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
//ModeAct : Multi-dimensional array for active space B : ModeAct[Idm(nn)+mm], (nn,mm) in [0,SizeActMax[ x [0,NMode[.
uint32_t Lin=0;
uint8_t *Tester=NULL; //Variable multi index during the loop
//Size of bytes to look for in linearsearch arrays 
size_t SizeBit=NMode*sizeof(uint8_t);      
//
int ss=0; //Index for surface
int go=1; //Index to continue if Tester is not negativ
int TestNeg=0; //Test if the sum is negative, here only check the max dimension
int Cpld;
double Freq0, Freq0T;
uint32_t SizeAv=0;
uint32_t SizeAvAv=0;
long CheckIn=0;
//
uint32_t SizeInit=1; //Zero state always included in position 0
allocate_TabUint8(Tester,NMode)
// 
int DegInc=1; //Excitation=one modal degree increment in each direction
//
//        
while(SizeInit-SizeAv) //Self consistency on sizes.
{
//
SizeAvAv=SizeAv;
SizeAv=SizeInit;
//
for (Lin=SizeAvAv ; Lin < SizeAv ; Lin++)
 {
//""
 Freq0=GetFreq0(&ModeAct[Idm(Lin)],NMode, KFC);
//One modal excitations
Cpld=1;
//
 for (int mm=0; mm< NMode; mm++)
  {
//  
//""
    memcpy(Tester,&ModeAct[Idm(Lin)], SizeBit);
//
//             
    ss=0;
    go=1;
   while ( (ss<Cpld) && go )
    {
//""
//ModeAct[Idm(Lin)+mm]
    TestNeg=DegInc+ModeAct[Idm(Lin)+mm];
    if(TestNeg >= Pid[mm]){go=0;}
    else
     {
      Tester[mm]=TestNeg; 
     }
    ss++;
    }   
if(go)
   {     
    Freq0T=Freq0;   
    for(ss=0;ss<Cpld;ss++)
     {
     Freq0T+=DegInc*KFC.KijNCpld[1][mm];
     }
//
if(Freq0T<=TargetMax)
   {
   if(SizeInit<SizeActMax)
   {
//Test if Tester already constitutes one of the neighboors of the configuration 
//It has to be negative if not.
//Tester=0 and Size0.DimAct-SizeAv : no problem because CheckIn >=0 
//Tester non equal zero no problem either because not in previous checking.
//Checking can be done starting from previous layer because before the gap
//will be bigger than 1 degree difference
//''
  CheckIn=LinearModeSearch(&ModeAct[Idm(SizeAvAv)], Tester, SizeBit,SizeInit-SizeAvAv+1);  
//
    if(CheckIn<0)
        {
//""
    memcpy(&ModeAct[Idm(SizeInit)],Tester,SizeBit); 
        SizeInit++;
        }
       }  
  else{printf("***Init subpsace oversized->Increase memory***\n");return 0;}
      }
     }
    }
   }
  }  
 delete [] Tester;
 return SizeInit; 
}
//
int ConvertX(uint8_t *Tab, uint8_t *Surf,int NMode) 
{
//Convert the multi-index Tab of size NMode
//on surface Surf={Positions({Tab[ii] # 0})}
   int mon_i=0;  
   int ss=0;
   for (mon_i=0; mon_i<NMode;mon_i++)
    { 
     if(Tab[mon_i]!=0)
     {Surf[ss]=mon_i;ss++;}
    }
    return ss;
}
//
void VoisinCpld (int IndexMonm, int NMode, KForce KFC, uint8_t *NVMode,uint8_t **ModeV)
{
/*Computes the neighboors of the zero state (0,0,0..0)=Mode0
considering the monom number IndexMonm.
H=sum_kk prod_mm H_{mm,kk}. Here kk=IndexMonm.
Return the Neighboor-indexes of the normal coordinate mm
available in the array ModeV[NVMode[mm]][mm]:
<ModeV(NVMode(mm),mm)|H_{mm,IndexMonm}|Mode0(mm)> # 0
*/
   InitTabInt(NVMode, NMode)      
//       
//Browse all the normal coordinates        
 for (int mm=0; mm< NMode; mm++)
   {     
      if (KFC.Monm[IndexMonm][mm]>0)
      {   
//First even degrees neighboors         
        int dmin=KFC.Monm[IndexMonm][mm]%2;
        int dmax=KFC.Monm[IndexMonm][mm];
//      
     if (dmax!=0)
        {
          for (int dd=dmin;dd<=dmax;dd+=2)
           {
            ModeV[NVMode[mm]][mm]=dd;
            NVMode[mm]++;           
           }            
        }      
      }
//Degree 0.
     else if( KFC.Monm[IndexMonm][mm]==0 )
     {
            ModeV[NVMode[mm]][mm]=0;
            NVMode[mm]++;
     }   
   }
}
//
void VoisinCpldRot (int NMode, uint8_t *Tester, uint8_t *NVMode,uint8_t **ModeV)
{
/*Computes the rotational terms neighboors of zero state Mode0=(0,0,0..0).
Return neighboor-indexes of the normal coordinate mm available in the array ModeV[NVMode[mm]][mm]:
<ModeV(NVMode(mm),mm)|H_{mm,IndexMonm}|Mode0(mm)>  # 0 .
Tester is the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors) 
*/
   InitTabInt(NVMode, NMode)     
//            
 for (int mm=0; mm< NMode; mm++)
   {     
      if (Tester[mm]>0)
      {      
//          
        int dmin=Tester[mm]%2;
        int dmax=Tester[mm];
//      
     if (dmax!=0)
        {
          for (int dd=dmin;dd<=dmax;dd+=2)
         {
            ModeV[NVMode[mm]][mm]=dd;
            NVMode[mm]++;           
           }            
        }       
      }
     //Degree 0.
     else if(Tester[mm]==0)
     {
            ModeV[NVMode[mm]][mm]=0;
            NVMode[mm]++;
     }  
   }
}
//
int LinearIntSearch(int *TabInt, int Key,  int SizeTab)
{
//Search for the int Key in array TabInt of size SizeTab.
 if (SizeTab)
 {
   int Cmpt=0;
 //memcmp return 0 if the configurations are identical else it returns a non nul value.
  while(TabInt[Cmpt]-Key)
    {
     Cmpt++;
     if(Cmpt >= SizeTab)
     {   
     return -1;
     } 
    }
    return Cmpt;
 }
else{return -1;}
}
//
 long LinearU8Search(uint8_t MatU8[MaxTarget][MaxNormal], uint8_t *ModeKey, size_t SizeBit,  uint32_t SizeTab)
{
//Search for the array ModeKey of size SizeBit in matrix MatU8.
//size_t SizeBit=NMode*sizeof(uint8_t);
if (SizeTab)
 {
  long Cmpt=0;
//memcmp return 0 if the configurations are identical else it returns a non nul value.
  while(memcmp(MatU8[Cmpt], ModeKey, SizeBit ))
    {
     Cmpt++;
     if(Cmpt >= SizeTab)
     {   
     return -1;
     } 
    }
    return Cmpt;
 }
else{return -1;}
}
//
void AfficheNu(uint8_t *Mode, int NMode)//Print Modal components
{//Print Modal components : dd*nu[mm]+, ... for a normal coordinate mm and
//a non null degree dd (corresponding to a degree of Hermite basis function)
int CmptMode=0;\
for (int mm = 0 ; mm < NMode ; mm++){\
if((Mode[mm]>1) && (!CmptMode)){\
printf(" %dnu[%d]", Mode[mm],mm+1);\
CmptMode++;}\
else if((Mode[mm]>1) && (CmptMode)){\
printf("+%dnu[%d]", Mode[mm],mm+1);}\
else if((Mode[mm]==1) && (!CmptMode)){\
printf(" nu[%d]",mm+1);\
CmptMode++;}\
else if((Mode[mm]==1) && (CmptMode)){\
printf("+nu[%d]",mm+1);}}\
if(!CmptMode){printf(" nu[0]");}
}
//
void AfficheNuTex(uint8_t *Mode, int NMode)//Print Modal components in Tex format
{//Print Modal components in Tex : dd*\nu{mm}+, ... for a normal coordinate mm and
//a non null degree dd (corresponding to a degree of Hermite basis function)
int CmptMode=0;\
for (int mm = 0 ; mm < NMode ; mm++){\
if((Mode[mm]>1) && (!CmptMode)){\
printf("%d\\nu_{%d}", Mode[mm],mm+1);\
CmptMode++;}\
else if((Mode[mm]>1) && (CmptMode)){\
printf("+%d\\nu_{%d}", Mode[mm],mm+1);}\
else if((Mode[mm]==1) && (!CmptMode)){\
printf("\\nu_{%d}",mm+1);\
CmptMode++;}\
else if((Mode[mm]==1) && (CmptMode)){\
printf("+\\nu_{%d}",mm+1);}}\
if(!CmptMode){printf("\\nu_{0}");}
}
//
void AfficheNuSupTex(uint8_t  *ModeAct, double *EigVec, int NMode, uint32_t DimAct, double ThrCoor)
{//Print Modal components in Tex format and vector components bigger than ThrCoor
     int Cpt=0;
     double Val;
     for (uint32_t ii=0;ii<DimAct;ii++)
           {
            Val=EigVec[ii]*SIGN<double>(EigVec[ii]);
            if(Val>ThrCoor && Cpt)
            {
           printf(", ");
           printf("$");
//""
     AfficheNuTex(&ModeAct[Idm(ii)], NMode); 
           printf("$");
           printf("(%.2f)",Val);
           Cpt++; 
            }
            else if(Val>ThrCoor && !Cpt)
            {
           printf("$");
//""
     AfficheNuTex(&ModeAct[Idm(ii)], NMode); 
           printf("$");
           printf("(%.2f)",Val);
           Cpt++;
            }
          }
}
//
void AfficheNuSup(uint8_t  *ModeAct, double *EigVec, int NMode, uint32_t DimAct, double ThrCoor)
{
//Print Modal components and vector components bigger than ThrCoor 
     int Cpt=0;
     double Val;
     for (uint32_t nn=0;nn<DimAct;nn++)
           {
           Val=EigVec[nn]*SIGN<double>(EigVec[nn]);
            if(Val>ThrCoor && Cpt)
            {
           printf(", ");
//""
         AfficheNu(&ModeAct[Idm(nn)], NMode); 
           printf("(%.2f)",Val);
           Cpt++; 
            }
            else if(Val>ThrCoor && !Cpt)
            {
//"" 
         AfficheNu(&ModeAct[Idm(nn)], NMode); 
           printf("(%.2f)",Val);
           Cpt++;
            }
          }
}
//
uint32_t CptLFFPerEx (int DegrePol, int NMode, int NPES, int DoRot, uint32_t NGenPoz,\
 KForce KFC, uint8_t *XPolSupPos, int *NLFFPerEx, double **ZetaXYZ, double ThrPES, int IncK2)
{
//Routine that counts the number of force constants per local force field.
//The routine return the number of positive excitations in H*.
//XPolSupPos contains the upper limit of positive excitations in H*.
//XPolSupPos[Idm(xx)+mm] is defined for (xx,mm) in [0,NGenPoz[ x [0,NMode[.
//Force constants:
//KFC.KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
//Local force field:
//LFF[ii(xx)] : local force field associated with positive excitation &DualHPos[Idm(xx)].
//LFF.Idx[ii(xx)] are defined for ii [LFF.Num[xx],LFF.Num[xx+1][
//LFF.Idx[ii(xx)] < 0 are indexes of non coupled force constants:
//KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
//0 <= LFF.Idx[ii(xx)] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF.Idx[ii(xx)]].
//LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
//with a unique corresponding (ni,nj,nk,nl)  
//NLFFPerEx[xx] : output array counting the number of force constant per excitation XPolSupPos[xx]
//ThrPES : threshold for acceptable force constants.
//NPES : Number of coupled PES coefficients stored in KFC.KijCpld
//ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}
//IncK2:Say if the the quadratic term should be included or not in the local force field
       uint8_t *NVMode; // Number of neighboors by dimension.
       uint8_t **ModeV; // Neighboors in each dimension. 
//Neighboor of index mm is given by given by ModeV[NVMode[mm]][mm].
       uint8_t *Tester; // Variable multi index during the loop
       uint8_t *MultiDegrees; //Degrees of the nested loop
//Allocation
       allocate_MatUint8(ModeV, DegrePol+1, NMode)
       allocate_TabUint8(NVMode,NMode)
       allocate_TabUint8(MultiDegrees, NMode)
       allocate_TabUint8(Tester,NMode)
//Initialisation of the degrees of the maximal mode
//int Cmpt=1; //First is zero degree
       int Deg=0;
       double ContribLFF=0;//Contribution of force constants summed up in absolute value
       long IndexEx;      
       size_t SizeBit=NMode;
//       
       uint32_t *Permuter=NULL;
       Permuter=new uint32_t [NGenPoz]; 
       uint32_t NGenPozUI=NGenPoz;
//
      for (uint32_t ii = 0; ii < NGenPozUI; ii++)
      {
      Permuter[ii]=ii;
      }
//
QuickSortMode(XPolSupPos, 0,NGenPoz-1, Permuter, SizeBit);
//
//Non coupled surfaces first for Local force constants identification
for (int mm=0; mm<NMode; mm++)
 { 
  InitTabInt(Tester, NMode)
 for (int kk=0; kk<DegrePol; kk++)
   { 
    if(kk!=1 || IncK2) //Avoid Harmonic
     {
    ContribLFF=SIGN<double>(KFC.KijNCpld[kk][mm])*KFC.KijNCpld[kk][mm];
     }
    else
     {
     ContribLFF=0;
     }
      Deg=kk+1;
     if(ContribLFF > ThrPES)
     {
      while(Deg >= 0 )
      {
      Tester[mm]=Deg;
      IndexEx=QsearchMode(Tester, XPolSupPos,Permuter,NGenPoz, SizeBit);
      NLFFPerEx[IndexEx]++;
      Deg-=2;
      }
     } 
    }
   }
//Local force constants identification for coupled terms
for (int kk=0; kk<NPES; kk++)
 {           
 ContribLFF=SIGN<double>(KFC.KijCpld[kk])*KFC.KijCpld[kk];
 if(ContribLFF>ThrPES)
  {
   VoisinCpld (kk, NMode, KFC, NVMode, ModeV);
//Mandatory intialisation of multi-degrees. 
    InitTabInt(MultiDegrees, NMode)
do {     
     for (int mm = 0 ; mm < NMode; mm++)
       {
    Tester[mm]=ModeV[MultiDegrees[mm]][mm];
       }
      IndexEx=QsearchMode(Tester, XPolSupPos,Permuter,NGenPoz, SizeBit);      
      NLFFPerEx[IndexEx]++;
    } while(nested_loop(MultiDegrees, NVMode, NMode));
   }
 }
//Identification of Local Coriolis Interactions                    
if(DoRot)
{
 for (int ni=0;ni<NMode;ni++)
  {
 for (int nj=ni+1;nj<NMode;nj++)
    {
//Avoid double count ijkl
  for (int nk=0;nk<NMode;nk++)
     {
    for (int nl=nk+1;nl<NMode;nl++)
      {
//COntribution is simply given by ZetaXYZ coefficients
       ContribLFF=ZetaXYZ[ni + ((nj-1)*nj)/2][nk + ((nl-1)*nl)/2];
       ContribLFF*=SIGN<double>(ContribLFF);
//
    if(ContribLFF > ThrPES)
      {
//Tester is the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors) 
     ContribRot(Tester, ni, nj, nk, nl, NMode);
     VoisinCpldRot (NMode, Tester, NVMode, ModeV);  
     InitTabInt(MultiDegrees, NMode)  
     do 
      {    
     for (int mm = 0 ; mm < NMode; mm++)
        {
    Tester[mm]=ModeV[MultiDegrees[mm]][mm];
        }
//
      IndexEx=QsearchMode(Tester, XPolSupPos,Permuter,NGenPoz, SizeBit);      
      NLFFPerEx[IndexEx]++;
//
       } while(nested_loop(MultiDegrees, NVMode, NMode));           
      }
     }
    }
   }
  }
}
      uint32_t NXDualHPlus=0;
//
      for (uint32_t ii = 0; ii < NGenPozUI; ii++)
      {
       if(NLFFPerEx[ii])
       {NXDualHPlus++;}
      }
//
         delete [] NVMode;
         delete [] MultiDegrees;
         delete [] Permuter;
         delete [] Tester;
//         
         FreeMat(ModeV, DegrePol+1)
//
         return NXDualHPlus;
}
//
uint32_t AssignLFFPerEx (int DegrePol, int NMode, int DoRot, int NPES, uint32_t NXDualHPlus,\
 KForce KFC, uint8_t *DualHPos, LocalFF LFF, double **ZetaXYZ, double ThrKX, double ThrPES,\
 uint32_t *Permuter, uint32_t *Corres, int IncK2)
{
//Routine that copies indexes of positive excitations of most contributing Local Force Fields
//into Corres.
//Corres : Correspondance array giving the indexes of the first NXDualHTruncPos most contributive positive excitations of H*.
//Select only the ones such that Sum|FC|>ThrKX
//DualHPos contains positive excitations in H*.
//DualHPos[Idm(xx)+mm] is defined for (xx,mm) in [0,NXDualHPlus[ x [0,NMode[.
//Force constants:
//KFC. KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode).
//KFC.KijCpld[kk] : coupled force constants, defined for kk in [0,NPES[.
//KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
//Local force field:
//LFF[ii(xx)] : local force field associated with positive excitation &DualHPos[Idm(xx)].
//LFF.Idx[ii(xx)] are defined for ii [LFF.Num[xx],LFF.Num[xx+1][
//LFF.Idx[ii(xx)] < 0 are indexes of non coupled force constants:
//KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode;
//0 <= LFF.Idx[ii(xx)] < NPES and are indexes of coupled force constants KFC.KijCpld[LFF.Idx[ii(xx)]].
//LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES
//with a unique corresponding (ni,nj,nk,nl)  
//Permuter : permutation index array of elements of DualHPos (memcmp order)
//ThrPES : threshold for acceptable force constants.
//NPES : Number of coupled PES coefficients stored in KFC.KijCpld
//ZetaXYZ[ijkl]=\sum_{(\alpha ,\beta)\in (x,y,z)} \mu_{\alpha,\beta}\zeta_{ij}^{\alpha}\zeta_{kl}^{\beta}
//NXDualHTruncPos : output = Total number of raising excitations for H* (after truncation through ThrKX)
//IncK2:Say if the the quadratic term should be included or not in the local force field
//
       uint8_t *NVMode; // Number of neighboors by dimension.
       uint8_t **ModeV; // Neighboors in each dimension.
       uint8_t *Tester; //Variable multi index during the loop
       uint8_t *MultiDegrees; //Degrees of the nested loop
//Allocation
       allocate_MatUint8(ModeV, DegrePol+1, NMode)
       allocate_TabUint8(NVMode,NMode)
       allocate_TabUint8(MultiDegrees, NMode)
       allocate_TabUint8(Tester,NMode)
//Initialisation of the degrees of the maximal mode
       int Deg=0;
       long IndexEx;      
       size_t SizeBit=NMode;
//
       uint32_t NXDualHTruncPos=0;
//
       double *Contrib=NULL; 
       Contrib=new double [NXDualHPlus]; 
       uint32_t *LFFNUM=NULL; 
       LFFNUM=new uint32_t [NXDualHPlus+1];
//"" 

// 
       double ContribLFF;//Contribution of force constants summed up in absolute value
//
      for (uint32_t ii = 0; ii < NXDualHPlus; ii++)
      {
      Permuter[ii]=ii;
      Contrib[ii]=0;   
      }
//""Starting points of element number xx
     for (uint32_t xx = 0; xx <= NXDualHPlus; xx++)
      {
       LFFNUM[xx]=LFF.Num[xx];
      }
//
//
//
 QuickSortMode(DualHPos,0,NXDualHPlus-1, Permuter, SizeBit);
//Non coupled surfaces first for Local Force Constants
for (int mm=0; mm<NMode; mm++)
 { 
  InitTabInt(Tester, NMode)
 for (int kk=0; kk<DegrePol; kk++)
   { 
   if(kk!=1 || IncK2) //Avoid Harmonic
       {
    ContribLFF=SIGN<double>(KFC.KijNCpld[kk][mm])*KFC.KijNCpld[kk][mm];
       }
    else
       {
     ContribLFF=0;
       }
    if(ContribLFF > ThrPES)
      {
      Deg=kk+1;
      while(Deg >= 0 )
       {
      Tester[mm]=Deg;
      IndexEx=QsearchMode(Tester, DualHPos,Permuter,NXDualHPlus, SizeBit);
      if(IndexEx>=0)
        {
//""       LFF[IndexEx].Idx[LFF[IndexEx].Num]=-((NMode)*(kk)+mm);
       LFF.Idx[LFFNUM[IndexEx]]=-((NMode)*(kk)+mm);
//""       LFF[IndexEx].Num++; Starting incremented up to LFF.Num[IndexEx+1];
       LFFNUM[IndexEx]++;
       Contrib[IndexEx]+=ContribLFF;
        }
       Deg-=2;
       }
      } 
    }
  }
//Coupled surfaces (Cpld> 1) for Local Force Constants
for (int kk=0; kk<NPES; kk++)
 {           
   VoisinCpld (kk, NMode, KFC, NVMode, ModeV);
   ContribLFF=SIGN<double>(KFC.KijCpld[kk])*KFC.KijCpld[kk];
if(ContribLFF > ThrPES)
  {
  //Mandatory intialisation of multi-degrees. 
    InitTabInt(MultiDegrees, NMode)
do {
//     
     for (int mm = 0 ; mm < NMode; mm++)
     {
    Tester[mm]=ModeV[MultiDegrees[mm]][mm];
     }
//
      IndexEx=QsearchMode(Tester,DualHPos,Permuter,NXDualHPlus, SizeBit);      
      if(IndexEx>=0)
       {
//""   LFF[IndexEx].Idx[LFF[IndexEx].Num]=kk;
//""   LFF[IndexEx].Num++;
       LFF.Idx[LFFNUM[IndexEx]]=kk;
       LFFNUM[IndexEx]++;//Starting incremented up to LFF.Num[IndexEx+1];
       Contrib[IndexEx]+=ContribLFF;
       }
//
    } while(nested_loop(MultiDegrees, NVMode, NMode));
   }
 }
uint64_t NumRot;
//Local Coriolis Interaction
if(DoRot)
 {
 for (int ni=0;ni<NMode;ni++)
  {
 for (int nj=ni+1;nj<NMode;nj++)
    {
//Avoid double count ijkl
  for (int nk=0;nk<NMode;nk++)
     {
    for (int nl=nk+1;nl<NMode;nl++)
      {
       ContribLFF=ZetaXYZ[ni + ((nj-1)*nj)/2][nk + ((nl-1)*nl)/2];
       ContribLFF*=SIGN<double>(ContribLFF);
    if(ContribLFF > ThrPES)
       {
       NumRot=AssignRot(Tester,  NMode,  ni, nj, nk, nl);
//Tester is the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical vectors) 
       VoisinCpldRot (NMode, Tester, NVMode, ModeV); 
       InitTabInt(MultiDegrees, NMode)
      do
        {     
     for (int mm = 0 ; mm < NMode; mm++)
          {
    Tester[mm]=ModeV[MultiDegrees[mm]][mm];
          }
      IndexEx=QsearchMode(Tester,DualHPos,Permuter,NXDualHPlus, SizeBit);      
      if(IndexEx>=0)
          {
//""       LFF[IndexEx].Idx[LFF[IndexEx].Num]=NPES+NumRot;
//""       LFF[IndexEx].Num++;
       LFF.Idx[LFFNUM[IndexEx]]=NPES+NumRot;
       LFFNUM[IndexEx]++;//Starting incremented up to LFF.Num[IndexEx+1];
       Contrib[IndexEx]+=ContribLFF;
          }
        } while(nested_loop(MultiDegrees, NVMode, NMode));
       }
      }
     }
    }
   } 
  }
//
Corres[0]=0; //First Index for 0 excitation always here (convention)
NXDualHTruncPos++;
//Corres will store the indexes of Sum|FC|>ThrKX
for (int unsigned xx=1; xx<NXDualHPlus; xx++)
 { 
  if(Contrib[xx]>ThrKX)
   {
   Corres[NXDualHTruncPos]=xx;
   NXDualHTruncPos++;
   }
 }
//
//""
         delete [] LFFNUM;
         delete [] NVMode;
         delete [] MultiDegrees;
         delete [] Tester;
         delete [] Contrib;
//
         FreeMat(ModeV, DegrePol+1)
//
         return NXDualHTruncPos;
//
}
