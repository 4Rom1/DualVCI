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
  

--------Additional source-------------------

The PrintBanner routine is taken from
www.cise.ufl.edu/class/cop4600/minix/src/commands/simple/banner.c
(By B.Wallis, 4 July 1988)
*/
//
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Shared.h"
#include "ReadData.h"
#include "Basic.h"
#include <cstring>
#include "stdint.h"
//
int CmptNonChar(char *Tab, char NC) 
{
/*
 Routine returning the number of counted elements
 after meeting the first character NC in Tab*/
//
  int mon_i=0;
  while((mon_i<MaxChar)) {
  if (Tab[mon_i] == NC || (Tab[mon_i]=='\0' || Tab[mon_i]=='\n'))
    { 
    return mon_i+1;
    }
     mon_i++;
   }
//
   return mon_i+1;
}
//
int CmptChar(char *Tab, char C, int limit) {
/*Routine returning the number of counted elements
 equal to character C in Tab*/
   int mon_i=0;
   int cmnt=0;
  while(mon_i<limit && !(Tab[mon_i]=='\n')) {
  if (Tab[mon_i] == C)
    { 
    cmnt++;
    }
     mon_i++;
   }
    return cmnt;
}
//
int DetectChar(char *Chaine, char C, int limit)
{
/*Routine returning 1 if character C is in Tab, 0 else.
Try to detect from the beginning of the line until limit=MinChar*/  
  int mon_i=0;
  do
  {
  if (Chaine[mon_i] == C)
    { 
    return 1;
    }
     mon_i++;
  }while(((mon_i<limit) || (Chaine[mon_i]!='\n'))\
 && (Chaine[mon_i-1]==' '));
     return 0;
}
//
int SubStringChar (const char *str, char *new_s, char C)
{//Read a string until C and copy it into news
   int jj=0;
         while ((str[jj] == ' '))
         { 
            jj++;
         };
   int ii=0;  
   if (str != NULL)
   {        
        do
         {  new_s[ii] = str[ii+jj];
            ii++;
         }while (str[ii+jj] != C);
         new_s[ii] = '\0';     
    }  
    return ii+jj;
}
//
int CmptEqKey(int *Tab, int Key, int taille) {
/*Routine returning the number of counted elements
 equal to integer Key in Tab until taille*/
   int mon_i=0;
   int NKey=0;
   
   for (mon_i=0; mon_i<taille;mon_i++)
    { 
    if((Tab[mon_i]==Key))
     {NKey++;}
    }
    return NKey;
}
//
int countlines(FILE *file)
{  
/*Routine counting the number of lines in file*/  
  char chaine[MaxChar]={0};
  int count=0;    
  rewind ( file ); 
  while ( fgets ( chaine , MaxChar , file ) != NULL )
  {
    if (!( (DetectChar(chaine, '/',MinChar)) || (DetectChar(chaine, '@',MinChar)) ))
    { count++; }  
  }  
   return count;
}
//
int FindComma (const char *str, int *PoCom)
{
/*Routine returning the number of comma in str,
and store positions of commas in PoCom*/
   int ii=0;  
   int NCom=1;
   if (str != NULL)
   {         
        do
         {  if(str[ii] == ',')
            {PoCom[NCom]=ii;
            NCom++;}
            ii++;
         }while (str[ii] != '\n');       
    }  
    return NCom;
}
//
int ShiftChar(char *Tab, char Spec) {
/*
 Routine returning the number of counted elements
 after meeting the first character Spec in Tab.
 Return -1 is no char Spec has been found before MaxChar.*/
   int mon_i=0;   
    while ((Tab[mon_i] != Spec) && mon_i < MaxChar )
    { 
    mon_i++;
    }
    if(mon_i==MaxChar)
    {return -1;}
    return mon_i+1;
}
//
int CmptBlank(char *Tab, int taille)
{
/*Routine returning the number of counted blank in Tab*/
  int mon_i=0;
  while ((mon_i < taille) && (Tab[mon_i] == ' '))
    { 
    mon_i++;
    }
    return mon_i;
}
//
int GetKeyWords(FILE *file1, int *NMode, int *DoRot, int *MaxEv, int *DeltaNev, int *AddTarget,int *PESType,\
    uint8_t TargetState[MaxTarget][MaxNormal], int *DoGraph, int *MAXNCV,\
    int *MaxQLevel, int *NAdd, int *MaxAdd, int *PrintOut,  double *EpsRez, double *KNREZ,\
    double *KNNZ, double *KNZREZ, double *EtaComp, double *Tol,double *Kappa, double *ThrMat,double *ThrPES, double *ThrCoor,\
    double *GroundState, double* MinFreq, double *MaxFreq, double *Freq0Max, float *ThrKX, double *Memory,\
    int *Verbose, char *OutName, char *PESName, char* RefName)
//
{
//
//Read the different key words as written in the manual. Case insensitive.
//
    char chaine[MaxChar]; 
    char chainedat[MaxChar];
    int BI; //Blank intervals
    int curs; //Position of the cursor
    int NTarget=0; //Number of target states
    int Pair[2]={0}; //Non zero entry for configs for exemple 1(2)
    char Pthesis='('; 
    int PoCom[MaxChar]={0}; //PoCom=integer localizing position of comma
    int NCom=0; //Number of commas
    int DetNC=0;
//
//Case Insensitive
//
    rewind ( file1 );
//    
    while (fgets(chaine, MaxChar, file1) != NULL) 
      {
//
        curs=SubString(chaine,chainedat);
        if (strcasecmp(chainedat,"NMode")==0)
        {sscanf(&chaine[curs], "%d", NMode);
         InitMat(TargetState, MaxTarget, *NMode)
         DetNC=1;
        }
        else if (strcasecmp(chainedat,"DoRot")==0)
        {sscanf(&chaine[curs], "%d", DoRot);}
        else if (strcasecmp(chainedat,"Verbose")==0)
        {sscanf(&chaine[curs], "%d", Verbose);}
        else if (strcasecmp(chainedat,"MaxEV")==0)
        {sscanf(&chaine[curs], "%d", MaxEv);}
        else if (strcasecmp(chainedat,"DeltaNev")==0)
        {sscanf(&chaine[curs], "%d", DeltaNev);}
        else if (strcasecmp(chainedat,"DoGraph")==0)
        {sscanf(&chaine[curs], "%d", DoGraph);}
        else if (strcasecmp(chainedat,"MAXNCV")==0)
        {sscanf(&chaine[curs], "%d", MAXNCV);}
        else if (strcasecmp(chainedat,"AddTarget")==0)
        {sscanf(&chaine[curs], "%d", AddTarget);}
        else if (strcasecmp(chainedat,"PESType")==0)
        {sscanf(&chaine[curs], "%d", PESType);}
        else if (strcasecmp(chainedat,"PESName")==0)
        {sscanf(&chaine[curs], "%s", PESName);}
        else if (strcasecmp(chainedat,"RefName")==0)
        {sscanf(&chaine[curs], "%s", RefName);}         
        else if (strcasecmp(chainedat,"PrintOut")==0)
        {sscanf(&chaine[curs], "%d", PrintOut);}
        else if (strcasecmp(chainedat,"MaxQLevel")==0)
        {sscanf(&chaine[curs], "%d", MaxQLevel);}
        else if (strcasecmp(chainedat,"NAdd")==0)
        {sscanf(&chaine[curs], "%d", NAdd);}
        else if (strcasecmp(chainedat,"MaxAdd")==0)
        {sscanf(&chaine[curs], "%d", MaxAdd);}
        else if (strcasecmp(chainedat,"Tol")==0)
        {sscanf(&chaine[curs], "%lf", Tol);}
        else if (strcasecmp(chainedat,"ThrKX")==0)
        {sscanf(&chaine[curs], "%f", ThrKX);}
        else if (strcasecmp(chainedat,"Kappa")==0)
        {sscanf(&chaine[curs], "%lf", Kappa);}
        else if (strcasecmp(chainedat,"KNREZ")==0)
        {sscanf(&chaine[curs], "%lf", KNREZ);}
        else if (strcasecmp(chainedat,"KNZREZ")==0)
        {sscanf(&chaine[curs], "%lf", KNZREZ);}
        else if (strcasecmp(chainedat,"KNNZ")==0)
        {sscanf(&chaine[curs], "%lf", KNNZ);}
        else if (strcasecmp(chainedat,"ThrMat")==0)
        {sscanf(&chaine[curs], "%lf", ThrMat);}
        else if (strcasecmp(chainedat,"ThrPES")==0)
        {sscanf(&chaine[curs], "%lf", ThrPES);}
        else if (strcasecmp(chainedat,"ThrCoor")==0)
        {sscanf(&chaine[curs], "%lf", ThrCoor);}
        else if (strcasecmp(chainedat,"MinFreq")==0)
        {sscanf(&chaine[curs], "%lf", MinFreq);}
        else if (strcasecmp(chainedat,"MaxFreq")==0)
        {sscanf(&chaine[curs], "%lf", MaxFreq);}
        else if (strcasecmp(chainedat,"Freq0Max")==0)
        {sscanf(&chaine[curs], "%lf", Freq0Max);}
        else if (strcasecmp(chainedat,"EtaComp")==0)
        {sscanf(&chaine[curs], "%lf", EtaComp);}
        else if (strcasecmp(chainedat,"EpsRez")==0)
        {sscanf(&chaine[curs], "%lf", EpsRez);}
        else if (strcasecmp(chainedat,"GroundState")==0)
        {sscanf(&chaine[curs], "%lf", GroundState);}
        else if (strcasecmp(chainedat,"Memory")==0)
        {sscanf(&chaine[curs], "%lf", Memory);} 
         else if (strcasecmp(chainedat,"OutName")==0)
        {sscanf(&chaine[curs], "%s", OutName);} 
//
      }
       rewind(file1);
//
      while (fgets(chaine, MaxChar, file1) != NULL) 
      {
        curs=SubString(chaine,chainedat);
        PoCom[0]=curs;
         if ((strcasecmp(chainedat,"TargetState")==0) && DetNC==1 )
        {
          BI=CmptBlank(&chaine[curs], MaxChar);
          if (strncasecmp(&chaine[curs+BI],"Fund",4)==0) 
          {        
          for (int cc=1;cc<= *NMode ;cc++)
           {
         TargetState[cc][cc-1]=1;
           }
          NTarget+=*(NMode)+1;
          }   
         else
         {
        NCom=FindComma (chaine, PoCom);
         for (int cc=0;cc<NCom;cc++)
             {
         BI=CmptBlank(&chaine[PoCom[cc]+1], MaxChar); 
         sscanf(&chaine[PoCom[cc]+1+BI], "%d%c%d",&Pair[0],&Pthesis,&Pair[1]);
         TargetState[NTarget][Max<int>(Pair[1]-1,0)]=Pair[0];         
             }
            NTarget++;
           }
         }
       }
          return NTarget;
}
//
void GetKijMonm(FILE *file, double KijCpld[],double **KijNCpld, int32_t **Monm, int NMode, double ThrPES)
{  
/*
This routine read the PES file and store the coupled terms in KijCpld, and the non coupled in KijNCpld.
KijNCpld[kk][mm]= non coupled force constants of degree kk+1 on normal coordinate mm < NMode. 
Monm[kk][mm] degree of monomial term kk on coordinate mm in PES.
ThrPES : threshold for acceptable force constants.
*/
    char chaine[MaxChar]={0};
    char subchaine[MaxChar]={0}; 
    int cnt=0;
    int curs1=0; //Position of the cursor
    int curs2=0,curs3=0,curs4=0;
    int Index[MaxDegrePI]={0};    
    int NCom=0;
    int NKeys=0;
    double TempVal=0;
    int Go=1;
    int cursf=0;
    int CtDeg=0;
    char Comma=',';
//
//
    while ( (fgets ( chaine , MaxChar , file ) != NULL) )
     {
        if (strncasecmp(chaine,"ENDFF",5)==0)
           {Go=0;}   
        if (!( (DetectChar(chaine, '/',MinChar)) || (DetectChar(chaine, '@',MinChar)) ) && Go )
      {
       curs1=ShiftChar(chaine,'[');
//
       curs2=SubString (&chaine[curs1], subchaine);
//       
       cursf=ShiftChar(chaine,']');
//       
       NCom=CmptChar(chaine,Comma,cursf);
       curs3=0;
       for (int ii=0;ii<NCom;ii++)
       {
        curs4=curs3;
        sscanf(&subchaine[curs4], "%d", &Index[ii]);
        curs3+=CmptNonChar(&subchaine[curs4],Comma);
       }
      sscanf(&chaine[curs1+curs2+2], "%lf", &TempVal);
//      
      if(SIGN<double>(TempVal)*TempVal > ThrPES)
      {
       int OKCpld=1;
       CtDeg=0;
       for (int ii=0;ii<NMode;ii++)
        {
          NKeys=CmptEqKey(Index,ii,NCom);
          if( (NKeys>0) &&  (NKeys<NCom) )
          {
           Monm[cnt][ii]=NKeys;
           CtDeg+=NKeys;
          }   
          else if(NKeys==NCom)
          {
           KijNCpld[NKeys-1][ii]=TempVal;  
           CtDeg=NKeys;
           OKCpld=0;
          }
        }
//            
        if(OKCpld){
                  KijCpld[cnt]=TempVal;
                  cnt++;
                  }    
       }
      }
     }
//    
 printf("\nLast derivative readed %e \n",TempVal); 
}
//
long GetNumKijCpld(FILE *file, int NMode, double ThrPES, int *NCPol, int *DegrePol, int DegreCoupl[MaxCpld])
{    
//Get the number of derivatives in PES file for PESType=1  
//Maximal degrees per couplings are saved in DegreCoupl
//DegreCoupl[0] corresponds to 1 mode coupling terms.
//Maximal degrees and couplings are respectively returned in  DegrePol and NCPol
//ThrPES : threshold for acceptable force constants.
    char chaine[MaxChar]={0};
    char subchaine[MaxChar]={0}; 
    int cnt=0;
    int curs1=0; //Position of the cursor
    int curs2=0,curs3=0,curs4=0;
    int Index[MaxDegrePI]={0};    
    int NCom=0;
    int NKeys=0;
    double TempVal=0;
    int Go=0;
    int OKCpld=0;
    int cursf=0;
    char Comma=',';
    int CtCpld=0;
    int CtDeg=0;
    int MaxDeg=0;
    int NTerms[MaxDegrePI]={0};
    int MaxCoup=0;
//
    while ( (fgets ( chaine , MaxChar , file ) != NULL) )
     {    
        if(Go==1)
        { 
         if (strncasecmp(chaine,"ENDFF",5)==0)
           {Go=0;}
//    
        if (!( (DetectChar(chaine, '/',MinChar)) || (DetectChar(chaine, '@',MinChar)) ) && Go )
        {
       curs1=ShiftChar(chaine,'[');

        if(curs1==-1)
          {printf("\n***Bad format or wrong value of PESType. \n Chaine lue %s \
            ***\n\n",chaine);return -1;}
//
       curs2=SubString (&chaine[curs1], subchaine);
//      
       cursf=ShiftChar(chaine,']');
//       
       NCom=CmptChar(chaine,Comma,cursf);
//      
       curs3=0;
       for (int ii=0;ii<NCom;ii++)
        {
        curs4=curs3;
        sscanf(&subchaine[curs4], "%d", &Index[ii]);
        curs3+=CmptNonChar(&subchaine[curs4],Comma);
        }
//
       sscanf(&chaine[curs1+curs2+2], "%lf", &TempVal);
       CtDeg=0;
      if(SIGN<double>(TempVal)*TempVal > ThrPES)
        {
        CtCpld=0;
        OKCpld=1;
       for (int ii=0;ii<NMode;ii++)
         {
          NKeys=CmptEqKey(Index,ii,NCom);
          CtDeg+=NKeys;
          //Couplings
          if(NKeys){CtCpld++;}
//
          if(NKeys==NCom)
          {
             OKCpld=0;
          }
         }
//
         if(CtDeg>MaxDeg){MaxDeg=CtDeg;}
         NTerms[CtDeg]++; 
          //Maximum couplings should equal NCPol       
         if(CtCpld>MaxCoup){MaxCoup=CtCpld;}
//         
         if(MaxCoup>MaxCpld)
          {
          printf("***Coupling polynomials only limited to %d***\n",MaxCpld);
          return 0;
          }
         //Maximum degrees per couplings
         if(DegreCoupl[CtCpld-1]<CtDeg){DegreCoupl[CtCpld-1]=CtDeg;}
         if(OKCpld)
          {
          cnt++;
          }              
        }

      }
    }
        if(strncasecmp(chaine,"FORCEFIELD",10)==0)
        {
        printf("\n***Force Field detected***\n");
        Go=1;
        }
  }
  rewind ( file );  
//
  *NCPol=MaxCoup;
  *DegrePol=MaxDeg;
//
    printf("\n Derivative orders details \n"); 
//  
    for (int kk=1;kk<MaxDegrePI;kk++)
        {
        if(NTerms[kk]!=0)
          {
           printf("%d Derivatives of order %d \n",NTerms[kk],kk);
          }
        }
  printf("Total number of terms = %d \n",cnt);
  return cnt;    
}
//
double AdjustCoeff(double **KijNCpld, double *KijCpld,\
int DegrePol, int NPES, int PESType, int32_t **Monm, int NMode)
{
/*
Routine that transformes all the force constant or derivative in cm-1 depending on the value of PEStype.
NPES : Number of coupled PES coefficients stored in KijCpld.
Monm[kk][mm] : degree of monomial term kk of coordinate mm in PES. 
KijNCpld[kk][mm] : Non coupled force constants of degree kk+1 on normal coordinate mm. 
For PESType = 0 the force constants are already in cm-1.
For PESType > 0 the derivatives are in atomic units.
The conversion in cm-1 is done by setting q_i=Q_i / Omega[i], Omega[i]=sqrt(\nu_i).
*/
int ii=0;
int jj=0;
int kk=0;
int mm=0;
int FF=0;
//int index;
double KFact=1;
double *Omega=NULL;
Omega=new double[NMode];
double KFMax=0;
//
 if (!PESType)//Force constants are already in cm-1
 {
 for (ii=0; ii<  NMode ; ii ++) 
   {
   KijNCpld[1][ii]=2*KijNCpld[1][ii];
   }
 }
//
 if (PESType)//Derivatives in a.u
 {
 for (jj=0; jj<  NMode ; jj ++) {
     Omega[jj]=pow(KijNCpld[1][jj],1.0/4.0);
   }
 for (jj=0; jj <  NMode ; jj++) 
  {     
	    KijNCpld[0][jj]=(KijNCpld[0][jj]*HA_TO_CM)/(Omega[jj]);
            KijNCpld[1][jj]=pow(Omega[jj],2)*HA_TO_CM;           
     for (ii=2; ii<  DegrePol ; ii ++) 
    {
  KFact=factorial(ii+1);
  KijNCpld[ii][jj]=(KijNCpld[ii][jj]*HA_TO_CM)/(KFact*pow(Omega[jj],ii+1));                                      
   }
 }
//   
   for (kk=0; kk <  NPES ; kk++) 
    {    
     for (mm=0; mm < NMode; mm++)       
      {
      FF=Monm[kk][mm];
      KFact=factorial(FF);
      KijCpld[kk]/=(KFact*pow(Omega[mm],FF));
      }
      KijCpld[kk]*=HA_TO_CM;  
     }
  }
//
  for (kk=0; kk <  NPES ; kk++) 
    { 
    if(KFMax < SIGN<double>(KijCpld[kk])*KijCpld[kk])
     {KFMax=SIGN<double>(KijCpld[kk])*KijCpld[kk];}
    }
//
 for (jj=0; jj <  NMode ; jj++) 
  {
    for (ii=0; ii<  DegrePol ; ii ++) 
    {
     if(ii-1)
     {     
      if(KFMax < SIGN<double>(KijNCpld[ii][jj])*KijNCpld[ii][jj])
      {KFMax=SIGN<double>(KijNCpld[ii][jj])*KijNCpld[ii][jj];}     
     }
    }
 }
//
delete[] Omega;
return KFMax;
}
//
int SubString (const char *str, char *new_s)
{
//Save str in new_s starting from non blank value
//until the next blank value.
//Returns the sum of first blank values plus the size of new_s
   int jj=0;
         while ((str[jj] == ' '))
         { 
            jj++;
         };
   int ii=0;  
   if (str != NULL)
   {        
        do
         {  new_s[ii] = str[ii+jj];
            ii++;
         }while (str[ii+jj] != ' ');
         new_s[ii] = '\0';     
    }  
    return ii+jj;
}


int ReadCoor(FILE *file, double *CoorCart,double *Mas, double **CoorMode, int NAtome, int NMode, int Verbose){
//
//Read cartesian coordinates of equilibrium geometry in Bohr, masses in me (electron rest mass) and normal coordinate vectors.
//The equilibrium geometry in cartesian coordinates reads:  
//         X                       Y                      Z
//CoorCart[3*(aa)],  CoorCart[3*(aa)+1],  CoorCart[3*(aa)+2], for aa (0<=aa<NAtome) 
//
//The normal coordinate vectors read :  
//         X                       Y                      Z
//CoorMode[mm][3*(aa)],  CoorMode[mm][3*(aa)+1],CoorMode[mm][3*(aa)+2], 
//For aa (0<=aa<NAtome) , 0<=mm<NMode
//
char chaine[MaxChar];
//
int CntLine=0;
//
int CntExclam=0; //Each section is delimited by an exclamation mark
//
int StartNC; //First line index per normal coordinate vector
//
int Go=1;
//
while ( (fgets ( chaine , MaxChar , file ) != NULL) && Go )
      {
      
      if(!( (DetectChar(chaine, '/',MinChar)) || (DetectChar(chaine, '@',MinChar)) )){//Used for comment lines

       if((DetectChar(chaine, '!',MinChar)))
         {
          if(Verbose)
          {
          printf("\n%s",chaine);
          }         
          CntExclam++;
         }
      else if(strncasecmp(chaine,"ENDCOOR",7)==0)
        {Go=0;}
      else
    {
if(CntExclam==1)
   {
   sscanf(chaine, "%lf %lf %lf", &CoorCart[3*(CntLine)],  &CoorCart[3*(CntLine)+1],  &CoorCart[3*(CntLine)+2] );
    if(Verbose)
     {
    printf("%s",chaine);
     }
   CntLine++;
    }
   else if(CntExclam==2)
    {
    sscanf(chaine, "%lf", &Mas[CntLine-NAtome] );
    if(Verbose)
     {
    printf("%s",chaine);
     }
    CntLine++;
    }
    else if(CntExclam>2)
    {
   StartNC=(CntExclam-1)*NAtome;
   sscanf(chaine, "%lf %lf %lf", &CoorMode[(CntExclam-3)][3*(CntLine-StartNC)],  &CoorMode[(CntExclam-3)][3*(CntLine-StartNC)+1],\
   &CoorMode[(CntExclam-3)][3*(CntLine-StartNC)+2] );
          if(Verbose)
          {
   printf("%e %e %e\n", CoorMode[(CntExclam-3)][3*(CntLine-StartNC)],  CoorMode[(CntExclam-3)][3*(CntLine-StartNC)+1],\
   CoorMode[(CntExclam-3)][3*(CntLine-StartNC)+2] );
          }
   CntLine++;
    }
   }  
  }
 }
     if(CntExclam==0)
      {
     printf("\n****Problem reading coordinates-> exclamation mark should precede a field****\n\n");
     return -1;
      }
//
     if(CntExclam-2 != NMode)
      {
       printf("\n****Problem reading coordinates-> exclamation mark should precede a field****\n\n");
       return -1;
      }
//
  return 0;
}
//
int GetDataPES(FILE *file, double KijCpld[],double **KijNCpld, int32_t **Monm, double *CoorCart,double *Mas, double **CoorMode, int NMode, int NAtome, int PESType, int Verbose, double ThrPES)
{   
//Detect Forcefield or derivatives in section FORCEFIELD of PES file. 
//Also detect equilibrium geometry and normal coordinates in section COORDINATE (if non void)
//Monm[kk][mm] degree du terme kk coordonnee mm in PES, 
//Everything is transformed in cm-1 depending on the value of PEStype.
//For PESType = 0 the force constants are already in cm-1.
//For PESType > 0 the derivatives are in atomic units.
//
//KijNCpld[kk][mm]=Non coupled force constants of degree kk+1 on normal coordinate mm. 
//
    char chaine[MaxChar]={0};
//
    int DoRot=0;       
//  
    int Check=1;
//
    while ( (fgets ( chaine , MaxChar , file ) != NULL) )
     {     
     if((strncasecmp(chaine,"COORDINATE",10)==0) && Check)
      {
      if(Verbose)
       {
      printf("\n **** Coordinates detected **** \n");
       }
      DoRot=1;
      Check=ReadCoor(file, CoorCart, Mas, CoorMode, NAtome, NMode, Verbose);
       if(Check<0){return -1;}
//
        rewind(file);     
      }
     else if(strncasecmp(chaine,"FORCEFIELD",10)==0)
      {
      if(Verbose)
       {printf("\n** Force Field detected for the second time **\n");}
      if(PESType)
       {
       GetKijMonm(file, KijCpld, KijNCpld, Monm, NMode, ThrPES);
       }
      else
       {
       GetKijMonm0(file, KijCpld, KijNCpld, Monm, NMode, ThrPES);
       }
      }
     }
   return DoRot;
}

void afficheTabDbFile(double *Tab, int lin, FILE *FF){
for (int i = 0 ; i < lin ; i++)
 {
 fprintf(FF,"%.15f \n", Tab[i]);
 } 
}

void GetKijMonm0(FILE *file, double KijCpld[],double **KijNCpld, int32_t **Monm, int NMode, double ThrPES)
{ 
/*   
This routine read the PES file and store the coupled terms in KijCpld, and the non coupled in KijNCpld.
Monm[kk][mm] degree of monomial term kk on normal coordinate mm in PES, 
KijNCpld[kk][mm]= non coupled force constants of degree kk+1 on normal coordinate mm < NMode.
ThrPES : threshold for acceptable force constants.*/
//
    char chaine[MaxChar]={0};
    int cnt=0;
    int curs1=0; //Position of the cursor
    int curs2=0,curs3=0,curs4=0;
    int Index[MaxDegrePI]={0};    
    int NKeys=0;
    double TempVal=0;
    int Go=1;
    int CtDeg=0;
    char Comma=',';
    char Blank=' ';
//
    while ( (fgets ( chaine , MaxChar , file ) != NULL) )
     {
        if (strncasecmp(chaine,"ENDFF",5)==0)
           {Go=0;}
//    
       if (!( (DetectChar(chaine, '/',MinChar)) || (DetectChar(chaine, '@',MinChar)) ) && Go )
      {
       curs1=CmptBlank(chaine,MaxChar);
       curs2=CmptNonChar (&chaine[curs1], Comma);
//       
       curs3=curs1;
//
       for (int ii=0;ii<NMode;ii++)
       {
        curs4=curs3;
        sscanf(&chaine[curs4], "%d", &Index[ii]);
        curs3+=CmptNonChar(&chaine[curs4],Blank);
        curs4=curs3;
        curs3+=CmptBlank(&chaine[curs4],MaxChar);       
       }
       sscanf(&chaine[curs1+curs2], "%lf", &TempVal);
      if(SIGN<double>(TempVal)*TempVal > ThrPES)
      {
       int OKCpld=1;
       CtDeg=0;
       NKeys=CmptEqKey(Index,0,NMode);
       for (int ii=0;ii<NMode;ii++)
        {
          if(NKeys<NMode-1)
          {
           Monm[cnt][ii]=Index[ii];
           CtDeg+=Index[ii];
          }   
          else if(NKeys==NMode-1)
          {
           if(Index[ii])
           {
           KijNCpld[Index[ii]-1][ii]=TempVal;  
           CtDeg=Index[ii];
           OKCpld=0;
           }
          }
        }
//                
        if(OKCpld){
                  KijCpld[cnt]=TempVal;
                  cnt++;
                  }    
       }
      }
     }
//        
 printf("\nLast force constant readed %e \n",TempVal); 
}
//
//
long GetNumKijCpld0(FILE *file, int NMode, double ThrPES, int *NCPol, int *DegrePol, int DegreCoupl[MaxCpld])
{
//Get the number of force constants in PES file for PESType=0  
//Maximal degrees per couplings are saved in DegreCoupl
//DegreCoupl[0] corresponds to 1 mode coupling terms.
//Maximal degrees and couplings are respectively returned in DegrePol and NCPol
//ThrPES : threshold for acceptable force constants.
    char chaine[MaxChar]={0};
    int cnt=0;
    int curs1=0; //Position of the cursor
    int curs2=0;
    int curs3=0; //Position of the cursor
    int curs4=0;
    int Index[MaxDegrePI]={0};    
    int NKeys=0;
    double TempVal=0;
    int Go=1;
    int NTerms[MaxCpld]={0};
    int CtDeg=0;
    char Comma=',';
    char Blank=' ';
    int CtCpld=0;
    int MaxDeg=0;
    int MaxCoup=0;
//
    while ( (fgets ( &chaine[0] , MaxChar , file ) != NULL) )
    {
//     
    if(Go==1)
     { 
      if (strncasecmp(chaine,"ENDFF",5)==0)
      {Go=0;}  
      if (!( (DetectChar(chaine,'/',MinChar)) || (DetectChar(chaine, '@',MinChar)) ) && Go )
      {
//        
       if(CmptChar(chaine, Comma,MaxChar)>1)
        {printf("\n***Bad format or wrong value of PESType \n Chaine lue %s \
            ***\n\n",chaine);return -1;}   
//    
       curs1=CmptBlank(chaine,MaxChar);
       curs2=CmptNonChar (&chaine[curs1], Comma);       
       curs3=curs1;
//
       for (int ii=0;ii<NMode;ii++)
       {
        curs4=curs3;
        sscanf(&chaine[curs4], "%d", &Index[ii]);
        curs3+=CmptNonChar(&chaine[curs4],Blank);
        curs4=curs3;
        curs3+=CmptBlank(&chaine[curs4],MaxChar);       
       }
       sscanf(&chaine[curs1+curs2], "%lf", &TempVal);
      if(SIGN<double>(TempVal)*TempVal > ThrPES)
      {
       int OKCpld=1;
       CtDeg=0;
       NKeys=CmptEqKey(Index,0,NMode);
      //Count max couplings
       CtCpld=(NMode-NKeys);
       if(MaxCoup<CtCpld){MaxCoup=CtCpld;}
//   
       if(MaxCoup>MaxCpld)
       {
       printf("***Coupling polynomials only limited to %d***\n",MaxCpld);
       return 0;
       }
//
       for (int ii=0;ii<NMode;ii++)
        {
          if(NKeys<NMode-1)
          {
           CtDeg+=Index[ii];
          }   
          else if(NKeys==NMode-1)
          {
           CtDeg+=Index[ii];
           OKCpld=0;
          }
        }
//        
         if(CtDeg>MaxDeg)
         {MaxDeg=CtDeg;}
//
         if(DegreCoupl[CtCpld-1]<CtDeg)
         {DegreCoupl[CtCpld-1]=CtDeg;}
//             
         NTerms[CtCpld-1]++; 
//       
        if(OKCpld){
                  cnt++;
                  }    
       }
      }
     }
//    
        if(strncasecmp(chaine,"FORCEFIELD",10)==0)
      {printf("\n****** Force Field detected ******\n");
        Go=1;}
  }
//  
  rewind ( file ); 
//  
  *NCPol=MaxCoup;
  *DegrePol=MaxDeg;
//
     printf("\n Force field details \n");    
//      
    for (int kk=1;kk<MaxCpld;kk++)
        {
        if(NTerms[kk-1]!=0)
          {
           printf("Total number of %d-mode couplings term = %d \n",kk,NTerms[kk-1]);
          }
        }
//
  printf("Total number of terms = %d \n",cnt);
  return cnt;
}
//
int PrintBanner (const char *Chaine)
{
//Taken from
//www.cise.ufl.edu/class/cop4600/minix/src/commands/simple/banner.c
//By B.Wallis, 4 July 1988
char glyyphs[150][150] = {
	  "         @@@  @@   @@  @ @   @@@@@          @@     @@@  ",
	  "         @@@  @@   @@  @ @  @  @  @@@   @  @  @    @@@  ",
	  "         @@@   @   @ @@@@@@@@  @   @@  @    @@      @   ",
	  "          @            @ @   @@@@@    @    @@@     @    ",
	  "                     @@@@@@@   @  @  @    @   @ @       ",
	  "         @@@           @ @  @  @  @ @  @@ @    @        ",
	  "         @@@           @ @   @@@@@ @   @@  @@@@ @       ",

	  "   @@    @@                                            @",
	  "  @        @   @   @    @                             @ ",
	  " @          @   @ @     @                            @  ",
	  " @          @ @@@@@@@ @@@@@   @@@   @@@@@           @   ",
	  " @          @   @ @     @     @@@                  @    ",
	  "  @        @   @   @    @      @            @@@   @     ",
	  "   @@    @@                   @             @@@  @      ",

	  "  @@@     @    @@@@@  @@@@@ @      @@@@@@@ @@@@@ @@@@@@@",
	  " @   @   @@   @     @@     @@    @ @      @     @@    @ ",
	  "@   @ @ @ @         @      @@    @ @      @          @  ",
	  "@  @  @   @    @@@@@  @@@@@ @@@@@@@ @@@@@ @@@@@@    @   ",
	  "@ @   @   @   @            @     @       @@     @  @    ",
	  " @   @    @   @      @     @     @ @     @@     @  @    ",
	  "  @@@   @@@@@ @@@@@@@ @@@@@      @  @@@@@  @@@@@   @    ",

	  " @@@@@  @@@@@          @@@      @           @     @@@@@ ",
	  "@     @@     @  @@@    @@@     @             @   @     @",
	  "@     @@     @  @@@           @     @@@@@     @        @",
	  " @@@@@  @@@@@@         @@@   @                 @     @@ ",
	  "@     @      @         @@@    @     @@@@@     @     @   ",
	  "@     @@     @  @@@     @      @             @          ",
	  " @@@@@  @@@@@   @@@    @        @           @       @   ",

	  " @@@@@    @   @@@@@@  @@@@@ @@@@@@ @@@@@@@@@@@@@@ @@@@@ ",
	  "@     @  @ @  @     @@     @@     @@      @      @     @",
	  "@ @@@ @ @   @ @     @@      @     @@      @      @      ",
	  "@ @ @ @@     @@@@@@@ @      @     @@@@@@  @@@@@  @  @@@@",
	  "@ @@@@ @@@@@@@@     @@      @     @@      @      @     @",
	  "@     @@     @@     @@     @@     @@      @      @     @",
	  " @@@@@ @     @@@@@@@  @@@@@ @@@@@@ @@@@@@@@       @@@@@ ",

	  "@     @  @*@        @@    @ @      @     @@     @@@@@@@@",
	  "@     @   @         @@   @  @      @@   @@@@    @@     @",
	  "@     @   @         @@  @   @      @ @ @ @@ @   @@     @",
	  "@@@@@@@   @         @@@@    @      @  @  @@  @  @@     @",
	  "@     @   @   @     @@  @   @      @     @@   @ @@     @",
	  "@     @   @   @     @@   @  @      @     @@    @@@     @",
	  "@     @  @@@   @@@@@ @    @ @@@@@@@@     @@     @@@@@@@@",

	  "@@@@@@  @@@@@ @@@@@@  @@@@@ @@@@@@@@     @@     @@     @",
	  "@     @@     @@     @@     @   @   @     @@     @@  @  @",
	  "@     @@     @@     @@         @   @     @@     @@  @  @",
	  "@@@@@@ @     @@@@@@@  @@@@@    @   @     @@     @@  @  @",
	  "@      @   @ @@   @        @   @   @     @ @   @ @  @  @",
	  "@      @    @ @    @ @     @   @   @     @  @ @  @  @  @",
	  "@       @@@@ @@     @ @@@@@    @    @@@@@    @    @@ @@ ",

	  "@     @@     @@@@@@@@ @@@@@ @       @@@@@    @          ",
	  " @   @  @   @      @  @      @          @   @ @         ",
	  "  @ @    @ @      @   @       @         @  @   @        ",
	  "   @      @      @    @        @        @               ",
	  "  @ @     @     @     @         @       @               ",
	  " @   @    @    @      @          @      @               ",
	  "@     @   @   @@@@@@@ @@@@@       @ @@@@@        @@@@@@@",

	  "  @@@                                                   ",
	  "  @@@     @@   @@@@@   @@@@  @@@@@  @@@@@@ @@@@@@  @@@@ ",
	  "   @     @  @  @    @ @    @ @    @ @      @      @    @",
	  "    @   @    @ @@@@@  @      @    @ @@@@@  @@@@@  @     ",
	  "        @@@@@@ @    @ @      @    @ @      @      @  @@@",
	  "        @    @ @    @ @    @ @    @ @      @      @    @",
	  "        @    @ @@@@@   @@@@  @@@@@  @@@@@@ @       @@@@ ",

	  "                                                        ",
	  " @    @    @        @ @    @ @      @    @ @    @  @@@@ ",
	  " @    @    @        @ @   @  @      @@  @@ @@   @ @    @",
	  " @@@@@@    @        @ @@@@   @      @ @@ @ @ @  @ @    @",
	  " @    @    @        @ @  @   @      @    @ @  @ @ @    @",
	  " @    @    @   @    @ @   @  @      @    @ @   @@ @    @",
	  " @    @    @    @@@@  @    @ @@@@@@ @    @ @    @  @@@@ ",

	  "                                                        ",
	  " @@@@@   @@@@  @@@@@   @@@@   @@@@@ @    @ @    @ @    @",
	  " @    @ @    @ @    @ @         @   @    @ @    @ @    @",
	  " @    @ @    @ @    @  @@@@     @   @    @ @    @ @    @",
	  " @@@@@  @  @ @ @@@@@       @    @   @    @ @    @ @ @@ @",
	  " @      @   @  @   @  @    @    @   @    @  @  @  @@  @@",
	  " @       @@@ @ @    @  @@@@     @    @@@@    @@   @    @",

	  "                       @@@     @     @@@   @@    @ @ @ @",
	  " @    @  @   @ @@@@@@ @        @        @ @  @  @ @ @ @ ",
	  "  @  @    @ @      @  @        @        @     @@ @ @ @ @",
	  "   @@      @      @  @@                 @@        @ @ @ ",
	  "   @@      @     @    @        @        @        @ @ @ @",
	  "  @  @     @    @     @        @        @         @ @ @ ",
	  " @    @    @   @@@@@@  @@@     @     @@@         @ @ @ @"
};
  int a, b, c, len, ind;
  char line[80];
	len = strlen(Chaine);
	if (len > 10) len = 10;
	for (a = 0; a < 7; a++) {
		for (b = 0; b < len; b++) {
			if ((ind = Chaine[b] - ' ') < 0) ind = 0;
			for (c = 0; c < 7; c++) {
				line[b * 8 + c] = glyyphs[(ind / 8 * 7) + a][(ind % 8 * 7) + c] == '@' ? ind + ' ' : ' ';
			}
			line[b * 8 + 7] = ' ';
		}
		for (b = len * 8 - 1; b >= 0; b--) {
			if (line[b] != ' ') break;
			line[b] = '\0';
		}
		printf("%s\n", line);
	}
	printf("\n");
  
  return(0);
}
//
int GetValRef(FILE *file, double TabRef[MaxRef])
{ 
//Store frequencies extracted from file in TabRef
    char chaine[MaxChar]={0};
    int Cnt=0;
    rewind(file);
    while ( (fgets ( chaine , MaxChar , file ) != NULL) && (Cnt < MaxRef) )
     {
      sscanf(chaine, "%lf", &TabRef[Cnt]);
      Cnt++;
     }
  return Cnt;
}
//
uint32_t GetConfsBin(FILE *FileBasis, ConfigId *&FinalBasis, int NMode, int* NScreen)
{    
//Read the final basis set file previously stored in binary format
//FinalBasis: where the configurations are stored
//NScreen number of Target screened
//FileBasis: name of the file of the final basis set
uint32_t FinalSize; 
//   
int NScreenTMP;
fread (&NScreenTMP, sizeof(int) , 1 , FileBasis); 
//
if(!NScreenTMP)
{
printf("\n*****DVCI should have exited with success*****\n");
return 0;
}
*NScreen=NScreenTMP; 
//read final size
fread (&FinalSize, sizeof(uint32_t) , 1 , FileBasis); 
//Allocate configs
FinalBasis = new ConfigId [FinalSize];
for (uint32_t kk = 0; kk < FinalSize; kk++)
{InitMode(FinalBasis[kk], NMode);}
//
//Copy configuration arrays in binary format
for (uint32_t ll=0;ll<FinalSize;ll++)
 {
fread (&FinalBasis[ll].Degrees[0], sizeof(uint8_t) , NMode , FileBasis);  
 }
return FinalSize;
} 
//
void PrintConfsBin(char *OutBasis, ConfigId *FinalBasis, int NMode, int NScreen, uint32_t FinalSize)
{ 
//Print configurations of final basis set stored in FinalBasis into files with extension name
//OutBasis, 
//Files to write Final basis set (FileBasis) 
     FILE *FileBasis=NULL;
//
     FileBasis=fopen(OutBasis,"wb");
//First the integer giving the number of screened states.
fwrite (&NScreen, sizeof(int) , 1 , FileBasis); 
// 
//Print final size
fwrite (&FinalSize, sizeof(uint32_t) , 1 , FileBasis); 
//Set up the ground state as the first binary eigenvector, no matter if amoung targets or not
//Copy configuration arrays in binary format
for (uint32_t ll=0;ll<FinalSize;ll++)
 {
fwrite (&FinalBasis[ll].Degrees[0], sizeof(uint8_t) , NMode , FileBasis);
//           
 }     
//
    fclose(FileBasis);
}
void PrintVecBin(char *OutVec, double *EigVec, double *EigVal, int NScreen,int *TabScreen, int PosZero, double GroundState, uint32_t FinalSize)
{ 
//OutBasis, binary eigenvector are also printed out into file having extension OutVec.
//Files to write eigenvectors in binary format (FileVec)
     FILE *FileVec=NULL;
//
     FileVec=fopen(OutVec,"wb");
//First the integer giving the number of screened states.
//Set up the ground state as the first binary eigenvector, no matter if amoung targets or not
fwrite (&EigVec[PosZero*(FinalSize)] , sizeof(double) , FinalSize , FileVec);
//Write corresponding eigenvalue in position FinalSize
fwrite (&EigVal[PosZero], sizeof(double) , 1 , FileVec);
//Copy configuration arrays in binary format    
//GroundState vector already copied no needs 2 times
int StartCpy=0;
if(GroundState <= 0)
 {
 StartCpy=1;
 }
//
int PosEig;
for (int pp=StartCpy;pp<NScreen;pp++)
 {
PosEig=TabScreen[pp];
fwrite (&EigVec[PosEig*(FinalSize)], sizeof(double), FinalSize, FileVec);
//Write corresponding eigenvalue in position FinalSize
fwrite (&EigVal[PosEig], sizeof(double) , 1 , FileVec );
 }
//
    fclose(FileVec);
}
void PrintMatCSCBin(char *OutMat, SizeArray *Size, int Iteration, int NEV, int NCV,
 uint32_t SizeInit, double Shift, double tol, double Ull, CSC IJ, double *ValAct)
{
//Print matrix Hb stored in CSC format (IJ,ValAct),
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//IJ.I=new uint32_t [NNZAct];
//IJ.NJ=new uint64_t [DimAct+1];
//Size[Iteration+1].NNZAct is the number of NNZ of active matrix Hb at current iteration.
//Size[Iteration+1].DimAct is the size of active space B at current iteration.
//ValAct is an array of size Size[Iteration+1].NNZAct containing NNZ of active matrix Hb 
//The matrix will be shifted as Hb=Hb-Shift*Id
//Ull Watson : -1/8*(Sum mu_ll) correction to add to groundstate value
//
FILE *FileMat=NULL;
//
FileMat=fopen(OutMat,"wb");
//Integers giving the number of ritz values and vectors
fwrite (&NEV, sizeof(int) , 1 , FileMat); 
fwrite (&NCV, sizeof(int) , 1 , FileMat); 
//Initial subspace size
fwrite (&SizeInit, sizeof(uint32_t) , 1 , FileMat); 
//The shift, tolerance and Watson correction
fwrite (&Shift, sizeof(double) , 1 , FileMat); 
fwrite (&tol, sizeof(double) , 1 , FileMat); 
fwrite (&Ull, sizeof(double) , 1 , FileMat);
//Integers giving the dimensions.
fwrite (&Size[Iteration+1].DimAct, sizeof(uint32_t) , 1 , FileMat); 
fwrite (&Size[Iteration+1].NNZAct, sizeof(uint64_t) , 1 , FileMat); 
//Matrix pointers and entries
fwrite (&IJ.NJ[0], sizeof(uint64_t) , Size[Iteration+1].DimAct+1 , FileMat); 
fwrite (&IJ.I[0], sizeof(uint32_t) , Size[Iteration+1].NNZAct , FileMat); 
fwrite (&ValAct[0], sizeof(double) , Size[Iteration+1].NNZAct , FileMat); 
//
fclose(FileMat);
}
int ReadMatCSCBin(FILE *FileMat,  int *NEV, int *NCV, uint32_t *SizeInit, uint32_t *DimAct, uint64_t *NNZAct,
 double *Shift, double *tol, double *Ull, CSC &IJ, double *&ValAct)
{
//Print matrix Hb stored in CSC format (IJ,ValAct),
//CSC : compressed sparse column format:
//IJ.NJ[jj] = number of the first nnz elmement of column jj
//IJ.I[nn] = line number of nnz element number nn 
//IJ.I=new uint32_t [NNZAct];
//IJ.NJ=new uint64_t [DimAct+1];
//Size[Iteration+1].NNZAct is the number of NNZ of active matrix Hb at current iteration.
//Size[Iteration+1].DimAct is the size of active space B at current iteration.
//ValAct is an array of size Size[Iteration+1].NNZAct containing NNZ of active matrix Hb 
//The matrix will be shifted as Hb=Hb-Shift*Id
//
//Start by dimensions
//Integers giving the number of ritz values and vectors
int NEVTmp,NCVTmp;
fread (&NEVTmp, sizeof(int) , 1 , FileMat); 
//Something bad happened previously
if(!NEVTmp)
{
printf("\n*****DVCI should have exited with success*****\n");
return 0;
}
//
fread (&NCVTmp, sizeof(int) , 1 , FileMat); 
*NEV=NEVTmp;
*NCV=NCVTmp;
//
uint32_t SizeTmp;
fread (&SizeTmp, sizeof(uint32_t) , 1 , FileMat); 
*SizeInit=SizeTmp;
//The shift and tolerance
double ShiftTmp,tolTmp,UllTmp;
fread (&ShiftTmp, sizeof(double) , 1 , FileMat); 
fread (&tolTmp, sizeof(double) , 1 , FileMat); 
fread (&UllTmp, sizeof(double) , 1 , FileMat); 
*Shift=ShiftTmp;
*tol=tolTmp;
*Ull=UllTmp;
//Integers giving the dimensions
uint32_t DimActTmp;
uint64_t NNZActTmp;
fread (&DimActTmp, sizeof(uint32_t) , 1 , FileMat); 
fread (&NNZActTmp, sizeof(uint64_t) , 1 , FileMat); 
*DimAct=DimActTmp;
*NNZAct=NNZActTmp;
//printf("\n *DimAct %u, *NNZAct %lu,\n",*DimAct,*NNZAct);
//Matrix pointers and entries
//Need to allocate first
   IJ.I=NULL;
   IJ.NJ=NULL;
   IJ.I=new uint32_t [*NNZAct];
   IJ.NJ=new uint64_t [*DimAct+1];
   ValAct=new double [*NNZAct];
//IJAct.NJ[jj],IJ.I[nn] in [0,DimAct[ x [0,NNZAct[
fread (&IJ.NJ[0], sizeof(uint64_t) , *DimAct+1 , FileMat); 
fread (&IJ.I[0], sizeof(uint32_t) , *NNZAct , FileMat); 
fread (&ValAct[0], sizeof(double) , *NNZAct , FileMat); 
//
return NEVTmp;
}
