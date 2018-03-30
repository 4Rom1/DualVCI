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
#ifndef ReadData_H
#define ReadData_H
#include "stdint.h"
//
int countlines(FILE *file);
/*Routine counting the number of lines in file*/  
//
int SubString (const char *str, char *new_s);
int ReadCoor(FILE *file, double *CoorCart,double *Mas, double **CoorMode,int NAtome, int NMode, int Verbose);
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
//
int CmptBlank(char *Tab, int taille);
/*Routine returning the number of counted blank in Tab*/
//
int FindComma (const char *str, int *PoCom);
/*Routine returning the number of comma in str,
and store positions of commas in PoCom*/
//
int GetKeyWords(FILE *file1, int *NMode, int *DoRot, int *MaxEv, int *DeltaNev, int *AddTarget,int *PESType,\
    uint8_t TargetState[MaxTarget][MaxNormal], int *DoGraph, int *MAXNCV,\
    int *MaxQLevel, int *NAdd, int *MaxAdd, int *PrintOut,  double *EpsRez, double *KNREZ,\
    double *KNNZ, double *KNZREZ, double *EtaComp, double *Tol,double *Kappa, double *ThrMat,double *ThrPES, double *ThrCoor,\
    double *GroundState, double* MinFreq, double *MaxFreq, double *Freq0Max, float *ThrKX, double *Memory,\
    int *Verbose, char *OutName, char *PESName, char* RefName);
//Read the different key words as written in the manual. Case insensitive.
//
int ShiftChar(char *Tab, char Spec);
/*
 Routine returning the number of counted elements
 after meeting the first character Spec in Tab.
 Return -1 is no char Spec has been found before MaxChar.*/
//
void GetKijMonm(FILE *file, double KijCpld[],double **KijNCpld, int32_t **Monm, int NMode, double ThrPES);
/*
This routine read the PES file and store the coupled terms in KijCpld, and the non coupled in KijNCpld.
KijNCpld[kk][mm]= non coupled force constants of degree kk+1 on normal coordinate mm < NMode. 
Monm[kk][mm] degree of monomial term kk on coordinate mm in PES.
ThrPES : threshold for acceptable force constants.
*/
//
int CmptEqKey(int *Tab, int Key, int taille);
/*Routine returning the number of counted elements
 equal to integer Key in Tab until taille*/
//
int CmptNonChar(char *Tab, char NC);
/*
 Routine returning the number of counted elements
 after meeting the first character NC in Tab*/
//
double AdjustCoeff(double **KijNCpld, double *KijCpld,\
int DegrePol, int NPES, int PESType, int32_t **Monm, int NMode);
/*
Routine that transformes all the force constant or derivative in cm-1 depending on the value of PEStype.
NPES : Number of coupled PES coefficients stored in KijCpld.
Monm[kk][mm] : degree of monomial term kk of coordinate mm in PES. 
KijNCpld[kk][mm] : Non coupled force constants of degree kk+1 on normal coordinate mm. 
For PESType = 0 the force constants are already in cm-1.
For PESType > 0 the derivatives are in atomic units.
The conversion in cm-1 is done by setting q_i=Q_i / Omega[i], Omega[i]=sqrt(\nu_i).
*/
//
int GetDataPES(FILE *file, double KijCpld[],double **KijNCpld, int32_t **Monm, double *CoorCart,double *Mas,\
 double **CoorMode, int NMode, int NAtome, int PESType, int Verbose, double ThrPES);
//Detect Forcefield or derivatives in section FORCEFIELD of PES file. 
//Also detect equilibrium geometry and normal coordinates in section COORDINATE (if non void)
//Monm[kk][mm] degree du terme kk coordonnee mm in PES, 
//Everything is transformed in cm-1 depending on the value of PEStype.
//For PESType = 0 the force constants are already in cm-1.
//For PESType > 0 the derivatives are in atomic units.
//
//KijNCpld[kk][mm]=Non coupled force constants of degree kk+1 on normal coordinate mm. 
//
//
int CmptChar(char *Tab, char C, int limit);
/*Routine returning the number of counted elements
 equal to character C in Tab*/
//
long GetNumKijCpld(FILE *file, int NMode,double ThrPES, int *NCPol, int *DegrePol, int DegreCoupl[MaxCpld]);
//Get the number of derivatives in PES file for PESType=1  
//Maximal degrees per couplings are saved in DegreCoupl
//DegreCoupl[0] corresponds to 1 mode coupling terms.
//Maximal degrees and couplings are respectively returned in  DegrePol and NCPol
//ThrPES : threshold for acceptable force constants.
//
void afficheTabDbFile(double *Tab, int lin, FILE *FF);
int DetectChar(char *Chaine, char C, int limit);
/*Routine returning 1 if character C is in Tab, 0 else.
Try to detect from the beginning of the line until limit=MinChar*/ 
//
int PrintBanner (const char *Chaine);
void GetKijMonm0(FILE *file, double KijCpld[],double **KijNCpld, int32_t **Monm, int NMode, double ThrPES);
/*   
This routine read the PES file and store the coupled terms in KijCpld, and the non coupled in KijNCpld.
Monm[kk][mm] degree of monomial term kk on normal coordinate mm in PES, 
KijNCpld[kk][mm]= non coupled force constants of degree kk+1 on normal coordinate mm < NMode.
ThrPES : threshold for acceptable force constants.*/
//
long GetNumKijCpld0(FILE *file, int NMode, double ThrPES, int *NCPol, int *DegrePol, int DegreCoupl[MaxCpld]);
//Get the number of force constants in PES file for PESType=0  
//Maximal degrees per couplings are saved in DegreCoupl
//DegreCoupl[0] corresponds to 1 mode coupling terms.
//Maximal degrees and couplings are respectively returned in DegrePol and NCPol
//ThrPES : threshold for acceptable force constants.
//
int GetValRef(FILE *file, double *TabRef);
//Store frequencies extracted from file in TabRef
//
unsigned int GetConfs(FILE *file, ConfigId *FinalBasis, int NMode, int* NScreens);
//Read the final basis set file
//FinalBasis: where the configurations are stored
//NScreens number of Target screened
//file:Name of the file of the final basis set
int SubStringChar (const char *str, char *new_s, char C);
//Read a string until C and copy it into news
#endif
