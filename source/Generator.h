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
//
#ifndef Generator_H
#define Generator_H
#include "Assemblage.h"
#include "Shared.h"
#include <vector>

int GeneratorOnePoz(int NMode, uint8_t *XPolSupPos, int MaxDeg);
// Compute 1-d excitations and store it in XPolSupPos[Idm(xx)+mm] <= MaxDeg
// (xx,mm) in [0,NMode*MaxDeg +1[ x [0,NMode[.
//
uint32_t GeneratorXPoz(int NMode, int Eccit, int **XPolSupPos,
                       uint32_t SizeGenOld, int MaxDeg);
// Generator of double triple quadruple etc... raising excitations
// The multi indexes of excitations are stored in XPolSupPos[Idm(xx)+n]
// (xx>=SizeGenOld, n<NMode) NEccit is the number of couplings The maximal
// degrees Sum_n XPolSupPos[Idm(xx)+n] <= MaxDeg
//
uint32_t GetSizeEccitPoz(int NMode, int NEccit, int Maxd);
// Return the number of positive excitation in function of NMode,
// NEccit(Couplings) and Maxd degree of excitation : Sum_{i in [1,NEccit]} e_i <=
// Maxd
//
uint32_t GetNXDualHTrunc(uint8_t *DualHPos, uint32_t NXDualHTruncPos,
                         uint32_t *Corres, int NMode);
// Return the total number of excitations (+ et -) for the most contributive
// positive excitations indexed by Corres[ee], with ee in [0,NXDualHPlus[.
// DualHPos[Idm(xx)+mm] is defined for (xx,mm) in [0,NXDualHPlus[ x [0,NMode[.
//
uint32_t GeneratorTotPoz(int NMode, uint8_t *XPolSupPos, int *DegreEx,
                         int NCPol);
// Return the overall number of positive excitations in function
// of DegreEx in each direction and the maximal number of couplings NCPol.
// The excitations are stored in XPolSupPos,
// defined for XPolSupPos[Idm(xx)+mm], (xx,mm) in [0,NGenPoz[ x [0,NMode[.
//
uint64_t PolyX(uint8_t *ModeRez, uint8_t *ModeAct, uint8_t *DualHPos,
               std::vector<MatrixElem> QQ, uint32_t *PermutRez,
               uint32_t *PermutAct, uint32_t *Corres, KForce KFC,
               SizeArray *Size, uint8_t *Pid, double Freq0Max, double Freq0Min,
               int NMode, uint32_t NXDualHTruncPos, int NCPol, int DegrePol,
               int NPES, int Iteration, int NScreen, int *TabScreen,
               SizeArray SizeMax, float *RezVect, LocalFF LFF, double *EigVec,
               double **ZetaXYZ, double *Omega, CSC IJRez);
/*Compute the secondary space stored in  ModeRez applying the NXDualHTruncPos
excitations of DualHPos on configurations stored in ModeAct. And at the same
time the incomplete matrix vector product Hsb*X is calculated on the fly.
Additionally it computes the graph of the residual matrix stored in IJRez. The
space H*(A) is fetched at each iteration where A is the set of added basis
functions.*/
//
uint64_t PolyXNoG(uint8_t *ModeRez, uint8_t *ModeAct, uint8_t *DualHPos,
                  std::vector<MatrixElem> QQ, uint32_t *PermutRez,
                  uint32_t *PermutAct, uint32_t *Corres, KForce KFC,
                  SizeArray *Size, uint8_t *Pid, double Freq0Max,
                  double Freq0Min, int NMode, uint32_t NXDualHTruncPos,
                  int NCPol, int DegrePol, int NPES, int Iteration, int NScreen,
                  int *TabScreen, SizeArray SizeMax, float *RezVect,
                  LocalFF LFF, double *EigVec, double **ZetaXYZ, double *Omega);
/*Compute the secondary space stored in  ModeRez applying the NXDualHTruncPos
excitations of DualHPos on configurations stored in ModeAct. And at the same
time the complete matrix vector product Hsb*X is calculated on the fly. No graph
of residual matrix is computed in here.
The whole space H*(B) is fetched at each iteration.*/
//
#endif
