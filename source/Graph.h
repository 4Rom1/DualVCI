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

/*--------Additional sources-------------------

Routines QuickSortMode, QuickSortModeShift are made from the Quicksort algorithm
[1].

[1] C. A. R. Hoare. 1961. Algorithm 64: Quicksort. Commun. ACM 4, 7 (July 1961),
321-. DOI=http://dx.doi.org/10.1145/366622.366644

+

Routines QsearchMode, QsearchModeShift are made from the half interval search
algorithm also known as binary search [2].

[2] Willams, Jr., Louis F. (1975). A modification to the half-interval search
(binary search) method. Proceedings of the 14th ACM Southeast Conference. pp.
95–101. doi:10.1145/503561.503582.

----------------------------------------------

The routines nested_loop and nested_loop1 are depth variable nested loops
sharing the same structure as the one presented in STACKOVERFLOW [3]

[3] How to set a variable for nested loops?
Posted on stackoverflow (by M Oehm) at the following address
http://stackoverflow.com/questions/25785863/how-to-set-a-variable-for-nested-loops
*/

#ifndef Graph_H
#define Graph_H
#include "Shared.h"
//
#define MinSequential 300000
//
// Maximal elements to sort
#define MaxQSortSeq 30000
//
#define LimQSort 5000000
//
#define MinSeqAss 30000
//
// Allocation for modes
//""
#define allocate_TabMode(Mode1, n)                                             \
  uint8_t *Mode1 = NULL;                                                       \
  Mode1 = new uint8_t[n * NMode];                                              \
  for (uint64_t k = 0; k < n * NMode; k++) {                                   \
    Mode1[k] = 0;                                                              \
  }
//
void DeleteMode(ConfigId &M1);
int nested_loop(uint8_t *ix, uint8_t *im, int depth);
// Simulation of a nested loop of a given depth to generate multi-indexes stored
// in ix with ix[Cpt] bounded by im[Cpt].
// ix must be initialized by 0.

int nested_loop1(uint8_t *ix, uint8_t *im, int depth);
// Same as nested_loop but ix is initialized to one instead

long LinearModeSearch(uint8_t *TabMode, uint8_t *ModeT, size_t NMode,
                      uint32_t SizeTab);
// Look if the multi-indexes ModeT is located in the array of multi-indexes
// TabMode. Size of TabMode is SizeTab*SizeBit. size_t
// SizeBit=NMode*sizeof(uint8_t)=NMode;
//
int LinearIntSearch(int *TabInt, int Key, int SizeTab);
// Search for the int Key in array TabInt of size SizeTab.
//
long LinearU8Search(uint8_t MatU8[MaxTarget][MaxNormal], uint8_t *ModeKey,
                    size_t SizeBit, uint32_t SizeTab);
// Search for the array ModeKey of size SizeBit in matrix MatU8.
//
long QsearchMode(uint8_t *key, uint8_t *base, uint32_t *Permuter, uint32_t size,
                 int NMode);
// Binary search of configurations "Key" in basis "base"
// already sorted in the memcmp order given by the indexes of Permuter array
//
void QuickSortMode(uint8_t *MIndex, int left, int right, uint32_t *Permuter,
                   int NMode);
// Sort multi-indexes in MIndex according to memcmp order and return ascending
// order indexes in Permuter
//
int32_t partition(uint8_t *MIndex, uint32_t Permuter[], int left, int right,
                  int NMode);
/*
Permuter[rr] is the pivot element.
All the smaller and bigger elements than the pivot are respectively
placed to the left and the right of the latter.
Comparison of elements is done with memcmp c subroutine.
 */
//
int ConvertX(uint8_t *Tab, uint8_t *Surf, int NMode);
// Convert the multi-index Tab of size NMode
// on surface Surf={Positions({Tab[ii] # 0})}
//
uint32_t InitB0(uint8_t *ModeAct, KForce KFC, uint8_t *Pid, int NMode,
                uint32_t SizeActMax, double TargetMax);
// Compute the initial subspace En<TargetMax, where En harmonic energy
// and n are configuration arrays stored in ModeAct.
// Use a greedy algorithm that sweeps all the region until TargetMax.
// Here is used only one modal nested excitations on the ground state (0,...,0).
//
void VoisinCpld(int IndexMonm, int NMode, KForce KFC, uint8_t *NVMode,
                uint8_t **ModeV);
/*Computes the neighboors of the zero state (0,0,0..0)=Mode0
considering the monom number IndexMonm.
 H=sum_kk prod_mm H_{mm,kk}. Here kk=IndexMonm.
Return the Neighboor-indexes of the normal coordinate mm
available in the array ModeV[NVMode[mm]][mm]:
<ModeV(NVMode(mm),mm)|H_{mm,IndexMonm}|Mode0(mm)> # 0
*/
//
void VoisinCpldRot(int NMode, uint8_t *Tester, uint8_t *NVMode,
                   uint8_t **ModeV);
/*Computes the rotational terms neighboors of zero state Mode0=(0,0,0..0).
Return the neighboor-indexes of the normal coordinate mm available in the array
ModeV[NVMode[mm]][mm]: <ModeV(NVMode(mm),mm)|H_{mm,IndexMonm}|Mode0(mm)>  # 0 .
Tester is the equivalent monomial c^(ijkl)= 1i+1j+1k+1l, (sum of canonical
vectors)
*/
//
void AfficheNu(uint8_t *Mode, int NMode);
// Print Modal components : dd*nu[mm]+, ... for a normal coordinate mm and
// a non null degree dd (corresponding to a degree of Hermite basis function)
//
void AfficheNuTex(uint8_t *Mode, int NMode);
// Print Modal components in Tex : dd*\nu{mm}+, ... for a normal coordinate mm
// and a non null degree dd (corresponding to a degree of Hermite basis function)
//
void AfficheNuSupTex(uint8_t *ModeAct, double *EigVec, int NMode,
                     uint32_t DimAct, double ThrCoor);
// Print Modal components in Tex format and vector components bigger than
// ThrCoor
//
void AfficheNuSup(uint8_t *ModeAct, double *EigVec, int NMode, uint32_t DimAct,
                  double ThrCoor);
// Print Modal components and vector components bigger than ThrCoor
//
void QuickSortModeShift(uint8_t *MIndex, int ll, int rr, uint32_t *Permuter,
                        int NMode, uint64_t Shift);
//
int32_t PartShift(uint8_t *MIndex, uint32_t Permuter[], int ll, int rr,
                  int NMode, uint64_t Shift);
/*
Permuter[rr]+Shift is the pivot element.
All the smaller and bigger elements than the pivot are respectively
placed to the left and the right of the latter.
There is an additional shift to avoid memory overflow of the residual indexes
than can be bigger than UINT32MAX. Comparison of elements is done with memcmp c
subroutine.
*/
//
long QsearchModeShift(uint8_t *Key, uint8_t *base, uint32_t *Permuter,
                      uint32_t size, int NMode, uint64_t Shift);
// Same as previous subroutine : Binary search of configurations "Key" in basis
// "base" already sorted in the memcmp order given by the indexes of Permuter
// array. But indexes are shifted by the number Shift in order to avoid troubles
// when they oversize MAXUINT32 Because Permuter is an array of UINT32_t
//
uint32_t CptLFFPerEx(int DegrePol, int NMode, int NPES, int DoCor,
                     uint32_t NGenPoz, KForce KFC, uint8_t *GenNeigh,
                     int *NLFFPerEx, double **ZetaXYZ, double ThrPES,
                     int IncK2);
// Routine that counts the number of force constants per local force field.
// The routine return the number of positive excitations in H*.
// XPolSupPos contains the upper limit of positive excitations in H*.
// XPolSupPos[Idm(xx)+mm] is defined for (xx,mm) in [0,NGenPoz[ x [0,NMode[.
// IncK2:Say if the the quadratic term should be included or not in the local
// force field
//
uint32_t AssignLFFPerEx(int DegrePol, int NMode, int DoRot, int NPES,
                        uint32_t NXDualHPlus, KForce KFC, uint8_t *DualHPos,
                        LocalFF LFF, double **ZetaXYZ, double ThrKX,
                        double ThrPES, uint32_t *Permuter, uint32_t *Corres,
                        int IncK2);
// Routine that copies indexes of positive excitations of most contributing
// Local Force Fields into Corres. Corres : Correspondance array giving the
// indexes of the first NXDualHTruncPos most contributive positive excitations of
// H*. Select only the ones such that Sum|FC|>ThrKX DualHPos contains positive
// excitations in H*. DualHPos[Idm(xx)+mm] is defined for (xx,mm) in
// [0,NXDualHPlus[ x [0,NMode[. IncK2:Say if the the quadratic term should be
// included or not in the local force field
//
#endif
