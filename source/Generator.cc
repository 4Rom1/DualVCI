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
#include "Generator.h"
#include "Assemblage.h"
#include "Basic.h"
#include "Graph.h"
#include "Shared.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
int GeneratorOnePoz(int NMode, uint8_t *XPolSupPos, int MaxDeg) {
  // Compute 1-d excitations and store it in XPolSupPos[Idm(xx)+mm] <= MaxDeg
  // (xx,mm) in [0,NMode*MaxDeg +1[ x [0,NMode[.
  int SizeGen = 0;
  int Cmpt = 0;
  SizeGen = (MaxDeg)*NMode + 1;
  int ss = 0;
  for (int mm = 0; mm < NMode; mm++) {
    XPolSupPos[mm] = 0;
  }
  Cmpt++;
  while (Cmpt < SizeGen) {
    for (int mm = 0; mm < NMode; mm++) {
      XPolSupPos[Idm(Cmpt) + mm] = 0;
    }
    for (int dd = 1; dd <= MaxDeg; dd++) {
      XPolSupPos[Idm(Cmpt) + ss] = dd;
      Cmpt++;
    }
    if (ss >= NMode) {
      printf("ss trop grand dans Generator one : problem %i \n", ss);
      printf("Cmpt  %i \n", Cmpt);
      printf("SizeGen  %i \n", SizeGen);
      return 0;
    }
    ss++;
  }
  if (SizeGen - Cmpt) {
    printf("SizeGen not equal to Cmpt %i %i \n", SizeGen, Cmpt);
    return 0;
  }
  return SizeGen;
}
uint32_t GetNXDualHTrunc(uint8_t *DualHPos, uint32_t NXDualHTruncPos,
                         uint32_t *Corres, int NMode) {
  // Return the total number of excitations (+ et -) for the most contributive
  // positive excitations indexed by Corres[ee], with ee in [0,NXDualHPlus[.
  // DualHPos[Idm(xx)+mm] is defined for (xx,mm) in [0,NXDualHPlus[ x [0,NMode[.
  int Cpld;
  uint32_t NXDualHTrunc = 0;
  uint32_t xx;
  for (uint32_t ee = 1; ee < NXDualHTruncPos;
       ee++) // First excitation = to zero
  {
    // Corres[ee] : index of the ee-th most contributive excitations
    xx = Corres[ee];
    Cpld = CmptNNULL(&DualHPos[Idm(xx)], NMode);
    NXDualHTrunc += Powo<int>(2, Cpld);
  }
  return NXDualHTrunc;
}
uint32_t GetSizeEccitPoz(int NMode, int NEccit, int Maxd) {
  // Return the number of positive excitation in function of NMode,
  // NEccit(Couplings) and Maxd degree of excitation : Sum_{i in [1,NEccit]} e_i
  // <= Maxd
  uint32_t CmptEct = 0;
  for (int ee = NEccit; ee <= Maxd; ee++) {
    CmptEct += (uint32_t)(factorial(ee - 1) /
                          (factorial(NEccit - 1) * factorial(ee - NEccit)));
  }
  // Number of NEccit-uplets = (NEccit among NMode)
  uint32_t NSurf = (uint32_t)(factorial(NMode) /
                              (factorial(NMode - NEccit) * factorial(NEccit)));
  uint32_t TotEct = NSurf * CmptEct;
  return TotEct;
}
uint32_t GeneratorXPoz(int NMode, int NEccit, uint8_t *XPolSupPos,
                       int SizeGenOld, int MaxDeg) {
  // Generator of double triple quadruple etc... raising excitations
  // The multi indexes of excitations are stored in XPolSupPos[Idm(xx)+n]
  // (xx>=SizeGenOld, n<NMode) NEccit is the number of couplings The maximal
  // degrees Sum_n XPolSupPos[Idm(xx)+n] <= MaxDeg
  if (NEccit < 2) {
    printf("Excitation must be equal to 2 at least\n");
    return 0;
  }
  uint32_t SizeGen = SizeGenOld;
  uint32_t NSurf;
  uint32_t Cmpt = 0;
  uint32_t TotEct = 0;
  uint32_t CmptEct = 0;
  for (int ee = NEccit; ee <= MaxDeg; ee++) {
    CmptEct += (uint32_t)(factorial(ee - 1) /
                          (factorial(NEccit - 1) * factorial(ee - NEccit)));
  }
  // Number of NEccit-uplets = (NEccit among NMode)
  NSurf = (uint32_t)(factorial(NMode) /
                     (factorial(NMode - NEccit) * factorial(NEccit)));
  TotEct = NSurf * CmptEct;
  // Define all NEccit-uplet combination
  uint8_t *MaxDegSurf;
  uint8_t *MultiDegrees;
  uint8_t *MaxDegrees;
  uint8_t *MultiSurf;
  allocate_TabUint8(MultiDegrees, NEccit) allocate_TabUint8(MaxDegrees, NEccit)
      allocate_TabUint8(MultiSurf, NEccit) MaxDegSurf = new uint8_t[NMode];
  // Initialsize degrees
  uint8_t Sum;
  for (int dd = 0; dd < NEccit; dd++) {
    MaxDegSurf[dd] = 1;
  }
  for (int dd = NEccit; dd < NMode; dd++) {
    MaxDegSurf[dd] = 0;
  }
  int ss;
  InitTab1(MultiDegrees, NEccit)
      // Combination i<j<k<l...
      // Surface i<j<l<m... until NEccit ={ MultiSurf[ss], 0<=ss<NEccit }
      // Given by the indexes xx of MaxDegSurf[xx]!=0
      // MultiSurf[ss]<NMode contains the indexes of surface for 0<=ss<NEccit
      // MaxDegrees[ss] will be the maximal degree of modal excitation on
      // coordinate MultiSurf[ss]
      do {
    ss = 0;
    for (int mm = 0; mm < NMode; mm++) {
      if (MaxDegSurf[mm]) {
        MultiSurf[ss] = mm;
        MaxDegrees[ss] =
            MaxDeg - NEccit +
            2; // Must take into account the others 1 (MaxDeg+1-(NEccit-1))
        MultiDegrees[ss] = 1;
        ss++;
      }
    }
    do {
      Sum = sum_array<uint8_t>(MultiDegrees, NEccit);
      // Sum_ee Multidegrees[ee] is all the possible sum of NEccit uplets lower
      // than MaxDegrees+1.
      if (Sum <= MaxDeg) {
        for (int ee = 0; ee < NEccit; ee++) {
          XPolSupPos[Idm(SizeGen) + MultiSurf[ee]] = MultiDegrees[ee];
        }
        SizeGen++;
      }
    } while (nested_loop1(MultiDegrees, MaxDegrees, NEccit));
    Cmpt++;
  }
  while (std::prev_permutation(MaxDegSurf, MaxDegSurf + NMode))
    ;
  if (SizeGen - SizeGenOld != TotEct) {
    printf("Problem total excitation evaluation,\n \
SizeGen-SizeGenOld,TotEct, MaxDeg, NEccit, Cmpt, NSurf\n \
%u %u %u %u %u %u \n",
           SizeGen - SizeGenOld, TotEct, MaxDeg, NEccit, Cmpt, NSurf);
    return 0;
  }
  FreeTab(MultiDegrees) FreeTab(MaxDegrees) delete[] MultiSurf;
  delete[] MaxDegSurf;
  return SizeGen;
}
uint32_t GeneratorTotPoz(int NMode, uint8_t *XPolSupPos, int *DegreEx,
                         int NCPol) {
  // Return the overall number of positive excitations in function
  // of DegreEx in each direction and the maximal number of couplings NCPol.
  // The excitations are stored in XPolSupPos,
  // defined for XPolSupPos[Idm(xx)+mm], (xx,mm) in [0,NGenPoz[ x [0,NMode[.
  uint32_t SizeGenOld = GeneratorOnePoz(NMode, XPolSupPos, DegreEx[0]);
  // double and Triple and superior excitation
  for (int ee = 2; ee <= NCPol; ee++) {
    SizeGenOld =
        GeneratorXPoz(NMode, ee, XPolSupPos, SizeGenOld, DegreEx[ee - 1]);
    if (!SizeGenOld) {
      return 0;
    }
  }
  return SizeGenOld;
}
uint64_t PolyX(uint8_t *ModeRez, uint8_t *ModeAct, uint8_t *DualHPos,
               std::vector<MatrixElem> QQ, uint32_t *PermutRez,
               uint32_t *PermutAct, uint32_t *Corres, KForce KFC,
               SizeArray *Size, uint8_t *Pid, double Freq0Max, double Freq0Min,
               int NMode, uint32_t NXDualHTruncPos, int NCPol, int DegrePol,
               int NPES, int Iteration, int NScreen, int *TabScreen,
               SizeArray SizeMax, float *RezVect, LocalFF LFF, double *EigVec,
               double **ZetaXYZ, double *Omega, CSC IJRez) {
  /*Compute the secondary space stored in  ModeRez applying the NXDualHTruncPos
  excitations of DualHPos on configurations stored in ModeAct. And at the same
  time the incomplete matrix vector product Hsb*X is calculated on the fly.
  Additionally it computes the graph of the residual matrix stored in IJRez. The
  space H*(A) is fetched at each iteration where A is the set of added basis
  functions. Pid : Maximal quantum level in each direction. DualHPos : Multi
  index arrays of raising excitations. Defined for DualHPos[Idm(xx)+mm], (xx,mm)
  in [0,NXDualHPlus[ x [0,NMode[. Corres : Correspondance array giving the
  indexes of the first NXDualHTruncPos most contributive excitations of
  DualHPos. Freq0Max : harmonic energy wall above the ground state. ModeAct :
  Multi-dimensional array for active space B : ModeAct[Idm(nn)+mm], (nn,mm) in
  [0,SizeActMax[ x [0,NMode[. ModeRez : Multi-dimensional array for residual
  space : ModeRez[Idm(nn)+mm], (nn,mm) in [0,SizeRezMax[ x [0,NMode[. PermutRez
  : indexes of permutations for sorted elements in ModeRez PermutAct : indexes
  of permutations for sorted elements in ModeAct TabScreen : Indexes of targeted
  eigen-pairs Force constants: KFC. KijNCpld[dd][mm] : non coupled force
  constants, defined for (dd,mm) in [0,DegrePol[x[0,NMode[ (degree dd+1, mode
  mm<NMode). KFC.KijCpld[kk] : coupled force constants, defined for kk in
  [0,NPES[. KFC.Monm[kk][mm] = degree of monomial kk for coordinate mm in PES.
  Local force field:
  LFF[ii(xx)] : local force field associated with positive excitation
  &DualHPos[Idm(xx)]. LFF.Idx[ii(xx)] are defined for ii
  [LFF.Num[xx],LFF.Num[xx+1][ LFF.Idx[ii(xx)] < 0 are indexes of non coupled
  force constants: KFC.KijNCpld[dd][mm] with dd=-LFF.Idx[ii]/NMode,
  mm=-LFF.Idx[ii]-dd*NMode; 0 <= LFF.Idx[ii(xx)] < NPES and are indexes of
  coupled force constants KFC.KijCpld[LFF.Idx[ii(xx)]]. LFF.Idx[ii] >= NPES are
  key numbers of the rotational coefficients
  nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES with a unique corresponding
  (ni,nj,nk,nl) Compressed sparse  column format for IJRez: IJRez.NJ[jj] =
  number of the first nnz elmement of column jj IJRez.I[nn] = line number of nnz
  element number nn NCPol: maximal number of couplings. DegrePol : Maximal
  degree of PES. Omega : harmonic frequencies in hartree for conversion of
  rotational elements. QQ[ii].Coeff[ii][jj] are matrix elements of operator Q^ii
  for ii in [0,DegrePol], QQ[DegrePol+1].Coeff[ii][jj] are matrix elements of
  operator D2Q (second order derivative) QQ[DegrePol+2].Coeff[ii][jj] are matrix
  elements of operator D1Q (first order derivative) QQ[DegrePol+3].Coeff[ii][jj]
  are matrix elements of operator QD1Q QQ[DegrePol+4].Coeff[ii][jj] are matrix
  elements of operator D1QQ Size[Iteration+1].DimRez: Size of residual space
  beeing incremented; Size[Iteration+1].DimAct: Current size of active space B
  already incremented at previous iteration; Size[Iteration].DimAct: Previous
  size of active space B used as a starting index in the general loop;
  Size[Iteration+1].NNZAct: Current number of NNZ of active matrix Hb already
  incremented in routine AssembleHarmCSC; Size[Iteration+1].NNZRez: NNZ of
  residual matrix beeing incremented because DoGraph # 0;
  */
  int DegrePolP1 = DegrePol + 1;
  uint32_t Lin = 0;
  long LinRez = -1; // Important to evaluate RestQSort at first
  // Size of bytes to look for in linearsearch arrays
  size_t SizeBit = NMode * sizeof(uint8_t);
  int ss = 0;      // Index for normal coordinate on a surface
  int go = 1;      // Index to continue if Tester is not negative
  int TestNeg = 0; // Test if the sum is negative
  int CptQSort = 0;
  long LinAct = 0;
  int RestQSort =
      0; // Rest of division of Size[Iteration+1].DimRez by MaxQsortSeq
  int KQSort = 0;
  int Cpld;
  double MatElem;
  uint32_t xx;
  double Freq0;
  double Freq0T;
  int AssignFirstCol = 1;
  uint64_t CheckOutputPolyX = 0;
  uint8_t *Tester = NULL; // Variable multi index during the loop
  uint8_t *Surf = NULL;
  uint8_t *MultiDegrees2 = NULL;
  uint8_t *MaxDegrees2 = NULL;
  int *Ect = NULL;
  Ect = new int[NCPol];
  allocate_TabUint8(Tester, NMode) allocate_TabUint8(Surf, NCPol)
      allocate_TabUint8(MultiDegrees2, NCPol)
          allocate_TabUint8(MaxDegrees2, NCPol)
      // Initialize the number of NNZ elements of the Active and residual matrix
      // At next iteration residual space has got minimum Size[Iteration].DimRez
      // elements without the nulified elements handled in routines AssmbleHR.
      Size[Iteration + 1]
          .DimRez = Size[Iteration].DimRez;
  Size[Iteration + 1].NNZRez = Size[Iteration].NNZRez;
  for (int dd = 0; dd < NCPol; dd++) {
    MaxDegrees2[dd] = 2; //-1 and +1
  }
  // Sort active space in order given by memcmp.
  QuickSortMode(ModeAct, 0, Size[Iteration + 1].DimAct - 1, PermutAct, SizeBit);
  uint32_t Start = Size[Iteration].DimAct;
  for (Lin = Start; Lin < Size[Iteration + 1].DimAct;
       Lin++) // Start from previous subspace construction
  {
    Freq0 = GetFreq0(&ModeAct[Idm(Lin)], NMode, KFC);
    AssignFirstCol = 1;
    // Browse all the non null excitations
    for (uint32_t ee = 1; ee < NXDualHTruncPos; ee++) {
      // Start with tester equal to active mode
      xx = Corres[ee];
      Cpld = ConvertX(&DualHPos[Idm(xx)], Surf, NMode);
      InitTabInt(MultiDegrees2, Cpld)
          memcpy(Tester, &ModeAct[Idm(Lin)], SizeBit);
      // Start with tester equal to active mode
      do {
        ss = 0;
        go = 1;
        while ((ss < Cpld) && go) { // Generate product of positive and negative
                                    // excitations from positive ones.
          Ect[ss] =
              DualHPos[Idm(xx) + Surf[ss]] * Powo<int>(-1, MultiDegrees2[ss]);
          TestNeg = Ect[ss] + ModeAct[Idm(Lin) + Surf[ss]];
          if (TestNeg < 0 || TestNeg >= Pid[Surf[ss]]) {
            go = 0;
          } else {
            Tester[Surf[ss]] = TestNeg;
          }
          ss++;
        }
        if (go) {
          Freq0T = Freq0;
          for (ss = 0; ss < Cpld; ss++) {
            Freq0T += Ect[ss] * KFC.KijNCpld[1][Surf[ss]];
          }
          if (Freq0T > Freq0Max || Freq0T < Freq0Min) {
            go = 0;
          }
        }
        if (go) {
          // Test if Tester already constitutes one of the neighboors of the
          // configuration having the number Lin. Search Tester in
          // ModeAct[PermutAct[1:Size[Iteration+1].DimAct]]
          LinAct = QsearchMode(Tester, ModeAct, PermutAct,
                               Size[Iteration + 1].DimAct, SizeBit);
          if (LinAct <
              0) // It has to be negative to built residual space and vectors
          {
            // If the corresponding neighboor has already been created.
            // Looking for already added residual modes
            // Useless if we cross new active space with residual one
            if ((LinRez < 0)) // From previous iteration :
                              // Size[Iteration+1].DimRez has been incremented
            {
              if (Size[Iteration + 1].DimRez < MaxQSortSeq) {
                KQSort = 0;
                RestQSort = Size[Iteration + 1].DimRez;
              } else {
                KQSort = Size[Iteration + 1].DimRez / MaxQSortSeq;
                RestQSort = Size[Iteration + 1].DimRez - KQSort * MaxQSortSeq;
              }
              if (RestQSort == 0 && KQSort) {
                QuickSortModeShift(ModeRez, 0, MaxQSortSeq - 1,
                                   &PermutRez[(KQSort - 1) * MaxQSortSeq],
                                   SizeBit, (KQSort - 1) * MaxQSortSeq);
                // MaxQSortSeq elements are sorted each time
                // Size[Iteration+1].DimRez is a multiple of the latter.
              }
            }
            LinRez = -1;
            if (KQSort) {
              CptQSort = 0;
              while (LinRez < 0 && CptQSort < KQSort) {
                // Search Tester in ModeRez[PermutRez[CptQSort*MaxQSortSeq :
                // (CptQSort+1)*MaxQSortSeq]+CptQSort*MaxQSortSeq]
                LinRez = QsearchModeShift(
                    Tester, ModeRez, &PermutRez[CptQSort * MaxQSortSeq],
                    MaxQSortSeq, SizeBit, CptQSort * MaxQSortSeq);
                CptQSort++;
              }
            }
            if (LinRez < 0) {
              {
                LinRez = LinearModeSearch(&ModeRez[Idm(KQSort * MaxQSortSeq)],
                                          Tester, SizeBit, RestQSort);
                if (LinRez >= 0) {
                  LinRez += KQSort * MaxQSortSeq;
                }
              }
            }
            MatElem = MatrixElement(NMode, LFF, xx, DegrePolP1, NPES,
                                    &ModeAct[Idm(Lin)], Tester, QQ, KFC,
                                    ZetaXYZ, Omega);
            if ((LinRez < 0)) {
              if (Size[Iteration + 1].DimRez < SizeMax.DimRez &&
                  Size[Iteration + 1].NNZRez < SizeMax.NNZRez) {
                // Last element will be copied in the the residual space ModeRez
                memcpy(&ModeRez[Idm(Size[Iteration + 1].DimRez)], Tester,
                       SizeBit);
                for (int ll = 0; ll < NScreen; ll++) {
                  RezVect[Size[Iteration + 1].DimRez + SizeMax.DimRez * ll] +=
                      MatElem *
                      (EigVec[TabScreen[ll] * Size[Iteration + 1].DimAct +
                              Lin]);
                }
                if (AssignFirstCol) {
                  IJRez.NJ[Lin] = Size[Iteration + 1].NNZRez;
                  AssignFirstCol = 0;
                }
                IJRez.I[Size[Iteration + 1].NNZRez] =
                    Size[Iteration + 1].DimRez;
                Size[Iteration + 1].DimRez++;
                Size[Iteration + 1].NNZRez++;
              } else if (Size[Iteration + 1].DimRez >= SizeMax.DimRez) {
                printf(
                    "**|Bs| Too Big : To poursue,increase Memory or KNREZ**\n");
                CheckOutputPolyX = 0;
                goto end;
              } // Keep
              else if (Size[Iteration + 1].NNZRez >= SizeMax.NNZRez) {
                printf("**NNZ(Hsb) Too Big : To poursue,increase Memory or "
                       "KNZREZ**\n");
                CheckOutputPolyX = 0;
                goto end;
              }
            } else if (LinRez >= 0) {
              for (int ll = 0; ll < NScreen; ll++) {
                RezVect[LinRez + SizeMax.DimRez * ll] +=
                    MatElem *
                    EigVec[TabScreen[ll] * Size[Iteration + 1].DimAct + Lin];
              }
              if (Size[Iteration + 1].NNZRez < SizeMax.NNZRez) {
                if (AssignFirstCol) // Check if elements Lin of IJRez.NJ has
                                    // already been set
                { // Remember that IJRez.NJ[Lin] give the number of NNZ for the
                  // first element of line Lin
                  IJRez.NJ[Lin] = Size[Iteration + 1].NNZRez;
                  AssignFirstCol = 0;
                }
                IJRez.I[Size[Iteration + 1].NNZRez] = LinRez;
                Size[Iteration + 1].NNZRez++;
              } else {
                printf("**NNZ(Hsb) Too Big : To poursue, increase Memory or "
                       "KNZREZ**\n");
                CheckOutputPolyX = 0;
                goto end;
              }
            }
          }
        }
      } while (nested_loop(MultiDegrees2, MaxDegrees2, Cpld));
    }
  }
  CheckOutputPolyX = Size[Iteration + 1].DimRez; // No problem of execution
end:
  IJRez.NJ[Size[Iteration + 1].DimAct] = Size[Iteration + 1].NNZRez;
  FreeTab(MultiDegrees2) FreeTab(MaxDegrees2) FreeTab(Surf) delete[] Tester;
  delete[] Ect;
  return CheckOutputPolyX;
}
uint64_t PolyXNoG(uint8_t *ModeRez, uint8_t *ModeAct, uint8_t *DualHPos,
                  std::vector<MatrixElem> QQ, uint32_t *PermutRez,
                  uint32_t *PermutAct, uint32_t *Corres, KForce KFC,
                  SizeArray *Size, uint8_t *Pid, double Freq0Max,
                  double Freq0Min, int NMode, uint32_t NXDualHTruncPos,
                  int NCPol, int DegrePol, int NPES, int Iteration, int NScreen,
                  int *TabScreen, SizeArray SizeMax, float *RezVect,
                  LocalFF LFF, double *EigVec, double **ZetaXYZ,
                  double *Omega) {
  /*Compute the secondary space stored in  ModeRez applying the NXDualHTruncPos
  excitations of DualHPos on configurations stored in ModeAct. And at the same
  time the complete matrix vector product Hsb*X is calculated on the fly. No
  graph of residual matrix is computed in here. The whole space H*(B) is fetched
  at each iteration. Pid : Maximal quantum level in each direction. DualHPos :
  Multi index arrays of raising excitations. Defined for DualHPos[Idm(xx)+mm],
  (xx,mm) in [0,NXDualHPlus[ x [0,NMode[. Changed to DualHPos[Idx(xx)+mm] Corres
  : Correspondance array giving the indexes of the first NXDualHTruncPos most
  contributive excitations of DualHPos. Freq0Max : harmonic energy wall above
  the ground state. ModeAct : Multi-dimensional array for active space. Defined
  for ModeAct[Idm(nn)+mm], (nn,mm) in [0,SizeActMax[ x [0,NMode[. ModeRez :
  Multi-dimensional array for residual space : ModeRez[Idm(nn)+mm], (nn,mm) in
  [0,SizeRezMax[ x [0,NMode[. PermutRez : indexes of permutations for sorted
  elements in ModeRez PermutAct : indexes of permutations for sorted elements in
  ModeAct TabScreen : Indexes of targeted eigen-pairs Force constants: KFC.
  KijNCpld[dd][mm] : non coupled force constants, defined for (dd,mm) in
  [0,DegrePol[x[0,NMode[ (degree dd+1, mode mm<NMode). KFC.KijCpld[kk] : coupled
  force constants, defined for kk in [0,NPES[. KFC.Monm[kk][mm] = degree of
  monomial kk for coordinate mm in PES. Local force field: LFF.Idx[ii(xx)] are
  defined for ii [LFF.Num[xx],LFF.Num[xx+1][ LFF[ii(xx)] : local force field
  associated with positive excitation &DualHPos[Idm(xx)]. LFF.Idx[ii(xx)] < 0
  are indexes of non coupled force constants: KFC.KijNCpld[dd][mm] with
  dd=-LFF.Idx[ii]/NMode, mm=-LFF.Idx[ii]-dd*NMode; 0 <= LFF.Idx[ii(xx)] < NPES
  and are indexes of coupled force constants KFC.KijCpld[LFF.Idx[ii(xx)]].
  LFF.Idx[ii] >= NPES are key numbers of the rotational coefficients
  nl*NMode^3+nk*NMode^2+nj*NMode+ni+NPES with a unique corresponding
  (ni,nj,nk,nl) NCPol : maximal number of couplings. DegrePol : maximal degree
  in the PES. Omega : harmonic frequencies in hartree for conversion of
  rotational elements. Matrix elements: QQ[dd].Coeff[ii][jj] are matrix elements
  of operator Q^dd for dd in [0,DegrePol], QQ[DegrePol+1].Coeff[ii][jj] are
  matrix elements of operator D2Q (second order derivative)
  QQ[DegrePol+2].Coeff[ii][jj] are matrix elements of operator D1Q (first order
  derivative) QQ[DegrePol+3].Coeff[ii][jj] are matrix elements of operator QD1Q
  QQ[DegrePol+4].Coeff[ii][jj] are matrix elements of operator D1QQ
  Size:
  Size[Iteration+1].DimRez: Size of residual space beeing incremented;
  Size[Iteration+1].DimAct: Current size of active space already incremented at
  previous iteration; Size[Iteration+1].NNZAct: Current number of NNZ of active
  matrix already incremented in routine AssembleHarmCSC;
  Size[Iteration+1].NNZRez: NNZ of residual matrix not beeing incremented
  because DoGraph=0;
  */
  int DegrePolP1 = DegrePol + 1;
  uint32_t Lin = 0;
  long LinRez = -1; // Important to evaluate RestQSort at first
  // Size of bytes to look for in linearsearch arrays
  size_t SizeBit = NMode * sizeof(uint8_t);
  int ss = 0;      // Index for normal coordinate on a surface
  int go = 1;      // Index to continue if Tester is not negative
  int TestNeg = 0; // Test if the sum is negative
  int CptQSort = 0;
  long LinAct = 0;
  int RestQSort =
      0; // Rest of division of Size[Iteration+1].DimRez by MaxQsortSeq
  int KQSort = 0;
  int Cpld;
  double MatElem;
  uint32_t xx;
  double Freq0;
  double Freq0T;
  uint64_t CheckOutputPolyX = 0;
  uint8_t *Tester = NULL; // Variable multi index during the loop
  uint8_t *Surf = NULL;
  uint8_t *MultiDegrees2 = NULL;
  uint8_t *MaxDegrees2 = NULL;
  int *Ect = NULL;
  Ect = new int[NCPol];
  allocate_TabUint8(Tester, NMode) allocate_TabUint8(Surf, NCPol)
      allocate_TabUint8(MultiDegrees2, NCPol)
          allocate_TabUint8(MaxDegrees2, NCPol)
      // Initialize the number of NNZ elements of the Active and residual matrix
      // At next iteration residual space has got minimum Size[Iteration].DimRez
      // elements
      Size[Iteration + 1]
          .DimRez = Size[Iteration].DimRez;
  for (int dd = 0; dd < NCPol; dd++) {
    MaxDegrees2[dd] = 2; //-1 and +1
  }
  // Sort Active space
  QuickSortMode(ModeAct, 0, Size[Iteration + 1].DimAct - 1, PermutAct, SizeBit);
  uint32_t Start = 0;
  for (Lin = Start; Lin < Size[Iteration + 1].DimAct;
       Lin++) // Start from 0 : takes more time compared to PolyX
  {
    Freq0 = GetFreq0(&ModeAct[Idm(Lin)], NMode, KFC);
    // Browse all the non null excitations
    for (uint32_t ee = 1; ee < NXDualHTruncPos; ee++) {
      // Start with tester equal to active mode
      xx = Corres[ee];
      Cpld = ConvertX(&DualHPos[Idm(xx)], Surf, NMode);
      InitTabInt(MultiDegrees2, Cpld)
          // Start with tester equal to active mode
          memcpy(Tester, &ModeAct[Idm(Lin)], SizeBit);
      do {
        ss = 0;
        go = 1;
        while ((ss < Cpld) && go) { // Generate product of positive and negative
                                    // excitations from positive ones.
          Ect[ss] =
              DualHPos[Idm(xx) + Surf[ss]] * Powo<int>(-1, MultiDegrees2[ss]);
          TestNeg = Ect[ss] + ModeAct[Idm(Lin) + Surf[ss]];
          if (TestNeg < 0 || TestNeg >= Pid[Surf[ss]]) {
            go = 0;
          } else {
            Tester[Surf[ss]] = TestNeg;
          }
          ss++;
        }
        if (go) {
          Freq0T = Freq0;
          for (ss = 0; ss < Cpld; ss++) {
            Freq0T += Ect[ss] * KFC.KijNCpld[1][Surf[ss]];
          }
          if (Freq0T > Freq0Max || Freq0T < Freq0Min) {
            go = 0;
          }
        }
        if (go) {
          // Test if Tester already constitutes one of the neighboors of the
          // configuration having the number Lin. Search Tester in
          // ModeAct[PermutAct[1:Size[Iteration+1].DimAct]]
          LinAct = QsearchMode(Tester, ModeAct, PermutAct,
                               Size[Iteration + 1].DimAct, SizeBit);
          if (LinAct <
              0) // It has to be negative to built residual space and vectors
          {
            // If the corresponding neighboor has already been created.
            // Looking for already added residual modes
            // Useless if we cross new active space with residual one
            if ((LinRez < 0)) // From previous iteration :
                              // Size[Iteration+1].DimRez has been incremented
            {
              if (Size[Iteration + 1].DimRez < MaxQSortSeq) {
                KQSort = 0;
                RestQSort = Size[Iteration + 1].DimRez;
              } else {
                KQSort = Size[Iteration + 1].DimRez / MaxQSortSeq;
                RestQSort = Size[Iteration + 1].DimRez - KQSort * MaxQSortSeq;
              }
              if (RestQSort == 0 && KQSort) {
                QuickSortModeShift(ModeRez, 0, MaxQSortSeq - 1,
                                   &PermutRez[(KQSort - 1) * MaxQSortSeq],
                                   SizeBit, (KQSort - 1) * MaxQSortSeq);
                // MaxQSortSeq elements are sorted each time
                // Size[Iteration+1].DimRez is a multiple of the latter.
              }
            }
            LinRez = -1;
            if (KQSort) {
              CptQSort = 0;
              while (LinRez < 0 && CptQSort < KQSort) {
                // Search Tester in ModeRez[PermutRez[CptQSort*MaxQSortSeq :
                // (CptQSort+1)*MaxQSortSeq]+CptQSort*MaxQSortSeq]
                LinRez = QsearchModeShift(
                    Tester, ModeRez, &PermutRez[CptQSort * MaxQSortSeq],
                    MaxQSortSeq, SizeBit, CptQSort * MaxQSortSeq);
                CptQSort++;
              }
            }
            if (LinRez < 0) {
              {
                LinRez = LinearModeSearch(&ModeRez[Idm(KQSort * MaxQSortSeq)],
                                          Tester, SizeBit, RestQSort);
                if (LinRez >= 0) {
                  LinRez += KQSort * MaxQSortSeq;
                }
              }
            }

            MatElem = MatrixElement(NMode, LFF, xx, DegrePolP1, NPES,
                                    &ModeAct[Idm(Lin)], Tester, QQ, KFC,
                                    ZetaXYZ, Omega);
            if ((LinRez < 0)) {
              if (Size[Iteration + 1].DimRez < SizeMax.DimRez) {
                // Last element will be copied in the the residual space mode
                // ModeRez
                memcpy(&ModeRez[Idm(Size[Iteration + 1].DimRez)], Tester,
                       SizeBit);
                for (int ll = 0; ll < NScreen; ll++) {
                  RezVect[Size[Iteration + 1].DimRez + SizeMax.DimRez * ll] +=
                      MatElem *
                      (EigVec[TabScreen[ll] * Size[Iteration + 1].DimAct +
                              Lin]);
                }
                Size[Iteration + 1].DimRez++;
              } else // Keep
              {
                printf("**|Bs| Too Big : To poursue, increase Memory or "
                       "KNREZ**\n");
                CheckOutputPolyX = 0;
                goto end;
              } // Keep
            } else if (LinRez >= 0) {
              for (int ll = 0; ll < NScreen; ll++) {
                RezVect[LinRez + SizeMax.DimRez * ll] +=
                    MatElem *
                    EigVec[TabScreen[ll] * Size[Iteration + 1].DimAct + Lin];
              }
            }
          }
        }
      } while (nested_loop(MultiDegrees2, MaxDegrees2, Cpld));
    }
  }
  CheckOutputPolyX = Size[Iteration + 1].DimRez; // No problem
end:
  FreeTab(MultiDegrees2) FreeTab(MaxDegrees2) FreeTab(Surf) delete[] Tester;
  delete[] Ect;
  return CheckOutputPolyX;
}
