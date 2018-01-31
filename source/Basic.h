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
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Updated sources can be found at 
    https://github.com/4Rom1/DualVCI 

--------Sources-------------------

This part contains only public domain routines: 
Sign of a number, Min, Max, Sum of arrays and bubble sort algorithm.

For bubble sort see for example https://en.wikipedia.org/wiki/Bubble_sort
*/

#ifndef Basic_H
#define Basic_H
//
#include <cmath>
#include <vector>
//
//using namespace std : no need
//
template<typename T>
int find_minimum(T *MonTab, int n) 
{
  int index, c;
  T min;
  min = MonTab[0];
  index = 0;
  for (c = 1; c < n; c++) {
    if (MonTab[c] < min) 
    {
       index = c;
       min = MonTab[c];
    }
  }
  return index;
}
//
template<typename T>
int find_maximum(T *MonTab, int n)
{
  int index, c;
  T max;
  max = MonTab[0];
  index = 0;
  for (c = 1; c < n; c++) 
  {
    if (MonTab[c] > max) 
    {
       index = c;
       max = MonTab[c];
    }
  }
  return index;
}
//
template<typename T>
T MinVal(T *MonTab, int n) 
{
  int c;
  T min;
  min = MonTab[0];
  for (c = 1; c < n; c++) 
   {
    if (MonTab[c] < min) 
    {
       min = MonTab[c];
    }
  }
  return min;
}
//
template<typename T>
T MaxVal(T *MonTab, int n) 
{
  int c;
  T max; 
  max = MonTab[0];
  for (c = 1; c < n; c++) 
   {
    if (MonTab[c] > max) 
    {
       max = MonTab[c];
    }
  }
  return max;
}
//
template<typename T>
T Min(T x, T y)
{
  return (x < y) ? x : y;
}

template<typename T>
T Max(T x, T y)
{
  return (x > y) ? x : y;
}
//
template<typename T>
T SIGN(T x){
if (x > 0) {return 1;}
if (x < 0) {return -1;}
return 0;
}
//
//Bubble sort+permutation array 
//Ascending order
template<typename T>
void SortAsc(T *MonTab, unsigned int *Permuter, unsigned int DimArray)
 {
 long ii, jj, tempi;
 T tempd;
for (ii = 0; ii < DimArray ; ii++)
{Permuter[ii]=ii;}

 for (ii = 0; ii < (DimArray-1) ; ii++)
 {
      for (jj = 0; jj < (DimArray - 1 - ii); jj++ )
      {
           if (MonTab[jj] > MonTab[jj+1])
          {
                tempd = MonTab[jj+1];
                MonTab[jj+1] = MonTab[jj];
                MonTab[jj] = tempd;
                tempi=Permuter[jj+1];
                Permuter[jj+1]=Permuter[jj];
                Permuter[jj]=tempi;                                
          }
      }
   }
}  

//Ascending order no permutation
template<typename T>
void SortAsc2(T *MonTab, long DimArray)
{
 long ii, jj;
 T tempd;

 for (ii = 0; ii < (DimArray-1) ; ii++)
 {
      for (jj = 0; jj < (DimArray - 1 - ii); jj++ )
      {
           if (MonTab[jj] > MonTab[jj+1])
           {
                tempd = MonTab[jj+1];
                MonTab[jj+1] = MonTab[jj];
                MonTab[jj] = tempd;                                
           }
      }
  }
}  

//Ascending order no permutation integer
template<typename T>
void SortAsc3(T *MonTab, int DimArray)
 {
 int ii, jj;
 T tempd;

 for (ii = 0; ii < (DimArray-1) ; ii++)
 {
      for (jj = 0; jj < (DimArray - 1 - ii); jj++ )
      {
           if (MonTab[jj] > MonTab[jj+1])
           {
                tempd = MonTab[jj+1];
                MonTab[jj+1] = MonTab[jj];
                MonTab[jj] = tempd;
                                
           }
      }
  }
}

//Sort in descending order and give the permutations
template<typename T>
void SortDsc(T *MonTab, int *Permuter, int DimArray)
 {
 int ii, jj, tempi;
 T tempd;

 for (ii = 0; ii < DimArray ; ii++){
   Permuter[ii]=ii;}

for (ii = 0; ii < (DimArray-1); ii++)
{
      for (jj = 0; jj < (DimArray - 1 - ii); jj++ )
      {
           if (MonTab[jj] < MonTab[jj+1])
           {
                tempd = MonTab[jj+1];
                MonTab[jj+1] = MonTab[jj];
                MonTab[jj] = tempd;
                tempi=Permuter[jj+1];
                Permuter[jj+1]=Permuter[jj];
                Permuter[jj]=tempi;                              
           }
      }
  }
} 

//Sort in descending order 
template<typename T>
void SortDsc2(T *MonTab, long DimArray)
{
 long ii, jj, tempi;
 T tempd;
 for (ii = 0; ii < (DimArray-1); ii++)
 {
      for (jj = 0; jj < (DimArray - 1 - ii); jj++ )
      {
           if (MonTab[jj] < MonTab[jj+1])
           {
                tempd = MonTab[jj+1];
                MonTab[jj+1] = MonTab[jj];
                MonTab[jj] = tempd;                              
           }
      }
  }
} 
//
template<typename T>
void SortAbsDsc(T *MonTab, long *Permuter, long DimArray)
 {
 long ii, jj, tempi;
 T tempd;
//
 for (ii = 0; ii < DimArray ; ii++){Permuter[ii]=ii;}
//
 for (ii = 0; ii < (DimArray-1); ii++)
 {
      for (jj = 0; jj < (DimArray - 1 - ii); jj++ )
      {
           if ( SIGN(MonTab[jj])*MonTab[jj] < SIGN(MonTab[jj+1])*(MonTab[jj+1]) )
           {
                tempd = MonTab[jj+1];
                MonTab[jj+1] = MonTab[jj];
                MonTab[jj] = tempd;
                tempi=Permuter[jj+1];
                Permuter[jj+1]=Permuter[jj];
                Permuter[jj]=tempi;                                
           }
      }
  }
} 
//
template<typename T>
long FindMaximum2(T *MonTab, long n) {
  long index, c;
  T max;
  max = pow(MonTab[0],2);
  index = 0;
// 
  for (c = 1; c < n; c++) {
    if (pow(MonTab[c],2) > max) {
       index = c;
       max = pow(MonTab[c],2);
    }
  }
  return index;
}
//
template<typename T>
long Find2Maximum2(T *MonTab, long n, unsigned int index1) {
  long index2, c;
  T max;
  max = -1;
  index2 = -1;
  for (c = 0; c < n; c++) {
    if ((pow(MonTab[c],2)) > max && (c !=index1)) {
       index2 = c;
       max = pow(MonTab[c],2);
    }
  }
// 
  return index2;
}
//
template<typename T>
T FindMaximumVal(T *MonTab, int n) {
  int index, c;
  T max; 
  max = pow(MonTab[0],2);
  index = 0;
// 
  for (c = 1; c < n; c++) {
    if (pow(MonTab[c],2) > max) {
       index = c;
       max = pow(MonTab[c],2);
    }
  }
  return MonTab[index];
}
//
template<typename T>
T Powo(T a, int n) {
T Val=1;
  for (int ii = 0; ii < n; ii++) {
  Val*=a;
  }
//
return Val;
}

template<typename T>
T FindMaximumValPoz(T *MonTab, int n) {
  int index, c;
  T max;
  max = pow(MonTab[0],2);
  index = 0;
// 
  for (c = 1; c < n; c++) {
    if (pow(MonTab[c],2) > max) {
       index = c;
       max = pow(MonTab[c],2);
    }
  }
// 
  return MonTab[index]*SIGN<T>(MonTab[index]);
}
//
template<typename T>
long FindMaximumAbs(T *MonTab, long n) {
  long index, c;
  T max; 
  max = SIGN<T>(MonTab[0])*MonTab[0];
  index = 0;
//
  for (c = 1; c < n; c++) {
    if (SIGN<T>(MonTab[c])*MonTab[c] > max) {
       index = c;
       max = SIGN<T>(MonTab[c])*MonTab[c];
    }
  }
// 
  return index;
}
//
template<typename T>
long FindMaximumInd(T *MonTab, long n) {
  long index, c;
  T max;
  max = MonTab[0];
  index = 0;
// 
  for (c = 1; c < n; c++) {
    if (MonTab[c] > max) {
       index = c;
       max = MonTab[c];
    }
  }
// 
  return index;
}
//
template<typename T>
int EqualArray(T *MonTab, T *b, int taille) {
//
  int mon_i=0;   
  while ((mon_i < taille) && (MonTab[mon_i] == b[mon_i]))
    { 
    mon_i++;
    }
//
 if(mon_i==taille)
{return 0;}
else{return 1;} 
//
}
//
template<typename T>
T sum_array(T MonTab[], int num_elements)
{
   int i; 
   T sum=0;
   for (i=0; i<num_elements; i++)
   {
	 sum = sum + MonTab[i];
   }
   return sum;
}
//
template<typename T>
T SumPM(T MonTab[], int num_elements)
{
   int i; 
   T sum=0;
   for (i=0; i<num_elements; i++)
   {
	 sum = sum + MonTab[i]*pow(-1,i);
   }
   return sum;
}
//
template<typename T>
T sum_array_abs(T MonTab[], int num_elements)
{
   int i; 
   T sum=0;
   for (i=0; i<num_elements; i++)
   {
	 sum = sum + SIGN<T>(MonTab[i])*MonTab[i];
   }
   return sum;
}
//
#endif
