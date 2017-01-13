#ifndef BITONIC_H
#define BITONIC_H
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

template <typename T>
void Bitonic_Compare(int low, int high, std::vector<T> &A, bool upSort)
{
  if((A[low] > A[high]) == upSort)
    swap(A[low], A[high]);
  return;
}

int Pow2_LT(int high)
{
  int k = 1;
  while(k < high)
    k = k<<1;
  return k>>1;
}

template <typename T>
void Bitonic_Merge(int low, int high, std::vector<T> &A, bool upSort)
{
  if(high > 1)
    {
      int m = Pow2_LT(high);
      for(int i = low; i < low + high - m; ++i)
	Bitonic_Compare(i, i + m, A, upSort);
      Bitonic_Merge(low, m, A, upSort);
      Bitonic_Merge(low + m, high - m, A, upSort);
    }
  return;
}

template <typename T>
void Bitonic_Sort(int low, int high, std::vector<T> &A, bool upSort)
{
  if(high > 1)
    {
      int midPoint = high/2;
      Bitonic_Sort(low, midPoint, A, !upSort);
      Bitonic_Sort(low + midPoint, high - midPoint, A, upSort);
      Bitonic_Merge(low, high, A, upSort);
    }
  return;
}

template <typename T>
void BitonicSort(std::vector<T> &A)
{
  Bitonic_Sort(0, A.size(), A, true);
  return;
}

#endif
