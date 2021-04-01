#include<iostream>
#include<vector>
#include<functional>
#include<cmath>
#include<random>
#include<algorithm>
#include <cstdlib>
#include "functions-DLA.h"

int main()
{
  std::vector<int> a(12100,0);
  a[6105]=1;
  const int m=sqrt(a.size());
  for(int particula=0;particula<10000;particula++)
    {
      int Na=rand() %12100;
      while(Na==60){  Na=rand() % 12100; }
      recorrido(Na,a);
    }
  print(a);
  return 0;
}
