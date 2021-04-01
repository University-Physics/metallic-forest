#include<iostream>
#include<cstring>
#include<vector>
#include<functional>
#include<cmath>
#include<fstream>
#include<random>
#include<algorithm>
#include <cstdlib>
#include "functions-DLA.h"

int main()
{
  std::vector<int> a(12100,0);
  a[6105]=1;
  const int m=sqrt(a.size());
  /*for(int particula=0;particula<10000;particula++)
    {
      int Na=rand() %12100;
      while(Na==60){  Na=rand() % 12100; }
      recorrido(Na,a);
    }
  std::string b="data.txt";
  print(a,b);
  */
  a.assign(12100,0);
  for(int i=0;i<110;i++) //En esta parte del código se inicializa, la frontera del electrodo
    {
      a[109*110+i]=1;
    }
  for(int particula=0;particula<172500;particula++)
    {
      int Na=(rand() %110); //Las partículas se inicializan al lado izquierdo deel arreglo
      recorrido(Na,a);
    }
  std::string b="data.txt";
  print(a,b);
  
  return 0;
}
