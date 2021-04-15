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
#include "constantes.h"
int main()
{
  std::vector<int> a(largo*ancho,0);
  a[largo*ancho/2]=1;
  /*for(int particula=0;particula<10000;particula++)
    {
      int Na=rand() %12100;
      while(Na==60){  Na=rand() % 12100; }
      recorrido(Na,a);
    }
  std::string b="data.txt";
  print(a,b);
  */
  a.assign(largo*ancho,0);
  for(int i=0;i<largo;i++) //En esta parte del código se inicializa, la frontera del electrodo
    {
      a[(ancho-1)*largo+i]=1;
    }
  for(int particula=0;particula<numberparticles;particula++)
    {
      
      int y=(rand() %(largo*ancho)); //Las partículas se inicializan al lado izquierdo deel arreglo.
      recorrido(y,a);
    }
  std::string b="data.txt";
  print(a,b);
  
  return 0;
}
