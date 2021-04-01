#include<iostream>
#include<vector>
#include<functional>
#include<cmath>
#include<random>
#include<algorithm>
#include <cstdlib>
std::vector<int> vecinos(int Na, std::vector <int> &a);
int casilla(int N);
void print(std::vector <int> &a);
int main()
{
  std::vector<int> a(12100,0);
  a[6105]=1;
  const int m=sqrt(a.size());
  for(int particula=0;particula<10000;particula++)
    {
   
 int Na=rand() %12100;
      while(Na==60){  Na=rand() % 12100; }
  while(true)
    {
    std::vector<int> b=vecinos(Na,a);
    if(find(begin(b),end(b),1) != std::end(b)){ a[Na]=1; break;}
    /*Numero aleatorio de 0 a 4*/
    int random=rand() % 4;
    if( Na+casilla(random)<a.size() && Na+casilla(random)>0 && b[random]!=2) Na=Na+casilla(random);
    else break;
    }
    }
  print(a);
    return 0;
  }


std::vector<int> vecinos(int Na, std::vector <int> &a)
{
  //Esta función da información acerca de los vecinos de la casilla Na
  int m=sqrt(a.size());
  std::vector<int> veci;
  if(Na+m<a.size()) veci.push_back(a[Na+m]); // Se refiere al vecino del sitio de abajo
  else veci.push_back(2);
  if(Na-m>=0) veci.push_back(a[Na-m]); // se refiere al vecino del sitio de arriba
  else veci.push_back(2);
  if(Na+1<a.size() && (Na+1)%110!=0) veci.push_back(a[Na+1]); // Se refiere al vecino del sitio de la derecha
  else veci.push_back(2);
  if(Na-1>=0 && (Na)%110!=0) veci.push_back(a[Na-1]); // Se refiere al vecino del sitio de la izquierda
  else veci.push_back(2);
  return veci;
}
int casilla(int N)
{
  //Esta función retorna las casillas que debe recorrer la partícula en el vector, para ir a los sitios vecinos
  int m=110;  
  if(N==0) return m; // Número de sitios de casillas para ir al sitio de abajo
  else if(N==1) return (-m); // Número de sitios de casillas para ir al sitio de arriba
  else if(N==2) return 1; // Número de sitios de casillas para ir al sitio de la derecha
  else return -1; // Número de sitios de casillas para ir al sitio de la izquierda
}
void print(std::vector <int> &a)
{
  for(int i=0;i<110;i++)
    {
      for(int ii=0;ii<110;ii++)
	{
	  std::cout<<i<<" "<<ii<<" "<<a[i*110+ii]<<" "<<std::endl;
	}
     
    }

}
