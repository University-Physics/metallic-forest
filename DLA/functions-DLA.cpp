#include<vector>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<random>
#include<fstream>
#include<string>
std::vector<int> vecinos(int Na, std::vector <int> &a)
{
  //Esta función da información acerca de los vecinos de la casilla Na
  /*Si el sitio vecino es un lugar prohibido se colocara el número 2
  Un lugar prohibido es aquel sitio que se sale de la grilla. 
  */
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

void recorrido( int Na, std::vector<int> &a) //Esta función indica el algoritmo seguido para mover aleatoriamente la partícula
{
while(true)
    {
     std::vector<int> b=vecinos(Na,a); // Se hallan sus vecinos
     if(find(begin(b),end(b),1) != std::end(b)){ a[Na]=1; break;} /* si uno de sus vecinos es una casilla ocupada (es decir si una
     casilla ocupada tiene un valor de 1), entonces se le asigna a la casilla a[Na] un valor de ocupada */ 
     int random=rand() % 4; //Se genera un número aleatorio
     if( Na+casilla(random)<a.size() && Na+casilla(random)>0 && b[random]!=2) {Na=Na+casilla(random);} /* La casilla a la que la partícula se movera debe cumplir unas condiciones, debe estar en el arreglo y no debe ser prohibida, es decir no se puede pasar del lado derecho del arreglo al izquierdo en un solo paso y viceversa.*/
     else {break;}
    }
}


void print(std::vector <int> &a, std::string &b)
{
  std::ofstream fout(b);
  for(int i=0;i<110;i++)
    {
      for(int ii=0;ii<110;ii++){fout<<i<<" "<<ii<<" "<<a[i*110+ii]<<" "<<std::endl;}
      //Esta función imprime la ubicación de la casilla [i][ii] (i es la fila, ii la columna) y su valor (si esta ocupada es 1 y vacia 0) 
    }
  fout.close();

}
