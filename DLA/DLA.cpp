#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include <cstdlib>
std::vector<int> vecinos(int Na, std::vector <int> &a);
int casilla(int N, std::vector <int> &b);
void print(std::vector <int> &a);
int main()
{
  std::vector<int> a(12100,0);
  a[60]=1;
  const int m=sqrt(a.size());
  for(int particula=0;particula<10000;particula++)
    {
      int Na=rand() % 12100;
      while(Na==60){  Na=rand() % 12100; }
  while(true)
    {
    std::vector<int> b=vecinos(Na,a);
    if(find(begin(b),end(b),1) != std::end(b)){ a[Na]=1; break;}
    /*Numero aleatorio de 0 a 4*/
    int random=rand() % 5; 
    if( Na+casilla(rand() % 5,b)<a.size() && Na+casilla(rand() % 5,b)>0 && a[Na+casilla(rand() % 5,b)]!=1) Na=Na+casilla(rand() % 5,b);
    else break;
    }
    }
  print(a);
    return 0;
  }


std::vector<int> vecinos(int Na, std::vector <int> &a)
{
  int m=sqrt(a.size());
  std::vector<int> veci;
  if(Na+m<a.size()) veci.push_back(a[Na+m]);
  else veci.push_back(2);
  if(Na-m>=0) veci.push_back(a[Na-m]);
  else veci.push_back(2);
  if(Na+1<a.size()) veci.push_back(a[Na+1]);
  else veci.push_back(2);
  if(Na-1>=0) veci.push_back(a[Na-1]);
  else veci.push_back(2);
  return veci;
}
int casilla(int N, std::vector <int> &b)
{
  int m=110;  
  if(N==0) return m;
  else if(N==1) return -m;
  else if(N==2) return 1;
  else return -1;
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
