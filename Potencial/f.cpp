#include<iostream>
int main(int argc, char **argv)
{
  double a=std::stod(argv[1]);
  std::cout<<a<<" "<<std::to_string(a)<<std::endl;
  return 0;
}
