#include"funciones.cpp"
#include"animate.h"
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include<random>
// alt-> M

std::string filename(int n);

int main(int argc, char **argv)
{
    Body Molecule[N];
    CRandom rand(42);
    double mu = 0, sigma = 0.001, sigma1 = 0.001*std::atoi(argv[1]); // The first argument by console is the seed.
    Vector3D move;
    double dt=0.01;
    double V=0.1*std::stod(argv[2]);  //The second is 10 times V where V is the voltage.
    double radio=0.001*std::atoi(argv[3]); //the third is 1000 times the radio of the particles.
    double x, y, z, vx, vy, vz, q0;
    //Declare potential and density array
    data_t potential(NX*NY);
    std:: random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0,1.0);
    int desbalance =  (1+std::atoi(argv[4])*0.01)*N/2;
    for (int i=0; i<N;i++)
      {
        // Initial positions in cubic lattice
        x = dis(gen);
        y = dis(gen);
        z = 0;

        // Initial nule velocities
        vx = rand.gauss(mu, sigma1);
        vy = rand.gauss(mu, sigma1);
        vz = 0;
	// k=sigma (Na-Nb)/N, Na -> + , Nb -> -.
	if(i<desbalance)
	  {
	    q0=-0.01;
	  }
	else
	  {
	    q0=0.01;
	  }
        Molecule[i].init(x, y, z, vx, vy, vz, m0, q0, false, radio);
      }

    //    start_animation(3);

    //Calculate initial potential
    initial_conditions(potential, NX, NY);
    /* 
	if(std::atoi(argv[6])==0){boundary_conditions(potential, NX, NY, Molecule, Lx, N, V);}
    if(std::atoi(argv[6])==1){boundary_conditions1(potential, NX, NY, Molecule, Lx, N, V);}
    if(std::atoi(argv[6])==2){boundary_conditions2(potential, NX, NY, Molecule, Lx, N, V);}
    if(std::atoi(argv[6])==3){boundary_conditions3(potential, NX, NY, Molecule, Lx, N, V);}
    */ 
    boundary_conditions_wnR(potential,NX,NY,Molecule,Lx,N,V,std::atoi(argv[6]),std::atoi(argv[7]));
    //evolve(potential, NX, NY, NSTEPS, NSTEPS);
    bool a=true;
    int count=0;
    while(a==true)
      { 
	a=relaxation_step(potential,NX,NY);
	  count+=1;	 
      }
     std::string filena="Fract_size"+std::to_string(std::atoi(argv[1]))+".txt";

    data_q distribution;
    for (int t = 0; t < 5000; t++)
    {
           
      /*if (t % 2 == 0)
      {
	  print_potential_size(NX, NY, potential, filena, t);
	  	 
	     begin_frame(argc);
	     print_potential_gnuplot(NX,NY,potential);
             for (int k = 0; k < N; k++)
	       Molecule[k].print();
	     end_frame(argc);
	  
        }
      */
      PEFRL(Molecule, potential, NX, NY, Lx, N, mu, sigma, dt, t+std::atoi(argv[5]), V,true);	
      /*
      if(check_fractal(Molecule,NX,N,std::atoi(argv[4])))
	{
	  distribution.push_back(Probability_distribution(Molecule,N));
	  break;
	}
      */
      distribution.push_back(Probability_distribution(Molecule,N));
    }
    std::string filename="data/Out"+std::to_string(std::atoi(argv[1]))+"T"+std::to_string(std::atoi(argv[2]))+"V"+std::to_string(std::atoi(argv[3]))+"R"+std::to_string(std::atoi(argv[4]))+"I"+std::to_string(std::atoi(argv[5]))+"S"+std::to_string(std::atoi(argv[6]))+"A"+std::to_string(std::atoi(argv[7]))+"F"+".txt";
    print_fractal(NX,NY, potential, filename);
    std::string namedis=std::to_string(std::atoi(argv[1]))+"T"+std::to_string(std::atoi(argv[2]))+"V"+std::to_string(std::atoi(argv[3]))+"R"+std::to_string(std::atoi(argv[4]))+"I"+std::to_string(std::atoi(argv[5]))+"S"+std::to_string(std::atoi(argv[6]))+"A"+std::to_string(std::atoi(argv[7]))+"F";
    print_fractal(NX,NY, potential, filename);
    print(distribution,namedis+"Probability_distribution.txt");
    //print_gnuplot(potential, NX, NY);
      
return 0;
}

std::string filename(int n)
{
    std::string name;
    std::stringstream n_s;
    n_s << n;
    name = "./" + n_s.str() + ".csv";
    return name;
}
