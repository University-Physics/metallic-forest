#include"funciones.cpp"
#include"animate.h"
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>


std::string filename(int n);

int main(int argc, char **argv)
{
    Body Molecule[N];
    CRandom rand(42);
    double mu = 0, sigma = 0.001, sigma1 = 0.001*std::atoi(argv[1]); // The first argument by console is the seed.
    double move_x, move_y;
    Vector3D move;
    double tdibujo = 0;
    double dt=0.01;
    double V=0.1*std::atoi(argv[2]);  //The second is 10 times V where V is the voltage.
    double radio=0.01*std::atoi(argv[3]); //the third is 100 times the radio of the particles.
    double x, y, z, vx, vy, vz, q0, x0 = 0.25, y0 = 0.25;

    double dx = Lx /(2*(Nx + 1));
    double dy = Ly /(2*(Ny + 1));
    //Declare potential and density array
    data_t potential(NX*NY);

    for (int i=0; i<N;i++)
      {
        // Initial positions in cubic lattice
        x = dx + (i % Nx) * dx + x0;
        y = dy + ((i / Nx) % Ny) * dy + y0;
        z = 0;

        // Initial nule velocities
        vx = rand.gauss(mu, sigma1);
        vy = rand.gauss(mu, sigma1);
        vz = 0;

	q0=-0.01;
	if(i%std::atoi(argv[4])==0) q0*=-1; // ratio of population between charges +q0 and -q0: (Number of q0)/(Number of -q0)= (argv[4]-1) if argv[4]|N


        Molecule[i].init(x, y, z, vx, vy, vz, m0, q0, false, radio);
	//Molecule[i].print();
      }
    //start_animation(argc);

    //Calculate initial potential
    initial_conditions(potential, NX, NY);
    boundary_conditions(potential, NX, NY, Molecule, Lx, N, V);
    evolve(potential, NX, NY, NSTEPS, NSTEPS);
    //std::string filename="Fract_size"+std::to_string(std::atoi(argv[1]))+".txt";
    
    for (int t = 0; t < 5000; t++)
    {
      /* 
      if (t % 200 == 0)
        {
	  //print_potential_size(NX, NY, potential, filename, t);
	 
	     begin_frame(argc);
             for (int k = 0; k < N; k++)
                Molecule[k].print();
	     end_frame(argc);
	 
        }
      */
      PEFRL(Molecule, potential, NX, NY, Lx, N, mu, sigma, dt, t+std::atoi(argv[5]), V);
    }
    std::string filename="data/Out"+std::to_string(std::atoi(argv[1]))+"T"+std::to_string(std::atoi(argv[2]))+"V"+std::to_string(std::atoi(argv[3]))+"R"+std::to_string(std::atoi(argv[4]))+"I"+std::to_string(std::atoi(argv[5]))+"S.txt";
    print_fractal(NX,NY, potential, filename);
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
