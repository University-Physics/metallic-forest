#include"funciones.cpp"
#include"animate.h"
#include <string>
#include <sstream>
#include <fstream>


std::string filename(int n);

int main(int argc, char **argv)
{
    Body Molecule[N];
    CRandom rand(42);
    double mu = 0, sigma = 0.1;
    double move_x, move_y;
    Vector3D move;
    double tdibujo = 0;
    double dt=0.1;

    double x, y, z, vx, vy, vz, q0, x0 = 0, y0 = 0;

    double dx = Lx / (Nx + 1);
    double dy = Ly / (Ny + 1);
    //Declare potential and density array
    data_t potential(NX*NY);
    data_q Q(NX*NY,{0});
    for (int i=0; i<Nmax;i++)
      {
        // Initial positions in cubic lattice
        x = dx + (i % Nx) * dx + x0;
        y = dy + ((i / Nx) % Ny) * dy + y0;
        z = 0;

        // Initial nule velocities
        vx = 0;
        vy = 0;
        vz = 0;

	q0=-1;
	if(i%2==0) q0=1;

        Molecule[i].init(x, y, z, vx, vy, vz, m0, q0);
	//Molecule[i].print();
      }
    //start_animation(argc);
    
    for (int t = 0; t < 1000; t++)
    {
      /*if (t % 2 == 0)
        {
	     begin_frame(argc);
             for (int k = 0; k < N; k++)
                Molecule[k].print();
	     end_frame(argc);
        }
      */
	// set initial and boundary conditions
	Get_Q(Molecule,Q,NX,NY,Lx, Nmax);
	initial_conditions(potential, NX, NY);
	boundary_conditions(potential,NX, NY, Molecule ,Lx, Nmax);

	// evolve
	evolve(potential, NX, NY, NSTEPS,Q);
        update_and_check_pos(Molecule, Nx, Ny, Lx, Nmax, potential, mu, sigma, dt);
	Update_boundary(Molecule, NX, NY, Lx, Nmax, potential);
    }
    print_fractal(NX,NY, potential);
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
