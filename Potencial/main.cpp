#include"funciones.cpp"
int main(int argc, char **argv)
{
    // declare data structures
    double x, y, z, vx, vy, vz, q0, m0=1;
    data_t potential(NX*NY);
    Body N[Nmax];
    for (int i=0; i<Nmax;i++)
      {
	if(i<Nmax/2)
	  {
          x=0.7;
          y=0.7;
          z=0;
          q0=1;
          vx=0;
          vy=0;
          vz= 0;
          N[i].init(x, y, z, vx, vy, vz, m0, q0);
	  }
	else
	  {
	x=0.5;
	y=0.5;
	q0=-1;
	vx=0;
	vy=0;
    N[i].init(x, y, z, vx, vy, vz, m0, q0);
	  }
      }
    data_q Q(NX*NY,{0});
    // set initial and boundary conditions
    Get_Q(N,Q,NX,NY,1.2, Nmax);
    initial_conditions(potential, NX, NY);
    boundary_conditions(potential,NX, NY, N ,1.2, Nmax);

    // evolve and print
    
    
    evolve(potential, NX, NY, NSTEPS,Q);
    return 0;
}
