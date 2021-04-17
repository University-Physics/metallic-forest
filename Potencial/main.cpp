
#include"constantes.h"
#include"encabezado.h"
int main(int argc, char **argv)
{
    // declare data structures
    data_t potential(NX*NY);
    data_p N(1000);
    for (int i=0; i<N.size();i++)
      {
	if(i<N.size()/2)
	  {
	N[i].x=0.7;
	N[i].y=0.7;
	N[i].q=1;
	N[i].Vx=0;
	N[i].Vy=0;
	  }
	else
	  {
	N[i].x=0.5;
	N[i].y=0.5;
	N[i].q=-1;
	N[i].Vx=0;
	N[i].Vy=0;
	  }
      }
    data_q Q(NX*NY,{0});
    // set initial and boundary conditions
    Get_Q(N,Q,NX,NY,1.2);
    initial_conditions(potential, NX, NY);
    boundary_conditions(potential,NX, NY, N ,1.2);

    // evolve and print
    
    
    evolve(potential, NX, NY, NSTEPS,Q);
    return 0;
}
